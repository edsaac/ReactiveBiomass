/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    unsatNutrientCycle

Description
    Microbial growth is substrate-limited and electron acceptor-limited via
    dual Monod.

    Generation and consumption of EPS, UAP and BAP are included.
    Only advection of dissolved and particulate species is considered

    Clogging is introduced some porosity-permeability relationship
    
    k/k0 = f(n/n0)

    Solves eight species:
    - Aerobic heterotrophs (XAR) {X}
    - Inert biomass  (XI) {Xi}
    - Extrapolymeric substances (EPS) {E}
    - Biomass-associated products (BAP) {B}
    - Substrate (DOC) {S}
    - Oxygen (O2) {O}

    Eqs:

IGNORE
    d(nS)/dt  = - rH*X - rDN*Xdn
    d(nN)/dt  = - rN*Xn + gammaN*fd*(b*X + bN*Xn + bDN*Xdn) - gammaN*Y'*rH*X

    dE/dt  = kE*(rH*X + rN*XN + rDN*Xdn) - khyd*E
    d(nU)/dt  = k1*(rH*X + rN*XN + rDN*Xdn) - rU*X - rUDN*Xdn
    d(nB)/dt  = khyd*E - rB*X - rBDN*Xdn

    dX/dt   = (Y'*rH + YP*(rU+rB) - b)*X
    dXn/dt  = (Yn'*rN - bN)*Xn
    dXdn/dt = (Ydn'*rDN + YP*(rUDN+rBDN)- bDN)*Xdn
    dXi/dt  = (1-fd)*(b*X + bN*XN + bDN*XDN)

    d(nO2)/dt = -alpha1*rH*X - alpha1P*(rU+rB)*X - alphaN*rN*Xn - 1.42*fd*(b*X+bN*Xn)
    d(nN3)/dt = (1 - gammaN*YN)*rN*Xn - beta1*rDN*Xdn - beta1P*(rUDN+rBDN)*Xdn - 0.69*fd*bDN*Xdn

    with:
    n : Porosity = (volVoids/volREV)
    M(C) = C/(KC+C)
    I(O2) = KI/(KI+O2)

    rH = q*M(S)*M(O2)
    rU = qUAP*M(U)*M(O2)
    rB = qBAP*M(B)*M(O2)

    rN = qN*M(N)*M(O2)

    rDN = qDN*M(S)*M(NO3)*I(O2)
    rUDN = qUAP*M(U)*M(NO3)*I(O2)
    rBDN = qBAP*M(B)*M(NO3)*I(O2)

    Y' = Y*(1 - k1 - kE)  ## Code assumes Y is Y'
ENDIGNORE
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

#include "cloggingModel.H"
#include "attachmentModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "createPhi.H"
    #include "CourantNo.H"
    #include "createFvOptions.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Foam::Info << "\nCalculating...\n" << endl;

    while (simple.loop(runTime))
    {
        Foam::Info<< "Time = " << runTime.timeName() << nl << endl;

        //Info << "\nUpdate clog space limitation" << endl;
        totalBiomass = XAR + XN + XDN + XI + EPS;
        clogLimiter = 1.0 - totalBiomass/XMAX;

        if (Foam::min(clogLimiter) < dimensionedScalar("zero",dimless,0.0))
        {
            Foam::SeriousError<< "Total biomass greater that XMAX" << endl;            
            runTime.write();
            break;
        }

        n  = clogging->nRef() - totalBiomass/rho_X;

        //Calculate hydraulic head using mass balance + Darcy's equation
        if (cloggingSwitch) { clogging->calcPerm();}
        hydraulicCond = perm * g * rho / mu;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix FlowEquation
            (
                fvm::laplacian(hydraulicCond, h)
            );
            fvOptions.constrain(FlowEquation);
            FlowEquation.solve();
            fvOptions.correct(h);
        }

        // Update flow field
        Foam::Info << "\n Update porosity field (n):" << endl;
        U   = - hydraulicCond * fvc::grad(h);
        phi = fvc::flux(U);
        #include "CourantNo.H"

        // Calculate attach/detach rates
        attachment->calcAttachment();
        detachment->calcAttachment();

        
        // Transport equations
        while (simple.correctNonOrthogonal())
        {
            // Info << "\nRates of utilization (...) " << endl;
            rH = 
                pXAR.qhat  
               * DOC/(pXAR.Kdonor + DOC)
               * O2/(pXAR.Kaccep + O2)
               * NH4/(K_minN + NH4)
               * clogLimiter;

            rN = 
                pXN.qhat 
               * NH4/(pXN.Kdonor + NH4) 
               * O2/(pXN.Kaccep + O2)
               * clogLimiter;

            rDN = 
                pXDN.qhat
                * DOC/(pXDN.Kdonor + DOC)
                * NO3/(pXDN.Kaccep + NO3)
                * K_I / (K_I + O2)
                * clogLimiter;
            
            //Info << "\nImmobile aerobic heterotrophs (XAR)" << endl;
            fvScalarMatrix ARGrowth
            (
                fvm::ddt(XAR)
                ==
                pXAR.yield * rH * XAR                   //Growth
                - pXAR.bDie * XAR                       //Decay XAR -> EPS + XI
                - kdet * XAR                            //Detachment XAR -> XARp
                + (Sw * n) * katt * clogLimiter * XARp         //Attachment XARp -> XAR
            );
            ARGrowth.relax();
            fvOptions.constrain(ARGrowth);
            ARGrowth.solve();
            fvOptions.correct(XAR);

            //Info << "\nMobile aerobic heterotrophs (XARp)" << endl;
            fvScalarMatrix ARGrowth_p
            (
                n * Sw * fvm::ddt(XARp)
                // + XARp * fvc::ddt(Sw,n)
                + fvm::div(phi, XARp)
                - fvm::laplacian(Sw * n * (mag(U)*DispTensor + molDiff), XARp)
                ==
                Sw * n * pXAR.yield * rH * XARp              //Growth
                - Sw * n * pXAR.bDie * XARp                  //Decay XARp -> BAP + POCr
                + kdet * XAR                                 //Detachment XAR -> XARp
                - Sw* n * katt * clogLimiter * XARp          //Attachment XARp -> XAR
            );
            ARGrowth_p.relax();
            fvOptions.constrain(ARGrowth_p);
            ARGrowth_p.solve();
            fvOptions.correct(XARp);

            // Info << "\nImmobile Nitrifiers (XN)" << endl;
            fvScalarMatrix NitrifiersGrowth
            (
                fvm::ddt(XN)
                ==
                  pXN.yield * rN * XN                        //Growth
                - pXN.bDie * XN                              //Decay XN -> EPS + XI
                - kdet * XN                                  //Detachment XN -> XNp
                + Sw * n * katt * clogLimiter * XNp          //Attachment XNp -> XN
            );
            NitrifiersGrowth.relax();
            fvOptions.constrain(NitrifiersGrowth);
            NitrifiersGrowth.solve();
            fvOptions.correct(XN);

            // Info << "\nMobile Nitrifiers (XN)" << endl;
            fvScalarMatrix NitrifiersGrowth_p
            (
                n * Sw * fvm::ddt(XNp)
                // + XNp * fvc::ddt(Sw,n)
                + fvm::div(phi, XNp)
                - fvm::laplacian(Sw * n * (mag(U)*DispTensor + molDiff), XNp)
                ==
                Sw * n * pXN.yield * rN * XNp                    //Growth
                - Sw * n * pXN.bDie * XNp                        //Decay XNp -> BAP + POCr
                + kdet * XN                                      //Detachment XN -> XNp
                - Sw * n * katt * clogLimiter * XNp              //Attachment XNp -> XN
            );
            NitrifiersGrowth_p.relax();
            fvOptions.constrain(NitrifiersGrowth_p);
            NitrifiersGrowth_p.solve();
            fvOptions.correct(XNp);

            // Info << "\nImmobile Denitrifiers (XDN)" << endl;
            fvScalarMatrix DenitrifiersGrowth
            (
                fvm::ddt(XDN)
                ==
                pXDN.yield * rDN * XDN                       //Growth
                - pXDN.bDie * XDN                            //Decay XDN -> EPS + XI
                - kdet * XDN                                 //Detachment XDN -> XDNp
                + Sw * n * katt * clogLimiter * XDNp         //Attachment XDNp -> XDN
            );
            DenitrifiersGrowth.relax();
            fvOptions.constrain(DenitrifiersGrowth);
            DenitrifiersGrowth.solve();
            fvOptions.correct(XDN);

            // Info << "\nMobile Denitrifiers (XDN)" << endl;
            fvScalarMatrix DenitrifiersGrowth_p
            (
                n * Sw * fvm::ddt(XDNp)
                // + XDNp * fvc::ddt(Sw,n)
                + fvm::div(phi, XDNp)
                - fvm::laplacian(Sw * n*(mag(U)*DispTensor + molDiff), XDNp)
                ==
                Sw * n * pXDN.yield * rDN * XDNp                 //Growth
                - Sw * n * pXDN.bDie * XDNp                      //Decay XDN -> BAP + POCr
                + kdet * XDN                                     //Detachment XDN -> XDNp
                - Sw * n * katt * clogLimiter * XDNp             //Attachment XDNp -> XDN
            );
            DenitrifiersGrowth_p.relax();
            fvOptions.constrain(DenitrifiersGrowth_p);
            DenitrifiersGrowth_p.solve();
            fvOptions.correct(XDNp);

	        // Minimum viable biomass
            XAR *= pos(XAR - X_min);
            XN  *= pos(XN  - X_min);
            XDN *= pos(XDN - X_min);

            //Info << "\nExtrapolymeric substances EPS transport" << endl;
            fvScalarMatrix EPSTransport
            (
                fvm::ddt(EPS)
                ==
                  (kEPS*rH  + fd*pXAR.bDie) * XAR       //Metabolism + decay XAR
                + (kEPS*rN  + fd*pXN.bDie)  * XN        //Metabolism + decay XN
                + (kEPS*rDN + fd*pXDN.bDie) * XDN       //Metabolism + decay XDN
                - kdet * EPS                            //Detachment EPS -> BAP
                + Sw * n * katt * clogLimiter * BAP     //Attachment BAP -> EPS
                - khyd_labil * EPS                      //Hydrolysis EPS -> DOC
            );
            EPSTransport.relax();
            fvOptions.constrain(EPSTransport);
            EPSTransport.solve();
            fvOptions.correct(EPS);

	        //Info << "\nInert biomass (XI)" << endl;
            fvScalarMatrix InertGeneration
            (
                fvm::ddt(XI)
                ==
                (1-fd) * (
                    pXAR.bDie*XAR                       //Decay XAR
                    + pXN.bDie*XN                       //Decay XN
                    + pXDN.bDie*XDN )                   //Decay XDN
                - kdet * XI                             //Detachment XI -> POCr
                + Sw * n * katt * clogLimiter * POCr    //Attachment POCr -> XI
                - khyd_recal * XI                       //Hydrolysis of XI -> DOC
            );
            InertGeneration.relax();
            fvOptions.constrain(InertGeneration);
            InertGeneration.solve();
            fvOptions.correct(XI);

            //Info << "\nDissolved Organic Carbon DOC transport" << endl;
            fvScalarMatrix DissolvedCarbonTransport
            (
                n * Sw * fvm::ddt(DOC)
                // + DOC * fvc::ddt(Sw,n)
                + fvm::div(phi, DOC)
                - fvm::laplacian(Sw * n*(mag(U)*DispTensor + molDiff), DOC)
                ==
                - rH * XAR                              //Metabolism XAR
                - rDN * XDN                             //Metabolism XDN
                + Sw * n * (
                    khyd_labil * BAP                    //Hydrolisis BAP -> DOC
                    + khyd_recal * POCr                 //Hydrolisis POCr -> DOC
                    + khyd_labil * EPS                  //Hydrolysis of EPS -> DOC
                    + khyd_recal * XI                   //Hydrolysis of XI -> DOC
                )
            );
            DissolvedCarbonTransport.relax();
            fvOptions.constrain(DissolvedCarbonTransport);
            DissolvedCarbonTransport.solve();
            fvOptions.correct(DOC);
            
            //Info << "\nDissolved oxygen (O2)" << endl;
            fvScalarMatrix OxygenTransport
            (
                Sw * n * fvm::ddt(O2)
                // + O2 * fvc::ddt(Sw,n)
                + fvm::div(phi, O2)
                - fvm::laplacian(Sw * n*(mag(U)*DispTensor + molDiff), O2)
                ==
                - alpha_1 * rH * (XAR + Sw*n*XARp)   //Metabolism XAR
                - alpha_N * rN * (XN + Sw*n*XNp)     //Metabolism XN
            );
            OxygenTransport.relax();
            fvOptions.constrain(OxygenTransport);
            OxygenTransport.solve();
            fvOptions.correct(O2);

            // Info << "\nNH4-N transport" << endl;
            fvScalarMatrix AmmoniumTransport
            (
                n * Sw * fvm::ddt(NH4)
                // + NH4 * fvc::ddt(Sw,n)
                + fvm::div(phi, NH4)
                - fvm::laplacian(Sw*n*(mag(U)*DispTensor + molDiff), NH4)
                ==
                - rN * XN                               //Metabolism XN
                + gamma_N * fd * (
                    pXAR.bDie*XAR                       //Decay XAR
                    + pXN.bDie*XN                       //Decay XN
                    + pXDN.bDie*XDN                     //Decay XDN
                    )
                + gamma_N * fd * Sw * n * (                      
                    pXAR.bDie*XARp                      //Decay XARp
                    + pXN.bDie*XNp                      //Decay XNp
                    + pXDN.bDie*XDNp                    //Decay XDNp
                    )                   
                - gamma_N * pXAR.yield * rH * XAR            //Metabolism XAR
                - Sw * n * gamma_N * pXAR.yield * rH * XARp  //Metabolism XARp
            );
            AmmoniumTransport.relax();
            fvOptions.constrain(AmmoniumTransport);
            AmmoniumTransport.solve();
            fvOptions.correct(NH4);

            // Info << "\nDissolved nitrate-nitrogen (NO3+N)" << endl;
            fvScalarMatrix NitrateTransport
            (
                n * Sw * fvm::ddt(NO3)
                // + NO3 * fvc::ddt(Sw,n) 
                + fvm::div(phi, NO3) 
                - fvm::laplacian(Sw*n*(mag(U)*DispTensor + molDiff), NO3)
                ==
                (1.0 - 0.31*pXN.yield - kEPS) * rN * (XN + Sw*n*XNp)   //Metabolism XN
                - beta_1 * rDN * (XDN + Sw*n*XDNp)                     //Metabolism XDN
            );
            NitrateTransport.relax();
            fvOptions.constrain(NitrateTransport);
            NitrateTransport.solve();
            fvOptions.correct(NO3);

	        //Info << "\nBiomass-associated products (BAP) transport" << endl;
            fvScalarMatrix BAPTransport
            (
                n * Sw * fvm::ddt(BAP)
                // + BAP * fvc::ddt(Sw,n)
                + fvm::div(phi, BAP)
                - fvm::laplacian(Sw*n*(mag(U)*DispTensor + molDiff), BAP)
                ==
                kdet * EPS                              //Detachment EPS -> BAP
                - Sw * n * katt * clogLimiter * BAP     //Attachment BAP -> EPS
                - Sw * n * khyd_labil * BAP             //Hydrolisis BAP -> DOC
                + Sw * n * fd * (
                      pXAR.bDie * XARp                  //Decay XARp
                    + pXN.bDie  * XNp                   //Decay XNp
                    + pXDN.bDie * XDNp                  //Decay XDNp
                )
            );
            BAPTransport.relax();
            fvOptions.constrain(BAPTransport);
            BAPTransport.solve();
	        fvOptions.correct(BAP);

            //Info << "\nInert biomass (XI)" << endl;
            fvScalarMatrix POCrGeneration
            (
                n * Sw * fvm::ddt(POCr)
                // + POCr * fvc::ddt(Sw,n)
                + fvm::div(phi, POCr)
                - fvm::laplacian(Sw * n*(mag(U)*DispTensor + molDiff), POCr)
                ==
                kdet * XI                       //Detachment XI -> POCr
                - Sw * n * POCr * (
                    katt * clogLimiter          //Attachment POCr -> XI
                    + khyd_recal                //Hydrolisis POCr -> DOC
                    )                 
                + Sw * n * (1-fd) * (
                      pXAR.bDie * XARp          //Decay XARp
                    + pXN.bDie  * XNp           //Decay XNp
                    + pXDN.bDie * XDNp          //Decay XDNp
                )
            );
            POCrGeneration.relax();
            fvOptions.constrain(POCrGeneration);
            POCrGeneration.solve();
            fvOptions.correct(POCr);
        }

        //End bits
        runTime.write();

        Foam::Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
    }

    Foam::Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
