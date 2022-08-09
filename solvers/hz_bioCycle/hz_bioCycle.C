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
    hz_bioFoam

Description
    Microbial growth is substrate-limited and electron acceptor-limited via
    dual Monod.
    Generation and consumption of EPS, UAP and BAP are included.
    Only advection of dissolved species is considered

    Clogging is introduced with Kozeny-Carman relationship affecting permeability
    k = max(k0*(n/no)^3*((1-n0)/(1-n))^2, kmin*k0)

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
        Info<< "Time = " << runTime.timeName() << nl << endl;

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
            
            //Info << "\nAerobic heterotrophs (XAR)" << endl;
            fvScalarMatrix ARGrowth
            (
                fvm::ddt(XAR)
                ==
                (pXAR.yield * rH - pXAR.bDie) * XAR
            );
            ARGrowth.relax();
            fvOptions.constrain(ARGrowth);
            ARGrowth.solve();
            fvOptions.correct(XAR);

            // Info << "\nNitrifiers (XN)" << endl;
            fvScalarMatrix NitrifiersGrowth
            (
                fvm::ddt(XN)
                ==
                (pXN.yield * rN - pXN.bDie) * XN
            );
            NitrifiersGrowth.relax();
            fvOptions.constrain(NitrifiersGrowth);
            NitrifiersGrowth.solve();
            fvOptions.correct(XN);

            // Info << "\nDenitrifiers (XDN)" << endl;
            fvScalarMatrix DenitrifiersGrowth
            (
                fvm::ddt(XDN)
                ==
                (pXDN.yield * rDN - pXDN.bDie) * XDN
            );
            DenitrifiersGrowth.relax();
            fvOptions.constrain(DenitrifiersGrowth);
            DenitrifiersGrowth.solve();
            fvOptions.correct(XDN);

	        // Minimum viable biomass
            XAR *= pos(XAR - X_min);
            XN  *= pos(XN  - X_min);
            XDN *= pos(XDN - X_min);

            //Info << "\nExtrapolymeric substances EPS transport" << endl;
            fvScalarMatrix EPSTransport
            (
                fvm::ddt(EPS)
                ==
                  (kEPS*rH  + fd*pXAR.bDie) * XAR
                + (kEPS*rN  + fd*pXN.bDie)  * XN
                + (kEPS*rDN + fd*pXDN.bDie) * XDN
                - kdet_EPS * EPS
                + n * katt_BAP * BAP
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
                (1-fd) * (pXAR.bDie*XAR + pXN.bDie*XN + pXDN.bDie*XDN)
                - kdet_XI * XI
                + n * katt_POCr * POCr
            );
            InertGeneration.relax();
            fvOptions.constrain(InertGeneration);
            InertGeneration.solve();
            fvOptions.correct(XI);

            //Info << "\nDissolved Organic Carbon DOC transport" << endl;
            fvScalarMatrix DissolvedCarbonTransport
            (
                fvm::ddt(n,DOC)
                + fvm::div(phi, DOC)
                - fvm::laplacian(n*(mag(U)*DispTensor + molDiff), DOC)
                ==
                - rH * XAR
                - rDN * XDN
                + n * khyd_BAP * BAP
                + n * khyd_POCr * POCr
            );
            DissolvedCarbonTransport.relax();
            fvOptions.constrain(DissolvedCarbonTransport);
            DissolvedCarbonTransport.solve();
            fvOptions.correct(DOC);
            
            //Info << "\nDissolved oxygen (O2)" << endl;
            fvScalarMatrix OxygenTransport
            (
                fvm::ddt(n,O2)
                + fvm::div(phi, O2)
                - fvm::laplacian(n*(mag(U)*DispTensor + molDiff), O2)
                ==
                - alpha_1 * rH * XAR
                - alpha_N * rN * XN
            );
            OxygenTransport.relax();
            fvOptions.constrain(OxygenTransport);
            OxygenTransport.solve();
            fvOptions.correct(O2);

            // Info << "\nNH4-N transport" << endl;
            fvScalarMatrix AmmoniumTransport
            (
                fvm::ddt(n,NH4)
                + fvm::div(phi, NH4)
                - fvm::laplacian(n*(mag(U)*DispTensor + molDiff), NH4)
                ==
                - rN * XN 
                + gamma_N * fd * (pXAR.bDie*XAR + pXN.bDie*XN + pXDN.bDie*XDN)
                - gamma_N * pXAR.yield * rH * XAR
            );
            AmmoniumTransport.relax();
            fvOptions.constrain(AmmoniumTransport);
            AmmoniumTransport.solve();
            fvOptions.correct(NH4);

            // Info << "\nDissolved nitrate-nitrogen (NO3+N)" << endl;
            fvScalarMatrix NitrateTransport
            (
                fvm::ddt(n,NO3) 
                + fvm::div(phi, NO3) 
                - fvm::laplacian(n*(mag(U)*DispTensor + molDiff), NO3)
                ==
                (1.0 - 0.31*pXN.yield - kEPS) * rN * XN
                - beta_1 * rDN * XDN
            );
            NitrateTransport.relax();
            fvOptions.constrain(NitrateTransport);
            NitrateTransport.solve();
            fvOptions.correct(NO3);

	        //Info << "\nBiomass-associated products (BAP) transport" << endl;
            fvScalarMatrix BAPTransport
            (
                fvm::ddt(n,BAP)
                + fvm::div(phi, BAP)
                - fvm::laplacian(n*(mag(U)*DispTensor + molDiff), BAP)
                ==
                kdet_EPS * EPS
                - n * katt_BAP * BAP
                - n * khyd_BAP * BAP
            );
            BAPTransport.relax();
            fvOptions.constrain(BAPTransport);
            BAPTransport.solve();
	        fvOptions.correct(BAP);

            //Info << "\nInert biomass (XI)" << endl;
            fvScalarMatrix POCrGeneration
            (
                fvm::ddt(n,POCr)
                + fvm::div(phi, POCr)
                - fvm::laplacian(mag(U)*DispTensor, POCr)
                ==
                kdet_XI * XI
                - n * katt_POCr * POCr
                - n * khyd_POCr * POCr
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
