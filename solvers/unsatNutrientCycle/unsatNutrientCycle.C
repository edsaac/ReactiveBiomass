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

    Generation and consumption of EPS and BAP are included.
    Only advection of dissolved and particulate species is considered

    Clogging is introduced by some porosity-permeability relationship
    
    k(n)/k0 = f(n/n0)

    Solves eight species:
    - Aerobic heterotrophs (XAR) {X}
    - Inert biomass  (XI) {Xi}
    - Extrapolymeric substances (EPS) {E}
    - Biomass-associated products (BAP) {B}
    - Substrate (DOC) {S}
    - Oxygen (O2) {O}
    - Ammonium (NH4) and nitrate (NO3)

    Eqs:

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

#define DEBUG false

#define updateHydCond hydraulicCond = K_0 * perm_clog * perm_satu
#define debug(message) if(DEBUG){Foam::Info << message << nl << endl;}

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

#include "timeStepper.H"
#include "declareClasses.H"
#include "cloggingModel.H"
#include "attachmentModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    unsigned short int RETURNCODE = 0;

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);
    volScalarField z (mesh.C().component(2));
    const label nbMesh = mesh.nCells();
    Foam::Info << "nCells: " << nbMesh;

    #include "readParameters.H"
    #include "createFields.H"
    #include "CourantNo.H"
    #include "createFvOptions.H"

    double convergeFlow = 1.0;
    int nCycles = 0;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Foam::Info << "\nInitialize derived fields...\n" << endl;
    
    debug("Water saturation...");
    soil.waterSaturationCalculator(h);
    Sw.write();
    
    debug("Clogging & porosity...");
    if (cloggingSwitch) { clogging->calcPerm();}
    soil.mualemCalculator(h);
    
    // debug("Hydraulic conductivity...");
    updateHydCond;
    hydraulicCond.write();

    // debug("Flow velocity...");
    U = - hydraulicCond * (fvc::grad(h) + fvc::grad(z));  
    U.write();

    Foam::Info << "\nCalculating...\n" << endl;

    while (simple.loop(runTime))
    {
        Foam::Info<< "Time = " << runTime.timeName() << nl << endl;

        
        /*
        - Calculate the total immobile biomass, the growth rate limiting function
          and check that the total biomass has not exceeded XMAX anywhere

          clogLimiter
            ^
            |
          1 -
            | *
            |    *
            |       *
          0 -----------|------> totalBiomass
                      XMAX
        */

        totalBiomass = XAR + XN + XDN + XI + EPS;
        clogLimiter = 1.0 - totalBiomass/XMAX;

        if (Foam::min(clogLimiter) < dimensionedScalar("zero",dimless,0.0))
        {
            Foam::SeriousError<< "Total biomass greater that XMAX" << endl;
            Info << "Exiting with code 102..." << endl;
            exit(102);
        }

        /*
        - Update the porosity field based on the void space filled with biomass
          and update the permeability perm_clog multiplier
        */
        porosity  = clogging->nRef() - totalBiomass/rho_X;
        
        if (cloggingSwitch) 
        { 
            clogging->calcPerm();
        }
        
        //- Start Richards' solver block
        nCycles = 0;
        h_before = h;
        h_after = h;

        while(true)
        {

            //- Calculate grad(K(h)) and extract z-component
            debug("Calculate perm_satu");
            soil.waterSaturationCalculator(h_before); // <- Updates Sw(h)
            soil.mualemCalculator(h_before); // <- Updates ksat(h)

            debug("Update K");
            updateHydCond; // <- Updates hydraulicCond 

            debug("Calculate K gradients");
            grad_k = fvc::grad(hydraulicCond);
            grad_kz = grad_k.component(vector::Z);

            //- Solve Richard's equation for undeformable porous media
            //--  To do: Add ddt(porosity) term for deformable media (done)
            fvScalarMatrix richardsEquation
            (
                porosity * soil.capillary(h_before) * fvm::ddt(h_after)
                // fvm::ddt(soil.capillary(h_before), h_after)
                + Sw * fvc::ddt(porosity)
                ==
                fvm::laplacian(hydraulicCond, h_after)
                + grad_kz
            );
            debug("Solve Richards");
            fvOptions.constrain(richardsEquation);
            richardsEquation.solve();
            fvOptions.correct(h_after);

            //- Check if solution converged
            debug("Check convergence");
            err = Foam::mag(h_after - h_before);
            convergeFlow = Foam::gSumMag(err)/nbMesh;
            Foam::Info << "nCycles: "    << nCycles      << "\t"
                       << "Converger: " << convergeFlow << endl;

            debug("Adjust time step if possible");
            
            if (adjustTimeStep) {
                #include "timeControl.H"
            }
            else { 
                break; 
            }
            //- If h converged, the h field is updated with the one found 
            //  from the iterations above and this loop is broken
        }

        //- Update Sw and perm_satu based on hydraulic head
        soil.waterSaturationCalculator(h);
        soil.mualemCalculator(h);
        updateHydCond;
        Sa = 1.0 - Sw;

        capillarity = soil.capillary(h);
        
        // Update flow field
        U = - hydraulicCond * (fvc::grad(h) + fvc::grad(z));
        phi = fvc::flux(U);
        #include "CourantNo.H"

        //- Calculate attach/detach rates
        attachment->calcAttachment();
        detachment->calcAttachment();

        //- Transport equations
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
                - fvm::laplacian(diffusiveGrowth, XAR)  //Pseudo-diffusion
                ==
                pXAR.yield * rH * XAR                   //Growth
                - pXAR.bDie * XAR                       //Decay XAR -> EPS + XI
                - kdet * XAR                            //Detachment XAR -> XARp
                + (Sw * porosity) * katt * clogLimiter * XARp         //Attachment XARp -> XAR
            );
            ARGrowth.relax();
            fvOptions.constrain(ARGrowth);
            ARGrowth.solve();
            fvOptions.correct(XAR);

            //Info << "\nMobile aerobic heterotrophs (XARp)" << endl;
            fvScalarMatrix ARGrowth_p
            (
                fvm::ddt(porosity, Sw, XARp)
                + fvm::div(phi, XARp)
                - fvm::laplacian(porosity * Sw * (mag(U)*DispTensor + molDiff), XARp)
                ==
                porosity * Sw * pXAR.yield * rH * XARp              //Growth
                - porosity * Sw * pXAR.bDie * XARp                  //Decay XARp -> BAP + POCr
                + kdet * XAR                                 //Detachment XAR -> XARp
                - porosity * Sw * katt * clogLimiter * XARp          //Attachment XARp -> XAR
            );
            ARGrowth_p.relax();
            fvOptions.constrain(ARGrowth_p);
            ARGrowth_p.solve();
            fvOptions.correct(XARp);

            // Info << "\nImmobile Nitrifiers (XN)" << endl;
            fvScalarMatrix NitrifiersGrowth
            (
                fvm::ddt(XN)
                - fvm::laplacian(diffusiveGrowth, XN)  //Pseudo-diffusion
                ==
                  pXN.yield * rN * XN                        //Growth
                - pXN.bDie * XN                              //Decay XN -> EPS + XI
                - kdet * XN                                  //Detachment XN -> XNp
                + porosity * Sw * katt * clogLimiter * XNp          //Attachment XNp -> XN
            );
            NitrifiersGrowth.relax();
            fvOptions.constrain(NitrifiersGrowth);
            NitrifiersGrowth.solve();
            fvOptions.correct(XN);

            // Info << "\nMobile Nitrifiers (XN)" << endl;
            fvScalarMatrix NitrifiersGrowth_p
            (
                fvm::ddt(porosity, Sw, XNp)
                + fvm::div(phi, XNp)
                - fvm::laplacian(porosity * Sw * (mag(U)*DispTensor + molDiff), XNp)
                ==
                porosity * Sw * pXN.yield * rN * XNp                    //Growth
                - porosity * Sw * pXN.bDie * XNp                        //Decay XNp -> BAP + POCr
                + kdet * XN                                      //Detachment XN -> XNp
                - porosity * Sw * katt * clogLimiter * XNp              //Attachment XNp -> XN
            );
            NitrifiersGrowth_p.relax();
            fvOptions.constrain(NitrifiersGrowth_p);
            NitrifiersGrowth_p.solve();
            fvOptions.correct(XNp);

            // Info << "\nImmobile Denitrifiers (XDN)" << endl;
            fvScalarMatrix DenitrifiersGrowth
            (
                fvm::ddt(XDN)
                - fvm::laplacian(diffusiveGrowth, XDN)       //Pseudo-diffusion
                ==
                pXDN.yield * rDN * XDN                       //Growth
                - pXDN.bDie * XDN                            //Decay XDN -> EPS + XI
                - kdet * XDN                                 //Detachment XDN -> XDNp
                + porosity * Sw * katt * clogLimiter * XDNp         //Attachment XDNp -> XDN
            );
            DenitrifiersGrowth.relax();
            fvOptions.constrain(DenitrifiersGrowth);
            DenitrifiersGrowth.solve();
            fvOptions.correct(XDN);

            // Info << "\nMobile Denitrifiers (XDN)" << endl;
            fvScalarMatrix DenitrifiersGrowth_p
            (
                porosity * Sw * fvm::ddt(porosity, Sw, XDNp)
                + fvm::div(phi, XDNp)
                - fvm::laplacian(porosity * Sw*(mag(U)*DispTensor + molDiff), XDNp)
                ==
                porosity * Sw * pXDN.yield * rDN * XDNp                 //Growth
                - porosity * Sw * pXDN.bDie * XDNp                      //Decay XDN -> BAP + POCr
                + kdet * XDN                                     //Detachment XDN -> XDNp
                - porosity * Sw * katt * clogLimiter * XDNp             //Attachment XDNp -> XDN
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
                + porosity * Sw * katt * clogLimiter * BAP     //Attachment BAP -> EPS
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
                + porosity * Sw * katt * clogLimiter * POCr    //Attachment POCr -> XI
                - khyd_recal * XI                       //Hydrolysis of XI -> DOC
            );
            InertGeneration.relax();
            fvOptions.constrain(InertGeneration);
            InertGeneration.solve();
            fvOptions.correct(XI);

            //Info << "\nDissolved Organic Carbon DOC transport" << endl;
            fvScalarMatrix DissolvedCarbonTransport
            (
                fvm::ddt(porosity, Sw, DOC)
                + fvm::div(phi, DOC)
                - fvm::laplacian(porosity * Sw*(mag(U)*DispTensor + molDiff), DOC)
                ==
                - rH * XAR                              //Metabolism XAR
                - rDN * XDN                             //Metabolism XDN
                + porosity * Sw * (
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
                fvm::ddt(porosity, Sw, O2)
                + fvm::div(phi, O2)
                // - fvm::laplacian(porosity * Sw*(mag(U)*DispTensor + molDiff), O2)
                ==
                - alpha_1 * rH * (XAR + (porosity * Sw * XARp))   //Metabolism XAR
                - alpha_N * rN * (XN + (porosity * Sw * XNp))     //Metabolism XN
                // - oxygen_mass_transfer * porosity * ((Sw * O2) - (Sa * Hacc * O2gas)) //Replenish to saturation if Sa?
                // + oxygen_mass_transfer * porosity * (Sa * Hacc * O2gas) //Replenish to saturation if Sa?

            );
            OxygenTransport.relax();
            fvOptions.constrain(OxygenTransport);
            OxygenTransport.solve();
            fvOptions.correct(O2);

            //Info << "\Gaseous oxygen (O2gas)" << endl;
            // fvScalarMatrix OxygenGasTransport
            // (
            //     porosity * Sa * fvm::ddt(O2gas)
            //     + O2gas * fvc::ddt(Sa, porosity)                   // Oxygen in the gaseous phase
            //     - fvm::laplacian(porosity * Sa * molDiff_air, O2gas)
            //     ==
            //     // oxygen_mass_transfer * porosity * ((Sw * O2) - (Sa * Hacc * O2gas)) 
            //     // - oxygen_mass_transfer * porosity * (Sa * Hacc * O2gas)
            // );
            // OxygenGasTransport.relax();
            // fvOptions.constrain(OxygenGasTransport);
            // OxygenGasTransport.solve();
            // fvOptions.correct(O2gas);

            if (Foam::max(O2) > O2_saturation)
            {
                Foam::SeriousError<< "Oxygen in the aqueous phase is getting over saturated" << endl;
                Info << "Exiting with code 103..." << endl;
                exit(103);
            }

            // Info << "\nNH4-N transport" << endl;
            fvScalarMatrix AmmoniumTransport
            (
                fvm::ddt(porosity, Sw, NH4)
                + fvm::div(phi, NH4)
                - fvm::laplacian(porosity * Sw*(mag(U)*DispTensor + molDiff), NH4)
                ==
                - rN * XN                               //Metabolism XN
                + gamma_N * fd * (
                    pXAR.bDie*XAR                       //Decay XAR
                    + pXN.bDie*XN                       //Decay XN
                    + pXDN.bDie*XDN                     //Decay XDN
                    )
                + gamma_N * fd * porosity * Sw * (                      
                    pXAR.bDie*XARp                      //Decay XARp
                    + pXN.bDie*XNp                      //Decay XNp
                    + pXDN.bDie*XDNp                    //Decay XDNp
                    )                   
                - gamma_N * pXAR.yield * rH * XAR            //Metabolism XAR
                - porosity * Sw * gamma_N * pXAR.yield * rH * XARp  //Metabolism XARp
            );
            AmmoniumTransport.relax();
            fvOptions.constrain(AmmoniumTransport);
            AmmoniumTransport.solve();
            fvOptions.correct(NH4);

            // Info << "\nDissolved nitrate-nitrogen (NO3+N)" << endl;
            fvScalarMatrix NitrateTransport
            (
                fvm::ddt(porosity, Sw, NO3)
                + fvm::div(phi, NO3) 
                - fvm::laplacian(porosity * Sw*(mag(U)*DispTensor + molDiff), NO3)
                ==
                (1.0 - 0.31*pXN.yield - kEPS) * rN * (XN + porosity * Sw*XNp)   //Metabolism XN
                - beta_1 * rDN * (XDN + porosity * Sw*XDNp)                     //Metabolism XDN
            );
            NitrateTransport.relax();
            fvOptions.constrain(NitrateTransport);
            NitrateTransport.solve();
            fvOptions.correct(NO3);

	        //Info << "\nBiomass-associated products (BAP) transport" << endl;
            fvScalarMatrix BAPTransport
            (
                fvm::ddt(porosity, Sw, BAP)
                + fvm::div(phi, BAP)
                - fvm::laplacian(porosity * Sw*(mag(U)*DispTensor + molDiff), BAP)
                ==
                kdet * EPS                              //Detachment EPS -> BAP
                - porosity * Sw * katt * clogLimiter * BAP     //Attachment BAP -> EPS
                - porosity * Sw * khyd_labil * BAP             //Hydrolisis BAP -> DOC
                + porosity * Sw * fd * (
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
                fvm::ddt(porosity, Sw, POCr)
                + fvm::div(phi, POCr)
                - fvm::laplacian(porosity * Sw*(mag(U)*DispTensor + molDiff), POCr)
                ==
                kdet * XI                       //Detachment XI -> POCr
                - porosity * Sw * POCr * (
                    katt * clogLimiter          //Attachment POCr -> XI
                    + khyd_recal                //Hydrolisis POCr -> DOC
                    )                 
                + porosity * Sw * (1-fd) * (
                      pXAR.bDie * XARp          //Decay XARp
                    + pXN.bDie  * XNp           //Decay XNp
                    + pXDN.bDie * XDNp          //Decay XDNp
                )
            );
            POCrGeneration.relax();
            fvOptions.constrain(POCrGeneration);
            POCrGeneration.solve();
            fvOptions.correct(POCr);

        //Info << "\nNonreactive tracer transport" << endl;
            fvScalarMatrix NonReactiveTracer
            (
                fvm::ddt(porosity, Sw, tracer)
                + fvm::div(phi, tracer)
                // - fvm::laplacian(porosity * Sw*(mag(U)*DispTensor + molDiff), tracer)
            );
            NonReactiveTracer.relax();
            fvOptions.constrain(NonReactiveTracer);
            NonReactiveTracer.solve();
            fvOptions.correct(tracer);

        }

        //End bits
        runTime.write();

        Foam::Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
    }

    runTime.writeNow();

    Foam::Info<< "End\n" << endl;

    return RETURNCODE;
}


// ************************************************************************* //
