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
    unsatFoam

Description
    Only Richards equation solver

\*---------------------------------------------------------------------------*/

#define DEBUG false

// #define updateHydCond hydraulicCond = K_0 * perm_satu
#define debug(message) if(DEBUG){Foam::Info << message << nl << endl;}

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

#include "timeStepper.H"
#include "declareClasses.H"

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
    soil.waterSaturationCalculator(h);      //<- Updates Sw(h)
    Sw.write();

    soil.mualemCalculator(h);               //<- Updates perm_satu(h)
    hydraulicCond = K_0 * perm_satu;        //<- Updates K(h)
    hydraulicCond.write();
    
    U = - hydraulicCond * (fvc::grad(h) + fvc::grad(z));  
    U.write();

    Foam::Info << "\nCalculating...\n" << endl;

    while (simple.loop(runTime))
    {
        Foam::Info<< "Time = " << runTime.timeName() << nl << endl;
        
        //- Start Richards' solver block
        nCycles = 0;
        h_before = h;
        h_after = h;

        while(true)
        {

            //- Calculate grad(K(h)) and extract z-component
            soil.waterSaturationCalculator(h_before);    // <- Updates Sw(h)
            soil.mualemCalculator(h_before);             // <- Updates perm_satu(h)
            hydraulicCond = K_0 * perm_satu;             // <- Updates hydraulicCond 

            grad_k = fvc::grad(hydraulicCond);
            grad_kz = grad_k.component(vector::Z);

            //- Solve Richard's equation for deformable porous media
            fvScalarMatrix richardsEquation
            (
                porosity * soil.capillary(h_before) * fvm::ddt(h_after)
                + Sw * fvc::ddt(porosity)
                ==
                fvm::laplacian(hydraulicCond, h_after)
                + grad_kz
            );

            fvOptions.constrain(richardsEquation);
            richardsEquation.solve();
            fvOptions.correct(h_after);

            //- Check if solution converged
            err = Foam::mag(h_after - h_before);
            convergeFlow = Foam::gSumMag(err)/nbMesh;
            Foam::Info << "nCycles: "    << nCycles      << "\t"
                       << "Converger: " << convergeFlow << endl;

            debug("Adjust time step if possible");
            if (adjustTimeStep) 
            {
                #include "timeControl.H"
            }
            else
            { 
                break; 
            }
            //- If h converged, the h field is updated with the one found 
            //  from the iterations above and this loop is broken
        }

        //- Update Sw and perm_satu based on hydraulic head
        soil.waterSaturationCalculator(h);
        soil.mualemCalculator(h);
        hydraulicCond = K_0 * perm_satu;      
        Sa = 1.0 - Sw;
        capillarity = soil.capillary(h);
        
        // Update flow field
        U = - hydraulicCond * (fvc::grad(h) + fvc::grad(z));
        phi = fvc::flux(U);
        #include "CourantNo.H"

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
