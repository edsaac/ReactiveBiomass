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
    richardsFoamE

Description
    Refactor RichardsFoam3 to adapt to a reactive transport solver

    Solves for saturation and head:
    - Water saturation (theta)
    - Head (h)
    Eqs:
    
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"

#include "UnsaturatedSoilParametersField.H"
#include "timeStepper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    volScalarField z (mesh.C().component(2));
    const label nbMesh = mesh.nCells();
    Foam::Info << "nCells: " << nbMesh;

    #include "createFields.H"
    #include "readParameters.H"
    #include "createFvOptions.H"

    double convergeFlow = 1.0;
    int nCycles = 0;
    // unsigned int nItersDebug = 0;

    Foam::Info << "\nCalculating...\n" << endl;

    while (runTime.loop())
    {
        Foam::Info<< "Time = "   << runTime.timeName() << endl; 
        
        nCycles = 0;
        h_before = h;
        h_after = h;

        // Solve saturation from Richards Equation    
        while(true)
        {
            // Solve Richard's equation    
            fvScalarMatrix richardsEquation
            (
                fvm::ddt(soil.capillary(h_before), h_after)
                ==
                fvm::laplacian(soil.K(h_before), h_after)
                + fvc::laplacian(soil.K(h_before), z)
            );
            fvOptions.constrain(richardsEquation);
            richardsEquation.solve();
            fvOptions.correct(h_after);

            // Check if solution converged
            err = Foam::mag(h_after - h_before);
            convergeFlow = Foam::gSumMag(err)/nbMesh;
            Foam::Info << "nCycles: "    << nCycles      << "\t"
                       << "Converger: " << convergeFlow << endl;

            if (adjustTimeStep)
            {
                #include "timeControl.H"
            }
            else 
            {
                break;
            }
        }
        
        // Calculate Darcy flow velocity
        hydrConduct = soil.K(h);
        U = - hydrConduct * (fvc::grad(h) + fvc::grad(z));
        phi = fvc::flux(U);

        // Calculate derived unsaturated stuff
        theta = soil.waterContentCalculator(h);

        //End bits
        runTime.write();

        Foam::Info << "ExecutionTime = " << runTime.elapsedCpuTime()   << " s"
                   << "    ClockTime = " << runTime.elapsedClockTime() << " s"
                   << "       deltaT = " << runTime.deltaTValue()      << " s"
                   << nl << endl;

        // nItersDebug++;
        // if (nItersDebug > 100) { break; }
    }

    Foam::Info<< "End\n" << endl;
    return 0;
}