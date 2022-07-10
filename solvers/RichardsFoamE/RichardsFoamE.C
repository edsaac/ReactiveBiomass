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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    volScalarField z (mesh.C().component(2));
      
    #include "readParameters.H"
    #include "createFields.H"
    #include "createFvOptions.H"
   
    const label nbMesh = mesh.nCells();
    double convergeFlow = 1.0;
    int nCycles = 0;
    // unsigned int nItersDebug = 0;

    Foam::Info << "\nCalculating...\n" << endl;
    
    while (runTime.loop())
    {
        Foam::Info<< "Time = "   << runTime.timeName() << nl 
                  << "deltaT = " << runTime.deltaTValue()  << endl;
        
        nCycles = 0;
        h_before = h;
        h_after = h;

        Foam::Info << "nCells: " << nbMesh;

        // Solve saturation from Richards Equation    
        while(true)
        {
            Foam::Info << "Enter the inner loop " << endl; 

            // Solve Richard's equation    
            fvScalarMatrix richardsEquation
            (
                fvm::ddt(soil.capillary(h_before), h_after)
                - fvm::laplacian(soil.K(h_before), h_after)
                ==
                fvc::laplacian(soil.K(h_before), z)
            );
            richardsEquation.relax();
            fvOptions.constrain(richardsEquation);
            richardsEquation.solve();
            fvOptions.correct(h_after);

            // Check if solution converged
            err = h_after - h_before;
            convergeFlow = gSumMag(err)/nbMesh;
            Foam::Info << "Converger: " << convergeFlow << "\n"
                       << "nCycles:"    << nCycles      << endl;

            if(convergeFlow < 1.0E-04) // converged
            {
                h = h_after;

                if(nCycles < 2)
                {
                    runTime.setDeltaT( 1.2 * runTime.deltaTValue() );
                    Foam::Info << "\n deltaT going UP ↑↑↑"
                               << "deltaT = " << runTime.deltaTValue()  << endl;
                }
                break;
            }

            else // did not converged, try again
            {
                h_before = h_after;
                h_after = h;
                nCycles++;
            }

            if(nCycles > 5) // has not converged yet
            {
                Foam::Info << "deltaT halved" << endl;
                runTime.setDeltaT
                (
                    min
                    (
                        0.5 * runTime.deltaTValue(),
                        100.0
                    )
                );

                Foam::Info << "\n deltaT going DOWN ↓↓↓"
                            << "deltaT = " << runTime.deltaTValue()  << endl;

                h_after = h;
                nCycles = 0;
            }
            
        }
        
        // Calculate Darcy flow velocity
        U = - soil.K(h) * (fvc::grad(h) + fvc::grad(z));
        phi = fvc::flux(U);

        //End bits
        runTime.write();

        Foam::Info << "ExecutionTime = " << runTime.elapsedCpuTime()   << " s"
                   << "    ClockTime = " << runTime.elapsedClockTime() << " s"
                   << nl << endl;

        // nItersDebug++;
        // if (nItersDebug > 4) { break; }
    }

    Foam::Info<< "End\n" << endl;
    return 0;
}