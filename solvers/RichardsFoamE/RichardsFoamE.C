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
   
    Foam::Info << "\nCalculating...\n" << endl;
    
    while (runTime.loop())
    {
        Foam::Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // Solve saturation from Richards Equation    
        for (size_t i = 0; i < 1; i++)
        {
            Foam::Info << "Enter the inner loop " << endl; 
                     
            // Solve Richard's equation    
            fvScalarMatrix richardsEquation
            (
                fvm::ddt(soil.capillary(h), h)
                - fvm::laplacian(soil.K(h), h)
                ==
                fvc::laplacian(soil.K(h), z)
            );
            richardsEquation.relax();
            fvOptions.constrain(richardsEquation);
            richardsEquation.solve();
            fvOptions.correct(h);
        }
        
        // Calculate Darcy flow velocity
        U = - soil.K(h) * (fvc::grad(h) + fvc::grad(z));
        phi = fvc::flux(U);

        //End bits
        runTime.write();

        Foam::Info << "ExecutionTime = " << runTime.elapsedCpuTime()   << " s"
                   << "    ClockTime = " << runTime.elapsedClockTime() << " s"
                   << nl << endl;
    }

    Foam::Info<< "End\n" << endl;
    return 0;
}