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
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    volScalarField z (mesh.C().component(2));
    simpleControl simple(mesh);
      
    #include "readParameters.H"
    #include "createFields.H"
    #include "createFvOptions.H"
   
    Foam::Info << "\nCalculating...\n" << endl;
    
    while (runTime.loop())
    {
        Foam::Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // Solve saturation from Richards Equation    
        while (simple.correctNonOrthogonal())
        {
            Foam::Info << "Enter the inner loop " << endl; 
            
            // Aux field of (-alpha*h)^n
            // if h < 0 
            innerAlphaPower = Foam::pow(- soil.alpha * h * Foam::neg(h), soil.n);
            Foam::Info << "theta_e calculated " << innerAlphaPower.dimensions() << endl; 
            
            // Calculate effective saturation from van Genutchen eq.
            theta_e = Foam::pow(1.0 + innerAlphaPower, - soil.m);
            Foam::Info << "theta_e calculated " << theta_e.dimensions() << endl;

            // Calculate water content
            theta = theta_e * (soil.theta_s - soil.theta_r) + soil.theta_r;
            Foam::Info << "theta calculated " << endl;
            
            // Calculate relative permebility
            perm = Foam::sqrt(theta_e) * Foam::sqr(1.0 - Foam::pow(1-Foam::pow(theta_e, 1.0/soil.m), soil.m));
            Foam::Info << "perm calculated " << endl;

            // Calculate calpillary capacity
            // Check sympy for the dtheta/dh derivation
            capillary = (soil.theta_s - soil.theta_r) 
                * soil.m * soil.n
                * innerAlphaPower
                * theta_e
                / (h * (innerAlphaPower + 1));
            Foam::Info << "capillary calculated" << endl;

            // Calculate hydraulic conductivity 
            hydrConduct = soil.K_s * perm;
            Foam::Info << "hydrConduct calculated" << endl;

            // Solve Richard's equation    
            fvScalarMatrix richardsEquation
            (
                fvm::ddt(capillary,h)
                - fvm::laplacian(hydrConduct, h)
                ==
                fvc::laplacian(hydrConduct, z)
            );
            Foam::Info << "FVM matrix built" << endl;
            richardsEquation.relax();
            fvOptions.constrain(richardsEquation);
            richardsEquation.solve();
            Foam::Info << "FVM matrix solved" << endl;
            fvOptions.correct(h);
        }
        
        // Calculate Darcy flow velocity
        U = - hydrConduct * (fvc::grad(h) + fvc::grad(z));
        Foam::Info << "Flow velocity calculated " << endl;
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