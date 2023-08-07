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

Description
    Set the time control

// \*---------------------------------------------------------------------------*/

// Was the non-linerity solved?
if (convergeFlow < timeControl.convergeThreshold)
{
    h = h_after;

    // Can we increase the timestep?
    if (nCycles <= timeControl.UpDeltaT_nIters)
    {
        runTime.setDeltaT
        (
            min
            (
                timeControl.UpDeltaT_Factor * runTime.deltaTValue(),
                timeControl.maxDeltaT
            )
        );

        Foam::Info << nl << "deltaT going UP ↑↑↑ = "
                   << runTime.deltaTValue() << endl;
    }
    break;
}

// The non-linerity was not solved, iterate over h.
else 
{
    h_before = (timeControl.hguess_coeff * h_after) + ((1.0 - timeControl.hguess_coeff) * h_before);
    // h_after = h;
    nCycles++;

    // Eq has not converged after the max nIterations
    // Decrease timeStep and start all over
    if (nCycles >= timeControl.maxItersToConverge)
    {
        // End simulation if deltaT is too small 
        if (runTime.deltaTValue() < timeControl.minDeltaT)
        {
            RETURNCODE = 1;
            
            Foam::Info << "\n deltaT is lower than minDeltaT: "
                       << runTime.deltaTValue() << " < " << timeControl.minDeltaT << nl
                       << "\n Failed to converge :( " << nl << endl;
            
            Foam::Info << "\n Final time: " << runTime.timeName()<< endl;
            runTime.writeAndEnd();
            break;
        }

        runTime.setDeltaT
        (
            timeControl.DownDeltaT_Factor * runTime.deltaTValue()
        );

        Foam::Info << nl << "deltaT going DOWN ↓↓↓: "
                   << runTime.deltaTValue() << endl;

        h_after = h;
        nCycles = 0;
    }
}