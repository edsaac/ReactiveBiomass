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

\*---------------------------------------------------------------------------*/

#include "cloggingModel.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cloggingModel, 0);
    defineRunTimeSelectionTable(cloggingModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cloggingModel::cloggingModel
(
    const word& name,
    const dictionary& cloggingProperties,
    const volScalarField& n,
    volScalarField* ptrperm
)
:
    name_(name),
    cloggingProperties_(cloggingProperties),
    n_(n),
    ptrperm_(ptrperm),
    permMin_(cloggingProperties.lookup("permMin")),
    permRef_(cloggingProperties.lookup("permRef")),
    nMin_(cloggingProperties.lookup("nMin")),
    nRef_(cloggingProperties.lookup("nRef"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::cloggingModel::read(const dictionary& cloggingProperties)
{
    cloggingProperties_ = cloggingProperties;

    return true;
}


// ************************************************************************* //
