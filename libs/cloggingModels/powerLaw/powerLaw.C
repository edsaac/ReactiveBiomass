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

#include "powerLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cloggingModels
{
    defineTypeNameAndDebug(powerLaw, 0);

    addToRunTimeSelectionTable
    (
        cloggingModel,
        powerLaw,
        dictionary
    );
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::cloggingModels::powerLaw::calcPerm()
{
    *ptrperm_ = ((permRef_ - permMin_) * Foam::pow((n_ - nMin_)/(nRef_ - nMin_),nExponent_.value()) * pos(n_ - nMin_)) + permMin_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cloggingModels::powerLaw::powerLaw
(
    const word& name,
    const dictionary& cloggingProperties,
    const volScalarField& n,
    volScalarField* ptrperm
)
:
    cloggingModel(name, cloggingProperties, n, ptrperm),
    powerLawCoeffs_(cloggingProperties.optionalSubDict(typeName + "Coeffs")),
    nExponent_("nExponent", dimless, powerLawCoeffs_)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::cloggingModels::powerLaw::read
(
    const dictionary& cloggingProperties
)
{
    cloggingModel::read(cloggingProperties);

    powerLawCoeffs_ = cloggingProperties.optionalSubDict(typeName + "Coeffs");
    powerLawCoeffs_.lookup("nExponent") >> nExponent_;

    return true;
}


// ************************************************************************* //
