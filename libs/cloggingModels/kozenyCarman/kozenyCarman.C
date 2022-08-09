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

#include "kozenyCarman.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cloggingModels
{
    defineTypeNameAndDebug(kozenyCarman, 0);

    addToRunTimeSelectionTable
    (
        cloggingModel,
        kozenyCarman,
        dictionary
    );
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::cloggingModels::kozenyCarman::calcPerm()
{
    *ptrperm_ = ((permRef_ - permMin_) 
                * Foam::pow((n_ - nMin_)/(nRef_ - nMin_),nExponent_.value()) 
                * Foam::pow((1.0 - nRef_)/(1.0 - n_),mExponent_.value())
                * pos(n_ - nMin_)) 
                + permMin_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cloggingModels::kozenyCarman::kozenyCarman
(
    const word& name,
    const dictionary& cloggingProperties,
    const volScalarField& n,
    volScalarField* ptrperm
)
:
    cloggingModel(name, cloggingProperties, n, ptrperm),
    kozenyCarmanCoeffs_(cloggingProperties.optionalSubDict(typeName + "Coeffs")),
    nExponent_("nExponent", dimless, kozenyCarmanCoeffs_),
    mExponent_("mExponent", dimless, kozenyCarmanCoeffs_)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::cloggingModels::kozenyCarman::read
(
    const dictionary& cloggingProperties
)
{
    cloggingModel::read(cloggingProperties);

    kozenyCarmanCoeffs_ = cloggingProperties.optionalSubDict(typeName + "Coeffs");
    kozenyCarmanCoeffs_.lookup("nExponent") >> nExponent_;
    kozenyCarmanCoeffs_.lookup("mExponent") >> mExponent_;

    return true;
}


// ************************************************************************* //
