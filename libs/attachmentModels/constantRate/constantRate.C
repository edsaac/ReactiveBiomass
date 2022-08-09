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

#include "constantRate.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace attachmentModels
{
    defineTypeNameAndDebug(constantRate, 0);

    addToRunTimeSelectionTable
    (
        attachmentModel,
        constantRate,
        dictionary
    );
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::attachmentModels::constantRate::constantRate
(
    const word& name,
    const dictionary& attachmentProperties,
    volScalarField* ptrkatt,
    const fvMesh& mesh,
    const Time& runTime,
    const volVectorField& U,
    const volScalarField& n
)
:
    attachmentModel(name, attachmentProperties, ptrkatt, mesh, runTime,U,n),
    constantRateCoeffs_(attachmentProperties.optionalSubDict(typeName + "Coeffs")),
    
    value_("value", dimensionSet(0,0,-1,0,0,0,0), constantRateCoeffs_)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::attachmentModels::constantRate::calcAttachment()
{   
    *ptrkatt_ = this->value_;
}


bool Foam::attachmentModels::constantRate::read
(
    const dictionary& attachmentProperties
)
{
    attachmentModel::read(attachmentProperties);

    constantRateCoeffs_ = attachmentProperties.optionalSubDict(typeName + "Coeffs");
    constantRateCoeffs_.lookup("value") >> value_;

    return true;
}


// ************************************************************************* //
