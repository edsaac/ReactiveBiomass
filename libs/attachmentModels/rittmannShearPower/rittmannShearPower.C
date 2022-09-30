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

#include "rittmannShearPower.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace attachmentModels
{
    defineTypeNameAndDebug(rittmannShearPower, 0);

    addToRunTimeSelectionTable
    (
        attachmentModel,
        rittmannShearPower,
        dictionary
    );
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::attachmentModels::rittmannShearPower::rittmannShearPower
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
    
    rittmannShearPowerCoeffs_(attachmentProperties.optionalSubDict(typeName + "Coeffs")),
    dc_("collectorSize", dimLength, rittmannShearPowerCoeffs_),
    mu_("fluidViscosity", dimDynamicViscosity, rittmannShearPowerCoeffs_),
    Cd_("rittmannCoeff", dimless, rittmannShearPowerCoeffs_),
    exponent_("rittmannExpon", dimless, rittmannShearPowerCoeffs_)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::attachmentModels::rittmannShearPower::calcAttachment()
{   
    *ptrkatt_ = 
        Cd_
        *
        pow
        (
            (100.0 * mu_ * mag(this->U_) * pow3(1.0 - this->n_)
            / (sqr(dc_) * pow3(this->n_) * (6.0/dc_)))
            / dimensionedScalar("dynecm2ToNm2",dimPressure,10.0),
            exponent_)
        /
        1.0/dimensionedScalar("unitTime",dimTime,1.0
        );
}


bool Foam::attachmentModels::rittmannShearPower::read
(
    const dictionary& attachmentProperties
)
{
    attachmentModel::read(attachmentProperties);

    rittmannShearPowerCoeffs_ = attachmentProperties.optionalSubDict(typeName + "Coeffs");
    
    rittmannShearPowerCoeffs_.lookup("collectorSize") >> dc_;
    rittmannShearPowerCoeffs_.lookup("fluidViscosity") >> mu_;
    rittmannShearPowerCoeffs_.lookup("rittmannCoeff") >> Cd_;
    rittmannShearPowerCoeffs_.lookup("rittmannExpon") >> exponent_;

    return true;
}


// ************************************************************************* //
