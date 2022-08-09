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

#include "colloidFiltrationTheory.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace attachmentModels
{
    defineTypeNameAndDebug(colloidFiltrationTheory, 0);

    addToRunTimeSelectionTable
    (
        attachmentModel,
        colloidFiltrationTheory,
        dictionary
    );
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //
void Foam::attachmentModels::colloidFiltrationTheory::updateLambda_()
{
    lambda_ = 
        min(
            3.0 * (1.0 - this->n_) * alpha_ * eta_ / (2.0 * dc_),
            dimensionedScalar(
                "maxLambda",
                dimensionSet(0,-1,0,0,0,0,0),
                1.0E+05)
            );
}

void Foam::attachmentModels::colloidFiltrationTheory::updateEta_()
{
    eta_ =   
        exp( // Diffusion mechanism
            0.875468737
            + (0.333 * logAs_)
            - (0.081 * logNR_)
            - (0.715 * logNP_)
            + (0.052 * logNv_) 
            )
        + 
        exp( // Interception mechanism
            - 0.597837001
            + logAs_ 
            + (1.550 * logNR_)
            - (0.125 * logNP_)
            + (0.125 * logNv_)
            )
        + 
        exp( // Gravitational deposition mechanism
            - 0.744440475
            - (1.350 * logNR_)
            - (1.110 * logNP_)
            + (0.053 * logNv_)
            + (1.110 * logNg_)
            );
}

void Foam::attachmentModels::colloidFiltrationTheory::updatelogAs_()
{
    logAs_ = 
        log( 
            8.99992117 * sqr(1.0/this->n_) 
            + 
            -7.49203318 * (1.0/this->n_) 
            + 
            0.4119361
            );
}

void Foam::attachmentModels::colloidFiltrationTheory::updatelogNP_()
{
    logNP_ = 
        log(
            ( mag(this->U_) * dc_ )
            /
            ( this->n_ * DiffusionCoef_ )
        );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::attachmentModels::colloidFiltrationTheory::colloidFiltrationTheory
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
    
    colloidFiltrationTheoryCoeffs_(attachmentProperties.optionalSubDict(typeName + "Coeffs")),
    dc_("collectorSize", dimLength, colloidFiltrationTheoryCoeffs_),
    dp_("particleSize", dimLength, colloidFiltrationTheoryCoeffs_),
    rhof_("fluidDensity", dimDensity, colloidFiltrationTheoryCoeffs_),
    rhop_("particleDensity", dimDensity, colloidFiltrationTheoryCoeffs_),
    mu_("fluidViscosity", dimDynamicViscosity, colloidFiltrationTheoryCoeffs_),
    alpha_("alphaEfficiency", dimless, colloidFiltrationTheoryCoeffs_),
    Hamaker_("hamakerConst", dimensionSet(1,2,-2,0,0,0,0), colloidFiltrationTheoryCoeffs_),
    Temp_("refTemperature", dimTemperature, colloidFiltrationTheoryCoeffs_),
    DiffusionCoef_(KBOLTZ * Temp_ / (3.0 * MATH_PI * mu_ * dp_)),

    logAs_(
        IOobject
        (
            "logAs",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("logAs", dimless,SMALL)),

    logNR_( log(dp_/dc_) ),

    logNP_(       
        IOobject
        (
            "logNP",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("logNP", dimless,SMALL)),
    
    logNv_( log(Hamaker_/(KBOLTZ*Temp_)) ),
    logNg_( 
        log(
            MATH_PI * pow(dp_,4) 
            * (rhop_ - rhof_) 
            * dimensionedScalar("g",dimAcceleration,9.81) 
            / ( 12.0 * KBOLTZ * Temp_ )
            )),

    lambda_(
        IOobject
        (
            "Lambda",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("Lambda", dimensionSet(0,-1,0,0,0,0,0) , SMALL)),

    eta_(
        IOobject
        (
            "eta",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("eta", dimless, SMALL))
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::attachmentModels::colloidFiltrationTheory::calcAttachment()
{   
    updatelogAs_();
    updatelogNP_();
    updateEta_();
    updateLambda_();
    *ptrkatt_ = lambda_ * mag(this->U_) / this->n_;

    // Info<< "Debugging CFT " << nl;
    // Info<< "********************** " << nl;
    // Info<< "Diffusion coeff" << DiffusionCoef_ <<endl;
    // Info<< "AS \n" << exp(logAs_) << endl;
    // Info<< "NR \n" << exp(logNR_) << endl;
    // Info<< "NP \n" << exp(logNP_) << endl;
    // Info<< "Nv \n" << exp(logNv_) << endl;
    // Info<< "Ng \n" << exp(logNg_) << endl;

    // Info<< "ETA \n" << eta_ << endl;
    // Info<< "LAMBDA \n" << lambda_ << endl;
    // Info<< "POROSITY \n" << n_ << endl;
    // Info<< "VELOCITY \n" << U_ << endl;
}


bool Foam::attachmentModels::colloidFiltrationTheory::read
(
    const dictionary& attachmentProperties
)
{
    attachmentModel::read(attachmentProperties);

    colloidFiltrationTheoryCoeffs_ = attachmentProperties.optionalSubDict(typeName + "Coeffs");
    
    colloidFiltrationTheoryCoeffs_.lookup("collectorSize") >> dc_;
    colloidFiltrationTheoryCoeffs_.lookup("particleSize") >> dp_;
    colloidFiltrationTheoryCoeffs_.lookup("fluidDensity") >> rhof_;
    colloidFiltrationTheoryCoeffs_.lookup("particleDensity") >> rhop_;
    colloidFiltrationTheoryCoeffs_.lookup("fluidViscosity") >> mu_;
    colloidFiltrationTheoryCoeffs_.lookup("alphaEfficiency") >> alpha_;
    colloidFiltrationTheoryCoeffs_.lookup("hamakerConst") >> Hamaker_;
    colloidFiltrationTheoryCoeffs_.lookup("refTemperature") >> Temp_;

    return true;
}


// ************************************************************************* //
