/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Class
    Foam::cloggingModel::colloidFiltrationTheory

Description
     Attachment rate calculated from colloid filtration theory

     katt = 3/2 (1-n)/d_c q/n /alpha /eta_0

        n:   Porosity [-]
        d_c: CollectorSize [L]
        q:   Darcy flow velocity [L/T]
        /alpha: Collision efficiency [-]
        /eta_0: Single collector contact efficiency [-]

    The collector efficiency is estimated from Tufenkji & Elimelech 
    (2004):

    /eta_0 = /eta_D + /eta_I + /eta_G

        /eta_D: contact due transport by diffusion
        /eta_I: contact due transport by interception
        /eta_G: contact due transport by gravity/settling

    Each transport component follows an empirical correlation with
    non-dimensional numbers:

    /eta_D = 2.400 A_s^{1/3} N_R^{-0.081} N_P^{-0.715} N_v^{0.052}
    /eta_I = 0.550 A_s N_R^{1.55} N_P^{-0.125} N_v^{0.125}
    /eta_G = 0.475 N_R^{-135} N_P^{-1.11} N_v^{0.053} N_g^{1.11}

    A_s is the Happel parameter (with $x = (1 - n)^{1/3}$), N_R is the 
    aspect ratio, N_P is the Peclet number, N_v is the van der Waals
    number amd N_g is the gravitational number.
    
    A_s = 2 (1 - x^5) / (2 - 3x + 3x^5 - 2x^6)
    N_R = d_p/d_c
    N_P = (q d_c)/(n D_{/infty})
    N_v = A/(kT)
    N_g = 4/3 (/pi d_p^4 (/rho_p - /rho_f) g)/(kT)

        d_p: particle size (diam.) [L]
        k:   Boltzmann constant (1.380649E-23 J/K)
        T:   Temperature in K [T]
        A:   Hakamer constant [J]
        /rho_f:   Fluid density [M/L³]
        /rho_p:   Particle density [M/L³]
        g:   Gravity acceleration [L/T²]
        D_{/infty}: Diffusion coefficient [L²/T]

    The bulk diffusion coefficient is calculated following the Stokes-
    Einstein equation

    D_{/infty} = kT/(3 /pi /mu d_p)

        /mu:    Fluid dynamic viscosity [Pa × s]


Usage
    In the cloggingProperties dictionary

    /verbatim
        attachmentModel colloidFiltrationTheory;

        colloidFiltrationTheoryCoeffs
        {
            collectorSize   <value>;
            particleSize    <value>;
            fluidDensity    <value>;
            particleDensity <value>;
            fluidViscosity  <value>;
            alphaEfficiency <value>;
            hamakerConst    <value>;
            refTemperature  <value>;
        }
    /endverbatim

SourceFiles
    colloidFiltrationTheory.C

\*---------------------------------------------------------------------------*/

#ifndef colloidFiltrationTheory_H
#define colloidFiltrationTheory_H

#include "attachmentModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

#include "constants.H"
#define MATH_PI Foam::constant::mathematical::pi
#define KBOLTZ Foam::constant::physicoChemical::k

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace attachmentModels
{

/*---------------------------------------------------------------------------*\
                           Class powerLaw Declaration
\*---------------------------------------------------------------------------*/

class colloidFiltrationTheory
:
    public attachmentModel
{
    // Private Data

        dictionary colloidFiltrationTheoryCoeffs_;
        dimensionedScalar dc_;
        dimensionedScalar dp_;
        dimensionedScalar rhof_;
        dimensionedScalar rhop_;
        dimensionedScalar mu_;
        dimensionedScalar alpha_;
        dimensionedScalar Hamaker_;
        dimensionedScalar Temp_;
        dimensionedScalar DiffusionCoef_;

        volScalarField logAs_;
        const dimensionedScalar logNR_;
        volScalarField logNP_;
        const dimensionedScalar logNv_;
        const dimensionedScalar logNg_;

        volScalarField lambda_;
        volScalarField eta_;

    // Private Member Functions
  
        virtual void updatelogAs_();
        virtual void updatelogNP_();
        virtual void updateEta_();       
        virtual void updateLambda_();

public:
   
    //- Runtime type information
    TypeName("colloidFiltrationTheory");

    // Constructors

        //- Construct from components
        colloidFiltrationTheory
        (
            const word& name,
            const dictionary& attachmentProperties,
            volScalarField* ptrkatt,
            const fvMesh& mesh,
            const Time& runTime,
            const volVectorField& U,
            const volScalarField& n
        );


    //- Destructor
    virtual ~colloidFiltrationTheory()
    {}


    // Member Functions

        //- Update the permeability field
        void calcAttachment();

        //- Read attachmentProperties dictionary
        virtual bool read(const dictionary& attachmentProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
