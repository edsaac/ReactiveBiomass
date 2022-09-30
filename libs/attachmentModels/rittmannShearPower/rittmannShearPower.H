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
    Foam::cloggingModel::rittmannShearPower

Description
     Detachment rate calculated from Rittman (1982)
     and misinterpreted by Brovelli (2009)

     kdet = coeff * (0.1 * shear)^nPower

        coeff:  Empirical coefficient [-]
                 (adjusted for SI units = 8.42E-2 × / 86400)
        nPower: Empirical exponent [-]
                 (0.58 from Rittman (1982))
        shear:  Shear stress [N/m²]
                 (the 0.1 factor comes from converting dyne/cm² to N/m²
                  which was used in the original paper)
        
    The characteristic shear is calculated as:

    shear = 100 * /mu * u * (1-n)³ / (dp² * n³ * a)

        /mu: Fluid dynamic viscosity [Pa × s]
        u:   Darcy velocity  [m/s]
        n:   Porosity [-]
        a:   Specific surface area [1/m]
        dp:  Collector size diameter [m]

    Specific surface area for spheres of the same size is given by:

    a = 6/dp

Usage
    In the cloggingProperties dictionary

    /verbatim
        attachmentModel rittmannShearPower;

        rittmannShearPowerCoeffs
        {
            collectorSize   <value>;
            fluidViscosity  <value>;
            rittmannCoeff   <value>;
            rittmannExpon   <value>;
        }
    /endverbatim

SourceFiles
    rittmannShearPower.C

\*---------------------------------------------------------------------------*/

#ifndef rittmannShearPower_H
#define rittmannShearPower_H

#include "attachmentModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace attachmentModels
{

/*---------------------------------------------------------------------------*\
                           Class powerLaw Declaration
\*---------------------------------------------------------------------------*/

class rittmannShearPower
:
    public attachmentModel
{
    // Private Data

        dictionary rittmannShearPowerCoeffs_;
        dimensionedScalar dc_;
        dimensionedScalar mu_;
        dimensionedScalar Cd_;
        dimensionedScalar exponent_;

    // Private Member Functions
  
public:
   
    //- Runtime type information
    TypeName("rittmannShearPower");

    // Constructors

        //- Construct from components
        rittmannShearPower
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
    virtual ~rittmannShearPower()
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
