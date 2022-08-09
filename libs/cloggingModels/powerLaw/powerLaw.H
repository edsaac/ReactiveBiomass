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
    Foam::cloggingModel::powerLaw

Description
     Standard power-law clogging model.

     perm = (permRef - permMin) * ((n-nMin)/(nRef-nMin))^nExponent
          + permMin

Usage
    In the cloggingProperties dictionary

    /verbatim
        cloggingModel powerLaw;

        powerLawCoeffs
        {
            nExponent dimless 3.0;
        }
    /endverbatim

SourceFiles
    powerLaw.C

\*---------------------------------------------------------------------------*/

#ifndef powerLaw_H
#define powerLaw_H

#include "cloggingModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cloggingModels
{

/*---------------------------------------------------------------------------*\
                           Class powerLaw Declaration
\*---------------------------------------------------------------------------*/

class powerLaw
:
    public cloggingModel
{
    // Private Data

        dictionary powerLawCoeffs_;
        dimensionedScalar nExponent_;

    // Private Member Functions

public:

    //- Runtime type information
    TypeName("powerLaw");

    // Constructors

        //- Construct from components
        powerLaw
        (
            const word& name,
            const dictionary& cloggingProperties,
            const volScalarField& n,
            volScalarField* ptrperm
        );


    //- Destructor
    virtual ~powerLaw()
    {}


    // Member Functions

        //- Update the permeability field
        virtual void calcPerm();

        //- Read cloggingProperties dictionary
        virtual bool read(const dictionary& cloggingProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
