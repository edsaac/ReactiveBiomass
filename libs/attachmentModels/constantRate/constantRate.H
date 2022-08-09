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
    Foam::cloggingModel::constantRate

Description
     Attachment rate is constant

     katt = Constant

Usage
    In the cloggingProperties dictionary

    /verbatim
        attachmentModel constantRate;

        constantRateCoeffs
        {
            value   <value>;
        }
    /endverbatim

SourceFiles
    constantRate.C

\*---------------------------------------------------------------------------*/

#ifndef constantRate_H
#define constantRate_H

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

class constantRate
:
    public attachmentModel
{
    // Private Data

        dictionary constantRateCoeffs_;
        dimensionedScalar value_;

    // Private Member Functions
  
public:
   
    //- Runtime type information
    TypeName("constantRate");

    // Constructors

        //- Construct from components
        constantRate
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
    virtual ~constantRate()
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
