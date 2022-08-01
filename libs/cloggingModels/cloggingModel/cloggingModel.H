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

Namespace
    Foam::cloggingModels

Description
    A namespace for various incompressible cloggingModel implementations.

Class
    Foam::cloggingModel

Description
    An abstract base class for incompressible cloggingModels.

    All clogging models must be bounded to a minimum and reference
    porosity and permeability, such that:

       k(n)
        |
   kRef |             ______
        |          .-'  
        |       .-'     
        |    .-'
   kMin |___/ 
        |___,_________,___ n
            nMin     nRef

    Implemented models in:
    https://www.desmos.com/calculator/sg2jg37j2n

SourceFiles
    cloggingModel.C
    cloggingModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef cloggingModel_H
#define cloggingModel_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class cloggingModel Declaration
\*---------------------------------------------------------------------------*/

class cloggingModel
{

protected:

    // Protected data

    word name_;
    dictionary cloggingProperties_;
    const volScalarField& n_;
    volScalarField* ptrperm_;
    
    dimensionedScalar permMin_;
    dimensionedScalar permRef_;
    dimensionedScalar nMin_;
    dimensionedScalar nRef_;

public:
    
    //- Runtime type information
    TypeName("cloggingModel");

    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        cloggingModel,
        dictionary,
        (
            const word& name,
            const dictionary& cloggingProperties,
            const volScalarField& n,
            volScalarField* ptrperm
        ),
        (name, cloggingProperties, n, ptrperm)
    );


    // Constructors

    //- Construct from components
    cloggingModel
    (
        const word& name,
        const dictionary& cloggingProperties,
        const volScalarField& n,
        volScalarField* ptrperm
    );

    //- Disallow default bitwise copy construction
    cloggingModel(const cloggingModel&);

    // Selectors

    //- Return a reference to the selected model
    static autoPtr<cloggingModel> New
    (
        const word& name,
        const dictionary& cloggingProperties,
        const volScalarField& n,
        volScalarField* ptrperm
    );


    //- Destructor
    virtual ~cloggingModel()
    {}

    // Member Functions

        //- Return the clogging properties dictionary
        const dictionary& cloggingProperties() const
        {
            return cloggingProperties_;
        }

        //- Update the permeability field
        virtual void calcPerm() = 0;

        //- Getters for reference values
        inline dimensionedScalar permRef()
        {
            return permRef_;
        }

        inline dimensionedScalar permMin()
        {
            return permMin_;
        }

        inline dimensionedScalar nRef()
        {
            return nRef_;
        }

        inline dimensionedScalar nMin()
        {
            return nMin_;
        }

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& cloggingProperties) = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cloggingModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
