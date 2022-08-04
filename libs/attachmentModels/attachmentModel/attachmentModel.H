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
    Foam::attachmentModels

Description
    A namespace for various attachment rate model implementations.

    Implemented models in:
    https://www.desmos.com/calculator/sg2jg37j2n

Class
    Foam::attachmentModels::attachmentModel

SourceFiles
    attachmentModel.C
    attachmentModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef attachmentModel_H
#define attachmentModel_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"
#include "Time.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class attachmentModel Declaration
\*---------------------------------------------------------------------------*/

class attachmentModel
{

protected:

    // Protected data

    word name_;
    dictionary attachmentProperties_;  
    volScalarField* ptrkatt_;
    const fvMesh& mesh_;
    const Time& runTime_;
    const volVectorField& U_;
    const volScalarField& n_;

public:
    
    //- Runtime type information
    TypeName("attachmentModel");

    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        attachmentModel,
        dictionary,
        (
            const word& name,
            const dictionary& attachmentProperties,
            volScalarField* ptrkatt,
            const fvMesh& mesh,
            const Time& runTime,
            const volVectorField& U,
            const volScalarField& n
        ),
        (name, attachmentProperties, ptrkatt, mesh, runTime, U, n)
    );


    // Constructors

    //- Construct from components
    attachmentModel
    (
        const word& name,
        const dictionary& attachmentProperties,
        volScalarField* ptrkatt,
        const fvMesh& mesh,
        const Time& runTime,
        const volVectorField& U,
        const volScalarField& n
    );

    //- Disallow default bitwise copy construction
    attachmentModel(const attachmentModel&);

    // Selectors

    //- Return a reference to the selected model
    static autoPtr<attachmentModel> New
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
    virtual ~attachmentModel()
    {}

    // Member Functions

        //- Return the clogging properties dictionary
        const dictionary& attachmentProperties() const
        {
            return attachmentProperties_;
        }

        //- Update the attchment rate field
        virtual void calcAttachment() = 0;

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& attachmentProperties) = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const attachmentModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
