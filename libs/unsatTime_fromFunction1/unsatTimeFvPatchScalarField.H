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
    Foam::unsatTimeFvPatchScalarField

Description
    This boundary condition provides a uniform fixed value condition.

Usage
    \table
        Property     | Description             | Required    | Default value
        uniformValue | uniform value           | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            unsatTime;
        uniformValue    constant 0.2;
    }
    \endverbatim

Note
    The uniformValue entry is a Function1 type, able to describe time
    varying functions.  The example above gives the usage for supplying a
    constant value.

See also
    Foam::Function1Types
    Foam::fixedValueFvPatchField

SourceFiles
    unsatTimeFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef unsatTimeFvPatchScalarField_H
#define unsatTimeFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Function1.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class unsatTimeFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class unsatTimeFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        autoPtr<Function1<scalar>> uniformValue_;


public:

    //- Runtime type information
    TypeName("unsatTime");


    // Constructors

        //- Construct from patch and internal field
        unsatTimeFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        // //- Construct from patch and internal field and patch field
        // unsatTimeFvPatchScalarField
        // (
        //     const fvPatch&,
        //     const DimensionedField<scalar, volMesh>&,
        //     const Field<scalar>& fld
        // );

        //- Construct from patch, internal field and dictionary
        unsatTimeFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given unsatTimeFvPatchScalarField
        //  onto a new patch
        unsatTimeFvPatchScalarField
        (
            const unsatTimeFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        unsatTimeFvPatchScalarField
        (
            const unsatTimeFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new unsatTimeFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        unsatTimeFvPatchScalarField
        (
            const unsatTimeFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new unsatTimeFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
//     #include "unsatTimeFvPatchScalarField.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
