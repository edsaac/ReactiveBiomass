/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::unsatFluxFvPatchVectorField

Group
    grpInletBoundaryConditions grpWallBoundaryConditions

Description
    Prescribes an inlet profile for a circular pipe inlet patch. Flow direction is determined
    automatically from face normals, only velocity magnitude needs to be provided.

    \heading Patch usage

    \table
        Property    		| Description             			| Required  | Default value
        targetFlow   		| Target flux [m/s]  				| yes     	|
    \endtable

    Example of the boundary condition specification:
    \verbatim
    inlet
	{
		// Calculates the head required to achieve that flux, assuming grad(h) was -1. 
	    type        unsatFlux;
	    // pipe radius
	    targetFlow  <scalar>;
	    // dummy value overwritten by the BC, used for initialisation only
	    value		$internalField;
	}
	\endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    unsatFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/
#ifndef unsatFluxFvPatchScalarField_H
#define unsatFluxFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Switch.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class prescribedPipeInletFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

// NOTE: the patch is derived from a fixedValue BC since it does prescribe
// a fixed value, although non-uniform across all of the faces, unlike the base
// class does.
class unsatFluxFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
private:
    // Private data

    	// flow velocity magnitude
    	scalar targetFlow_;

public:

    //- Runtime type information
    // NOTE: this gets used by the runtimeSelector. In other words, this is the
    // name under which OpenFOAM knows this BC.
    TypeName("unsatFlux");

    // Constructors

        //- Construct from patch and internal field
        unsatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct as copy setting internal field reference
        unsatFluxFvPatchScalarField
        (
            const unsatFluxFvPatchScalarField&,
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
                new unsatFluxFvPatchScalarField(*this, iF)
            );
        }

        //- Construct from patch, internal field and dictionary
        unsatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping given prescribedPipeInletFvPatchVectorField
        //  onto a new patch
        unsatFluxFvPatchScalarField
        (
            const unsatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );

        //- Disallow copy without setting internal field reference
        unsatFluxFvPatchScalarField(const unsatFluxFvPatchScalarField&) = delete;

    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
