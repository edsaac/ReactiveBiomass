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

\*---------------------------------------------------------------------------*/

#include "unsatTimeFvPatchScalarField.H"
#include "fvPatch.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::unsatTimeFvPatchScalarField::unsatTimeFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    uniformValue_()
{}

// Foam::unsatTimeFvPatchScalarField::unsatTimeFvPatchScalarField
// (
//     const fvPatch& p,
//     const DimensionedField<scalar, volMesh>& iF,
//     const Field<scalar>& fld
// )
// :
//     fixedValueFvPatchScalarField(p, iF, fld),
//     uniformValue_()
// {}

Foam::unsatTimeFvPatchScalarField::unsatTimeFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    uniformValue_(Function1<scalar>::New("uniformValue", dict))
{
    this->evaluate();
    updateCoeffs();
}

Foam::unsatTimeFvPatchScalarField::unsatTimeFvPatchScalarField
(
    const unsatTimeFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper, false), // Don't map
    uniformValue_(ptf.uniformValue_, false)
{
    // Evaluate since value not mapped
    this->evaluate();
    updateCoeffs();
}


Foam::unsatTimeFvPatchScalarField::unsatTimeFvPatchScalarField
(
    const unsatTimeFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    uniformValue_(ptf.uniformValue_, false)
{}


Foam::unsatTimeFvPatchScalarField::unsatTimeFvPatchScalarField
(
    const unsatTimeFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    uniformValue_(ptf.uniformValue_, false)
{
    // Evaluate the profile if defined
    if (ptf.uniformValue_.valid())
    {
        this->evaluate();
        updateCoeffs();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::unsatTimeFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    // Initialize a field for the calculated h
    scalar value = uniformValue_->value(t);
	
    scalarField value_as_field(this->patch().size(), value);
    

    // fvPatchField<Type>::operator==(uniformValue_->value(t));
    // fixedValueFvPatchField<Type>::operator==(uniformValue_->value(t));   
    // fixedValueFvPatchField<Type>::updateCoeffs();

    // set the value_ of this patch to the newly computed flow speed
    this->operator==(value_as_field);

    Foam::Info << "updateCoeffs: Is this always zero?: " << t << endl;
    Foam::Info << "updateCoeffs: what does this evaluate?: " << value << endl;
    // call the base class method to make sure all the other bits and pieces get updated
    // fixedValueFvPatchScalarField::updateCoeffs();
    fixedValueFvPatchScalarField::updateCoeffs();

}


void Foam::unsatTimeFvPatchScalarField::write(Ostream& os) const
{
    const scalar t = this->db().time().timeOutputValue();
    Foam::Info << "Is this always zero?: " << t << endl;
    Foam::Info << "what does this evaluate?: " << uniformValue_->value(t) << endl;

    fvPatchScalarField::write(os);
    // os.writeKeyword("value") << "uniform " << uniformValue_->value(t) << token::END_STATEMENT << nl;
    writeEntry(os, uniformValue_());
    writeEntry(os, "value", *this);
    // writeEntry(os, "value", value_text << uniformValue_->value(t));
}

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        unsatTimeFvPatchScalarField
    );
}
// ************************************************************************* //
