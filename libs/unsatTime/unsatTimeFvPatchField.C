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

#include "unsatTimeFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::unsatTimeFvPatchField<Type>::unsatTimeFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    uniformValue_()
{}


template<class Type>
Foam::unsatTimeFvPatchField<Type>::unsatTimeFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& fld
)
:
    fixedValueFvPatchField<Type>(p, iF, fld),
    uniformValue_()
{}


template<class Type>
Foam::unsatTimeFvPatchField<Type>::unsatTimeFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
    uniformValue_(Function1<Type>::New("uniformValue", dict))
{
    // Added
    // fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
	// updateCoeffs();
    // End-Added

    this->evaluate();
}


template<class Type>
Foam::unsatTimeFvPatchField<Type>::unsatTimeFvPatchField
(
    const unsatTimeFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper, false), // Don't map
    uniformValue_(ptf.uniformValue_, false)
{
    // Evaluate since value not mapped
    this->evaluate();
}


template<class Type>
Foam::unsatTimeFvPatchField<Type>::unsatTimeFvPatchField
(
    const unsatTimeFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    uniformValue_(ptf.uniformValue_, false)
{}


template<class Type>
Foam::unsatTimeFvPatchField<Type>::unsatTimeFvPatchField
(
    const unsatTimeFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    uniformValue_(ptf.uniformValue_, false)
{
    // Evaluate the profile if defined
    if (ptf.uniformValue_.valid())
    {
        this->evaluate();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::unsatTimeFvPatchField<Type>::updateCoeffs()
{
    // if (this->updated())
    // {
    //     return;
    // }

    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<Type>::operator==(uniformValue_->value(t));
    fixedValueFvPatchField<Type>::operator==(uniformValue_->value(t));   
    
    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::unsatTimeFvPatchField<Type>::write(Ostream& os) const
{
    const scalar t = this->db().time().timeOutputValue();

    fvPatchField<Type>::write(os);
    os.writeKeyword("value") << "uniform " << uniformValue_->value(t) << token::END_STATEMENT << nl;
    writeEntry(os, uniformValue_());
    // writeEntry(os, "value", *this);
    // writeEntry(os, "value", value_text << uniformValue_->value(t));
}


// ************************************************************************* //
