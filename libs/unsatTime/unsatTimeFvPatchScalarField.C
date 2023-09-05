/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    // NOTE: call the default constructor to make sure everything gets initialised properly
    fixedValueFvPatchScalarField(p, iF),
   
    // NOTE: assign default values to the members using an initialiser list
    initialHead_(0.0),
    initialTime_(0.0),
    finalHead_(0.0),    
    finalTime_(0.0)
{}

Foam::unsatTimeFvPatchScalarField::unsatTimeFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    // NOTE: this constructor reads all of the control parameters from the boundary
    // condition definition specified in the time folder h file, imported here
    // as a dictionary reference.
    fixedValueFvPatchScalarField(p, iF),

    initialHead_(0.0),
    initialTime_(0.0),
    finalHead_(0.0),    
    finalTime_(0.0)
{
    // NOTE: calls the = operator to assign the value to the faces held by this BC
    // fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator==
        (
            scalarField(this->patch().size(), this->result_)
        );
    }
    

    // NOTE: looks up the necessary parameters
    Foam::Info << "Look up necessary parameters in BC" << endl;
    dict.lookup("initialHead") >> initialHead_;
    dict.lookup("initialTime") >> initialTime_;
    dict.lookup("finalHead") >> finalHead_;
    dict.lookup("finalTime") >> finalTime_;
    
    // NOTE: calls the .updateCoeffs() method to calculate in
    // accordance with the controls which have just been read.
	// updateCoeffs();
    fixedValueFvPatchScalarField::evaluate();
}

Foam::unsatTimeFvPatchScalarField::unsatTimeFvPatchScalarField
(
    const unsatTimeFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    // NOTE: this constructor, and the two subsequent ones, transfer data to the
    // instance being created from another one.
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),

    initialHead_(ptf.initialHead_),
    initialTime_(ptf.initialTime_),    
    finalHead_(ptf.finalHead_),  
    finalTime_(ptf.finalTime_)

{}

Foam::unsatTimeFvPatchScalarField::unsatTimeFvPatchScalarField
(
    const unsatTimeFvPatchScalarField& rifvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(rifvpvf, iF),

    initialHead_(rifvpvf.initialHead_),
    initialTime_(rifvpvf.initialTime_),    
    finalHead_(rifvpvf.finalHead_),  
    finalTime_(rifvpvf.finalTime_)

{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// NOTE: this is the key method which implements the actual maths for calculating
// the inlet profiles.
void Foam::unsatTimeFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const scalar t = this->db().time().timeOutputValue();
    this->slope_ = (finalHead_ - initialHead_)/(finalTime_ - initialTime_);
     // const fvMesh& mesh = patch().boundaryMesh().mesh();
    // label target_patch = mesh.boundaryMesh().findPatchID(this->patch().name());
    
    if (t < initialTime_){
        this->result_ = initialHead_;
    }
    else if(t > finalTime_){
        this->result_ = finalHead_;
    }
    else{
        this->result_ = slope_ * (t - initialTime_) + initialHead_;
    }

    Info << "Head boundary state:" << nl;
    Info << "Slope: " << this->slope_ << endl;
    Info << "h(t): " << this->result_ << endl;

	// set the value_ of this patch to the newly computed flow speed

    fixedValueFvPatchScalarField::operator==(
        scalarField(this->patch().size(), this->result_)
    );

    fvPatchScalarField::operator=(
        scalarField(this->patch().size(), this->result_)
    );

    // call the base class method to make sure all the other bits and pieces get updated
    fixedValueFvPatchScalarField::updateCoeffs();
    
    
    Info << "*this "  << nl << *this << endl;

}


void Foam::unsatTimeFvPatchScalarField::write(Ostream& os) const
{

    fvPatchScalarField::write(os);
    // fixedValueFvPatchScalarField::write(os);

    os.writeKeyword("initialHead") << initialHead_ << token::END_STATEMENT << nl;
    os.writeKeyword("initialTime") << initialTime_ << token::END_STATEMENT << nl;
    os.writeKeyword("finalHead") << finalHead_ << token::END_STATEMENT << nl;
    os.writeKeyword("finalTime") << finalTime_ << token::END_STATEMENT << nl;

    writeEntry(os, "value", *this);
    // os.writeKeyword("value") << result_ << token::END_STATEMENT << nl;
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        unsatTimeFvPatchScalarField
    );
}

// ************************************************************************* //
