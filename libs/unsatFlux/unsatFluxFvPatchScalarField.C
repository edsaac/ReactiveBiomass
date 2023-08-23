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

#include "unsatFluxFvPatchScalarField.H"
#include "unsaturatedCalculator.H"
#include "fvPatch.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::unsatFluxFvPatchScalarField::unsatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    // NOTE: call the default constructor to make sure everything gets initialised properly
    fixedValueFvPatchScalarField(p, iF),
   
    // NOTE: assign default values to the members using an initialiser list
    targetFlow_(0.)
{}

Foam::unsatFluxFvPatchScalarField::unsatFluxFvPatchScalarField
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
    targetFlow_(0.)
{
    // NOTE: calls the = operator to assign the value to the faces held by this BC
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    // NOTE: looks up the necessary paramters
    // approximationType_ = dict.lookupOrDefault<word>("approximationType","exponential");
    dict.lookup("targetFlow") >> targetFlow_;
	// dict.lookup("deltaByR") >> deltaByR_;
	// centrepoint_ = dict.lookupOrDefault<vector>("centrepoint",vector::zero);
	// dict.lookup("R") >> R_;
	// lambda_ = dict.lookupOrDefault<scalar>("lambda",0.);

    // NOTE: calls the .updateCoeffs() method to calculate the inlet profile in
    // accordance with the controls which have just been read.
	updateCoeffs();
}

Foam::unsatFluxFvPatchScalarField::unsatFluxFvPatchScalarField
(
    const unsatFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    // NOTE: this constructor, and the two subsequent ones, transfer data to the
    // instance being created from another one.
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    targetFlow_(ptf.targetFlow_)
    // approximationType_(ptf.approximationType_),
    // flowSpeed_(ptf.flowSpeed_),
	// deltaByR_(ptf.deltaByR_),
	// centrepoint_(ptf.centrepoint_),
	// R_(ptf.R_),
	// lambda_(ptf.lambda_)
{}

Foam::unsatFluxFvPatchScalarField::unsatFluxFvPatchScalarField
(
    const unsatFluxFvPatchScalarField& rifvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(rifvpvf, iF),
    targetFlow_(rifvpvf.targetFlow_)
    // approximationType_(rifvpvf.approximationType_),
    // flowSpeed_(rifvpvf.flowSpeed_),
    // deltaByR_(rifvpvf.deltaByR_),
    // centrepoint_(rifvpvf.centrepoint_),
    // R_(rifvpvf.R_),
    // lambda_(rifvpvf.lambda_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// NOTE: this is the key method which implements the actual maths for calculating
// the inlet profiles.
void Foam::unsatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    label target_patch = mesh.boundaryMesh().findPatchID("top");
    
    // Get the other fields
    const volScalarField& Ks = this->db().objectRegistry::lookupObject<volScalarField>("K_0");    
    const volScalarField& alpha = this->db().objectRegistry::lookupObject<volScalarField>("alpha");
    const volScalarField& n_vang = this->db().objectRegistry::lookupObject<volScalarField>("n_vangenucthen");
    const volScalarField& Swr = this->db().objectRegistry::lookupObject<volScalarField>("Sw_r");
    const volScalarField& Sws = this->db().objectRegistry::lookupObject<volScalarField>("Sw_s");
    const volScalarField& perm_clogging = this->db().objectRegistry::lookupObject<volScalarField>("perm_clog");

    // Get their boundaries
    const fvPatchScalarField& Ks_bc = Ks.boundaryField()[target_patch];
    const fvPatchScalarField& alpha_bc = alpha.boundaryField()[target_patch];
    const fvPatchScalarField& n_vang_bc = n_vang.boundaryField()[target_patch];
    const fvPatchScalarField& Swr_bc = Swr.boundaryField()[target_patch];
    const fvPatchScalarField& Sws_bc = Sws.boundaryField()[target_patch];
    const fvPatchScalarField& perm_clogging_bc = perm_clogging.boundaryField()[target_patch];

	scalarField h_solution(this->patch().size(), 1.0);

    // go over each face and add the BL profile for faces close to the wall
	forAll(h_solution, faceI)
	{
        // Build soil top
        UnsaturatedSoilTop soiltop(
            Sws_bc[faceI],
            Swr_bc[faceI],
            alpha_bc[faceI],
            n_vang_bc[faceI],
            Ks_bc[faceI] * perm_clogging_bc[faceI], // Ks but penalized by clogging
            targetFlow_ 
        );
        
        Info << "Soil state:" << nl;
        Info << "Ks [m/s]: " << soiltop.Ks << endl;
        Info << "k(n): " << perm_clogging_bc[faceI] << endl;
        
        h_solution[faceI] = soiltop.bisection_h_from_target(-10.0, 0.0);
        Info << "Bisection solution h: " << h_solution[faceI] << " m" << endl;
        // h_solution = Ks_bc[faceI];
	}

	// set the value_ of this patch to the newly computed flow speed
    this->operator==(h_solution);

    // call the base class method to make sure all the other bits and pieces get updated
    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::unsatFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("targetFlow") << targetFlow_ << token::END_STATEMENT << nl;
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        unsatFluxFvPatchScalarField
    );
}

// ************************************************************************* //
