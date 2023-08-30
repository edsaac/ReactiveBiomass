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

    // NOTE: looks up the necessary parameters
    dict.lookup("targetFlow") >> targetFlow_;
    
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
{}

Foam::unsatFluxFvPatchScalarField::unsatFluxFvPatchScalarField
(
    const unsatFluxFvPatchScalarField& rifvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(rifvpvf, iF),
    targetFlow_(rifvpvf.targetFlow_)
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
    label target_patch = mesh.boundaryMesh().findPatchID(this->patch().name());
    
    // Get the other fields
    const volScalarField& Ks = this->db().objectRegistry::lookupObject<volScalarField>("K_0");    
    const volScalarField& alpha = this->db().objectRegistry::lookupObject<volScalarField>("alpha");
    const volScalarField& n_vang = this->db().objectRegistry::lookupObject<volScalarField>("n_vangenucthen");
    const volScalarField& Swr = this->db().objectRegistry::lookupObject<volScalarField>("Sw_r");
    const volScalarField& Sws = this->db().objectRegistry::lookupObject<volScalarField>("Sw_s");
    
    const volScalarField& perm_clogging = this->db().objectRegistry::lookupObject<volScalarField>("perm_clog");
    const volScalarField& head = this->db().objectRegistry::lookupObject<volScalarField>("h");

    // Get their boundaries
    const fvPatchScalarField& Ks_bc = Ks.boundaryField()[target_patch];
    const fvPatchScalarField& alpha_bc = alpha.boundaryField()[target_patch];
    const fvPatchScalarField& n_vang_bc = n_vang.boundaryField()[target_patch];
    const fvPatchScalarField& Swr_bc = Swr.boundaryField()[target_patch];
    const fvPatchScalarField& Sws_bc = Sws.boundaryField()[target_patch];
    const fvPatchScalarField& perm_clogging_bc = perm_clogging.boundaryField()[target_patch];

    scalar headCell = 0.0;
    scalar distanceBoundaryCell = 1.0;
    
    // Initialize a field for the calculated h
	scalarField h_solution(this->patch().size(), 1.0);

    // Go over each face and add the BL profile for faces close to the wall
	forAll(h_solution, faceI)
	{
        // Head value in the cell next to the boundary
        label cellI = patch().faceCells()[faceI];
        Info << "head[cellI]: "<< head[cellI] << endl;
        headCell = head[cellI];
        
        Info << "Distance: " << 1.0/(mesh.deltaCoeffs().boundaryField()[target_patch][faceI]) << endl;
        distanceBoundaryCell = 1.0/(2 * mesh.deltaCoeffs().boundaryField()[target_patch][faceI]);

        // Build soil top
        UnsaturatedSoilTop soiltop(
            Sws_bc[faceI],
            Swr_bc[faceI],
            alpha_bc[faceI],
            n_vang_bc[faceI],
            Ks_bc[faceI] * perm_clogging_bc[faceI], // Ks but penalized by clogging
            targetFlow_,
            headCell,
            distanceBoundaryCell
        );
        
        Info << "Soil state:" << nl;
        Info << "Ks [m/s]: " << soiltop.Ks << endl;
        Info << "k(n): " << perm_clogging_bc[faceI] << endl;
        
        h_solution[faceI] = soiltop.bisection_h_from_target(-20.0, 5.0);
        Info << "Bisection solution h: " << h_solution[faceI] << " m" << endl;
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
