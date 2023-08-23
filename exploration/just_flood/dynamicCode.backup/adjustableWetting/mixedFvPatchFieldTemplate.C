/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "mixedFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude
#line 139 "/home/edsaac/Repos/ReactiveBiomass/exploration/just_flood/coded_twenty/0.000/h.boundaryField.top"
#include "fvCFD.H"
            #include <iostream>
            #include <math.h>
            #include "unsaturatedCalculator.H"
//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = fcd8aab833ad6398d0abc01ad576bf6322ca495a
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void adjustableWetting_fcd8aab833ad6398d0abc01ad576bf6322ca495a(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    adjustableWettingMixedValueFvPatchScalarField
);


const char* const adjustableWettingMixedValueFvPatchScalarField::SHA1sum =
    "fcd8aab833ad6398d0abc01ad576bf6322ca495a";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjustableWettingMixedValueFvPatchScalarField::
adjustableWettingMixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct adjustableWetting sha1: fcd8aab833ad6398d0abc01ad576bf6322ca495a"
            " from patch/DimensionedField\n";
    }
}


adjustableWettingMixedValueFvPatchScalarField::
adjustableWettingMixedValueFvPatchScalarField
(
    const adjustableWettingMixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct adjustableWetting sha1: fcd8aab833ad6398d0abc01ad576bf6322ca495a"
            " from patch/DimensionedField/mapper\n";
    }
}


adjustableWettingMixedValueFvPatchScalarField::
adjustableWettingMixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct adjustableWetting sha1: fcd8aab833ad6398d0abc01ad576bf6322ca495a"
            " from patch/dictionary\n";
    }
}


adjustableWettingMixedValueFvPatchScalarField::
adjustableWettingMixedValueFvPatchScalarField
(
    const adjustableWettingMixedValueFvPatchScalarField& ptf
)
:
    mixedFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct adjustableWetting sha1: fcd8aab833ad6398d0abc01ad576bf6322ca495a"
            " as copy\n";
    }
}


adjustableWettingMixedValueFvPatchScalarField::
adjustableWettingMixedValueFvPatchScalarField
(
    const adjustableWettingMixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct adjustableWetting sha1: fcd8aab833ad6398d0abc01ad576bf6322ca495a "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

adjustableWettingMixedValueFvPatchScalarField::
~adjustableWettingMixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy adjustableWetting sha1: fcd8aab833ad6398d0abc01ad576bf6322ca495a\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void adjustableWettingMixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs adjustableWetting sha1: fcd8aab833ad6398d0abc01ad576bf6322ca495a\n";
    }

//{{{ begin code
    #line 74 "/home/edsaac/Repos/ReactiveBiomass/exploration/just_flood/coded_twenty/0.000/h.boundaryField.top"
// Read the soil parameters
            const volScalarField& Ks = this->db().objectRegistry::lookupObject<volScalarField>("K_0");
            const volScalarField& alpha = this->db().objectRegistry::lookupObject<volScalarField>("alpha");
            const volScalarField& n_vang = this->db().objectRegistry::lookupObject<volScalarField>("n_vangenucthen");
            const volScalarField& Swr = this->db().objectRegistry::lookupObject<volScalarField>("Sw_r");
            const volScalarField& Sws = this->db().objectRegistry::lookupObject<volScalarField>("Sw_s");

            // Read current clogging state
            const volScalarField& perm_clogging = this->db().objectRegistry::lookupObject<volScalarField>("perm_clog");
  
            // Call geometry
            const fvPatch& boundaryPatch = patch(); 

            // Initialize field
            const vectorField& Cf = boundaryPatch.Cf(); 
            scalarField& field = *this; 

            // Call time
            // scalar t = this->db().time().value();
            const scalar qtarget = 3.36E-6;
            scalar solution;

            // Loop through all faces in the boundary
            forAll(Cf, faceI)
            {
                // Build soil top
                UnsaturatedSoilTop soiltop(
                    Sws[faceI],
                    Swr[faceI],
                    alpha[faceI],
                    n_vang[faceI],
                    Ks[faceI],
                    qtarget * perm_clogging[faceI] // Ks but penalized by clogging
                );

                // Info << "Location: "<< nl << endl;
                // x = Cf[faceI].x();
                
                Info << "Soil state:" << nl;
                Info << "Ks [m/s]: " << soiltop.Ks << endl;
                Info << "k(n): " << perm_clogging[faceI] << endl;
                
                solution = soiltop.bisection_h_from_target(-10.0, 0.0);
                Info << "Bisection solution h: " << solution << " m" << endl;
                
                // Assign value to h
                field[faceI] = solution;
            }

            // Assign uniform value to value
            // operator==(field);
            // this->updateCoeffs();
            this->refValue() = field;
//}}} end code

    this->mixedFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

