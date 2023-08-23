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

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude
#line 131 "/home/edsaac/Repos/ReactiveBiomass/exploration/just_flood/coded_twenty/0.000/h.boundaryField.top"
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
    // SHA1 = 275742c1b689f944bb554834bffc7e16d65ff3e8
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void adjustableWetting_275742c1b689f944bb554834bffc7e16d65ff3e8(bool load)
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
    adjustableWettingFixedValueFvPatchScalarField
);


const char* const adjustableWettingFixedValueFvPatchScalarField::SHA1sum =
    "275742c1b689f944bb554834bffc7e16d65ff3e8";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjustableWettingFixedValueFvPatchScalarField::
adjustableWettingFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct adjustableWetting sha1: 275742c1b689f944bb554834bffc7e16d65ff3e8"
            " from patch/DimensionedField\n";
    }
}


adjustableWettingFixedValueFvPatchScalarField::
adjustableWettingFixedValueFvPatchScalarField
(
    const adjustableWettingFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct adjustableWetting sha1: 275742c1b689f944bb554834bffc7e16d65ff3e8"
            " from patch/DimensionedField/mapper\n";
    }
}


adjustableWettingFixedValueFvPatchScalarField::
adjustableWettingFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct adjustableWetting sha1: 275742c1b689f944bb554834bffc7e16d65ff3e8"
            " from patch/dictionary\n";
    }
}


adjustableWettingFixedValueFvPatchScalarField::
adjustableWettingFixedValueFvPatchScalarField
(
    const adjustableWettingFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct adjustableWetting sha1: 275742c1b689f944bb554834bffc7e16d65ff3e8"
            " as copy\n";
    }
}


adjustableWettingFixedValueFvPatchScalarField::
adjustableWettingFixedValueFvPatchScalarField
(
    const adjustableWettingFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct adjustableWetting sha1: 275742c1b689f944bb554834bffc7e16d65ff3e8 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

adjustableWettingFixedValueFvPatchScalarField::
~adjustableWettingFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy adjustableWetting sha1: 275742c1b689f944bb554834bffc7e16d65ff3e8\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void adjustableWettingFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs adjustableWetting sha1: 275742c1b689f944bb554834bffc7e16d65ff3e8\n";
    }

//{{{ begin code
    #line 71 "/home/edsaac/Repos/ReactiveBiomass/exploration/just_flood/coded_twenty/0.000/h.boundaryField.top"
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
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

