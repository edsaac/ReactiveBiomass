/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  7
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0.000";
    object      h;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*
    For the low flowrate experiment
        Q_in = 1 mL/min = 1.667E-8 m³/s
    
    The column's cross sectional area is
        A = π D²/4 = π (0.080 m)²/4 = 5.03E-3 m²
    
    Targeted constant flow is 
        q_in = Q_in/A = 3.316E-6 m/s
*/

dimensions      [0 1 0 0 0 0 0];

internalField   uniform <<h_INITIAL_CONDITION>>;

boundaryField
{
    top
    {
        // //- ♣ Fixed head value
        type            fixedValue;
        value           uniform <<h_TOP_BOUNDARY>>;

        // //- ♣ Fixed gradient
        // type            fixedGradient;
        // gradient        uniform -1.0001;

        // //- ♣ Fixed mixed
        // type            mixed;
        // refValue        uniform -0.69;
        // refGradient     uniform -1.0001;
        // valueFraction   uniform 0.5;          // If 0: fixedGradient, if 1: fixedValue


        // //- ♣ Fixed Darcy Flux

        // type            codedMixed;
        // refValue        uniform 0;
        // refGradient     uniform 0;
        // valueFraction   uniform 0;          // If 0: fixedGradient, if 1: fixedValue
        // name            DarcyInflux;        // name of generated BC
    
        // code
        // #{
        //     // Read the hydraulic conductivity field
        //     const volScalarField& myK = this->db().objectRegistry::lookupObject<volScalarField>("hydrConduct");
            
        //     // Specify the target Darcy flow velocity (Pay attention to the direction)
        //     dimensionedScalar targetDarcyFlow ("targetDarcyFlow", dimVelocity, scalar(-3.32E-6));
            
        //     // Assign the corresponding gradient to each boundary face 
        //     forAll(patch(), faceI)
        //     { 
        //         Info << myK[faceI] << endl;
        //         this->refGrad()[faceI] = - (1.0 + targetDarcyFlow.value()/ myK[faceI]);
        //     }
            
        //     //this->refGrad() = Zero;
        //     //this->valueFraction() = 1.0;
        // #};
    }
    
    bottom
    {
        //- ♣ Fixed head value
        // type            fixedValue;
        // value           uniform 0.0;

        //- ♣ Fixed gradient
        /*- No flow boundary: set equal to 1.0
            Free flow boundary: set equal to 0.0
        */
        type            fixedGradient;
        gradient        uniform 0.0;
    }

    "(left|right|front|back)"
    {
        type            empty;
    }
}


// ************************************************************************* //
