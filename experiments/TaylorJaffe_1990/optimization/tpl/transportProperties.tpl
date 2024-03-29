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
    class       dictionary;
    location    "constant";
    object      transportProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// [kg m s K mol * *]

// Densities
rho       rho       [1 -3 0 0 0 0 0]    1000.0;
rho_X     rho_X     [1 -3 0 0 0 0 0]    <<RHO_X>>;

//Other physical stuff
g         g         [0 1 -2 0 0 0 0]    9.81;
mu        mu        [1 -1 -1 0 0 0 0]   0.001;

//Molecular diffusion and dispersion
/*
    Molecular diffusion values at 
         https://www.aqion.de/site/diffusion-coefficients

    For methanol in water check:
        Derlacki et al. (1985) https://doi.org/10.1021/j100270a039
            Table 1 : 1.563E-9 m²/s @(298.15 K and 0.1MPa)
*/
mDif       1.563E-9;   //
LongDisp   0.001;
TransDisp  0.0001;

molDiff    molDiff    [0 2 -1 0 0 0 0]
    ($mDif 0.0 0.0 
     0.0 $mDif 0.0 
     0.0 0.0 $mDif);

DispTensor DispTensor [0 1 0 0 0 0 0]
    ($LongDisp $TransDisp $TransDisp 
     $TransDisp $LongDisp $TransDisp 
     $TransDisp $TransDisp $LongDisp);


// ************************************************************************* //
