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
    object      attachmentProperties_BAP;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// [kg m s K mol * *]

attachmentModel colloidFiltrationTheory;

colloidFiltrationTheoryCoeffs
{
    collectorSize   2.0E-3;
    particleSize    10.0E-6;
    fluidDensity    999.79;
    particleDensity 1050.0;
    fluidViscosity  0.0008891;
    alphaEfficiency <<ATT_BAP>>;
    hamakerConst    5.0E-21;
    refTemperature  283.0;
}

// attachmentModel constantRate;

// constantRateCoeffs
// {
//     value   1.0E-6;
// }

// ************************************************************************* //
