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
    object      detachmentProperties_EPS;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// [kg m s K mol * *]

attachmentModel rittmannShearPower;

rittmannShearPowerCoeffs
{
    collectorSize   2.0E-3;
    fluidViscosity  0.0008891;
    rittmannCoeff   <<DET_EPS>>;
    rittmannExpon   0.58;
}

// attachmentModel constantRate;

// constantRateCoeffs
// {
//     value   0; //0.6E-7;
// }

// ************************************************************************* //
