/*---------------------------------------------------------------------------*\
            Declarations of the fields of primary variables
\*---------------------------------------------------------------------------*/

//- Total water pressure (pressure head/matrix potential) [m] at cell centers
//  primary variable of resolution for Richards equation
Info<< "Reading field h\n" << endl;
volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField h_before 
(
    IOobject
    (
        "h_before",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("h_before",dimensionSet(0,1,0,0,0,0,0),scalar(0.0))
);

volScalarField h_after = h;

// volScalarField h_tmp
// (
//     IOobject
//     (
//         "h_tmp",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("h_tmp",dimensionSet(0,1,0,0,0,0,0),scalar(0.0))
// );

volScalarField err
(
    IOobject
    (
        "h_tmp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("err",dimensionSet(0,1,0,0,0,0,0),scalar(0.0))
);

//- Velocity vector field at the center of cells [m/s]
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("U", dimensionSet(0,1,-1,0,0,0,0) , vector(0.0,0.0,0.0))
);

#include "createPhi.H"

//- Water contents and effective water content [-]
Info<< "Reading field saturation\n" << endl;
volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField theta_e
(
    IOobject
    (
        "theta_e",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("theta_e",dimless,scalar(0.0))
);

volScalarField innerAlphaPower
(
    IOobject
    (
        "innerAlphaPower",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("innerAlphaPower",dimless,scalar(0.0))
);

volScalarField capillary
(
    IOobject
    (
        "capillary",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("capillary",dimensionSet(0,-1,0,0,0),scalar(0.0))
);

//- Permeability and hydraulic conductivity [-]
Info<< "Reading relative permeability\n" << endl;
volScalarField perm
(
    IOobject
    (
        "perm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("perm",dimless,scalar(0.0))
);

volScalarField hydrConduct
(
    IOobject
    (
        "hydrConduct",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("hydrConduct",dimensionSet(0,1,-1,0,0),scalar(0.0))
);