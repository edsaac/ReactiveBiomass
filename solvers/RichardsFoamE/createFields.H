//  Hydraulic head H = h + z [m]
//  We only use h and explicitely handle z

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

//  For Picard loops
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

//  Velocity vector field [m/s]
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

//- Face fluxes [m³/s]
#include "createPhi.H"

//- Water contents and effective water content [-]
Info<< "Initialize water content field\n" << endl;
volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("theta", dimensionSet(0,0,0,0,0,0,0) , scalar(0.0))
);

//- Permeability and hydraulic conductivity [m/s]
Info<< "Initialize hydraulic conductivity\n" << endl;
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
    dimensionedScalar("hydrConduct",dimensionSet(0,1,-1,0,0),scalar(1.0))
);

////////////////// Constant volScalarFields ////////////////
volScalarField theta_s
(
    IOobject
    (
        "theta_s",
        runTime.constantPath()/"soilParameters",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField theta_r
(
    IOobject
    (
        "theta_r",
        runTime.constantPath()/"soilParameters",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField K_s
(
    IOobject
    (
        "K_s",
        runTime.constantPath()/"soilParameters",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.constantPath()/"soilParameters",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField n
(
    IOobject
    (
        "n",
        runTime.constantPath()/"soilParameters",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField m
(
    1.0 - 1.0/n
);
