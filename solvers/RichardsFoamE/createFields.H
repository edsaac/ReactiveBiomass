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
        "waterContent",
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
        "waterContent_effective",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("theta_e",dimless,scalar(0.0))
);

volScalarField capillary
(
    IOobject
    (
        "saturation",
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
        "relativePermeability",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("perm",dimensionSet(0,2,0,0,0),scalar(0.0))
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