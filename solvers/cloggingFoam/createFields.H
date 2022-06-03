// Read the transport properties dictionary
#define PI 3.14169264
#include "readParameters.H"

Foam::Info << "Reading field: Hydraulic head (h) \n" << endl;
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

Foam::Info << "Reading field: porosity (n) \n" << endl;
volScalarField n
(
    IOobject
    (
        "n",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField perm
(
    IOobject
    (
        "perm",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField hydraulicCond (perm * rho * g / mu);

// U in this code refers to the Darcy velocity (a.k.a q)
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
    dimensionedVector("U", dimensionSet(0,1,-1,0,0,0,0) , vector(0.0, 0.0, 0.0))
);


Foam::Info << "Reading field: Suspended clay \n" << endl;
volScalarField suspendedClay
(
    IOobject
    (
        "suspendedClay",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::Info << "Reading field: deposited clay \n" << endl;
volScalarField depositedClay
(
    IOobject
    (
        "depositedClay",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// Calculate deposited mass and clogLimiter
Foam::Info << "Calculate total biomass and clogLimiter \n" << endl;
dimensionedScalar XMAX ((n_0 - n_c) * rho_X);
volScalarField clogLimiter (1.0 - depositedClay/XMAX);
clogLimiter *= Foam::pos(clogLimiter);

#include "initializeCFT.H"
#include "createPhi.H"
#include "createFvOptions.H"
