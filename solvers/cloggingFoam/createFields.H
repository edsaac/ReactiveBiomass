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
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField hydraulicCond (perm_0 * rho * g / mu);

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
    - hydraulicCond * fvc::grad(h) 
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


// Calculate diffusion coefficient from Einstein-Stokes equation
//   >> use this in the Péclet number calculation
dimensionedScalar D_diffusion ( (kBoltz*Temp)/(3.0*PI*mu*dp) );

Foam::Info << "D_diff: ("<< D_diffusion.value() <<")\n" << endl;


// Calculate deposited mass and clogLimiter
Foam::Info << "Calculate total biomass and clogLimiter \n" << endl;
volScalarField XMAX ((n_0 - n_c) * rho_X);
volScalarField clogLimiter (1.0 - depositedClay/XMAX);
clogLimiter *= Foam::pos(clogLimiter);


#include "calculateCFT.H"

// Field of filtration coefficient
volScalarField Lambda
(
IOobject
(
    "Lambda",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
),
(3.0 * (1.0-n) * alphaCFT * eta)/(2.0*ds)
);

volScalarField katt (Lambda * mag(U));

#include "createPhi.H"
#include "CourantNo.H"
#include "createFvOptions.H"
