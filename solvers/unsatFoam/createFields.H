//- Unsaturated flow constants

volScalarField Sw_s
// Maximum water saturation
(
    IOobject
    (
        "Sw_s",
        runTime.constantPath()/"soilParameters",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField Sw_r
// Residual water saturation
(
    IOobject
    (
        "Sw_r",
        runTime.constantPath()/"soilParameters",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField K_0
// Clean & saturated hydraulic conductivity
(
    IOobject
    (
        "K_0",
        runTime.constantPath()/"soilParameters",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField alpha
// van Genuchten - Capillary head scale parameter 
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

volScalarField n_vangenucthen
// van Genuchten - fitting exponent n
(
    IOobject
    (
        "n_vangenucthen",
        runTime.constantPath()/"soilParameters",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField m_vangenucthen
// van Genuchten - fitting exponent m
(
    1.0 - 1.0/n_vangenucthen
);


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

Foam::Info << "Reading field: porosity \n" << endl;
volScalarField porosity
(
    IOobject
    (
        "porosity",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::Info << "Reading field: water saturation (Sw) \n" << endl;
volScalarField Sw
(
    IOobject
    (
        "Sw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Sw", dimless, 1.0)
);

Foam::Info << "Reading field: air saturation (Sa) \n" << endl;
volScalarField Sa
(
    IOobject
    (
        "Sa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Sa", dimless, 0.0)
);

Foam::Info << "Creating field: relative permeability by saturation \n" << endl;
volScalarField perm_satu
(
    IOobject
    (
        "perm_satu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("perm_satu", dimless, 1.0)
);

Foam::Info << "Creating field: hydraulic conductivity [m/s] \n" << endl;
volScalarField hydraulicCond
(
    IOobject
    (
        "hydraulicCond",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("hydraulicCond0", dimVelocity, 1.0)
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
    dimensionedScalar("h_before", dimLength, scalar(0.0))
);

volScalarField h_after = h;

volScalarField err
(
    IOobject
    (
        "err",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("err", dimLength, scalar(0.0))
);

//- gradient fields for dK/dz 
volVectorField grad_k
(
    fvc::grad(hydraulicCond)
);

volScalarField grad_kz
( 
    grad_k.component(vector::Z) 
);

Foam::Info << "Creating field: hydraulic capacity \n" << endl;
volScalarField capillarity
(
    IOobject
    (
        "capillarity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("capillarity", dimensionSet(0,-1,0,0,0,0,0), 1.0)
);

Foam::Info << "Creating field: Darcy specific discharge [m/s] \n" << endl;
volVectorField U
(
    // U refers to the Darcy velocity (a.k.a q), not to the pore water velocity
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

//- Face fluxes [m³/s]
#include "createPhi.H"

//- Transport

Foam::Info << "Reading field: Non-reactive tracer (tracer) \n" << endl;
volScalarField tracer
(
    IOobject
    (
        "tracer",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

UnsaturatedSoilParametersField soil
// Instantiate soil
(
    &Sw_s,
    &Sw_r,
    &alpha,
    &n_vangenucthen,
    &m_vangenucthen,
    &perm_satu, // <-- It will modify this field
    &Sw // <-- It will modify this field
);