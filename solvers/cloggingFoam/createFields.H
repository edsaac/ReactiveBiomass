const dimensionedScalar GRAVITY ("g",dimAcceleration,9.81);
const dimensionedScalar WATERDENSITY ("rhoWater",dimDensity,1000.0);
const dimensionedScalar WATERDYNAMICVISCOSITY ("muWater",dimDynamicViscosity,0.001);

Foam::Info << "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Diffusion/Dispersion
Foam::Info << "Molecular diffusion\n" << endl;
dimensionedTensor molDiff(transportProperties.lookup("molDiff"));

Foam::Info << "Dispersion tensor\n" << endl;
dimensionedTensor DispTensor(transportProperties.lookup("DispTensor"));

Foam::Info << "Detachment rate (first order rate)\n" << endl;
dimensionedScalar kdet(transportProperties.lookup("kdet"));

Foam::Info << "Reading clog density (rho_clog)\n" << endl;
dimensionedScalar rho_clog(transportProperties.lookup("rho_clog"));

///////////////////
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

volScalarField katt
(
    IOobject
    (
        "katt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("katt0", dimensionSet(0,0,-1,0,0,0,0) , scalar(0.0))
);

Info<< "katt field initialized" <<"\n";

////////////////////
IOdictionary cloggingProperties
(
    IOobject
    (
        "cloggingProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Foam::Info << "Switch for clogging?? (cloggingSwitch)\n" << endl;
Foam::Switch cloggingSwitch(cloggingProperties.lookup("clog?"));

autoPtr<cloggingModel> clogging
(
    cloggingModel::New("clog",cloggingProperties,n,&perm)
);

clogging->read(cloggingProperties);

IOdictionary attachmentProperties
(
    IOobject
    (
        "attachmentProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

autoPtr<attachmentModel> attachment
(
    attachmentModel::New("attach",attachmentProperties,&katt,mesh,runTime,U,n)
);

attachment->read(attachmentProperties);

// Calculate deposited mass and clogLimiter
dimensionedScalar XMAX ((clogging->nRef() - clogging->nMin()) * rho_clog);
Foam::Info << "Calculate max possible deposited mass: " << XMAX.value() << "\n";

volScalarField clogLimiter (1.0 - depositedClay/XMAX);
//clogLimiter *= Foam::pos(clogLimiter);

volScalarField hydraulicCond (perm * WATERDENSITY * GRAVITY / WATERDYNAMICVISCOSITY);

#include "createPhi.H"
#include "createFvOptions.H"


