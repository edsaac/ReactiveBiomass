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

// Densities
Foam::Info << "Reading water density (rho)\n" << endl;
dimensionedScalar rho(transportProperties.lookup("rho"));

Foam::Info << "Reading biomass density (rho_X)\n" << endl;
dimensionedScalar rho_X(transportProperties.lookup("rho_X"));

//Other physical stuff
Foam::Info << "Reading gravity (g)\n" << endl;
dimensionedScalar g(transportProperties.lookup("g"));

Foam::Info << "Reading dynamic viscosity (mu)\n" << endl;
dimensionedScalar mu(transportProperties.lookup("mu"));

// Diffusion and dispersion
Foam::Info << "Molecular diffusion\n" << endl;
dimensionedTensor molDiff(transportProperties.lookup("molDiff"));

Foam::Info << "Dispersion tensor\n" << endl;
dimensionedTensor DispTensor(transportProperties.lookup("DispTensor"));

Foam::Info << "Diffusive growth\n" << endl;
dimensionedScalar diffusiveGrowth(transportProperties.lookup("diffusiveGrowth"));

Foam::Info << "Henry's law for oxygen partition\n" << endl;
dimensionedScalar Hacc(transportProperties.lookup("Hacc"));

Foam::Info << "Oxygen saturation concentration\n" << endl;
dimensionedScalar O2_saturation(transportProperties.lookup("O2_saturation"));


//-----------

Foam::Info << "Reading cloggingProperties\n" << endl;

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

//-----------

Foam::Info << "Reading <microbe>Properties\n" << endl;


IOdictionary aerobicHeterotrophKinetics
(
    IOobject
    (
        "kineticsAerobicHeterotrophs",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOdictionary nitrifierKinetics
(
    IOobject
    (
        "kineticsNitrifiers",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOdictionary denitrifierKinetics
(
    IOobject
    (
        "kineticsDenitrifiers",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//Biochemical parameters
Kinetics pXAR{aerobicHeterotrophKinetics};
Kinetics pXN{nitrifierKinetics};
Kinetics pXDN{denitrifierKinetics};

IOdictionary biochemicalProperties
(
    IOobject
    (
        "biochemicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Oxygen kinetics
Foam::Info << "Inhibition concentration KI (O2)\n" << endl;
dimensionedScalar K_I(biochemicalProperties.lookup("K_I"));

Foam::Info << "Minimum concentration of NH4 as N-source (K_minN)\n" << endl;
dimensionedScalar K_minN(biochemicalProperties.lookup("K_minN"));

// Hydrolisis of particulates into dissolved carbon
Foam::Info << "Labile species hydrolisis rate into DOC (khyd_labil)\n" << endl;
dimensionedScalar khyd_labil(biochemicalProperties.lookup("khyd_labil"));

Foam::Info << "Recalcitrant species hydrolisis rate into DOC (khyd_recal)\n" << endl;
dimensionedScalar khyd_recal(biochemicalProperties.lookup("khyd_recal"));

Foam::Info << "EPS formation coefficient (kEPS)\n" << endl;
dimensionedScalar kEPS(biochemicalProperties.lookup("kEPS"));

// Other stuff?
Foam::Info << "Biodegradable fraction of active biomass (fd)\n" << endl;
dimensionedScalar fd(biochemicalProperties.lookup("fd"));

Foam::Info << "Ratio of N to VS in biomass mg_N/mg_VSS (gamma_N)\n" << endl;
dimensionedScalar gamma_N(biochemicalProperties.lookup("gamma_N"));

Foam::Info << "Minimum viable biomass concentraton (X_min)\n" << endl;
dimensionedScalar X_min(biochemicalProperties.lookup("X_min"));

//----

//  Time stepping properties 
//  Currently it only supports and unique homogeneous material

Foam::Info<< "Reading time stepping properties\n" << endl;
dictionary timeProperties
( 
    mesh.solutionDict().subDict("timeStepping")
);

TimeStepper timeControl {timeProperties};

bool adjustTimeStep = runTime.controlDict().lookupOrDefault("adjustTimeStep", true);
