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

Foam::Info << "Detachment of EPS to BAP (first order rate)\n" << endl;
dimensionedScalar kdet_EPS(transportProperties.lookup("kdet_EPS"));

Foam::Info << "Detachment of XI to POCr (first order rate)\n" << endl;
dimensionedScalar kdet_XI(transportProperties.lookup("kdet_XI"));


////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////

class Kinetics{
    public:
        dimensionedScalar qhat;
        dimensionedScalar yield;
        dimensionedScalar bDie;
        dimensionedScalar Kdonor;
        dimensionedScalar Kaccep;

    Kinetics(IOdictionary& ParameterDict){
        Foam::Info << "Maximum reaction rate qhat\n" << endl;
        qhat.read(ParameterDict.lookup("qhat"));

        Foam::Info << "Yield\n" << endl;
        yield.read(ParameterDict.lookup("yield"));

        Foam::Info << "Half-saturation e-donor Kdonor\n" << endl;
        Kdonor.read(ParameterDict.lookup("Kdonor"));

        Foam::Info << "Half-saturation e-acceptor Kaccep \n" << endl;
        Kaccep.read(ParameterDict.lookup("Kaccep"));

        Foam::Info << "Die-off rate bdie\n" << endl;
        bDie.read(ParameterDict.lookup("bdie"));
    }
};

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
Foam::Info << "BAP hydrolisis rate into DOC (khyd_BAP)\n" << endl;
dimensionedScalar khyd_BAP(biochemicalProperties.lookup("khyd_BAP"));

Foam::Info << "POCr hydrolisis rate into DOC (khyd_POCr)\n" << endl;
dimensionedScalar khyd_POCr(biochemicalProperties.lookup("khyd_POCr"));

Foam::Info << "EPS formation coefficient (kEPS)\n" << endl;
dimensionedScalar kEPS(biochemicalProperties.lookup("kEPS"));

// Other stuff?
Foam::Info << "Biodegradable fraction of active biomass (fd)\n" << endl;
dimensionedScalar fd(biochemicalProperties.lookup("fd"));

Foam::Info << "Ratio of N to VS in biomass mg_N/mg_VSS (gamma_N)\n" << endl;
dimensionedScalar gamma_N(biochemicalProperties.lookup("gamma_N"));

Foam::Info << "Minimum viable biomass concentraton (X_min)\n" << endl;
dimensionedScalar X_min(biochemicalProperties.lookup("X_min"));


