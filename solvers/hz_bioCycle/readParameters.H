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

// Diameters
Foam::Info << "Reading bed grain size (ds)\n" << endl;
dimensionedScalar ds(transportProperties.lookup("ds"));

Foam::Info << "Reading clay size (dp)\n" << endl;
dimensionedScalar dp(transportProperties.lookup("dp"));

// Densities
Foam::Info << "Reading water density (rho)\n" << endl;
dimensionedScalar rho(transportProperties.lookup("rho"));

Foam::Info << "Reading clay density (rho_clay)\n" << endl;
dimensionedScalar rho_clay(transportProperties.lookup("rho_clay"));

Foam::Info << "Reading biomass density (rho_X)\n" << endl;
dimensionedScalar rho_X(transportProperties.lookup("rho_X"));

Foam::Info << "Reading sand density (rho_sand)\n" << endl;
dimensionedScalar rho_sand(transportProperties.lookup("rho_sand"));

//Other physical stuff
Foam::Info << "Reading gravity (g)\n" << endl;
dimensionedScalar g(transportProperties.lookup("g"));

Foam::Info << "Reading dynamic viscosity (mu)\n" << endl;
dimensionedScalar mu(transportProperties.lookup("mu"));

//Foam::Info << "Reading Slope\n" << endl;
//dimensionedVector Slope(transportProperties.lookup("Slope"));

//Dispersion + Diffusion
// Foam::Info << "Longitudinal dispersion \n" << endl;
// dimensionedScalar LongDisp(transportProperties.lookup("LongDisp"));

// Foam::Info << "Transversal dispersion \n" << endl;
// dimensionedScalar TransDisp(transportProperties.lookup("TransDisp"));

// Foam::Info << "Molecular diffusion\n" << endl;
// dimensionedScalar molDiff(transportProperties.lookup("molDiff"));

Foam::Info << "Molecular diffusion\n" << endl;
dimensionedTensor molDiff(transportProperties.lookup("molDiff"));

Foam::Info << "Dispersion tensor\n" << endl;
dimensionedTensor DispTensor(transportProperties.lookup("DispTensor"));

//Dynamic-var Filtration
Foam::Info << "Reading attachement efficiency (alpha)\n" << endl;
dimensionedScalar alpha(transportProperties.lookup("alpha"));

//Chemical parameters
Foam::Info << "Reading Hamaker constant (Hamaker)\n" << endl;
dimensionedScalar Hamaker(transportProperties.lookup("Hamaker"));

Foam::Info << "Reading Boltzmann constan (kBoltz)\n" << endl;
dimensionedScalar kBoltz(transportProperties.lookup("kBoltz"));

Foam::Info << "Reading Temperature (Temp)\n" << endl;
dimensionedScalar Temp(transportProperties.lookup("Temp"));

//Static-vars for kinetic filtration

Foam::Info << "Attachment of BAP to EPS (first order rate)\n" << endl;
dimensionedScalar katt_BAP(transportProperties.lookup("katt_BAP"));

Foam::Info << "Attachment of POCr to XI (first order rate)\n" << endl;
dimensionedScalar katt_POCr(transportProperties.lookup("katt_POCr"));

Foam::Info << "Detachment of EPS to BAP (first order rate)\n" << endl;
dimensionedScalar kdet_EPS(transportProperties.lookup("kdet_EPS"));

Foam::Info << "Detachment of XI to POCr (first order rate)\n" << endl;
dimensionedScalar kdet_XI(transportProperties.lookup("kdet_XI"));

//Max growth due space
//Foam::Info << "Max biomass concentration\n" << endl;
//dimensionedScalar X_max(transportProperties.lookup("X_max"));

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

Foam::Info << "Switch for clogging?? (cloggingSwitch)\n" << endl;
Foam::Switch cloggingSwitch(cloggingProperties.lookup("clog?"));

Foam::Info << "Initial porosity (n_0)\n" << endl;
dimensionedScalar n_0 (cloggingProperties.lookup("n_0"));

Foam::Info << "Initial relative permeability (perm_0)\n" << endl;
dimensionedScalar perm_0 (cloggingProperties.lookup("perm_0"));

Foam::Info << "Critical porosity (n_c)\n" << endl;
dimensionedScalar n_c(cloggingProperties.lookup("n_c"));

Foam::Info << "Minimum relative permeability (perm_c)\n" << endl;
dimensionedScalar perm_c(cloggingProperties.lookup("perm_c"));

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


