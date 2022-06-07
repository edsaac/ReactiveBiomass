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

Foam::Info << "Reading sand density (rho_sand)\n" << endl;
dimensionedScalar rho_sand(transportProperties.lookup("rho_sand"));

//Other physical stuff
Foam::Info << "Reading gravity (g)\n" << endl;
dimensionedScalar g(transportProperties.lookup("g"));

Foam::Info << "Reading dynamic viscosity (mu)\n" << endl;
dimensionedScalar mu(transportProperties.lookup("mu"));

//Foam::Info << "Reading Slope\n" << endl;
//dimensionedVector Slope(transportProperties.lookup("Slope"));

Foam::Info << "Molecular diffusion\n" << endl;
dimensionedTensor molDiff(transportProperties.lookup("molDiff"));

Foam::Info << "Dispersion tensor\n" << endl;
dimensionedTensor DispTensor(transportProperties.lookup("DispTensor"));

//Dynamic-var Filtration
Foam::Info << "Reading attachment efficiency (alphaCFT)\n" << endl;
dimensionedScalar alphaCFT(transportProperties.lookup("alphaCFT"));

//Chemical parameters
Foam::Info << "Reading Hamaker constant (Hamaker)\n" << endl;
dimensionedScalar Hamaker(transportProperties.lookup("Hamaker"));

Foam::Info << "Reading Boltzmann constan (kBoltz)\n" << endl;
dimensionedScalar kBoltz(transportProperties.lookup("kBoltz"));

Foam::Info << "Reading Temperature (Temp)\n" << endl;
dimensionedScalar Temp(transportProperties.lookup("Temp"));

//Static-vars for kinetic filtration

// Foam::Info << "Attachment rate (first order rate)\n" << endl;
// dimensionedScalar katt(transportProperties.lookup("katt"));

// Foam::Info << "Attachment rate (first order rate)\n" << endl;
// dimensionedScalar alphaCFT(transportProperties.lookup("alphaCFT"));

Foam::Info << "Detachment rate (first order rate)\n" << endl;
dimensionedScalar kdet(transportProperties.lookup("kdet"));

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

