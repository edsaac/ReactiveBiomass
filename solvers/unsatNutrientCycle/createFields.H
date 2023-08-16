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

Foam::Info << "Creating field: relative permeability by clogging \n" << endl;
volScalarField perm_clog
(
    IOobject
    (
        "perm_clog",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("perm_clog", dimless, 1.0)
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

Foam::Info << "Reading field: Dissolved organic carbon (DOC) \n" << endl;
volScalarField DOC
(
    IOobject
    (
        "DOC",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::Info << "Reading field: Extrapolymeric substances (EPS) \n" << endl;
volScalarField EPS
(
    IOobject
    (
        "EPS",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::Info << "Reading field: Biomass-associated products (BAP) \n" << endl;
volScalarField BAP
(
    IOobject
    (
        "BAP",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::Info << "Reading field: Ammonium (NH4+) \n" << endl;
volScalarField NH4
(
    IOobject
    (
        "NH4",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::Info << "Reading field: Oxygen (O2) \n" << endl;
volScalarField O2
(
    IOobject
    (
        "O2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::Info << "Reading field: Oxygen gas (O2gas) \n" << endl;
volScalarField O2gas
(
    IOobject
    (
        "O2gas",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::Info << "Reading field: Nitrate (NO3) \n" << endl;
volScalarField NO3
(
    IOobject
    (
        "NO3",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::Info << "Reading field: Biomass Aerobic Respirators (XAR) \n" << endl;
volScalarField XAR
(
    IOobject
    (
        "XAR",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::Info << "Reading field: Biomass Aerobic Respirators (XARp) \n" << endl;
volScalarField XARp
(
    IOobject
    (
        "XARp",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::Info << "Reading field: Biomass Nitrifiers (XN) \n" << endl;
volScalarField XN
(
    IOobject
    (
        "XN",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::Info << "Reading field: Biomass Nitrifiers (XNp) \n" << endl;
volScalarField XNp
(
    IOobject
    (
        "XNp",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::Info << "Reading field: Biomass Denitrifiers (XDN) \n" << endl;
volScalarField XDN
(
    IOobject
    (
        "XDN",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::Info << "Reading field: Biomass Denitrifiers (XDNp) \n" << endl;
volScalarField XDNp
(
    IOobject
    (
        "XDNp",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::Info << "Reading field: Inert Biomass (XI) \n" << endl;
volScalarField XI
(
    IOobject
    (
        "XI",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Foam::Info << "Reading field: recalcitrant particulate carbon (POCr) \n" << endl;
volScalarField POCr
(
    IOobject
    (
        "POCr",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Calculate diffusion coefficient from Einstein-Stokes equation
//   >> use this in the Péclet number calculation
//dimensionedScalar D_diffusion
//(
//    (kBoltz*Temp)/(3.0*ppi*mu*dp)
//);
//Foam::Info << "D_diff: ("<< D_diffusion.value() <<")\n" << endl;
////////////////////


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

Foam::Info << "Switch for clogging?? (cloggingSwitch)\n" << endl;
Foam::Switch cloggingSwitch(cloggingProperties.lookup("clog?"));

autoPtr<cloggingModel> clogging
(
    cloggingModel::New(
        "clog",
        cloggingProperties,
        porosity,
        &perm_clog // <-- It will modify this field
        )
);

clogging->read(cloggingProperties);

soil.mualemCalculator(h); // This will modify perm_satu
if (cloggingSwitch) { clogging->calcPerm();} // This will modify perm_clog

// updateHydCond;
// U = - hydraulicCond * (fvc::grad(h) + fvc::grad(z));
// phi = fvc::flux(U);

// Calculate deposited mass and clogLimiter
dimensionedScalar XMAX 
(
    (clogging->nRef() - clogging->nMin()) * rho_X
);

Foam::Info << "Calculate max possible deposited mass: " << XMAX << "\n";

volScalarField totalBiomass (XAR + XN + XDN + XI + EPS);
volScalarField clogLimiter (1.0 - totalBiomass/XMAX);

// Reaction utilization rates
//// Aerobic heterotrophs
Foam::Info << "Calculate initial rates of utilization \n" << endl;

volScalarField rH
(
    IOobject
    (
        "rH",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pXAR.qhat  
        * DOC/(pXAR.Kdonor + DOC)
        * O2/(pXAR.Kaccep + O2)
        * NH4/(K_minN + NH4)
        * clogLimiter
);

Foam::Info << "Calculate initial rates of utilization rN \n" << endl;
//// Aerobic nitrifiers
volScalarField rN
(
    IOobject
    (
        "rN",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pXN.qhat 
        * NH4/(pXN.Kdonor + NH4) 
        * O2/(pXN.Kaccep + O2)
        * clogLimiter
);

//// Anoxic denitrifiers
volScalarField rDN
(
    IOobject
    (
        "rDN",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pXDN.qhat
        * DOC/(pXDN.Kdonor + DOC)
        * NO3/(pXDN.Kaccep + NO3)
        * K_I / (K_I + O2)
        * clogLimiter
);


// Stoichiometric relations
// alpha1  = 1.00 * (1. - 1.42*pXAR.Y - pXAR.kEPS) #gO2/gBOD
// alphaN  = 4.57 * (1. - 0.31*pXN.Y  -  pXN.kEPS) #gO2/gNH4  //Need to remove NH4 used by heterotrophs?
// beta1   = 0.35 * (1. - 1.42*pXDN.Y - pXDN.kEPS) #gNO3/gBOD

Foam::Info << "Ratio of mg O2 / mg VSS-Xa from DOC (alpha_1)\n" << endl;
dimensionedScalar alpha_1(1.0 * (1.0 - 1.42*pXAR.yield - kEPS));

Foam::Info << "Ratio of mg O2 / mg VSS-Xn (alpha_N)\n" << endl;
dimensionedScalar alpha_N(4.57 * (1.0 - 0.31*pXN.yield - kEPS));   

Foam::Info << "Ratio of mg NO3 / mg VSS-Xdn from DOC (beta_1)\n" << endl;
dimensionedScalar beta_1(0.35 * (1.0 - 1.42*pXDN.yield - kEPS));


//////////// Attachment rate calculators
//Attachment kinetics 

Foam::Info << "Attachment of mobile particulates to immobile species(first order rate)\n" << endl;

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
    dimensionedScalar("katt_", dimensionSet(0,0,-1,0,0,0,0) , 0.0)
);

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
    attachmentModel::New(
        "attach",
        attachmentProperties,
        &katt,
        mesh,
        runTime,
        U,
        porosity)
);

attachment->read(attachmentProperties);
Info<< attachment->attachmentProperties() << endl;


Foam::Info << "Detachment of immobile species to mobile particulates\n" << endl;
volScalarField kdet
(
    IOobject
    (
        "kdet",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kdet", dimensionSet(0,0,-1,0,0,0,0) , 0.0)
);

IOdictionary detachmentProperties
(
    IOobject
    (
        "detachmentProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

autoPtr<attachmentModel> detachment
(
    attachmentModel::New(
        "detach",
        detachmentProperties,
        &kdet,
        mesh,
        runTime,
        U,
        porosity)
);

detachment->read(detachmentProperties);
Info<< detachment->attachmentProperties() << endl;
