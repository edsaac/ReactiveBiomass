// Read the transport properties dictionary
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

volScalarField hydraulicCond
(
    IOobject
    (
        "hydraulicCond",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("hydraulicCond0",dimVelocity,1.0)
);

hydraulicCond = perm * rho * g / mu;

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

Foam::Info << "Switch for clogging?? (cloggingSwitch)\n" << endl;
Foam::Switch cloggingSwitch(cloggingProperties.lookup("clog?"));

autoPtr<cloggingModel> clogging
(
    cloggingModel::New("clog",cloggingProperties,n,&perm)
);

clogging->read(cloggingProperties);

// Calculate deposited mass and clogLimiter
dimensionedScalar XMAX ((clogging->nRef() - clogging->nMin()) * rho_X);
Foam::Info << "Calculate max possible deposited mass: " << XMAX << "\n";

volScalarField totalBiomass (XAR + XN + XDN + XI + EPS);
volScalarField clogLimiter (1.0 - totalBiomass/XMAX);

// Reaction utilization rates
//// Aerobic heterotrophs
Foam::Info << "Calculate initial rates of utilization \n" << endl;
Foam::Info << pXAR.Kdonor << endl;

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
    pXAR.qhat * DOC/(pXAR.Kdonor + DOC) * O2/(pXAR.Kaccep + O2) * NH4/(K_minN + NH4) * clogLimiter
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
dimensionedScalar alpha_1(1.0 - 1.42*pXAR.yield - kEPS);

Foam::Info << "Ratio of mg O2 / mg VSS-Xn (alpha_N)\n" << endl;
dimensionedScalar alpha_N(4.57 * (1.0 - 0.31*pXN.yield - kEPS));   

Foam::Info << "Ratio of mg NO3 / mg VSS-Xdn from DOC (beta_1)\n" << endl;
dimensionedScalar beta_1(0.35 * (1.0 - 1.42*pXDN.yield - kEPS));


//////////// Attachment rate calculators
//Attachment kinetics 

Foam::Info << "Attachment of BAP to EPS (first order rate)\n" << endl;

volScalarField katt_BAP
(
    IOobject
    (
        "katt_BAP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("katt_BAP", dimensionSet(0,0,-1,0,0,0,0) , 0.0)
);

IOdictionary attachmentProperties_BAP
(
    IOobject
    (
        "attachmentProperties_BAP",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

autoPtr<attachmentModel> attachment_BAP
(
    attachmentModel::New(
        "attachBAP",
        attachmentProperties_BAP,
        &katt_BAP,
        mesh,
        runTime,
        U,
        n)
);

attachment_BAP->read(attachmentProperties_BAP);
Info<< attachment_BAP->attachmentProperties() << endl;

/////////

Foam::Info << "Attachment of POCr to XI (first order rate)\n" << endl;

volScalarField katt_POCr
(
    IOobject
    (
        "katt_POCr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("katt_POCr", dimensionSet(0,0,-1,0,0,0,0) , 0.0)
);

IOdictionary attachmentProperties_POCr
(
    IOobject
    (
        "attachmentProperties_POCr",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

autoPtr<attachmentModel> attachment_POCr
(
    attachmentModel::New(
        "attachPOCr",
        attachmentProperties_POCr,
        &katt_POCr,
        mesh,
        runTime,
        U,
        n)
);

attachment_POCr->read(attachmentProperties_POCr);
Info<< attachment_POCr->attachmentProperties() << endl;

Foam::Info << "Detachment of EPS to BAP (first order rate)\n" << endl;
volScalarField kdet_EPS
(
    IOobject
    (
        "kdet_EPS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kdet_EPS", dimensionSet(0,0,-1,0,0,0,0) , 0.0)
);

IOdictionary detachmentProperties_EPS
(
    IOobject
    (
        "detachmentProperties_EPS",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

autoPtr<attachmentModel> detachment_EPS
(
    attachmentModel::New(
        "detachEPS",
        detachmentProperties_EPS,
        &kdet_EPS,
        mesh,
        runTime,
        U,
        n)
);

detachment_EPS->read(detachmentProperties_EPS);
Info<< detachment_EPS->attachmentProperties() << endl;

Foam::Info << "Detachment of XI to POCr (first order rate)\n" << endl;
volScalarField kdet_XI
(
    IOobject
    (
        "kdet_XI",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("kdet_XI", dimensionSet(0,0,-1,0,0,0,0) , 0.0)
);

IOdictionary detachmentProperties_XI
(
    IOobject
    (
        "detachmentProperties_XI",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

autoPtr<attachmentModel> detachment_XI
(
    attachmentModel::New(
        "detachXI",
        detachmentProperties_XI,
        &kdet_XI,
        mesh,
        runTime,
        U,
        n)
);

detachment_XI->read(detachmentProperties_XI);
Info<< detachment_XI->attachmentProperties() << endl;