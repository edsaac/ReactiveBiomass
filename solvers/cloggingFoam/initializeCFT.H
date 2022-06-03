// Happel model for filtration
//   Check HappelPolyFit.ipynb on how this replaces the
//   original Happel expression. 
const double As_Coeffs[3] = { 8.99992117, -7.49203318, 0.4119361 };

// Calculate diffusion coefficient from Einstein-Stokes equation
//   >> use this in the Péclet number calculation
dimensionedScalar D_diffusion ( (kBoltz*Temp)/(3.0*PI*mu*dp) );
Foam::Info << "D_diff: ("<< D_diffusion.value() <<")\n" << endl;

// Static dimensionless numbers for filtration model
dimensionedScalar log_N_Ratio ( log(dp/ds) );
Info<< "Aspect Ratio:"<< log_N_Ratio.value() <<"\n";

dimensionedScalar log_N_vdW ( log(Hamaker/(kBoltz*Temp)) );
Info<< "van der Waals:"<< log_N_vdW.value() <<"\n";

dimensionedScalar log_N_gravit ( log(PI * pow(dp,4) * (rho_clay - rho) * g / ( 12.0 * kBoltz * Temp )) );
Info<< "Gravitational:"<< log_N_gravit.value() <<"\n";

volScalarField log_As
(
    IOobject
    (
        "log_As",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("log_As", dimensionSet(0,0,0,0,0,0,0) , scalar(1.0))
);
Info<< "Happel field initialized" <<"\n";

volScalarField log_N_Peclet
(
    IOobject
    (
        "log_N_Peclet",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("log_N_Peclet", dimensionSet(0,0,0,0,0,0,0) , scalar(1.0))
);
Info<< "Peclet field initialized" <<"\n";

volScalarField eta
(
    IOobject
    (
        "log_N_Peclet",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("log_N_Peclet", dimensionSet(0,0,0,0,0,0,0) , scalar(1.0))
);
Info<< "eta0 field initialized" <<"\n";

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
    mesh,
    dimensionedScalar("Lambda0", dimensionSet(0,-1,0,0,0,0,0) , scalar(0.0))
);
Info<< "Lambda field initialized" <<"\n";

volScalarField katt
(
    IOobject
    (
        "katt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("katt0", dimensionSet(0,0,-1,0,0,0,0) , scalar(0.0))
);
Info<< "katt field initialized" <<"\n";