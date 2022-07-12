/*---------------------------------------------------------------------------*\
                    Reading of scalar input data
\*---------------------------------------------------------------------------*/

// Info<< "Reading transportProperties\n" << endl;
// IOdictionary transportProperties
// (
//     IOobject
//     (
//         "transportProperties",
//         runTime.constant(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::NO_WRITE
//     )
// );

Info<< "Reading soilProperties\n" << endl;
IOdictionary soilProperties
(
    IOobject
    (
        "soilProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

class UnsaturatedSoilParameters{
    public:
        // Foam::dimensionedScalar sat_s ("sat_s", dimless , scalar(0.0));
        // Foam::dimensionedScalar sat_r ("sat_r", dimless , scalar(0.0));
        // Foam::dimensionedScalar K_s ("K_s", dimensionSet(0,1,-1,0,0,0,0) , scalar(0.0));
        // Foam::dimensionedScalar alpha ("alpha", dimensionSet(0,-1,0,0,0,0,0) , scalar(0.0));
        // Foam::dimensionedScalar n ("n", dimless , scalar(0.0));
        // Foam::dimensionedScalar m ("m", dimless , scalar(0.0));
        
        Foam::dimensionedScalar theta_s;
        Foam::dimensionedScalar theta_r;
        Foam::dimensionedScalar K_s;
        Foam::dimensionedScalar alpha;
        Foam::dimensionedScalar n;
        Foam::dimensionedScalar m;

    
    UnsaturatedSoilParameters(IOdictionary& ParameterDict)
    {
        Foam::Info << "Saturation water content [-]\n" << endl;
        theta_s.read(ParameterDict.lookup("theta_s"));
        
        Foam::Info << "Residual water content [-]\n" << endl;
        theta_r.read(ParameterDict.lookup("theta_r"));

        Foam::Info << "Saturated hydraulic conductivity field [m/s]\n" << endl;
        K_s.read(ParameterDict.lookup("K_s"));

        Foam::Info << "Inverse of the capillary length [1/m]\n" << endl;
        alpha.read(ParameterDict.lookup("alpha"));

        Foam::Info << "van Genuchten exponent n [-]\n" << endl;
        n.read(ParameterDict.lookup("n"));

        m = dimensionedScalar("m", dimless, scalar(1.0) - scalar(1.0)/n.value());
        Foam::Info << "van Genuchten exponent m [-]\n" << m << endl;
    }

    // void innerVanGenuchtenCalculator(const volScalarField& h, volScalarField& innerAlphaPower)
    // {
    //     innerAlphaPower = Foam::pow(- alpha * h * Foam::neg(h), n);
    // }

    // void vanGenuchtenCalculator(const volScalarField& h, const volScalarField& innerAlphaPower, volScalarField& theta_e)
    // {
    //     theta_e = Foam::pow( 1.0 + innerAlphaPower, -m);
    // }
    
    // void waterContentCalculator(const volScalarField& theta_e, volScalarField& theta)
    // {
    //     theta = theta_e * (theta_s - theta_r) + theta_r;
    // }
    
    // void mualemCalculator(const volScalarField& theta_e, volScalarField& perm)
    // {
    //     perm = Foam::sqrt(theta_e) * Foam::sqr(1.0 - Foam::pow(1-Foam::pow(theta_e, 1.0/m), m));
    // }

    // void capillaryCapacityCalculator(const volScalarField& innerAlphaPower, const volScalarField& theta_e, const volScalarField& h, volScalarField& capillary)
    // {
    //     capillary = (theta_s - theta_r) * m * n * innerAlphaPower * theta_e / (h * (innerAlphaPower + 1.0));
    // }

    volScalarField innerVanGenuchtenCalculator(const volScalarField& h)
    {
        // innerAlphaPower
        return Foam::pow(-alpha * h * Foam::neg(h), n);
    }

    volScalarField vanGenuchtenCalculator(const volScalarField& h)
    {
        // theta_e
        return  Foam::pow( 1.0 + innerVanGenuchtenCalculator(h), -m);
    }

    volScalarField waterContentCalculator(const volScalarField& h)
    {
        // theta
        return vanGenuchtenCalculator(h) * (theta_s - theta_r) + theta_r;
    }

    volScalarField mualemCalculator(const volScalarField& h)
    {
        // perm
        return Foam::sqrt(vanGenuchtenCalculator(h)) * Foam::sqr(1.0 - Foam::pow(1-Foam::pow(vanGenuchtenCalculator(h), 1.0/m), m));
    }

    volScalarField capillary(const volScalarField& h)
    {
        // C(h)
        return (theta_s - theta_r) * m * n * innerVanGenuchtenCalculator(h) * vanGenuchtenCalculator(h) / (h * (innerVanGenuchtenCalculator(h) + 1.0));
    }

    volScalarField K(const volScalarField& h)
    {
        // Hydraulic conductivity
        return K_s * mualemCalculator(h);
    }
};


//Biochemical parameters
UnsaturatedSoilParameters soil {soilProperties};

/////////////////// Time control parameters //////////////////
Info<< "Reading time stepping properties\n" << endl;
dictionary timeProperties = mesh.solutionDict().subDict("timeStepping");

class TimeStepper{
    public:
     
        scalar convergeThreshold;
        int maxItersToConverge;
        
        int UpDeltaT_nIters;
        scalar UpDeltaT_Factor;
        scalar maxDeltaT;

        scalar DownDeltaT_Factor;
        scalar minDeltaT;
           
    TimeStepper(dictionary& ParameterDict)
    {
        Foam::Info << "Convergence threshold [-]\n" << endl;
        convergeThreshold = ParameterDict.lookupOrDefault<scalar>("convergeThreshold", 1.0E-04);

        Foam::Info << "Give up after iterations [-]\n" << endl;
        maxItersToConverge = ParameterDict.lookupOrDefault<int>("maxItersToConverge", 10);
        
        Foam::Info << "n Iters to increase deltaT [-]\n" << endl;
        UpDeltaT_nIters = ParameterDict.lookupOrDefault<int>("UpDeltaT_nIters", 3);

        Foam::Info << "Factor to increase deltaT [-]\n" << endl;
        UpDeltaT_Factor = ParameterDict.lookupOrDefault<scalar>("UpDeltaT_Factor", 1.3);

        Foam::Info << "Factor to decrease deltaT [-]\n" << endl;
        DownDeltaT_Factor = ParameterDict.lookupOrDefault<scalar>("DownDeltaT_Factor", 0.7);

        Foam::Info << "maxDeltaT [s] \n" << endl;
        maxDeltaT = ParameterDict.lookupOrDefault<scalar>("maxDeltaT", 100);

        Foam::Info << "minDeltaT [s] \n" << endl;
        minDeltaT = ParameterDict.lookupOrDefault<scalar>("minDeltaT", 0.01);   

    }
};

//Biochemical parameters
TimeStepper timeControl {timeProperties};

// TimeStep adjustable?
bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);
