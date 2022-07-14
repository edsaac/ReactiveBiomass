class UnsaturatedSoilParameters
{
    public:   
        Foam::dimensionedScalar theta_s;
        Foam::dimensionedScalar theta_r;
        Foam::dimensionedScalar K_s;
        Foam::dimensionedScalar alpha;
        Foam::dimensionedScalar n;
        Foam::dimensionedScalar m;

    // Constructors
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

    // Mualem / Van Genucthen functions
    // All these function calls should be cached :S
    
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
        return - (theta_s - theta_r) * m * n * innerVanGenuchtenCalculator(h) * vanGenuchtenCalculator(h) / (h * (innerVanGenuchtenCalculator(h) + 1.0));
    }

    volScalarField K(const volScalarField& h)
    {
        // Hydraulic conductivity
        return K_s * mualemCalculator(h);
    }
};