class UnsaturatedSoilParametersField
{
    public:   
        volScalarField* ptr_theta_s = nullptr;
        volScalarField* ptr_theta_r = nullptr;
        volScalarField* ptr_K_s = nullptr;
        volScalarField* ptr_alpha = nullptr;
        volScalarField* ptr_n = nullptr;
        volScalarField* ptr_m = nullptr;

    // Constructors
    UnsaturatedSoilParametersField
    (
        volScalarField* theta_s,
        volScalarField* theta_r,
        volScalarField* K_s,
        volScalarField* alpha,
        volScalarField* n,
        volScalarField* m
    )
        :
        ptr_theta_s(theta_s),
        ptr_theta_r(theta_r),
        ptr_K_s(K_s),
        ptr_alpha(alpha),
        ptr_n(n),
        ptr_m(m)
    {   
    }

    // Mualem / Van Genucthen functions
    // All these function calls should be cached :S
    
    volScalarField innerVanGenuchtenCalculator(const volScalarField& h)
    {
        // innerAlphaPower
        return Foam::pow(-(*ptr_alpha) * h * Foam::neg(h), (*ptr_n));
    }

    volScalarField vanGenuchtenCalculator(const volScalarField& h)
    {
        // theta_e
        return  Foam::pow( 1.0 + innerVanGenuchtenCalculator(h), -(*ptr_m));
    }

    volScalarField waterContentCalculator(const volScalarField& h)
    {
        // theta
        return vanGenuchtenCalculator(h) * ((*ptr_theta_s) - (*ptr_theta_r)) + (*ptr_theta_r);
    }

    volScalarField mualemCalculator(const volScalarField& h)
    {
        // perm
        return Foam::sqrt(vanGenuchtenCalculator(h)) * Foam::sqr(1.0 - Foam::pow(1-Foam::pow(vanGenuchtenCalculator(h), 1.0/(*ptr_m)), (*ptr_m)));
    }

    volScalarField capillary(const volScalarField& h)
    {
        // C(h)
        return ((*ptr_theta_s) - (*ptr_theta_r)) 
            * (*ptr_alpha) * (*ptr_m) * (*ptr_n)
            * Foam::pow(-(*ptr_alpha) * h * Foam::neg(h), *ptr_n - 1)
            * Foam::pow(1.0 + innerVanGenuchtenCalculator(h), -(*ptr_m)-1)  ;
    }

    volScalarField K(const volScalarField& h)
    {
        // Hydraulic conductivity
        return (*ptr_K_s) * mualemCalculator(h);
    }
};