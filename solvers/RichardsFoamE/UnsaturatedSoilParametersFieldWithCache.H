class UnsaturatedSoilParametersField
{
    public:
        volScalarField* ptr_theta_s = nullptr;
        volScalarField* ptr_theta_r = nullptr;
        volScalarField* ptr_K_s = nullptr;
        volScalarField* ptr_alpha = nullptr;
        volScalarField* ptr_n = nullptr;
        volScalarField* ptr_m = nullptr;

        volScalarField* ptr_innerVanGenuchten = nullptr;
        volScalarField* ptr_theta_e = nullptr;
        volScalarField* ptr_theta = nullptr;
        volScalarField* ptr_permRel = nullptr;
        volScalarField* ptr_capillary = nullptr;
        volScalarField* ptr_hydrConduct = nullptr;

    // Constructors
    UnsaturatedSoilParametersField
    (
        volScalarField* theta_s,
        volScalarField* theta_r,
        volScalarField* K_s,
        volScalarField* alpha,
        volScalarField* n,
        volScalarField* m,

        volScalarField* innerVanGenuchten,
        volScalarField* theta_e,
        volScalarField* theta,
        volScalarField* permRel,
        volScalarField* capillary,
        volScalarField* hydrConduct,

        const volScalarField& h
    )
        :
        ptr_theta_s(theta_s),
        ptr_theta_r(theta_r),
        ptr_K_s(K_s),
        ptr_alpha(alpha),
        ptr_n(n),
        ptr_m(m),

        ptr_innerVanGenuchten(innerVanGenuchten),
        ptr_theta_e(theta_e),
        ptr_theta(theta),
        ptr_permRel(permRel),
        ptr_capillary(capillary),
        ptr_hydrConduct(hydrConduct)
    {
        update(h);
    }

    // Mualem / Van Genucthen functions
    // All these function calls should be cached :S
    
    void innerVanGenuchtenCalculator(const volScalarField& h)
    {
        // innerAlphaPower
        *ptr_innerVanGenuchten = Foam::pow(-(*ptr_alpha) * h * Foam::neg(h), (*ptr_n));
    }

    void vanGenuchtenCalculator()
    {
        // theta_e
        *ptr_theta_e = Foam::pow( 1.0 + *ptr_innerVanGenuchten, -(*ptr_m));
    }

    void waterContentCalculator()
    {
        // theta
        *ptr_theta = vanGenuchtenCalculator(h) * ((*ptr_theta_s) - (*ptr_theta_r)) + (*ptr_theta_r);
    }

    void mualemCalculator()
    {
        // perm
        *ptr_permRel = Foam::sqrt((*ptr_theta_e)) * Foam::sqr(1.0 - Foam::pow(1-Foam::pow(*ptr_theta, 1.0/(*ptr_m)), (*ptr_m)));
    }

    void capillary(const volScalarField& h)
    {
        // C(h)
        *ptr_capillary = - ((*ptr_theta_s) - (*ptr_theta_r)) 
            * (*ptr_m) * (*ptr_n) * (*ptr_innerVanGenuchten)
            * (*ptr_theta_e)
            / (h * ( (*ptr_innerVanGenuchten) + 1.0));
    }

    void Kcalculator()
    {
        // Hydraulic conductivity
        return (*ptr_K_s) * (*ptr_permRel);
    }

    void update(const volScalarField& h)
    {
        innerVanGenuchtenCalculator(h);
        vanGenuchtenCalculator();
        waterContentCalculator();
        mualemCalculator();
        capillary(h);
        Kcalculator();
    }
};