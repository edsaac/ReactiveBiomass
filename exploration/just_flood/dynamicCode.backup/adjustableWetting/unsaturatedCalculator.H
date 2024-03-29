class UnsaturatedSoilTop
{
public:
    scalar Sw_s;
    scalar Sw_r;
    scalar alpha;
    scalar n;
    scalar m;
    scalar Ks;
    scalar qtarget;

    // Constructors
    UnsaturatedSoilTop(
        scalar Sw_s,
        scalar Sw_r,
        scalar alpha,
        scalar n,
        scalar Ks,
        scalar qtarget)
        : Sw_s(Sw_s),
          Sw_r(Sw_r),
          alpha(alpha),
          n(n),
          m(1.0 - 1.0/n),
          Ks(Ks),
          qtarget(qtarget)
    {
    }

    // Mualem / Van Genucthen functions

    scalar innerVanGenuchtenCalculator(const scalar h)
    {
        // innerAlphaPower
        // (-αh)^n
        return Foam::pow(-(alpha) * h * Foam::neg(h), (n));
    }

    scalar vanGenuchtenCalculator(const scalar h)
    {
        // Se = (1 + (-αh)^n)^-m
        return Foam::pow(1.0 + innerVanGenuchtenCalculator(h), -m);
    }

    scalar waterSaturationCalculator(const scalar h)
    {
        // Sw = Se * (Sw_s - Sw_r) + Sw_r
        return (vanGenuchtenCalculator(h) * (Sw_s - Sw_r)) + Sw_r;
    }

    scalar mualemCalculator(const scalar h)
    {
        // perm
        // k = √Se (1 - (1 - (Se)^(1/m))^m)²
        // return Foam::sqrt(vanGenuchtenCalculator(h)) * Foam::sqr(1.0 - Foam::pow(1-Foam::pow(vanGenuchtenCalculator(h), 1.0/(*ptr_m)), (*ptr_m)));
        return
            Foam::sqrt(vanGenuchtenCalculator(h)) 
            * Foam::sqr(
                1.0 - Foam::pow(1.0 - Foam::pow(vanGenuchtenCalculator(h), 1.0 / m), m)
                );
    }

    scalar capillary(const scalar h)
    {
        // C(h)
        return (Sw_s - Sw_r) * alpha * m * n * Foam::pow(-alpha * h * Foam::neg(h), n - 1) * Foam::pow(1.0 + innerVanGenuchtenCalculator(h), -m - 1);
    }

    scalar target_flow(const scalar h)
    {
        return qtarget - (Ks * mualemCalculator(h));
    }
    scalar bisection_h_from_target(scalar lower_end, scalar upper_end)
    {
        // q = -K(h) grad(h)
        // Assume grad(h) = -1.0

        /* Bisection method, adapted from: 
           https://pythonnumericalmethods.berkeley.edu/notebooks/chapter19.03-Bisection-Method.html
        */
        const scalar tol = 1e-15;
        scalar midpoint = 0;

        // check if lower_end and upper_end bound lower_end root
        if ( Foam::sign(target_flow(lower_end)) == Foam::sign(target_flow(upper_end)) ){
            return 0.0;
        }
        
        midpoint = (lower_end + upper_end)/2;    // get midpoint

        // Info << "  current midpoint: " << midpoint << endl;
        // Info << "target_flow(" << midpoint << "): " << target_flow(midpoint) << endl;
    
        if ( Foam::mag(target_flow(midpoint)) < tol )
        {
            return midpoint;  // stopping condition, report m as root
        } 
        else if ( Foam::sign(target_flow(lower_end)) == Foam::sign(target_flow(midpoint)) ) 
        {
            // case where midpoint is an improvement on lower_end. 
            // Make recursive call with lower_end = midpoint
            return bisection_h_from_target(midpoint, upper_end);
        } 
        else if ( Foam::sign(target_flow(upper_end)) == Foam::sign(target_flow(midpoint)) )
        {
            // # case where midpoint is an improvement on upper_end. 
            // # Make recursive call with upper_end = midpoint
            return bisection_h_from_target(lower_end, midpoint);
        }

        return -1000;
    }
};