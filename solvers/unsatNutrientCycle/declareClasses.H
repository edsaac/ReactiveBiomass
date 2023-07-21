// Unsaturated flow
class UnsaturatedSoilParametersField
{
public:
    volScalarField *ptr_Sw_s = nullptr;
    volScalarField *ptr_Sw_r = nullptr;
    volScalarField *ptr_alpha = nullptr;
    volScalarField *ptr_n = nullptr;
    volScalarField *ptr_m = nullptr;
    volScalarField *ptr_perm_satu = nullptr;
    volScalarField *ptr_Sw = nullptr;

    // Constructors
    UnsaturatedSoilParametersField(
        volScalarField *Sw_s,
        volScalarField *Sw_r,
        volScalarField *alpha,
        volScalarField *n,
        volScalarField *m,
        volScalarField *perm_satu,
        volScalarField *Sw)
        : ptr_Sw_s(Sw_s),
          ptr_Sw_r(Sw_r),
          ptr_alpha(alpha),
          ptr_n(n),
          ptr_m(m),
          ptr_perm_satu(perm_satu),
          ptr_Sw(Sw)
    {
    }

    // Mualem / Van Genucthen functions
    // All these function calls should be cached :S

    volScalarField innerVanGenuchtenCalculator(const volScalarField &h)
    {
        // innerAlphaPower
        // (-αh)^n
        return Foam::pow(-(*ptr_alpha) * h * Foam::neg(h), (*ptr_n));
    }

    volScalarField vanGenuchtenCalculator(const volScalarField &h)
    {
        // Se = (1 + (-αh)^n)^-m
        return Foam::pow(1.0 + innerVanGenuchtenCalculator(h), -(*ptr_m));
    }

    void waterSaturationCalculator(const volScalarField &h)
    {
        // Sw
        // Sw = Se * (Sw_s - Sw_r) + Sw_r
        *(this->ptr_Sw) = vanGenuchtenCalculator(h) * ((*ptr_Sw_s) - (*ptr_Sw_r)) + (*ptr_Sw_r);
    }

    void mualemCalculator(const volScalarField &h)
    {
        // perm
        // k = √Se (1 - (1 - (Se)^(1/m))^m)²
        // return Foam::sqrt(vanGenuchtenCalculator(h)) * Foam::sqr(1.0 - Foam::pow(1-Foam::pow(vanGenuchtenCalculator(h), 1.0/(*ptr_m)), (*ptr_m)));
        *(this->ptr_perm_satu) =
            Foam::sqrt(vanGenuchtenCalculator(h)) * Foam::sqr(
                                                        1.0 - Foam::pow(
                                                                  1.0 - Foam::pow(
                                                                            vanGenuchtenCalculator(h), 1.0 / (*ptr_m)),
                                                                  (*ptr_m)));
    }

    volScalarField capillary(const volScalarField &h)
    {
        // C(h)
        return ((*ptr_Sw_s) - (*ptr_Sw_r)) * (*ptr_alpha) * (*ptr_m) * (*ptr_n) * Foam::pow(-(*ptr_alpha) * h * Foam::neg(h), *ptr_n - 1) * Foam::pow(1.0 + innerVanGenuchtenCalculator(h), -(*ptr_m) - 1);
    }
};

// Kinetics of microbes

class Kinetics
{
public:
    dimensionedScalar qhat;
    dimensionedScalar yield;
    dimensionedScalar bDie;
    dimensionedScalar Kdonor;
    dimensionedScalar Kaccep;

    Kinetics(IOdictionary &ParameterDict)
    {
        Foam::Info << "Maximum reaction rate qhat\n"
                   << endl;
        qhat.read(ParameterDict.lookup("qhat"));

        Foam::Info << "Yield\n"
                   << endl;
        yield.read(ParameterDict.lookup("yield"));

        Foam::Info << "Half-saturation e-donor Kdonor\n"
                   << endl;
        Kdonor.read(ParameterDict.lookup("Kdonor"));

        Foam::Info << "Half-saturation e-acceptor Kaccep \n"
                   << endl;
        Kaccep.read(ParameterDict.lookup("Kaccep"));

        Foam::Info << "Die-off rate bdie\n"
                   << endl;
        bDie.read(ParameterDict.lookup("bdie"));
    }
};
