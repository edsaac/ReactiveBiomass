/*---------------------------------------------------------------------------*\
                    Reading of scalar input data
\*---------------------------------------------------------------------------*/

Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Reading unsaturatedSoilParametersDict\n" << endl;
IOdictionary unsaturatedSoilParametersDict
(
    IOobject
    (
        "unsaturatedSoilParametersDict",
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
        Foam::dimensionedScalar sat_s;
        Foam::dimensionedScalar sat_r;
        Foam::dimensionedScalar K_s;
        Foam::dimensionedScalar alpha;
        Foam::dimensionedScalar n;
        Foam::dimensionedScalar m;

    UnsaturatedSoilParameters(IOdictionary& ParameterDict)
    {
        Foam::Info << "Saturation water content [-]\n" << endl;
        sat_s.read(ParameterDict.lookup("sat_s"));
        
        Foam::Info << "Residual water content [-]\n" << endl;
        sat_r.read(ParameterDict.lookup("sat_r"));

        Foam::Info << "Saturated hydraulic conductivity field [m/s]\n" << endl;
        K_s.read(ParameterDict.lookup("K_s"));

        Foam::Info << "Inverse of the capillary length [1/m]\n" << endl;
        alpha.read(ParameterDict.lookup("alpha"));

        Foam::Info << "van Genuchten exponent n [-]\n" << endl;
        n.read(ParameterDict.lookup("n"));

        m = dimensionedScalar("m", dimless, scalar(1.0) - scalar(1.0)/n.value());
        Foam::Info << "van Genuchten exponent m [-]\n" << m << endl;
    }
};

//Biochemical parameters
UnsaturatedSoilParameters soil {unsaturatedSoilParametersDict};
