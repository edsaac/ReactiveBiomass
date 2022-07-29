#include "fvCFD.H"

class cloggingModel
{

protected: 

    // Protected data
    word name_;
    dictionary cloggingProperties_;
    const volScalarField& n_;

public:

    TypeName("cloggingModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        cloggingModel,
        dictionary,
        (
            const word& name,
            const dictionary& cloggingProperties,
            const volScalarField& n
        ),
        (name, cloggingProperties,n)
    );

    // Constructors

    cloggingModel
    (
        const word& name,
        const dictionary& cloggingProperties,
        const volScalarField& n
    );

    cloggingModel(const cloggingModel&);

    // Selectors
    static autoPtr<cloggingModel> New
    (
        const word& name,
        const dictionary& cloggingProperties,
        const volScalarField& n
    );

    // Destructor
    virtual ~cloggingModel()
    {}

    // Member Functions 

    // Return the clogging properties dictionary
    const dictionary& cloggingProperties() const
    {
        return cloggingProperties_;
    }

    // Return the calculated permeability
    virtual void updatePerm(const volScalarField& n, volScalarField* perm) = 0;

    // Read cloggingProperties dictionary
    virtual bool read(const dictionary& cloggingProperties) = 0;

    void operator=(const cloggingModel&) = delete;
};

/////////////// cloggingModel.C

namespace Foam
{
    defineTypeNameAndDebug(cloggingModel,0);
    defineRunTimeSelectionTable(cloggingModel,dictionary);
}

// ********************* Constructors ************************

cloggingModel::cloggingModel
(
    const word& name,
    const dictionary& cloggingProperties,
    const volScalarField& n
)
:
    name_(name),
    cloggingProperties_(cloggingProperties),
    n_(n)
{}

// ********************** Member Functions *********************

bool cloggingModel::read(const dictionary& cloggingProperties)
{
    cloggingProperties_ = cloggingProperties;
    return true;
} 

///////////////////////////////////////////////////////////////

class powerLawClog
:
    public cloggingModel
{
    // Private Data
    dictionary powerLawCoeffs_;
    dimensionedScalar nMin_;
    dimensionedScalar nRef_;
    dimensionedScalar permMin_;
    dimensionedScalar permRef_;
    scalar powExponent_;

public: 
    TypeName("powerLawClog");

    // Constructors

    powerLawClog
    (
        const word& name,
        const dictionary& cloggingProperties,
        const volScalarField& n
    );

    // Destructor

    virtual ~powerLawClog()
    {}

    // Member Functions
    virtual void updatePerm(const volScalarField& n, volScalarField* perm);
    virtual bool read(const dictionary& cloggingProperties);
};

/////////////// powerLawClog.C

namespace Foam
{
    namespace cloggingModels
    {
        defineTypeNameAndDebug(powerLawClog,0);
        addToRunTimeSelectionTable
        (
            cloggingModel,
            powerLawClog,
            dictionary
        );
    }
} 

Foam::cloggingModels::powerLawClog::updatePerm(const volScalarField& n, volScalarField* perm) const
{
    perm = ((permRef_ - permMin_) * Foam::pow((n - nMin_)/(nRef_ - nMin_),powExponent_) * pos(n - nMin_)) + permin_;
}


