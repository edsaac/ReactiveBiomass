class TimeStepper
{
    public:
     
        scalar convergeThreshold;
        int maxItersToConverge;
        
        int UpDeltaT_nIters;
        scalar UpDeltaT_Factor;
        scalar maxDeltaT;

        scalar DownDeltaT_Factor;
        scalar minDeltaT;

        scalar hguess_coeff;
    
    //  Constructor
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

        Foam::Info << "hguess_coeff [-] \n" << endl;
        hguess_coeff = ParameterDict.lookupOrDefault<scalar>("hguess_coeff", 0.30);

        
    }
};