//  Soil properties 
//  Currently it only supports and unique homogeneous material

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

#include "UnsaturatedSoilParameters.H"
UnsaturatedSoilParameters soil {soilProperties};


//  Time stepping properties 
//  Currently it only supports and unique homogeneous material

Info<< "Reading time stepping properties\n" << endl;
dictionary timeProperties
( 
    mesh.solutionDict().subDict("timeStepping")
);

# include "timeStepper.H"
TimeStepper timeControl {timeProperties};

bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);
