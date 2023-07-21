Foam::Info << "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


//  Time stepping properties 
//  Currently it only supports and unique homogeneous material

Foam::Info<< "Reading time stepping properties\n" << endl;
dictionary timeProperties
( 
    mesh.solutionDict().subDict("timeStepping")
);

TimeStepper timeControl {timeProperties};

bool adjustTimeStep = runTime.controlDict().lookupOrDefault("adjustTimeStep", true);
