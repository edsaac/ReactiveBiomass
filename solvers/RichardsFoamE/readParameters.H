//  Soil properties 
//  Currently it only supports and unique homogeneous material

UnsaturatedSoilParametersField soil
(
    &theta_s,
    &theta_r,
    &K_s,
    &alpha,
    &n,
    &m
);


//  Time stepping properties 
//  Currently it only supports and unique homogeneous material

Info<< "Reading time stepping properties\n" << endl;
dictionary timeProperties
( 
    mesh.solutionDict().subDict("timeStepping")
);

TimeStepper timeControl {timeProperties};

bool adjustTimeStep = runTime.controlDict().lookupOrDefault("adjustTimeStep", true);
