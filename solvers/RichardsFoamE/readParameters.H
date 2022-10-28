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

// Initialize all derived fields
// Calculate Darcy flow velocity
hydrConduct = soil.K(h);
U = - hydrConduct * (fvc::grad(h) + fvc::grad(z));
phi = fvc::flux(U);

////////////////// gradient Fields for dK/dz ////////////////
volVectorField grad_k( fvc::grad(soil.K(h)) );
volScalarField grad_kz( grad_k.component(vector::Z) );

// Calculate derived unsaturated stuff
theta = soil.waterContentCalculator(h);

//  Time stepping properties 
//  Currently it only supports and unique homogeneous material

Info<< "Reading time stepping properties\n" << endl;
dictionary timeProperties
( 
    mesh.solutionDict().subDict("timeStepping")
);

TimeStepper timeControl {timeProperties};

bool adjustTimeStep = runTime.controlDict().lookupOrDefault("adjustTimeStep", true);


