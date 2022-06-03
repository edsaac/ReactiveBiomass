// Field of Happel parameter As
log_As = log( As_Coeffs[0]*sqr(1/n) + As_Coeffs[1]*(1/n) + As_Coeffs[2] );
Info<< "Happel:"<< max(log_As) <<"\n";

// Peclét number field 
log_N_Peclet = log( (mag(U)*ds)/D_diffusion ) ;
Info<< "Peclet:"<< max(log_N_Peclet) <<"\n";

// Field of filtration efficiency
eta = 
// Diffusion mechanism
exp(
    0.875468737
+ (0.333 * log_As)
- (0.081 * log_N_Ratio)
- (0.715 * log_N_Peclet)
+ (0.052 * log_N_vdW) )

+ 

// Interception mechanism
exp(
- 0.597837001
+ log_As
+ (1.550 * log_N_Ratio)
- (0.125 * log_N_Peclet)
+ (0.125 * log_N_vdW) )

+ 

// Gravitational deposition mechanism
exp(
- 0.744440475
- (1.350 * log_N_Ratio)
- (1.110 * log_N_Peclet)
+ (0.053 * log_N_vdW)
+ (1.110 * log_N_gravit) );

// Info<< "eta_0:"<< max(eta) <<"\n" << endl;

// Calculate attachment rate
Lambda = (3.0 * (1.0 - n) * alphaCFT * eta)/(2.0*ds);
katt   = Lambda * mag(U) / n;