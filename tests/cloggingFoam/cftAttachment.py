import numpy as np

''' GLOBAL CONSTANTS '''
PI = 3.141592
BOLTZMANN = 1.380649E-23 #N·m/K
g = 9.81 #m/s²


def dispCoef(D_m,alpha_L,U,n=1.0): 
    '''
    
    Returns
    ----------
    float 
        Dispersion coefficient $D_T$ [m²/s]
        .. math:: D_T = D_m + \\alpha_LU^n

    Parameters
    ----------
    D_m : float
        Molecular difussion coefficient [m²/s]
    alpha_L : float
        Longitudinal dispersion coefficient [m]
    U : float
        Pore-water intersticial velocity [m/s]
    n : float
        empirical fitting exponent [-]

    Notes
    ----------
    It neglects transversal dispersion

    '''
    return D_m + alpha_L*(U**n)


def poreVel(q,theta): 
    '''
    Returns
    ----------
    U : float
        Returns the intersticial flow velocity $U$, aka the pore-water velocity from the darcy velocity
        .. math:: \vec{U} = \\dfrac{\vec{q}}{\theta}
    
    Parameters
    ----------
    q : float
        Darcy velocity [m/s]
    theta : float (0.,1.)
        Porosity [-]
    '''
    return q/theta

def molecularDiff(visco,dp,T):
    '''
    Returns
    ----------
    Dm : float
        Molecular diffusion coefficient $D_m$ calculated from the Stokes-Einstein equation:
        .. math:: D_m = \\dfrac{k_BT}{3\\pi\\eta d_p}
    
    Parameters
    ----------
    visco: float
        Fluid dynamic viscosity :math:`\\eta` [N·s/m²]
    dp: float
        Particle diameter [m]        
    T : float
        Temperature [K]
    
    Notes
    ----------
    - $k_B$ : Boltzmann constant
    ''' 
    return (BOLTZMANN*T)/(3*PI*visco*dp)


def collectorEff(etaD,etaI,etaG):
    '''
    Returns
    ----------
    eta0 : float
        Return the single collector efficiency: :math:`\\eta_0`
        .. math:: \\eta_0 = \\eta_{\rm D} + \\eta_{\rm I} + \\eta_{\rm G}
    
    Parameters
    ----------
    etaD : float
        Collector efficiency due diffusion [-]
    etaI : float
        Collector efficiency due direct interception [-]
    etaG : float
        Collector efficiency due gravitational deposition [-]
    '''
    return etaD + etaI + etaG


def collectorEfficiency_Diffusion(A_s,N_R,N_Pe,N_vdW):
    '''
    Returns
    ----------
    etaD : float
        Collector efficiency due diffusion $eta_D$ following the approximation by Tufenkji & Elimelech (2004)[1]
        .. math:: \\eta_{\rm D} =& 2.4 A_s^{1/3}N_{\rm R}^{-0.081}N_{\rm Pe}^{-0.715}N_{\rm vdW}^{0.052}
    
    Parameters
    ----------
    A_s : float
        Happel parameter for a collection of spheres [-]
    N_R : float
        Size ratio [-]
    N_Pe : float
        Péclet number [-]
    N_vdW : float    
        van der Waals nondimensional number [-]

    References: 
    ----------
    .. [1] https://pubs.acs.org/doi/10.1021/es034049r
    '''
    return 2.40 * (A_s**(1./3.)) * (N_R**-0.081) * (N_Pe**-0.715) * (N_vdW**0.052)

def collectorEfficiency_Interception(A_s,N_R,N_Pe,N_vdW):
    '''
    Returns
    ----------
    etaI : float
        Collector efficiency due direct interception $eta_I$ following the approximation by Tufenkji & Elimelech (2004)[1]
        .. math:: \\eta_{\rm I} =& 0.55 A_sN_{\rm R}^{1.55}N_{\rm Pe}^{-0.125}N_{\rm vdW}^{0.125}

    Parameters
    ----------
    A_s : float
        Happel parameter for a collection of spheres [-]
    N_R : float
        Size ratio [-]
    N_Pe : float
        Péclet number [-]
    N_vdW : float    
        van der Waals nondimensional number [-]
    
    References: 
    ----------
    .. [1] https://pubs.acs.org/doi/10.1021/es034049r
    '''
    return 0.55 * A_s * (N_R**1.55) * (N_Pe**-0.125) * (N_vdW**0.125)

def collectorEfficiency_GDeposition(N_gr,N_R,N_Pe,N_vdW):
    '''
    Returns
    ----------
    etaG : float
        Collector efficiency due gravitational deposition $eta_G$ following the approximation by Tufenkji & Elimelech (2004) [1]
        .. math:: \\eta_{\rm G} =& 0.475 N_{\rm gr}^{1.11} N_{\rm R}^{-1.35}N_{\rm Pe}^{-1.11}N_{\rm vdW}^{0.053}

    Parameters
    ----------
    N_gr : float
        Gravitational number [-]
    N_R : float
        Size ratio [-]
    N_Pe : float
        Péclet number [-]
    N_vdW : float    
        van der Waals nondimensional number [-]
    
    References: 
    ----------
    .. [1] https://pubs.acs.org/doi/10.1021/es034049r
    '''
    return 0.475 * (N_gr**1.11)   * (N_R**-1.35)  * (N_Pe**-1.11)  * (N_vdW**0.053)


def happelParameter(theta):
    '''
    Returns 
    ----------
    As : float
        Happel parameter for packed spheres
        .. math:: A_s = \\dfrac{2(1-s^{5/3})}{2-3s^{1/3}+3s^{5/3}-2s^2} &\\quad& s = 1-\theta
    
    Parameters
    ----------
    theta : float (0.,1.)
        Porosity [-]
    '''
    s = 1-theta
    s53 = s**(5./3.)
    s13 = s**(1./3.)
    s21 = s**2
    return (2*(1-s53))/(2 - (3*s13) + (3*s53) - (2*s21))


def noDim_SizeRatio(dp,dc):
    '''
    Returns 
    ----------
    NR : float
        Size ratio:
        $N_{\rm R} = \\dfrac{d_p}{d}$

    Parameters
    ----------
    dp : float
        Particle diameter [m]
    dc : float
        Collector diameter $d$ [m]
    '''
    return dp/dc


def noDim_Péclet(q,dc,Dm):
    '''
    Returns 
    ----------
    NPe : float
        Péclet number
        .. math:: N_{\rm Pe} = \\dfrac{qd}{D_m}
    
    Parameters
    ----------
    q : float
        Darcy velocity [m/s]
    dc : float
        Collector diameter $d$ [m]
    Dm : float
        Molecular diffusion coefficient [m²/s]
    '''
    return q*dc/Dm


def noDim_vanderWaals(A,T):
    '''
    Returns 
    ----------
    NvdW : float
        Van der Waals nondimensional number:
        .. math:: N_{\rm vdW} = \\dfrac{A}{k_BT}
    
    Parameters
    ----------    
    A : float
        Hamaker constant between particle and collector [J]
    T : float
        Temperature [K]
    '''
    return A/(BOLTZMANN*T)

def noDim_Gravitational(dp,rho_f,rho_p,T):
    '''
    Returns 
    ----------
    NGr : float
        Gravitational number:
        .. math:: N_{\rm gr} = \\dfrac{4\\pi r_p^4 (\rho_p - \rho_f)g}{3k_BT} = \\dfrac{\\pi d_p^4 (\rho_p - \rho_f)g}{12k_BT}
    
    Parameters
    ----------    
    dp : float
        Particle diameter [m]
    rho_f : float
        Fluid mass density [kg/m³]
    rho_p : float
        Particle mass density [kg/m³]
    T : float
        Temperature [K]
    '''
    return (PI*(dp**4)*(rho_p-rho_f)*g)/(12.*BOLTZMANN*T)

def attachmentRate_CFT(dc,theta,alpha,U,eta0): 
    '''
    Notes
    ----------
    Just the definition, check attachmentRate for a complete workflow

    Returns 
    ----------
    kAtt : float
        Returns the attachment rate coefficient: $k_{\rm att}$ calculated via colloid filtration theory
        .. math:: k_{\rm att} = \\dfrac{3 }{2d}(1-\theta)\alpha||\vec{U}||\\eta_0
    
    Parameters
    ----------
    dc : float
        Collector diameter (soil grain size) $d$ [m]
    theta : float (0.,1.)
        Porosity [-]
    alpha : float (0.,1.)
        Collision/attachment efficiency [-], i.e., the rate at which particles attach to the collector over the rate at which particles collide with the collector
        alpha = 1.0 for favorable attachment conditions, \alpha < 1.0 otherwise.  
    U : float
        Interstitial velocity [m/s]
    eta0 : float
        Collector efficiency $\\eta_0$ [-]
    '''
    return (3*(1-theta)*alpha*U*eta0)/(2*dc)

def attachmentRate(q,theta,dp,dc,visco,rho_f,rho_p,A,T,alpha=1.0,debug=False):
    '''
    Returns 
    ----------
    kAtt : float
        Using particle/medium/fluid parameters, it returns the attachment rate coefficient calculated via colloid filtration theory.
        ..math:: k_{\\rm att} = \\dfrac{3}{2d}(1-\\theta)\\alpha||\\vec{U}||\\eta_0
    
    Parameters
    ----------
    dp : float
        Particle diameter [m]
    dc : float
        Collector diameter $d$ [m]
    q : float
        Darcy velocity [m/s]
    visco: float
        Fluid dynamic viscosity $\\mu$ [N·s/m²]
    rho_f : float
        Fluid mass density [kg/m³]
    rho_p : float
        Particle mass density [kg/m³]
    A : float
        Hamaker constant between particle and collector [J]
    T : float
        Temperature [K]
    alpha : float (0.,1.)
        Collision/attachment efficiency [-] 
    debug : bool
        Prints a list of all calculations done
    '''

    #Molecular diffusion
    Dm  = molecularDiff(visco,dp,T)
    
    #Pore water velocity
    U   = poreVel(q,theta)
    
    #Non-dimensional numbers
    As  = happelParameter(theta)
    NR  = noDim_SizeRatio(dp,dc)
    NPe = noDim_Péclet(q,dc,Dm)
    NvW = noDim_vanderWaals(A,T)
    NGr = noDim_Gravitational(dp,rho_f,rho_p,T)
    
    #Collector efficiency
    etaD = collectorEfficiency_Diffusion(As,NR,NPe,NvW)
    etaI = collectorEfficiency_Interception(As,NR,NPe,NvW)
    etaG = collectorEfficiency_GDeposition(NGr,NR,NPe,NvW)
    eta0 = collectorEff(etaD,etaI,etaG)
    
    #Attachment rate
    kAtt = attachmentRate_CFT(dc,theta,alpha,U,eta0)
     
    #Print report
    if(debug):
        print("Diffusion coeff:  {0:.4E}".format(Dm))
        print("Darcy velocity:   {0:.4E}".format(q))
        print("Pore-water vel:   {0:.4E}".format(U))
        print("---")
        print("Happel parameter: {0:.4E}".format(As))
        print("NR number:        {0:.4E}".format(NR))
        print("NPe number:       {0:.4E}".format(NPe))
        print("NvW number:       {0:.4E}".format(NvW))
        print("NGr number:       {0:.4E}".format(NGr))
        print("---")
        print("etaD collector:   {0:.4E}".format(etaD))
        print("etaI collector:   {0:.4E}".format(etaI))
        print("etaG collector:   {0:.4E}".format(etaG))
        print("eta0 collector:   {0:.4E}".format(eta0))
        print("---")
        print("Attach rate   :   {0:.4E}".format(kAtt))
        
        # htmlOut = """
        # <b>Diffusion coeff:</b>  {0:.4E}</br>
        # <b>Darcy velocity</b>   {1:.4E}</mainbr>
        # <b>NR number</b>        {4:.4E}</br>
        # <b>NPe number</b>       {5:.4E}</br>
        # <b>NvW number</b>       {6:.4E}</br>
        # <b>NGr number</b>       {7:.4E}</br>
        # </br>
        # <b>etaD collector</b>   {8:.4E}</br>
        # <b>etaI collector</b>   {9:.4E}</br>
        # <b>etaG collector</b>   {10:.4E}</br>
        # <b>eta0 collector</b>   {11:.4E}</br>
        # </br>
        # <b>Attach rate   </b>   {12:.4E}
        # """.format(Dm,q,U,As,NR,NPe,NvW,NGr,etaD,etaI,etaG,eta0,kAtt)
        
    return kAtt
    

class CFT():
    def __init__(self,q,n,**kwargs) -> None:
        #Porosity
        self.dp = kwargs.get("dp",10.0E-6)     # m
        self.dc = kwargs.get("dc",2.0E-3)      # m
        self.mu = kwargs.get("mu",0.0008891)   # N·s/m²
        self.rhof = kwargs.get("rhof",999.79)  # kg/m³
        self.rhop = kwargs.get("rhop",1050.0)  # kg/m³
        self.Aham = kwargs.get("Aham",5.0E-21) # J
        self.Temp = kwargs.get("Temp",283.)    # K
        self.alpha = kwargs.get("alpha",1.0)  # -
        self._q = q
        self._n = n
        self._katt = attachmentRate(q,n,self.dp,self.dc,self.mu,self.rhof,self.rhop,self.Aham,self.Temp,self.alpha,debug=kwargs.get("debug",False))
    
    def __str__(self) -> str:
        return f"{self.katt:.3E}"

    @property
    def q(self):
        return self._q
    
    @q.setter
    def q(self,value):
        self._q = value
        self._katt = attachmentRate(self.q,self.n,self.dp,self.dc,self.mu,self.rhof,self.rhop,self.Aham,self.Temp,self.alpha)

    @property
    def n(self):
        return self._n
    
    @n.setter
    def n(self,value):
        self._n = value
        self._katt = attachmentRate(self.q,self.n,self.dp,self.dc,self.mu,self.rhof,self.rhop,self.Aham,self.Temp,self.alpha)

    @property
    def katt(self):
        return self._katt


if __name__ == "__main__":
    main()

