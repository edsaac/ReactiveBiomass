from scipy.optimize import root
import numpy as np

def _aux_wetfrondetpth(        
        L: float,
        t:  float, 
        Ks: float,
        h0: float, 
        hs: float, 
        theta_s: float,
        theta_i: float):

    dtheta = theta_s - theta_i
    dhead = h0 - hs
    error = L - dhead * np.log(1.0 + L/dhead) - Ks*t/dtheta
    return error

def calc_wetfront_depth(
        t:  float, 
        Ks: float,
        h0: float, 
        hs: float, 
        theta_s: float,
        theta_i: float
        ):
    """
    Calculates the wet front depth at a given time 

    Parameters:
    t: float
        Time [T]
    Ks: float
        Saturated hydraulic conductivity [L/T]
    h0: float
        Pondind head [L]
    hs: float
        Suction head [L]
    theta_s: float
        Water content at saturation [-]
    theta_i: float
        Water content at the soil [-]
        This has to be consistent with the suction head given that
        theta_i = vanGenuchten(hs)

    Returns:
    depth: [L]

    """ 
    solve = root(
        _aux_wetfrondetpth,
        x0 = 10.0,
        args = (t, Ks, h0, hs, theta_s, theta_i),
        method = 'lm'
    )

    return solve

def main():
    pass

if __name__ == '__main__':
    main()