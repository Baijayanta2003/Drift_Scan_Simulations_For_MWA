r"""
**This can generate Gaussian Random Field For a given APS.**

**Functions:**

- APS_func

- GRF_hpmap_gen

"""

#The aim of this file is to produce the Gaussian Random Field from a given input Angular Power Spectrum
#Import Necessary Libraries
import numpy as np
import scipy as sp 
import healpy as hp
import builtins as blt
from functools import partial

import Params_MWA as psim
psim.paramsim()
psim.params("MWA")

def APS_func(l):  
    r"""
    This is the input model Angular Power Spectrum(APS) Function.
    
    .. math::
        
            C_{\ell}=Amp \ \ell^{\beta}
    
    Amp here defines the Amplitude of APS.beta(:math:`\beta`) here is the power index.

    Parameters
    ----------
    l : float or array 
        The Angular Multipole value.
        
    Returns
    -------
    float or array 
        The Angular Power Spectrum Value for that multipole(s).
        
    """
    return Amp * np.power(l, beta)


def GRF_hpmap_gen(nside, l_max, nu, func, iseed = -1):
    r"""
    Generates The GRF(Gaussian Random Field for all the pixels) for the input APS.

    Parameters
    ----------
    nside : int  
        Healpix parameter Usually power of 2.
    l_max : int 
        Maximum Value of Angular Multipole.
    nu :  float 
        The Observing frequency :math:`\nu` in MHz.
    func : function 
        The input APS.
    iseed : int or float 
        A seed value.In order to get different realizations for same APS.
        
    Returns
    -------
    array 
        Contains the GRF value for each pixels.Shape :math:`1 \times12\ {nside}^2`.
   
    """
    kk= np.zeros(l_max)
    kk[1:]=func(np.arange(1,l_max, dtype=np.float64))
    if (iseed != -1):
        np.random.seed(iseed)
        #print(f"iseed={iseed}")
    map_grf = hp.synfast(kk,nside=nside, lmax=l_max)
    return map_grf
