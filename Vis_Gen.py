import numpy as np
import scipy as sp
import healpy as hp
import builtins as blt
import time
from functools import partial,lru_cache


def visgen_mwa_multi(nside, in_sky_map, ra_ptg, dec_ptg, bl_file, nu):
    r"""
    This is the main module which calculates the visibilities for all the baselines
    given the RA,DEC,Frequency,GRF map and the baseline files.(Also does Multiple Realizations)

    Parameters
    ----------
    nside : int 
        Healpix parameter.Usually Power of 2.
    in_sky_map : array 
        Contains the Sky Brightness Temperature simulated from model APS for all the pixels.(Can contain Multiple Realizations).Shape :math:`N_{Realizations}\times 12\ {nside}^2`.
    ra_ptg : int or float 
        RA of the pointing direction.(In Degrees)
    dec_ptg : int or float 
        DEC of the pointing direction.(In Degrees)
    bl_file : array 
        The array conatins the components of the baselines.(defined by basis vectors :math:`\hat{e_1},\hat{e_2},\hat{e_3}`).Shape :math:`N_{Baselines}\times 3`
    nu : float 
        Observing Frequency :math:`\nu` in MHz unit.

    Returns
    -------
    vis : Complex array 
        Contains the simulated visibility for all the baselines given by bl_file.(Along axis=0 is the different Realizations,axis=1 is the baselines).Shape :math:`N_{realizations}\ \times\ N_{Baselines}`
    
    """
    start = time.time()
    Npix = hp.nside2npix(nside)     #Number of pixels given N_side
    res = hp.nside2resol(nside, arcmin = True)      #resolution given N_side
    lmax = nside * 3 - 1        #maximum number of \ell multipoles used
    print("Nside=",nside, 'resolution=', res, 'min', Npix, 'pixels' )
    dOmega = hp.nside2resol(nside)**2.0  # resolution in s.rad
    Pb_map=Primary_Beam_generate(nside,ra_ptg,dec_ptg,nu_c,beam_mwa)
    dot_product=dot_cal_superfast(nside,ra_ptg,dec_ptg)
    lam = C/nu  # wavelength in m
    print(lam)
    Q_nu = 2.0 * KB / lam**2.0  # in Jy/mK

    ipix_mask= needful_pixels(ra_ptg,dec_ptg,nside)
    GRF_PB_Product =map_grf[:,ipix_mask]*Pb_map
    print(GRF_PB_Product.shape)
    nbl, nn = bln.shape  # number of baselines
    chunk_size=400 #adjust the chunk size as per system requirements.
    chunks = np.array_split(np.arange(nbl), np.arange(chunk_size, nbl, chunk_size))
    #print(chunks)
    vis = np.zeros((Nrea, nbl), dtype=np.complex64)
    for ii in chunks:
        phase=calculate_phase(dot_product,bln[ii,:]) #phase calculation
        vis[:, ii]=GRF_PB_Product@phase
    vis =  Q_nu * dOmega * vis 
    #compensate for -u dot p term for all baselines
    ppi=np.pi
    w=bln[:,2]
    vis=ne.evaluate('vis*(exp(-2.0*ppi*1j*w))')
    end = time.time()
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds))  
    return vis
