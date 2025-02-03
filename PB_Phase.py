r"""
**This can Primary Beam Pattern of MWA**.

**Functions:**

- needful_pixels

- hat_n

- beam_mwa

- Primary_Beam_generate

- dot_cal_superfast

- calculate_phase

"""

#import the necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import builtins as blt
import time
from functools import partial,lru_cache
import numexpr as ne
#import the parameters
import Params_MWA as psim
psim.paramsim()
psim.params("MWA")

@lru_cache(maxsize=2)
def needful_pixels(ra_ptg,dec_ptg,nside,radius=np.radians(90)):
    r"""
    This calculates the needful pixels which makes less than radius radian
    w.r.t the pixel in the direction given by ra_ptg,dec_ptg.

    Parameters
    ----------
    ra_ptg : int or float 
        RA of the pointing direction.(In Degrees).
    dec_ptg :  int or float 
        DEC of the pointing direction.(In Degrees).
    nside : int 
        Healpix parameter.Usually Power of 2.
    radius  : int or float 
        The maximum angle subtended by a pixel given the pixel center by ra_ptg,dec_ptg.(In Radians)
        By Default 90 degrees only upper hemisphere.
    Returns
    -------
    pix_mask : array 
        The pixels lie within the given disk radius.Shape :math:`6\ {nside}^2 \times 1`.
    """
    MWA_vec = hp.ang2vec(ra_ptg, dec_ptg, lonlat=True)
    ipix_mask = hp.query_disc(nside=nside, vec=MWA_vec, radius=np.radians(90))
    return ipix_mask

@lru_cache(maxsize=2)
def hat_n(nside,ra_ptg,dec_ptg):
    r"""
    This can calculate the :math:`xyz` components of :math:`\hat{n}` for all the pixels which makes angle less
    than 90 deg by default w.r.t. the pixel in the direction given by ra_ptg,dec_ptg.

    Parameters
    ----------
    nside : int 
        Healpix parameter.Usually Power of 2.
    ra_ptg : int or float 
        RA of the pointing direction.(In Degrees).
    dec_ptg : int or float 
        DEC of the pointing direction.(In Degrees).

    Returns
    -------
    Array 
        Containg the :math:`xyz` comps of :math:`\hat{n}`. :math:`6\ {nside}^2 \times 1`

    """
    theta, phi = hp.pix2ang(nside,needful_pixels(ra_ptg,dec_ptg,nside))
    hat_n=np.array([np.sin(theta) * np.cos(phi),np.sin(theta) * np.sin(phi),np.cos(theta)]).T
    return hat_n

def beam_mwa(nu,ne1,ne2):
    r"""
    This is the primary beam function of the MWA Telescope.
    
    .. math::

       A(\hat{n},\nu) = sinc^2\left(\frac{{\pi}b{\nu}\ \hat{n}\cdot\hat{e}_1 (\alpha_{p})}{c}\right)sinc^2\left(\frac{{\pi}b{\nu}\ \hat{n}\cdot\hat{e}_2 (\alpha_{p})}{c}\right)

    
    Parameters
    ----------
    nu : float 
        Observing Frequency :math:`\nu` in MHz unit.
    ne1 : float or array 
        The dot product by :math:`\hat{n}\cdot\hat{e}_ 1`.
    ne2 : float or array 
        The dot product by :math:`\hat{n}\cdot\hat{e}_ 2`.
    Returns
    -------
    float or array 
        The primary beam pattern value for square aperture(MWA).
                
    """
    return (np.sinc((b*nu/C)*ne1)**2.0)*(np.sinc((b*nu/C)*ne2)**2.0)
def Primary_Beam_generate(nside,a0,d0,nu,beam_mwa):
    r"""
    This is the main function which can generate the Primary Beam (for all the pixels which makes angle less
    than 90 deg by default w.r.t. the pixel in the direction given by a0,d0) for MWA given the RA,DEC and frequency.

    Parameters
    ----------
    nside : int 
        Healpix parameter.Usually Power of 2.
    a0 : int or float 
        RA of the pointing direction.(In Degrees)
    d0 : int or float 
        DEC of the pointing direction.(In Degrees)
    nu : float 
        Observing Frequency :math:`\nu` in MHz unit.
    beam_mwa : function 
        The PB pattern of MWA.
    Returns
    -------
    PB :array 
        The value of Primary beam of The Telescope.Shape :math:`6\ {nside}^2 \times 1`.
                
    """
    start=time.time()
    hatn=hat_n(nside,a0,d0)
    #convert RA,Dec into radians
    a0=np.deg2rad(a0)   #RA of the phase center in radians
    d0 = np.deg2rad(d0) #Phase center dec of the observation  in radians
    #create the basis vectors
    e1 = np.array([-np.sin(a0), np.cos(a0), 0])
    e2 = np.array([-np.cos(a0)*np.sin(d0), -np.sin(d0)*np.sin(a0), np.cos(d0) ])
    e3 = np.array([np.cos(a0)*np.cos(d0), np.sin(a0)*np.cos(d0), np.sin(d0) ])
    # Calculate the dot product e1 . hat_n for all pixels
    ne1 = np.dot(hatn,e1) #n⋅e1
    # Calculate the dot product e2 . hat_n for all pixels
    ne2 = np.dot(hatn,e2) #n⋅e2
    PB=beam_mwa(nu,ne1,ne2)
    end=time.time()
    
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    print(f"Time taken To calculate Beam Pattern for RA {np.degrees(a0)} , Dec {np.degrees(d0)}")
    print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds))
    
    return PB 
def dot_cal_superfast(nside,ra_ptg,dec_ptg):
    r"""
    This function can calculate the :math:`xyz` component of :math:`\hat{n}` for all the pixels which makes angle less
    than 90 deg by default w.r.t. the pixel in the direction given by ra_ptg,dec_ptg 
    along the basis vectors :math:`\hat{e_1},\hat{e_2},\hat{e_3}` given by the RA,DEC.

    Parameters
    ----------
    nside : int 
        Healpix parameter.Usually Power of 2.
    ra_ptg : int or float 
        RA of the pointing direction.(In Degrees)
    dec_ptg :int or float 
        DEC of the pointing direction.(In Degrees)
    Returns
    -------
    dot_products : array 
        Contains the :math:`xyz` comps for those pixels.Shape :math:`6\ {nside}^2 \times 3`.
                
    """
    
    #convert the RA,Dec into radians
    a0=np.deg2rad(ra_ptg)   #RA of the phase center in radians
    d0 = np.deg2rad(dec_ptg) #Phase center dec of the observation  in radians
    #create the basis functions e1,e2,e3
    e1 = np.array([-np.sin(a0), np.cos(a0), 0])
    e2 = np.array([-np.cos(a0)*np.sin(d0), -np.sin(d0)*np.sin(a0), np.cos(d0)])
    e3 = np.array([np.cos(a0)*np.cos(d0), np.sin(a0)*np.cos(d0), np.sin(d0)])
    e=np.stack((e1,e2,e3),axis=1)
    dot_products = np.dot(hat_n(nside,ra_ptg,dec_ptg),e)
    return dot_products
def calculate_phase(dot_product,bl):
    r"""
    Given the baselines it can calculate the phase factor for
    for all the pixels which makes angle less than 90 deg by default w.r.t. the pixel in 
    the direction given by ra_ptg,dec_ptg of the Telescope.

    Parameters
    ----------
    dot_product : array  
        Contains the :math:`xyz` comps for those pixels that makes less than 90 degree angle. 
    bl : array 
        The array conatins the components of the baselines.(defined by basis vectors :math:`\hat{e_1},\hat{e_2},\hat{e_3}`).
        Shape :math:`N_{Baselines}\times 3`
    Returns
    -------
    array 
        The phase factor . Shape :math:`6\ {nside}^2 \times  N_{Baselines}`.
        Excluding the factor :math:`e^{-2 \pi i \vec{U} \cdot \hat{p}}` 
        Which will be multiplied later on.
    """
    dot=dot_product@bl.T
    ppi=np.pi
    return ne.evaluate('exp(2*ppi*1j*dot)') #np.exp(2j*np.pi*dot).astype(np.complex64)

