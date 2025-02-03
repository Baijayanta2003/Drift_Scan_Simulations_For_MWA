"""
**This contains useful parameters of MWA.**

`The Important things are listed here:`

- KB - Boltzmann Constant in :math:`Jy\ m^2/mK`

- nside - Healpix parameter.

- Amp - Amplitude of APS

- beta - Power index.

- Nrea - No of Realizations.

- ra_ptg - Pointing direction of telescope.

- b - Dimensions of Square Aperture in m.

- nu_c - Central Frequency.

- dec_mwa - Declination of the telescope.

"""
import builtins as blt
import numpy as np
import scipy as sp 
import healpy as hp
import time    
from functools import partial

# ======== constants ========#
blt.C = 299.792  # velocity of light in Km/s

def paramsim():
    """
    This Function contains the simulation parameters.
    """
    blt.KB = 1.38		 	# Boltzman constant in Jy.m^2/mK
    #blt.pi=np.pi
    # Healpy parameters for simulation
    blt.nside =32  # Nside for Healpix maps which sets the resolution
    
    # Input < C_{\ell}= Amp * \ell^{\beta} >parameters for simulation
    blt.Amp = 100  # amplitude of the input angular power spectrum in mK^2
    blt.beta = -2
    
     # No of realizations
    blt.Nrea = 10 # number of realizations of the simulations
    blt.ra_ptg=0
    
    # # observation input parameters
    blt.ra_ini = 0  # initial pointing RA of MWA during observation in degree
    # initial phase center of observation w.r.t. dec of MWA (-26.7 deg)
    blt.dec_ini = 0.0
    blt.obs_time = 0.0  # total observation time in mins
    blt.ptg_time = 0.0  # pointing time in secs

def params(tele):
	if (tele == 'MWA'):
		#MWA system params:
		blt.b = 4.0				#antenna dimensions in meters
		blt.nu_c = 154.25 		#central frequency in MHz
		#blt.Nchan = 154.24			# No. of channel
		blt.B_bw = 30.72			#in MHz units bandwidth
		#blt.dnu_c = blt.B_bw/blt.Nchan 		# channel width in MHz units
		#blt.lam_c = C/nu_c 		# Wavelength corresponding to central frequency in meter
		blt.NA = 126			#no of antennas
		blt.Nbl = NA*(NA-1)/2	#no of baselines
		blt.A_dBdT = 1			#(A*(dB/dT)^2)-1 (mk/Jy)^2 # 4 factor due to modified pbeam
		blt.dec_mwa = -26.79522222222223		#over head declination of MWA in degree
#		blt.Tsys = 150 			# system temperature in K
#		blt.eta = 0.6		 	# efficiency of the telescope	
	
	

			
		#Cosmology of MWA	
#		blt.r = 6845.5			#comoving distance in Mpc
#		blt.rp = 11.5			#dr/dnu in Mpc/MHz
#		blt.dBdT = 3.27			#(dB/dT)_326.5 in Jy/mK 
	
		
	else :
		print ("Please enter telescope name")
	return (0)

