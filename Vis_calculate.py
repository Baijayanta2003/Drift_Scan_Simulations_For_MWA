r"""
This is the file which calls the other modules
This can simulate the visibilities for a given Baseline file.
The Visibility is given by

.. math::

    \mathcal{V}(\vec{U},\nu)=Q_\nu \int_{UH}\ d\Omega_{\hat{n}}\ T(\hat{n},\nu) \ A(\Delta\hat{n},\nu)\ e^{2\pi i\vec{U}\cdot\Delta\hat{n}}

In Healpix we discretize the sky and the integral becomes summation.But the integral is restricted to only upper hemisphere.
So we find out those pixel indexing and sum only over them.

:strong:`The discritized version is given as`:

.. math::

    \mathcal{V}(\alpha_{p},\vec{U},\nu) = Q_{\nu}\ \Delta\Omega_{pix}\ \sum_{q} \ T(\hat{n}_{q},\nu) \ A(\Delta\hat{n}_{q},\nu) \ e^{2\pi i\vec{U}\cdot\Delta\hat{n}}
    
    
where

.. math::
    
    Q_{\nu}  &=2k_B /\lambda^2
    
    \Delta\hat{n}&= \hat{n}- \hat{p}
    
    \vec{\text{U}} &= u\ {\hat{e}_ 1 (\alpha_{p})} + v \ {\hat{e}_ 2 (\alpha_{p})} + w \ {\hat{e}_3(\alpha_{p})}
    

Where :math:`\Delta\Omega_{pix}` refers to the solid angle subtended by each simulation pixel and :math:`\hat{p}` refers to the pointing direction of the antenna.Here :math:`\hat{p}=\hat{e}_3`.

:strong:`Overview`

- First we generate the **Gaussian Random Field** from given **Angular Power Spectrum**.

- Then We find out Which Pixels are in the Upper Hemisphere given by the pointing direction.Here Pointing Direction is vertically overhead.

- Calculate the **Primary Beam Pattern** for the Telescope.

- Calculate the components of :math:`\hat{n}` along the basis vectors :math:`\hat{e}_1,\hat{e}_2,\hat{e}_3`.This will help us to determine the dot product :math:`\vec{\text{U}}\cdot\hat{n}`.

- Now We calculate the **Phase factor** :math:`e^{2\pi i\vec{\text{U}}\cdot\hat{n}}`.

- Multiply the **GRF and PB** first and then multiply by phase factor ,sum over the pixels.The whole thing is `done by matrix multiplication.`

- Finally multiply by :math:`e^{-2\pi i\vec{\text{U}}\cdot\hat{n}}=e^{-2\pi i w}`,so that we get the correct Phase Factor.

- The Visibility is Simulated.

The **Basis** vectors :math:`\hat{e_1},\hat{e_2},\hat{e_3}` are dependent on the R.A. :math:`\alpha` and DEC :math:`\delta` \
and the components of them in the cartesian system :math:`xyz` is given by:

.. math::

    \hat{e_1}&=[-\sin{\alpha},\cos{\alpha},0]\\
    \hat{e_2}&=[-\cos{\alpha}\sin{\delta},-\sin{\alpha}\sin{\delta},\cos{\delta}]\\
    \hat{e_3}&=[\cos{\alpha}\cos{\delta},\sin{\alpha}\cos{\delta},\sin{\delta}]


"""

#Import Necessary Libraries
import numpy as np
import scipy as sp
import healpy as hp
import builtins as blt
import time
from functools import partial,lru_cache
from sys import argv
import os
from astropy.table import Table
from astropy.io import fits


from GRF_Generate import *
from Params_MWA import *
from PB_Phase import *
from Vis_Gen import *

#call the paramertes
paramsim()
params("MWA")

#check UAPS or APS
if (Amp == 1) and (beta == 0):
    uaps = True
else:
    uaps = False
if uaps:
    print('UAPS simulation')
else:
    print(f'APS simulation for Amp = {Amp} and beta = {beta}')
    
Npix = hp.nside2npix(nside)
start=time.time()

#Generate the GRF
map_grf = np.zeros((Nrea, Npix))
for ii in range(Nrea):
    iseed = 10*ii+1 # 0, 11, ... etc. iseed must chane with realizations
    map_grf[ii] =GRF_hpmap_gen(nside = nside, l_max =nside * 3 - 1 , 
                                   nu = nu_c, func =APS_func, iseed = iseed) # -1 for random
print(f"{Nrea} Sky maps generated.")

end1 = time.time()
hours, rem = divmod(end1-start, 3600)
minutes, seconds = divmod(rem, 60)
print("Took time for generating GRF",end= '\t ')
print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds))
print(map_grf.shape)


#into fits file you can modify it as per you want

#Get the path and name of fits file
# in_fits = str(argv[1])
# print(in_fits)
# head, tail = os.path.split(in_fits)
# fitsName = os.path.splitext(tail)[0]
# print('head : ', head) #path
# print('\ntail : ', tail) #name total
# print('\nfitsName : ', fitsName) #only name


# # working with the input file
# hdulist = fits.open(in_fits,mode='readonly')
# # print (hdulist.info()) 				#to check different headers
# data_tmp = hdulist[0].data 				#to save visibilities in data
# dataT = Table(data_tmp)                 #to convert data from array form to table form, readable

# nu_c = hdulist[0].header['CRVAL4']*1.e-6
# pol = hdulist[0].header['NAXIS3']
# Nchan = hdulist[0].header['NAXIS4']
# dnu_c = hdulist[0].header['CDELT4']*1.e-6

# #EXTRACT the baselines
# u = dataT['UU']*3.e8 #*nu_c for baseline unit   #*3.e8 (for meter unit)
# v = dataT['VV']*3.e8
# w = dataT['WW']*3.e8
# bln = np.stack((u,v,w),axis = 1)
# print(bln.shape)

print('RA = %s'%ra_ptg)
print('DEC = %s'%dec_mwa)

#for loading from a file
bln=np.load("MWA_BASELINE_7875.npy")
print(bln.shape)


#generate the visibilities 
vis=visgen_mwa_multi(nside, map_grf, ra_ptg, dec_mwa, bln,nu_c)
print(vis.shape)



#save file in a directory
output_dir = "./simulated_visibility"
#check if directory exists,if not make one
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

    
#If you want to save it in a .npy file uncomment this
if uaps:
    op_filename = f"{output_dir}/Nside_{nside}_RA_{ra_ptg}_UAPS"
else:
    op_filename = f"{output_dir}/Nside_{nside}_RA_{ra_ptg}_APS_{Amp}_beta{beta}"

np.save(op_filename,vis)


#data output for multi channel multiple realization different file
#Freq_chan=dataT['DATA'].shape[4] #no of frequency channels #np.tile to repeat the same value
# for i in range(Nrea):
#     for k in range(pol):
#         hdulist[0].data["DATA"][:, 0, 0, 0,:, k, 0] =np.tile(vis[i].T.real.astype(np.float64)[:,np.newaxis],Freq_chan) 
#         hdulist[0].data["DATA"][:, 0, 0, 0,:, k, 1] =np.tile(vis[i].T.imag.astype(np.float64)[:,np.newaxis],Freq_chan)
#         if uaps:
#             op_filename = f"{output_dir}/Nside_{nside}_RA_{ra_ptg}_UAPS_{i+1}.fits"
#         else:
#             op_filename = f"{output_dir}/Nside_{nside}_RA_{ra_ptg}_APS_{Amp}_beta{beta}_{i+1}.fits"
#         hdulist.writeto(op_filename, overwrite=True)



#data output for multiple realization in same file along the frequency channels
# for j in range(pol):
#     hdulist[0].data['DATA'][:,0,0,0,:Nrea,j, 0] = vis.T.real.astype(np.float64)
#     hdulist[0].data['DATA'][:,0,0,0,:Nrea,j, 1] = vis.T.imag.astype(np.float64)

# if uaps:
#     op_filename = f"{output_dir}/Nside_{nside}_RA_{ra_ptg}_UAPS_{i+1}.fits"
# else:
#     op_filename = f"{output_dir}/Nside_{nside}_RA_{ra_ptg}_APS_{Amp}_beta{beta}_{i+1}.fits"

# hdulist.writeto(op_filename,overwrite = True) # note that the files will be overwritten

#hdulist.close()
end2 = time.time()
hours, rem = divmod(end2-start, 3600)
minutes, seconds = divmod(rem, 60)
print("Took time all total",end='\t')
print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds))
