## Drift Scan Simulations For MWA

How to Simulate Drift scan data for the MWA Telescope.

## Setting Up A Virtual environment

Activate the Virtual Environment.Go there for more details : https://github.com/Baijayanta2003/Python-Documentation-With-Sphinx#here-are-the-steps .

**Once activated Install the required packages**.

## Required Packages

```bash
pip install numpy healpy matplotlib astropy 

```

- Change your directory path to where you have downloaded the scripts.
- In the Terminal Write the following command:
```bash
python3 Vis_calculate.py 1000003720-unSub.fits
```
Followed by the name of fits file.
This will generate a output fits file in a directory named  `simulated_visibility`.

The structure of the fits file (`1000003720-unSub.fits`)I have used is given as follows:

```python

Filename: 1000003720-unSub.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 GroupsHDU       88   (3, 4, 1, 1, 1, 1)   float32   86625 Groups  9 Parameters
  1  AIPS FQ       1 BinTableHDU     28   1R x 5C   [1J, 1D, 1E, 1E, 1J]   
  2  AIPS AN       1 BinTableHDU     62   128R x 13C   [8A, 3D, 0D, 1J, 1J, 1E, 1A, 1E, 1E, 1A, 1E, 1E, 1E]   
  3  AIPS SU       1 BinTableHDU     69   1R x 19C   [1J, 20A, 1J, 4A, 1E, 1E, 1E, 1E, 1D, 1D, 1D, 1D, 1D, 1D, 1D, 1D, 1D, 1D, 1D]   
```
The header is given as:
```python
SIMPLE  =                    T /Standard FITS format                            
BITPIX  =                  -32 /Floating point values                           
NAXIS   =                    7                                                  
NAXIS1  =                    0 /Random groups, NOT image                        
NAXIS2  =                    3                                                  
NAXIS3  =                    4                                                  
NAXIS4  =                    1                                                  
NAXIS5  =                    1                                                  
NAXIS6  =                    1                                                  
NAXIS7  =                    1                                                  
EXTEND  =                    T /Tables may follow                               
BLOCKED =                    T /File may be blocked                             
GROUPS  =                    T /Random Group UV data                            
PCOUNT  =                    9 /Number of random parameters                     
GCOUNT  =                86625 /Number of groups (rows) in the file             
                                                                                
EPOCH   =   2.000000000000E+03                                                  
BSCALE  =   1.000000000000E+00                                                  
BZERO   =   0.000000000000E+00                                                  
BUNIT   = 'UNCALIB '                                                            
CTYPE2  = 'COMPLEX '                                                            
CRVAL2  =   1.000000000000E+00                                                  
CDELT2  =   1.000000000000E+00                                                  
CRPIX2  =   1.000000000000E+00                                                  
CROTA2  =   0.000000000000E+00                                                  
CTYPE3  = 'STOKES  '                                                            
CRVAL3  =  -5.000000000000E+00                                                  
CDELT3  =  -1.000000000000E+00                                                  
CRPIX3  =   1.000000000000E+00                                                  
CROTA3  =   0.000000000000E+00                                                  
CTYPE4  = 'FREQ    '                                                            
CRVAL4  =   1.542750000000E+08                                                  
CDELT4  =   7.200000000000E+05                                                  
CRPIX4  =   1.000000000000E+00                                                  
CROTA4  =   0.000000000000E+00                                                  
CTYPE5  = 'IF      '                                                            
CRVAL5  =   1.000000000000E+00                                                  
CDELT5  =   1.000000000000E+00                                                  
CRPIX5  =   1.000000000000E+00                                                  
CROTA5  =   0.000000000000E+00                                                  
CTYPE6  = 'RA      '                                                            
CRVAL6  =   0.000000000000E+00                                                  
CDELT6  =   1.000000000000E+00                                                  
CRPIX6  =   1.000000000000E+00                                                  
CROTA6  =   0.000000000000E+00                                                  
CTYPE7  = 'DEC     '                                                            
CRVAL7  =   0.000000000000E+00                                                  
CDELT7  =   1.000000000000E+00                                                  
CRPIX7  =   1.000000000000E+00                                                  
CROTA7  =   0.000000000000E+00                                                  
PTYPE1  = 'UU      '                                                            
PSCAL1  =   1.000000000000E+00                                                  
PZERO1  =   0.000000000000E+00                                                  
PTYPE2  = 'VV      '                                                            
PSCAL2  =   1.000000000000E+00                                                  
PZERO2  =   0.000000000000E+00                                                  
PTYPE3  = 'WW      '                                                            
PSCAL3  =   1.000000000000E+00                                                  
PZERO3  =   0.000000000000E+00                                                  
PTYPE4  = 'DATE    '           /Day number                                      
PSCAL4  =   1.000000000000E+00                                                  
PZERO4  =   0.000000000000E+00                                                  
PTYPE5  = 'DATE    '           /Day fraction                                    
PSCAL5  =   1.000000000000E+00                                                  
PZERO5  =   0.000000000000E+00                                                  
PTYPE6  = 'BASELINE'                                                            
PSCAL6  =   1.000000000000E+00                                                  
PZERO6  =   0.000000000000E+00                                                  
PTYPE7  = 'FREQSEL '                                                            
PSCAL7  =   1.000000000000E+00                                                  
PZERO7  =   0.000000000000E+00                                                  
PTYPE8  = 'SOURCE  '                                                            
PSCAL8  =   1.000000000000E+00                                                  
PZERO8  =   0.000000000000E+00                                                  
PTYPE9  = 'INTTIM  '                                                            
PSCAL9  =   1.000000000000E+00                                                  
PZERO9  =   0.000000000000E+00                                                  
OBJECT  = 'MULTI   '                                                            
DATE-OBS= '2016-10-12T12:00:00.000000'                                          
TELESCOP= 'MWA     '                                                            
INSTRUME= 'MWA     '                                                            
OBSERVER= '        '                                                            
SORTORD = 'TB      '                                                            
SPECSYS = 'TOPOCENT'                                                            
HISTORY AIPS WTSCAL = 1.0                                                       
                                                                                
ORIGIN  = 'casacore '
```
For more details Please read the `documentation` carefully.

## Important Note:

You can save the visibility in a binary file named `Vis_Sim_{Amp}_{bln.shape[0]}.npy` format.

If you have provided a `fits` file depending on the fits file you can write the following styles to save it in different ways.

- If your file has multiple frequency channels and you just want to put the different realizations data in there write or uncomment the following part in the file `Vis_calculate.py`.(line no 188)

```python

for j in range(pol):
    hdulist[0].data['DATA'][:,0,0,0,:Nrea,j, 0] = vis.T.real.astype(np.float64)
    hdulist[0].data['DATA'][:,0,0,0,:Nrea,j, 1] = vis.T.imag.astype(np.float64)

if uaps:
    op_filename = f"{output_dir}/Nside_{nside}_RA_{ra_ptg}_UAPS_{i+1}.fits"
else:
    op_filename = f"{output_dir}/Nside_{nside}_RA_{ra_ptg}_APS_{Amp}_beta{beta}_{i+1}.fits"

hdulist.writeto(op_filename,overwrite = True) # note that the files will be overwritten


hdulist.close()
```
- If you want to put the same visibility across all the frequency channels and make different fits file for different realizations do this.(line no 174)
```python
#data output for multi channel multiple realization different file
Freq_chan=dataT['DATA'].shape[4] #no of frequency channels
 for i in range(Nrea):
     for k in range(pol):
         hdulist[0].data["DATA"][:, 0, 0, 0,:, k, 0] =np.tile(vis[i].T.real.astype(np.float64)[:,np.newaxis],Freq_chan) 
         hdulist[0].data["DATA"][:, 0, 0, 0,:, k, 1] =np.tile(vis[i].T.imag.astype(np.float64)[:,np.newaxis],Freq_chan)
         if uaps:
             op_filename = f"{output_dir}/Nside_{nside}_RA_{ra_ptg}_UAPS_{i+1}.fits"
         else:
             op_filename = f"{output_dir}/Nside_{nside}_RA_{ra_ptg}_APS_{Amp}_beta{beta}_{i+1}.fits"
         hdulist.writeto(op_filename, overwrite=True)
```
In my case the fits file has only one frequency channel. 


