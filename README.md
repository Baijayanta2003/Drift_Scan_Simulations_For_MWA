## Drift Scan Simulations For MWA

How to simulate drift scan data for the MWA Telescope.

## Setting Up A Virtual environment

Activate the Virtual Environment.Go there for more https://github.com/Baijayanta2003/Python-Documentation-With-Sphinx#here-are-the-steps .

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

```bash

Filename: 1000003720-unSub.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 GroupsHDU       88   (3, 4, 1, 1, 1, 1)   float32   86625 Groups  9 Parameters
  1  AIPS FQ       1 BinTableHDU     28   1R x 5C   [1J, 1D, 1E, 1E, 1J]   
  2  AIPS AN       1 BinTableHDU     62   128R x 13C   [8A, 3D, 0D, 1J, 1J, 1E, 1A, 1E, 1E, 1A, 1E, 1E, 1E]   
  3  AIPS SU       1 BinTableHDU     69   1R x 19C   [1J, 20A, 1J, 4A, 1E, 1E, 1E, 1E, 1D, 1D, 1D, 1D, 1D, 1D, 1D, 1D, 1D, 1D, 1D]   
```
