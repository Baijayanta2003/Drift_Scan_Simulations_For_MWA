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
