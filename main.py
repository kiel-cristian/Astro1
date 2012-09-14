import pyfits as pf
import pylab as pl
import scipy.signal as sig
import scipy as sp
import numpy as np
from fits_lib import *
from metodos_t1 import *

hdulist = pf.open('blank.fits')
hdu= hdulist[0].data
catalog1="stellar.dat"
catalog2="galaxy.dat"
background = 1000.0
sigma_psf = 5.0
sigma_noise = 20.0

addStellarCatalog(hdu,catalog1)
addGalaxyCatalog(hdu,catalog2)
# addGalaxy(hdu, 10.0, 22.1, 0.11, 4.0, 10.0, 1.2, 0.5)
addBackground(hdu,background)
convolvePSF (hdu, sigma_psf)
addNoise(hdu,sigma_noise)
plot_image(hdu,log_scale=True)

hdulist.close()