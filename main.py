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
background = 0

addStellarCatalog(hdu,catalog1)
addGalaxyCatalog(hdu,catalog2)
addBackground(hdu,background)

# plot_image(hdu,log_scale=True)
convolvePSF (hdu, 5)

hdulist.close()