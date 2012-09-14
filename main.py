import pyfits as pf
import pylab as pl
import scipy.signal as sig
import scipy as sp
import numpy as np
from fits_lib import *

import metodos_t1 as m
hdulist = pf.open('blank.fits')
hdu= hdulist[0].data
catalog1="stellar.dat"
catalog2="galaxy.dat"

m.addStellarCatalog(hdu,catalog1)
m.addGalaxyCatalog(hdu,catalog2)
hdulist.close()

hdr1 = get_fits_header('blank.fits')
img1 = get_fits_matrix('blank.fits')
print "\n>> ARCHIVO 1:\n"
print_header(hdr1)
plot_image(img1,log_scale=True)