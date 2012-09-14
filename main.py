import pyfits as pf
import pylab as pl
import scipy.signal as sig
import scipy as sp
import numpy as np

import metodos_t1 as m
hdulist = pf.open('blank.fits')
hdu= hdulist[0].data
catalog1="stellar.dat"
catalog2="galaxy.dat"

m.addStellarCatalog(hdu,catalog1)
m.addGalaxyCatalog(hdu,catalog2)
hdulist.close()
