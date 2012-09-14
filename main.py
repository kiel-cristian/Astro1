#!/usr/local/bin/python
# -*- coding: utf-8 -*-
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
f_cut = 9
params = [f_cut]

print "Añadiendo Estrellas"
addStellarCatalog(hdu,catalog1)
print "Añadiendo Galaxias"
addGalaxyCatalog(hdu,catalog2)
# addGalaxy(hdu, 10.0, 22.1, 0.11, 4.0, 10.0, 1.2, 0.5)
print "Añadiendo Background"
addBackground(hdu,background)
print "Convolucionando"
convolvePSF (hdu, sigma_psf)
print "Agregando Ruido"
addNoise(hdu,sigma_noise)
print "Filtrando Imagen"
filterImage(hdu,params)
# plot_image(hdu,log_scale=True)

hdulist.close()