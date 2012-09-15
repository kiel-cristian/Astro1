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
header = hdulist[0].header

catalog1="stellar.dat"
catalog2="galaxy.dat"
background = 1000.0
sigma_psf = 6.0
sigma_noise = 10.0
f_cut = 20
params = [f_cut,header]

print "Añadiendo Estrellas"
addStellarCatalog(hdu,catalog1)
pf.append('stars.fits',hdu,header)

print "Añadiendo Galaxias"
addGalaxyCatalog(hdu,catalog2)
pf.append('galaxys.fits',hdu,header)

print "Añadiendo Background"
addBackground(hdu,background)
pf.append('background.fits',hdu,header)

print "Convolucionando"
convolvePSF (hdu, sigma_psf)
pf.append('convolve.fits',hdu,header)

print "Agregando Ruido"
addNoise(hdu,sigma_noise)
pf.append('noise.fits',hdu,header)

print "Filtrando Imagen"
filterImage(hdu,params)
pf.append('filtered.fits',hdu,header)

hdulist.close()