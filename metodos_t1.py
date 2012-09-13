#!/usr/local/bin/python
# -*- coding: utf-8 -*-
import pyfits as pf
import pylab as pl
import scipy.signal as sig
import scipy as sp
import numpy as np
from math import pow

def mToCounts (m, m0, F0):
    return F0*pow(10,-2.5*(m-m0))

def addStar (hdu, m, RA, DEC):
    # Definir rutina que agrega una estrella de magnitud m en la
    # coordenada RA, DEC
    return

def addStellarCatalog (hdu, catalog):
    # Definir rutina que agregue un catalogo de estrellas a hdu.
    return

def addGalaxy (hdu, m, RA, DEC, n, Re, el, theta):
    # Definir rutina que agrega una galaxia de magnitud m en la
    # coordenada RA, DEC, representada por un perfil de Sersic
    # definido por n, Re, el y rotado en un ángulo theta.
    return

def addGalaxyCatalog (hdu, catalog):
    # Definir rutina que agregue un catalogo de galaxias a hdu.
    return

def addBackground (hdu, background):
    # Agrega un background constante a hdu.
    return

def convolvePSF (hdu, sigma):
    # Convoluciona hdu con un PSF Gaussiano de desviación estandard
    # sigma.
    return
    
def addNoise (hdu, sigma):
    # Agrega ruido Gaussiano de desviacion estandard sigma a hdu.
    return

def filterImage (hdu, params):
    # Filtra hdu utilizando los parametros params.
    return
