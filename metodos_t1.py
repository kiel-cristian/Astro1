#!/usr/local/bin/python
# -*- coding: utf-8 -*-
import pyfits as pf
import pylab as pl
import scipy.signal as sig
import scipy as sp
import numpy as np
import math

def mToCounts (m, m0, F0):
	counts=exptime*F0*10**(-2/5*(m-m0))
	return (counts)

hdulist = pf.open('blank.fits')
hdulist.close()
flux20=hdulist[0].header['FLUX20']
exptime=hdulist[0].header['EXPTIME']
crpix1=hdulist[0].header['CRPIX1']
crpix2=hdulist[0].header['CRPIX2']
crval1=hdulist[0].header['CRVAL1']
crval2=hdulist[0].header['CRVAL2']
cd1_1=hdulist[0].header['CD1_1']
cd1_2=hdulist[0].header['CD1_2']
cd2_1=hdulist[0].header['CD2_1']
cd2_2=hdulist[0].header['CD2_2']

def RADECtoRowCol(RA,DEC):
	row =(cd2_1*(crval1-RA)-cd1_1*(crval2-DEC))/(cd1_2*cd2_1-cd1_1*cd2_2)+crpix2
	col =(-cd2_2*(crval1-RA)+cd1_2*(crval2-DEC))/(cd1_2*cd2_1-cd1_1*cd2_2)+crpix1
	return (row,col)#revisar bien esto despues

def addStar (hdu, m, RA, DEC):
	(ROW,COL)=RADECtoRowCol(RA,DEC)
	if 0<=ROW<=4096 and 0<=COL<=4096:	
		hdu[ROW,COL]=mToCounts(m,20,flux20)
		print ROW,COL,mToCounts(m,20,flux20)
	return

def addStellarCatalog(hdu, catalog):
	for linea in open(catalog):
		linea = linea.strip()
		obj, ra, dec,mag,sed,index,tipo = linea.split()
		ra=float(ra)
		dec=float(dec)
		mag=float(mag)
		#print ra,dec,RADECtoRowCol(ra,dec)
		addStar(hdu,mag,ra,dec)
	return

def E(xc,yc,x,y,el,theta):
	e=math.sqrt(((x-xc)*math.cos(theta)+(y-yc)*math.sin(theta))**2+((x-xc)*math.sin(theta)-(y-yc)*math.cos(theta))**2/(1-el)**2)
	return (e)

def psersic (e,Re,n,i0,m):
	bn=2*n-0.324
	ln=mToCounts(m,20,flux20)
	I0=ln*bn**2/((Re**2)*2*math.pi*n*math.gamma(2*n))
	I=I0*math.exp(-bn*(e/RE)**(1/n))
	return (I)
	
def addGalaxy (hdu, m, RA, DEC, n, Re, el, theta):
	(ROW,COL)=RADECtoRowCol(RA,DEC)
    # coordenada RA, DEC, representada por un perfil de Sersic
    # definido por n, Re, el y rotado en un angulo theta.

def addGalaxyCatalog (hdu, catalog):
	for linea in open(catalog):
		linea = linea.strip()
		obj, ra, dec,mag,sed,redshift,tipo,n,re,el,theta = linea.split()
		ra=float(ra)
		dec=float(dec)
		mag=float(mag)
		n=float(n)
		re=float(re)
		el=float(el)
		theta=float(theta)
		addGalaxy(hdu,mag,ra,dec,n,re.el,theta)
	return

def addBackground (hdu, background):
    # Agrega un background constante a hdu.
    return

def convolvePSF (hdu, sigma):
    # Convoluciona hdu con un PSF Gaussiano de desviaci�n estandard
    # sigma.
    return
    
def addNoise (hdu, sigma):
    # Agrega ruido Gaussiano de desviacion estandard sigma a hdu.
    return

def filterImage (hdu, params):
    # Filtra hdu utilizando los parametros params.
    return
