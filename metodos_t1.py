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
	row =1/(cd1_2*cd2_1-cd1_1*cd2_2)*(cd2_1*(RA-crval1)-cd1_1*(DEC-crval2))+crpix2
	col =1/(cd1_2*cd2_1-cd1_1*cd2_2)*(-cd2_2*(RA-crval1)+cd1_2*(DEC-crval2))+crpix1
	return (int(row),int(col))#revisar bien esto despues

def addStar (hdu, m, RA, DEC):
	(ROW,COL)=RADECtoRowCol(RA,DEC)
	if 0<=ROW<4096 and 0<=COL<4096:	
		hdu[ROW,COL]+=mToCounts(m,20,flux20)
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

def psersic (Re,n,m,xc,yc,x,y,el,theta):
	bn=2*n-0.324
	ln=mToCounts(m,20,flux20)
	I0=ln*bn**2/((Re**2)*2*math.pi*n*math.gamma(2*n))
	E=math.sqrt(((x-xc)*math.cos(theta)+(y-yc)*math.sin(theta))**2+((x-xc)*math.sin(theta)-(y-yc)*math.cos(theta))**2/(1-el)**2)
	I=I0*math.exp(-bn*(E/Re)**(1/n))
	return (I)
	
def addGalaxy (hdu, m, RA, DEC, n, Re, el, theta):
	(ROW,COL)=RADECtoRowCol(RA,DEC)
	if 0<=ROW<4096 and 0<=COL<4096:	
		a1=int(ROW-4*Re)
		b1=int(ROW+4*Re)
		a2=int(COL-4*Re)
		b2=int(COL+4*Re)
		for (y) in range(a1,b1):
			for(x) in range(a2,b2):
				if  0<=y<4096 and 0<=x<4096:
					hdu[y,x]+=psersic(Re,n,m,COL,ROW,x,y,el,theta)
		#hdu[ROW,COL]+=mToCounts(m,20,flux20)*(2*n-0.324)**2/((Re**2)*2*math.pi*n*math.gamma(2*n))
		print ROW,COL, mToCounts(m,20,flux20)*(2*n-0.324)**2/((Re**2)*2*math.pi*n*math.gamma(2*n))
	return

def addGalaxyCatalog (hdu, catalog):
	i=0
	for linea in open(catalog):
		linea = linea.strip()
		obj, ra, dec,mag,sed,redshift,tipo,n,re,elip,o = linea.split()
		ra=float(ra)
		dec=float(dec)
		mag=float(mag)
		n=float(n)
		re=float(re)
		elip=float(elip)
		o=float(o)
		addGalaxy(hdu,mag,ra,dec,n,re,elip,o)
		i+=1
		if i>100000: break
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
