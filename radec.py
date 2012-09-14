#!/usr/local/bin/python
# -*- coding: utf-8 -*-
import pyfits as pf

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

if __name__=="__main__":
	ra = 22.1142135471
	dec = -0.112283531445
	print "RA: "+ str(ra) +"\t\t\tDEC: "+str(dec)
	(x,y) = RADECtoRowCol(ra,dec)
	print "X: "+ str(x) +"\t\t\tY: "+str(y)