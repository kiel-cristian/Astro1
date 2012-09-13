#!/usr/local/bin/python
# -*- coding: utf-8 -*-
import pyfits

hdulist = pyfits.open('blank.fits')
hdulist.close()
crpix1=hdulist[0].header['CRPIX1']
crpix2=hdulist[0].header['CRPIX2']
crval1=hdulist[0].header['CRVAL1']
crval2=hdulist[0].header['CRVAL2']
cd1_1=hdulist[0].header['CD1_1']
cd1_2=hdulist[0].header['CD1_2']
cd2_1=hdulist[0].header['CD2_1']
cd2_2=hdulist[0].header['CD2_2']

def RADECtoXY(RA,DEC):
    row =(cd2_1*(crval1-RA)-cd1_1*(crval2-DEC))/(cd1_2*cd2_1-cd1_1*cd2_2)+crpix2
    col =(-cd2_2*(crval1-RA)+cd1_2*(crval2-DEC))/(cd1_2*cd2_1-cd1_1*cd2_2)+crpix1
    return (row,col)


ra = 22.1142135471
dec = -0.112283531445

print "RA: "+ str(ra) +"\t\t\tDEC: "+str(dec)
(x,y) = RADECtoXY(ra,dec)
print "X: "+ str(x) +"\t\t\tY: "+str(y)