#!/usr/local/bin/python
# -*- coding: utf-8 -*-
import pyfits as pf
from fits_lib import *

def plot_image_fit(name):
	hdulist = pf.open(name)
	hdu= hdulist[0].data
	plot_image(hdu,log_scale=True)
	hdulist.close()


if __name__ == "__main__":
	plot_image_fit('stars.fits')
	plot_image_fit('galaxys.fits')
	plot_image_fit('background.fits')
	plot_image_fit('convolve.fits')
	plot_image_fit('noise.fits')
	plot_image_fit('filtered.fits')