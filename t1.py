import pyfits as pf
import pylab as pl
import scipy.signal as sig
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import metodos_t1 as m

hdulist = pf.open('blank.fits')
hdu= hdulist[0].data
catalogo1="stellar.dat"
catalogo2="galaxy.dat"
m.addStellarCatalog(hdu,catalogo1)
m.addGalaxyCatalog(hdu,catalogo2)
print hdu [10,10]
m.addBackground (hdu, 1000)
print hdu[10,10]
def plot_image(image,interpolation="nearest",log_scale=False,title=None):
	if log_scale:
		min_value = np.min(image)		
		if min_value <= 0: image2 = image + np.abs(min_value) + 10.0
		else: image2 = image - min_value + 10.0
		plt.imshow(np.log(image2),interpolation=interpolation)
	else:
		plt.imshow(image,interpolation=interpolation)
	if title != None:
		plt.title(str(title))
	plt.show()
	return
m.convolvePSF (hdu, 5)

plot_image(hdu,log_scale=True)
hdulist.close()
