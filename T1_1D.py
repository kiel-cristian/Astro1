import pyfits as pf
import pylab as pl
import numpy as np
import scipy.signal
import numpy.random as random
import numpy.fft as fft

def addStar (data, I, i):
   data [i] += I

def addGalaxy (data, I0, i, n, Re):
    x = np.arange (len (data))
    bn = 2. * n - 0.324
    data += I0 * np.exp (-bn * (np.abs(x-i)/Re)**(1./n))

def addBackground (data, background):
    data += background

def convolvePSF (data, sigma):
    k = 4.
    h = np.ceil (k * sigma)
    x_f = np.arange (-h,h+1)
    filter =  np.exp(-(x_f*x_f)/(2.*sigma**2)) / np.sqrt(2.*np.pi*sigma**2)
    filter /= filter.sum()
    C = scipy.signal.convolve (data, filter, 'same')
    data[:] = C[:]

def addNoise (data, sigma):
    data += random.standard_normal(len(data)) * sigma
    
def filterImage (data, fft_cut):
    Y = fft.fft(data)
    plot (np.log(np.abs(Y)), "FFT", "06_FFT")

    N = len(data)
    Ymin = np.abs(Y).min()
    Y[fft_cut : N - fft_cut] = 0.
    Yabs = np.abs(Y)
    plot (np.log(Yabs * (Yabs > 0) + (Yabs <= 0)*Ymin*1e-1), "FFT", "06_FFT_cut")

    y1 = np.abs((fft.ifft (Y)))
    data [:] = y1[:]

def plot (x, title, filename, ymin = -1e16):
    pl.clf()
    pl.plot (x, "-")
    pl.ylim ([max([ymin, x.min()]), x.max()])
    pl.title(title)
    pl.savefig(filename)

N = 1000
data = np.zeros (N)

stars_i = np.array([100, 500, 750])
stars_I = np.array([100., 300., 200.])

galaxies_i = np.array([250, 850])
galaxies_I0 = np.array([50., 150.])
galaxies_n = np.array([1., 4.])
galaxies_Re = np.array([20., 50.])

backg = 1000
psf = 20.
noise = 2.

for i in range (len(stars_i)):
    addStar (data, stars_I[i], stars_i[i])

plot (data, "Stars", "01_stars")

for i in range (len(galaxies_i)):
    addGalaxy (data, galaxies_I0[i], galaxies_i[i], galaxies_n[i], galaxies_Re[i])

plot (data, "Stars and Galaxies", "02_stars_gals")

addBackground (data, backg)

plot (data, "Sources + Background", "03_backg")

convolvePSF (data, psf)
data[:] = (data * (data > backg) + backg * (data < backg))[:]

plot (data, "PSF Convolution", "04_PSF", ymin = backg)

addNoise (data, noise)

plot (data, "Noise added", "05_noise", ymin = backg - noise)

filterImage (data, 9)

plot (data, "Filtered", "06_denoised")
