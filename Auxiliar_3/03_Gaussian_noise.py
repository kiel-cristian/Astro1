import pylab as pl
import numpy as np
import numpy.fft as fft
from scipy import random

def plot(x, y, filename = False, clear = True):
    if clear:
        pl.clf()
    pl.plot (x, y, "-")
    if filename:
        pl.savefig (filename)
    else:
        pl.show()

N = 1024
x = np.arange(N)
f = fft.fftfreq (N)

# Gaussiana
sigma = 50.
y = np.exp (-(x-N/2)**2/(2*sigma**2))/np.sqrt(2*np.pi*sigma**2)
plot (x, y, "Gaussiana")

Y = abs(fft.fft(y))
plot (f, Y, "f_Gaussiana")

# Ruido
s = 0.1 / np.sqrt(2*np.pi*sigma**2)
y = random.standard_normal(N) * s

plot (x, y, "ruido")

Y = abs(fft.fft(y))
plot (f, Y, "f_ruido")


# Gaussiana + ruido
y = np.exp (-(x-N/2)**2/(2*sigma**2))/np.sqrt(2*np.pi*sigma**2)
y += random.standard_normal(N) * s

plot (x, y, "Gaussiana_ruido")

Y = abs(fft.fft(y))
plot (f, Y, "f_Gaussiana_ruido")
