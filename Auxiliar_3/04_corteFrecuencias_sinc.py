import numpy as np
import numpy.fft as fft
import pylab as pl
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

# Corte de frecuencias: funcion rectangular
Y = np.zeros(N)
Y[N/4:3*N/4] = 1
Y = fft.ifftshift(Y)
plot (f, Y, "f_rectangulo")

y = (fft.fftshift(fft.ifft(Y))).real
plot (x, y, "sinc")
plot (x[3*N/8:5*N/8], y[3*N/8:5*N/8], "sinc_zoom")

# Gaussiana + ruido + recorte en frecuencias
sigma = 50.
s = 0.1 / np.sqrt(2*np.pi*sigma**2)
y = np.exp (-(x-N/2)**2/(2*sigma**2))/np.sqrt(2*np.pi*sigma**2)
y += random.standard_normal(N) * s

plot (x, y, "Gaussiana_ruido")

Y = abs(fft.fft(y))
plot (f, Y, "f_Gaussiana_ruido")

N_corte = 10.
Y[N_corte : N - N_corte] = 0.
plot (f, Y, "f_Gaussiana_ruido_recorte")

y1 = fft.fftshift(fft.ifft (Y)).real
plot (x, y1, "Gaussiana_ruido_recorte")
plot (x, np.exp (-(x-N/2)**2/(2*sigma**2))/np.sqrt(2*np.pi*sigma**2), "Gaussiana_ruido_recorte", clear = False)
