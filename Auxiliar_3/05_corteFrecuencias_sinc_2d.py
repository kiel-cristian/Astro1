import pylab as pl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import numpy.fft as fft
from scipy import random

def plot2D (x, y, z, filename = False):

    pl.clf()
    fig  = pl.figure()
    ax   = fig.gca(projection='3d')
    surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=False)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)

    if filename:
        pl.savefig (filename)
    else:
        pl.show()

N  = 128
x  = np.zeros((N,N)) + np.arange(N)
y  = x.transpose()
fx = np.zeros((N,N)) + fft.fftfreq(N)
fy = fx.transpose()

# Corte de frecuencias: funcion rectangular
R = N/4
Z = np.zeros((N,N))
for i in range(N):
    for j in range(N):
	if (i - N/2)**2 + (j - N/2)**2 <= R**2:
	   Z[i][j] = 1
Z = fft.ifftshift(Z)
plot2D(fx,fy,Z,"f_rectangulo")

z = (fft.fftshift(fft.ifft2(Z))).real
zoom1 = 3*N/8
zoom2 = 5*N/8
plot2D (x[zoom1:zoom2,zoom1:zoom2], y[zoom1:zoom2,zoom1:zoom2],z[zoom1:zoom2,zoom1:zoom2], "sinc_zoom")

# Gaussiana + ruido + recorte en frecuencias
sigma_x = 5.
sigma_y = 5.
s = 0.01/(2*np.pi*sigma_x*sigma_y)
x_zero = N/2
y_zero = N/2

gaussian = np.exp(-((x-x_zero)**2.0/(2*sigma_x**2.0)+(y-y_zero)**2.0/(2*sigma_y**2.0)))/(2*np.pi*sigma_x*sigma_y)
z = gaussian + random.standard_normal((N,N)) * s

plot2D (x, y, z, "Gaussiana_ruido")

Z = fft.fft2(z)
plot2D (fx, fy, np.abs(Z), "f_Gaussiana_ruido")

N_corte  = 10
for i in range(N):
    for j in range(N):
	if ((i+N/2)%N - N/2)**2 + ((j+N/2)%N - N/2)**2 > (N_corte)**2:
	   Z[i][j] = 0.0

plot2D (fx, fy, np.abs(Z), "f_Gaussiana_ruido_recorte")
z1 = fft.ifft2 (Z).real
plot2D (x, y, z1, "Gaussiana_ruido_recorte")
plot2D (x, y, z1 * (z1 < 0.001) + 0.001 * (z1 >= 0.001), "corte_Gaussiana_ruido_recorte")
plot2D (x, y, gaussian, "Gaussiana_original")
