import pylab as pl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

N = 200
xmax = 10.

# 1D ---------
x  = np.arange(0, xmax, xmax/N)
y1 = np.exp(-x) * np.sin(x*2*np.pi)
y2 = np.exp(-x)

pl.clf ()
pl.plot (x, y1, "-b", label = r'$e^{-x}\sin(2\pi x)$')
pl.plot (x, y2, "-r", label = r'$e^{-x}$')
pl.xlabel("x")
pl.ylabel("y")
pl.legend(loc = 4)
pl.savefig ("1D")

# 2D ---------
def plot2D (x, y, z, filename = False):

    pl.clf()
    #pl.close()
    fig  = pl.figure()
    ax   = fig.gca(projection='3d')
    surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.jet,linewidth=0, antialiased=False)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)

    if filename:
        pl.savefig (filename)
    else:
        pl.show()
    

N = 100
xmax = 2.
x = np.zeros((N, N)) + np.arange (-xmax, xmax, 2*xmax/N)
y = x.transpose()
r = np.sqrt(x**2 + y**2)

plot2D (x, y, r, "radius")
z = np.exp(-r)*np.sin(r*2*np.pi)
plot2D (x, y, z, "function")
plot2D (x, y, z)
