#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import pyfits
import pylab as pl

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

def get_fits_matrix(path,index=0):

    arch = pyfits.open(path)
    mtx  = arch[index].data
    return mtx

def get_fits_header(path,index=0):

    arch = pyfits.open(path)
    hdr  = arch[index].header
    resp = []
    for key in list(hdr):
        resp.append((key,hdr[key]))
    return resp

def print_header(list_header):

    for key,value in list_header:
        key_str = key + " "*(8-len(key))
        print ">> %s\t= %s" % (key_str,str(value))

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
