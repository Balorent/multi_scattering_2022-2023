#                           COMPLEXTOOLS.PY
# ------------------------------------------------------------------------
# This file contains functions found on the GitHub repository of the user
# empet : https://github.com/empet/Math/blob/master/DomainColoring.ipynb
# ------------------------------------------------------------------------


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import hsv_to_rgb
import math
import view


def Hcomplex(z):# computes the hue corresponding to the complex number z
    H = np.angle(z) / (2*np.pi) + 1
    return np.mod(H, 1)


def func_vals(f, re, im, N):  # evaluates the complex function at the nodes of the grid
    # re and im are  tuples, re=(a, b) and im=(c, d), defining the rectangular region
    # N is the number of discrete points per unit interval

    l = re[1] - re[0]
    h = im[1] - im[0]
    resL = N * l  # horizontal resolution
    resH = N * h  # vertical resolution
    x = np.linspace(re[0], re[1], int(resL))
    y = np.linspace(im[0], im[1], int(resH))
    x, y = np.meshgrid(x, y)
    z = x + 1j * y
    return f(z)


def domaincol_c(w, s):  # Classical domain coloring
    # w is the  array of values f(z)
    # s is the constant saturation

    H = Hcomplex(w)
    S = s * np.ones(H.shape)
    modul = np.absolute(w)
    g = 0.04
    ctr = 0.1
    if view.scale_type.get() == 0:
        V = 1/math.pi * np.arctan2(g/2, (ctr - modul))
    elif view.scale_type.get() == 1:
        V = view.xy_plot.scale_min + modul * (view.xy_plot.scale_max - view.xy_plot.scale_min) / modul.max()
    elif view.scale_type.get() == 2:
        V = 1 / math.pi * np.arctan2(g / 2, (ctr - modul))
    # V = (1.0 - 1.0 / (1 + modul ** 2)) ** 0.2
    # the points mapped to infinity are colored with white; hsv_to_rgb(0, 0, 1)=(1, 1, 1)=white

    HSV = np.dstack((H, S, V))
    RGB = hsv_to_rgb(HSV)
    return RGB


def plot_domain(color_func, f, re=[-1, 1], im=[-1, 1], title='',
                s=0.9, N=200, daxis=None):
    w = func_vals(f, re, im, N)
    domc = color_func(w, s)
    plt.xlabel("$\Re(z)$")
    plt.ylabel("$\Im(z)$")
    plt.title(title)
    if (daxis):
        plt.imshow(domc, origin="lower", extent=[re[0], re[1], im[0], im[1]])

    else:
        plt.imshow(domc, origin="lower")
        plt.axis('off')