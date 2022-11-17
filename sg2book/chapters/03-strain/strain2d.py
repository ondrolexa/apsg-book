# Library of routines for 2D deformation visualization
#
# Ondrej Lexa 2017, 2022

import numpy as np
import matplotlib.pyplot as plt

def plot_defgrad(F, **kwargs):
    # default args
    title = kwargs.pop('title', 'F = ' + str(list(F)))
    xmin = kwargs.pop('xmin', -3)
    xmax = kwargs.pop('xmax', 3)
    xsteps = kwargs.pop('xsteps', 21)
    ymin = kwargs.pop('ymin', -2)
    ymax = kwargs.pop('ymax', 2)
    ysteps = kwargs.pop('ysteps', 17)
    xpoly = kwargs.pop('xpoly', [-1, -1, 1, 1, -1])
    ypoly = kwargs.pop('xpoly', [-1, 1, 1, -1, -1])
    if 'ax' in kwargs:
        ax = kwargs.pop('ax')
    else:
        fig, ax = plt.subplots(**kwargs)
    # go
    X, Y = np.meshgrid(np.linspace(xmin, xmax, xsteps),
                       np.linspace(ymin, ymax, ysteps))
    u, v = np.tensordot(F - np.eye(2), [X, Y], axes=1)
    # field
    ax.quiver(X, Y, u, v, angles='xy')
    # polygon
    xp, yp = np.dot(F, [xpoly, ypoly])
    ax.plot(xpoly, ypoly, 'r', xp, yp, 'g')
    # deformation ellipse
    theta = np.linspace(0, 2*np.pi, 180)
    xc, yc = np.cos(theta), np.sin(theta)
    xe, ye = np.dot(F, [xc,yc])
    ax.plot(xc, yc, 'r', xe, ye, 'g')
    # principal axes
    u, s, v = np.linalg.svd(F)
    ax.quiver(np.zeros(4), np.zeros(4),
               np.hstack((s*u[0], -s*u[0])), np.hstack((s*u[1], -s*u[1])),
               scale=1, units='xy')
    # title
    ax.set_title(title)
    ax.set_aspect(1)

