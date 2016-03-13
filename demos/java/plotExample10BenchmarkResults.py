#!/usr/bin/python

"""
plotExample10BenchmarkResults.py: Show 3D plots with the result of benchmarking
the CPU and GPU implementations of the KSG algorithm for mutual information.
"""

__author__   = "Pedro AM Mediano"
__email__    = "pmediano@imperial.ac.uk"
__license__  = "GPL"

import sys
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D

## Load files
filename1 = sys.argv[1]
filename2 = sys.argv[2]
res1 = np.loadtxt(filename1)
res2 = np.loadtxt(filename2)

## There is a bug with the legends of 3D scatterplots, so we have to save
# the plot style for the legend separately
img1 = pl.Line2D([0],[0], linestyle="none", c='r', marker='o')
img2 = pl.Line2D([0],[0], linestyle="none", c='b', marker='^')

## Plot times for first benchmark
fig1 = pl.figure()
ax1 = fig1.add_subplot(111, projection='3d')
ax1.scatter(np.log10(res1[:,0]), res1[:,1], res1[:,2], c='r', marker='o')
ax1.scatter(np.log10(res1[:,0]), res1[:,1], res1[:,3], c='b', marker='^')
ax1.legend([img1, img2], ['CPU', 'GPU'], numpoints = 1)
ax1.set_xlabel(r'$\log_{10} N$'); ax1.set_ylabel('D'); ax1.set_zlabel('Time [ms]');
ax1.set_title(filename1.partition(".")[0])

## Plot times for second benchmark
fig2 = pl.figure()
ax2 = fig2.add_subplot(111, projection='3d')
ax2.scatter(np.log10(res2[:,0]), res2[:,1], res2[:,2], c='r', marker='o')
ax2.scatter(np.log10(res2[:,0]), res2[:,1], res2[:,3], c='b', marker='^')
ax2.legend([img1, img2], ['CPU', 'GPU'], numpoints = 1)
ax2.set_xlabel(r'$\log_{10} N$'); ax2.set_ylabel(r'$\rho$'); ax2.set_zlabel('Time [ms]');
ax2.set_title(filename2.partition(".")[0])

pl.show()

