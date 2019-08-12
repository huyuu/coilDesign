import math as ma
import numpy as nu
from matplotlib import pyplot as pl
from matplotlib import cm
from mpl_toolkits import mplot3d
import csv
import pandas as pd

# Variables

I = 100.0
N = 500
h = 0.05
R = h

Y0Coeff = 0.1
Z0Coeff = 0.1
Y0 = Y0Coeff * h
Z0 = Z0Coeff * h

dsLowerCoeff = 0.2
dsHigherCoeff = 1.0
lsLowerCoeff = 1.0
lsHigherCoeff = 5.0

# y = 0.1  # y=0.1h
# x = 0.1  # x=0.1l
dirName = f'I={I}_N={N}_h={h*100}cm_Y0=Z0={Y0Coeff}h'


# Functions

def readData():
    maxX0s = pd.read_csv(f'{dirName}/maxX0s_D={dsLowerCoeff}hTo{dsHigherCoeff}h_L={lsLowerCoeff}hTo{lsHigherCoeff}h.csv', header=None)
    meanBs = pd.read_csv(f'{dirName}/meanBs_D={dsLowerCoeff}hTo{dsHigherCoeff}h_L={lsLowerCoeff}hTo{lsHigherCoeff}h.csv', header=None)
    return (maxX0s, meanBs)


def readSamplePoints():
    dSamples = pd.read_csv(f'{dirName}/dSamplesFrom{dsLowerCoeff}hTo{dsHigherCoeff}h.csv')
    lSamples = pd.read_csv(f'{dirName}/lSamplesFrom{lsLowerCoeff}hTo{lsHigherCoeff}h.csv')
    return (dSamples, lSamples)


# Main

maxX0s, meanBs = readData()
dSamples, lSamples = readSamplePoints()


fig = pl.figure()
ax = pl.axes(projection='3d')
ax.set_xlabel('d [h]', fontsize=22, labelpad=24)
ax.set_ylabel('l [h]', fontsize=22, labelpad=24)
ax.set_zlabel('maxX0s [h]', fontsize=22, labelpad=24)
ax.set_title(f'maxX0s_D={dsLowerCoeff}hTo{dsHigherCoeff}h_L={lsLowerCoeff}hTo{lsHigherCoeff}h', fontsize=28)
ax.tick_params(labelsize=22)

ax.plot_surface(dSamples.values, lSamples.values.reshape(1, len(lSamples)), maxX0s.values, rstride=20, cstride=20, cmap='Reds')
pl.show()


ax = pl.axes(projection='3d')
ax.set_xlabel('d [h]', fontsize=22, labelpad=24)
ax.set_ylabel('l [h]', fontsize=22, labelpad=24)
ax.set_zlabel('meanBs [h]', fontsize=22, labelpad=24)
ax.set_title(f'meanBs_D={dsLowerCoeff}hTo{dsHigherCoeff}h_L={lsLowerCoeff}hTo{lsHigherCoeff}h', fontsize=28)
ax.tick_params(labelsize=22)

ax.plot_surface(dSamples.values, lSamples.values.reshape(1, len(lSamples)), meanBs.values, rstride=20, cstride=20, cmap='Reds')
pl.show()
