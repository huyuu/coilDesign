import math as ma
import numpy as nu
from matplotlib import pyplot as pl
from matplotlib import cm
from mpl_toolkits import mplot3d
import csv
import pandas as pd

# Variables

I = 3.0
N = 400
h = 0.05
X0Coeff = 0.1
Y0Coeff = 0.1
Z0Coeff = 0.1
thickness = 2.0  # cm
conductorPhi = 1.12 # mm
dirName = f'precise_I={I}_N={N}_h={h*100}cm_X0={X0Coeff}h_Y0={Y0Coeff}h_Z0={Z0Coeff}h_thickness={thickness}cm_conductorPhi={conductorPhi}mm'


# Functions

def readData():
    resultsInPlane = pd.read_csv(f'{dirName}/resultsInZ0Plane.csv', header=None)
    result = pd.read_csv(f'{dirName}/result.csv', header=None)
    return (resultsInPlane, result)


def readSamplePoints():
    xSamples = pd.read_csv(f'{dirName}/xSamples.csv')
    ySamples = pd.read_csv(f'{dirName}/ySamples.csv')
    return (xSamples, ySamples)


# Main

resultsInPlane, result = readData()
xSamples, ySamples = readSamplePoints()


fig = pl.figure()
ax = pl.axes(projection='3d')
ax.set_xlabel('x [m]', fontsize=22, labelpad=24)
ax.set_ylabel('y [m]', fontsize=22, labelpad=24)
ax.set_zlabel('B [mT]', fontsize=22, labelpad=24)
ax.set_title(f'h = {h*100}cm, l = 2h, d = 0.5h', fontsize=28)
ax.tick_params(labelsize=22)

ax.plot_surface(xSamples.values, ySamples.values.reshape(1, len(ySamples)), (resultsInPlane*1000).values, rstride=20, cstride=20, cmap='Reds')
pl.show()


# ax = pl.axes(projection='3d')
# ax.set_xlabel('d [h]', fontsize=22, labelpad=24)
# ax.set_ylabel('l [h]', fontsize=22, labelpad=24)
# ax.set_zlabel('meanBs [h]', fontsize=22, labelpad=24)
# ax.set_title(f'meanBs_D={dsLowerCoeff}hTo{dsHigherCoeff}h_L={lsLowerCoeff}hTo{lsHigherCoeff}h', fontsize=28)
# ax.tick_params(labelsize=22)
#
# ax.plot_surface(dSamples.values, lSamples.values.reshape(1, len(lSamples)), meanBs.values, rstride=20, cstride=20, cmap='Reds')
# pl.show()
