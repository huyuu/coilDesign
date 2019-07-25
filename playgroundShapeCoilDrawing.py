import math as ma
import numpy as nu
from matplotlib import pyplot as pl
from matplotlib import cm
from mpl_toolkits import mplot3d
import csv
import pandas as pd

# Variables

preferredZ = 0.1
I = 100.0
N = 500
d = 1.0  # d=1.0h
y = 0.1  # y=0.1h
x = 0.1  # x=0.1l
dirName = f'I={I}_N={N}_d={d}h_y={y}h_x={x}l'


# Functions

def readDataOfPlane(elementAxis, planeZValue):
    xPoints = pd.read_csv(f'{dirName}/xSamplePointsAtZ={planeZValue}.csv', header=None)
    yPoints = pd.read_csv(f'{dirName}/ySamplePointsAtZ={planeZValue}.csv', header=None)
    df = pd.read_csv(f'{dirName}/{elementAxis}ElementsOfBAtZ={planeZValue}.csv', header=None)
    return df.dropna(axis=1)


def readSamplePointsOf(var, planeZValue):
    return pd.read_csv(f'{dirName}/{var}SamplePointsAtZ={planeZValue}.csv', header=None)


# Main

xElementsOfB = abs(readDataOfPlane(elementAxis='x', planeZValue=preferredZ))
yElementsOfB = abs(readDataOfPlane(elementAxis='y', planeZValue=preferredZ))
zElementsOfB = readDataOfPlane(elementAxis='z', planeZValue=preferredZ) * 1e3

xSamplePoints = readSamplePointsOf(var='x', planeZValue=preferredZ) * 1e2
ySamplePoints = readSamplePointsOf(var='y', planeZValue=preferredZ) * 1e2


fig = pl.figure()
ax = pl.axes(projection='3d')
ax.set_xlabel('X Distance [cm]', fontsize=22, labelpad=25)
ax.set_ylabel('Y Distance [cm]', fontsize=22, labelpad=25)
ax.set_zlabel('Magnetic Field B [mT]', fontsize=22, labelpad=25)
ax.set_title(f'At z = {preferredZ}d', fontsize=26)
ax.tick_params(labelsize=16)

ax.plot_surface(xSamplePoints.values, ySamplePoints.values.reshape(1, len(ySamplePoints)), zElementsOfB.values, rstride=20, cstride=20, cmap='Reds')
pl.show()
