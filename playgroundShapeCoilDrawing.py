import math as ma
import numpy as nu
from matplotlib import pyplot as pl
from matplotlib import cm
from mpl_toolkits import mplot3d
import csv
import pandas as pd

# Variables

preferredZ = 0.0


# Functions

def readDataOfPlane(elementAxis, expectedZValue):
    df = pd.read_csv(f'{elementAxis}ElementsOfB.csv', names=('x', 'y', 'z', 'B'))

    temp = df.copy()
    temp['z'] = abs(temp['z'] - expectedZValue)
    print(temp['z'])

    nearestZIndex = temp['z'].idxmin()
    print(f'Nearest Z index: {nearestZIndex}')

    nearestZValue = df.loc[nearestZIndex, 'z']
    print(f'Nearest Z Value: {nearestZValue}')

    return df[ df['z']==nearestZValue ]


# Main

xDataInPlane = readDataOfPlane(elementAxis='x', expectedZValue=preferredZ)
# yDataInPlane = readDataOfPlane(elementAxis='y', expectedZValue=preferredZ)
# zDataInPlane = readDataOfPlane(elementAxis='z', expectedZValue=preferredZ)

xAxisPoints = xDataInPlane['x']
yAxisPoints = xDataInPlane['y']

xElementsOfB = xDataInPlane['B']

fig = pl.figure()
ax = pl.axes(projection='3d')

ax.plot_surface(xAxisPoints, yAxisPoints, xElementsOfB, rstride=20, cstride=20, cmap='Reds')
pl.show()
