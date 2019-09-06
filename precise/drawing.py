import math as ma
import numpy as nu
from matplotlib import pyplot as pl
from matplotlib import cm
from mpl_toolkits import mplot3d
import csv
import pandas as pd

# Variables

data = {}
sampleFilePaths =
resultFilePaths =


# Functions

def readData(filePath):
    return pd.read_csv(f'{filePath}', header=None)


def readSamplePoints(filePath):
    return pd.read_csv(f'{filePath}')


# Main

for filePath in sampleFilePaths:
    data[f'{filePath}'] = readSamplePoints(filePath)
for filePath in resultFilePaths:
    data[f'{filePath}'] = readData(filePath)

xData = data
yData = data
filePathsToPlot =
zDatas = [ (filePath, data[filePath]) for filePath in filePathsToPlot ]

fig = pl.figure()
for (name, zData) in zDatas:
    ax = pl.axes(projection='3d')
    ax.set_xlabel('x [m]', fontsize=22, labelpad=24)
    ax.set_ylabel('y [m]', fontsize=22, labelpad=24)
    name = name.split('.')[0]
    ax.set_zlabel(f'{name}', fontsize=22, labelpad=24)
    ax.tick_params(labelsize=22)
    # for surface plot
    ax.plot_surface(xData.values, yData.values.reshape(1, len(yData)), zData.values, rstride=20, cstride=20, cmap='Reds')
    # for scatter plot
    xDataForScatter, yDataForScatter = nu.meshgrid(xData.values*100, yData.values*100)
    ax.scatter3D(xDataForScatter.ravel(), yDataForScatter.ravel(), zData.values.ravel(), c=zDataForScatter, cmap='Reds')
    pl.show()
