import math as ma
import numpy as nu
from matplotlib import pyplot as pl
from matplotlib import cm
from mpl_toolkits import mplot3d
import csv
import pandas as pd

# Variables

data = {}
sampleFilePaths = ["hSamples.csv", "lSamples.csv"]
resultFilePaths = ['meanBx.csv', 'meanBy.csv', 'meanBz.csv', 'varriationRateX.csv', 'varriationRateY.csv', 'varriationRateZ.csv']


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

xData = data['hSamples.csv']
yData = data['lSamples.csv']
filePathsToPlot = ['meanBz.csv', 'varriationRateZ.csv']
zDatas = [ (filePath, data[filePath]) for filePath in filePathsToPlot ]

fig = pl.figure()
ax = pl.axes(projection='3d')
ax.set_xlabel('h [cm]', fontsize=22, labelpad=24)
ax.set_ylabel('l [cm]', fontsize=22, labelpad=24)
ax.set_zlabel('meanBz [mT]', fontsize=22, labelpad=24)
ax.tick_params(labelsize=22)
xDataForScatter, yDataForScatter = nu.meshgrid(xData.values*100, yData.values*100)
ax.scatter3D(xDataForScatter, yDataForScatter, data['meanBz.csv'].values*1000)
pl.show()

ax = pl.axes(projection='3d')
ax.set_xlabel('h [cm]', fontsize=22, labelpad=24)
ax.set_ylabel('l [cm]', fontsize=22, labelpad=24)
ax.set_zlabel('varRateZ [%]', fontsize=22, labelpad=24)
ax.tick_params(labelsize=22)
xDataForScatter, yDataForScatter = nu.meshgrid(xData.values*100, yData.values*100)
ax.scatter3D(xDataForScatter, yDataForScatter, data['varriationRateZ.csv'].values*100)
pl.show()
