import math as ma
import numpy as nu
from matplotlib import pyplot as pl
from matplotlib import cm
from mpl_toolkits import mplot3d
import csv
import pandas as pd

# Variables

data = {}
sampleFilePaths = ["dSamples.csv", "lSamples.csv"]
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

xData = data['dSamples.csv']
yData = data['lSamples.csv']
filePathsToPlot = ['meanBx.csv', 'meanBz.csv', 'varriationRateX.csv', 'varriationRateZ.csv']
zDatas = [ (filePath, data[filePath]) for filePath in filePathsToPlot ]

fig = pl.figure()
ax = pl.axes(projection='3d')
ax.set_xlabel('d [h]', fontsize=22, labelpad=24)
ax.set_ylabel('l [h]', fontsize=22, labelpad=24)
ax.set_zlabel('meanBx [mT]', fontsize=22, labelpad=24)
ax.tick_params(labelsize=22)
zData = data['meanBx.csv']
ax.plot_surface(xData.values, yData.values.reshape(1, len(yData)), zData.values*1000, cmap='Blues')
# xDataForScatter, yDataForScatter = nu.meshgrid(xData.values, yData.values)
# zDataForScatter = data['meanBx.csv'].values.ravel()*1000
# ax.scatter3D(xDataForScatter, yDataForScatter, zDataForScatter, c=zDataForScatter, cmap='viridis')
pl.show()

ax = pl.axes(projection='3d')
ax.set_xlabel('d [h]', fontsize=22, labelpad=24)
ax.set_ylabel('l [h]', fontsize=22, labelpad=24)
ax.set_zlabel('meanBz [mT]', fontsize=22, labelpad=24)
ax.tick_params(labelsize=22)
zData = data['meanBz.csv']
ax.plot_surface(xData.values, yData.values.reshape(1, len(yData)), zData.values*1000, cmap='Reds')
# xDataForScatter, yDataForScatter = nu.meshgrid(xData.values, yData.values)
# zDataForScatter = data['meanBz.csv'].values.ravel()*1000
# ax.scatter3D(xDataForScatter, yDataForScatter, zDataForScatter, c=zDataForScatter, cmap='viridis')
pl.show()

ax = pl.axes(projection='3d')
ax.set_xlabel('d [h]', fontsize=22, labelpad=24)
ax.set_ylabel('l [h]', fontsize=22, labelpad=24)
ax.set_zlabel('varRateX [%]', fontsize=22, labelpad=24)
ax.tick_params(labelsize=22)
zData = data['varriationRateX.csv']
ax.plot_surface(xData.values, yData.values.reshape(1, len(yData)), zData.values*100, cmap='Blues')
# xDataForScatter, yDataForScatter = nu.meshgrid(xData.values, yData.values)
# zDataForScatter = data['varriationRateX.csv'].values.ravel()*100
# ax.scatter3D(xDataForScatter, yDataForScatter, zDataForScatter, c=zDataForScatter, cmap='viridis')
pl.show()

ax = pl.axes(projection='3d')
ax.set_xlabel('d [h]', fontsize=22, labelpad=24)
ax.set_ylabel('l [h]', fontsize=22, labelpad=24)
ax.set_zlabel('varRateZ [%]', fontsize=22, labelpad=24)
ax.tick_params(labelsize=22)
zData = data['varriationRateZ.csv']
ax.plot_surface(xData.values, yData.values.reshape(1, len(yData)), zData.values*100, cmap='Reds')
# xDataForScatter, yDataForScatter = nu.meshgrid(xData.values, yData.values)
# zDataForScatter = data['varriationRateZ.csv'].values.ravel()*100
# ax.scatter3D(xDataForScatter, yDataForScatter, zDataForScatter, c=zDataForScatter, cmap='viridis')
pl.show()
