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
resultFilePaths = ['meanBx.csv', 'meanBz.csv', 'maxBx.csv', 'maxBz.csv', 'minBx.csv', 'minBz.csv', 'varriationRateX.csv', 'varriationRateZ.csv']


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
filePathsToPlot = ['meanBx.csv', 'meanBz.csv', 'maxBx.csv', 'maxBz.csv', 'minBx.csv', 'minBz.csv', 'varriationRateX.csv', 'varriationRateZ.csv']
zDatas = [ (filePath, data[filePath]) for filePath in filePathsToPlot ]

fig = pl.figure()
ax = fig.add_subplot(2, 2, 1, projection='3d')
ax.set_xlabel('h [cm]', fontsize=22, labelpad=24)
ax.set_ylabel('l [cm]', fontsize=22, labelpad=24)
ax.set_zlabel('meanBx [mT]', fontsize=22, labelpad=24)
ax.tick_params(labelsize=22)
zData = data['meanBx.csv']
ax.plot_surface(xData.values*100, yData.values.reshape(1, len(yData))*100, zData.values*1000, cmap='Blues')
# xDataForScatter, yDataForScatter = nu.meshgrid(xData.values*100, yData.values*100)
# zDataForScatter = data['meanBz.csv'].values.ravel()*1000
# ax.scatter3D(xDataForScatter.ravel(), yDataForScatter.ravel(), zDataForScatter, c=zDataForScatter, cmap='Reds')

ax = fig.add_subplot(2, 2, 2, projection='3d')
ax.set_xlabel('h [cm]', fontsize=22, labelpad=24)
ax.set_ylabel('l [cm]', fontsize=22, labelpad=24)
ax.set_zlabel('varRateX [%]', fontsize=22, labelpad=24)
ax.tick_params(labelsize=22)
zData = data['varriationRateX.csv']
ax.plot_surface(xData.values*100, yData.values.reshape(1, len(yData))*100, zData.values*100, cmap='Blues')
# xDataForScatter, yDataForScatter = nu.meshgrid(xData.values*100, yData.values*100)
# zDataForScatter = data['varriationRateZ.csv'].values.ravel()*100
# ax.scatter3D(xDataForScatter.ravel(), yDataForScatter.ravel(), zDataForScatter, c=zDataForScatter, cmap='Reds')

ax = fig.add_subplot(2, 2, 3, projection='3d')
ax.set_xlabel('h [cm]', fontsize=22, labelpad=24)
ax.set_ylabel('l [cm]', fontsize=22, labelpad=24)
ax.set_zlabel('meanBz [mT]', fontsize=22, labelpad=24)
ax.tick_params(labelsize=22)
zData = data['meanBz.csv']
ax.plot_surface(xData.values*100, yData.values.reshape(1, len(yData))*100, zData.values*1000, cmap='Reds')
# xDataForScatter, yDataForScatter = nu.meshgrid(xData.values*100, yData.values*100)
# zDataForScatter = data['varriationRateZ.csv'].values.ravel()*100
# ax.scatter3D(xDataForScatter.ravel(), yDataForScatter.ravel(), zDataForScatter, c=zDataForScatter, cmap='Reds')

ax = fig.add_subplot(2, 2, 4, projection='3d')
ax.set_xlabel('h [cm]', fontsize=22, labelpad=24)
ax.set_ylabel('l [cm]', fontsize=22, labelpad=24)
ax.set_zlabel('varRateZ [%]', fontsize=22, labelpad=24)
ax.tick_params(labelsize=22)
zData = data['varriationRateZ.csv']
ax.plot_surface(xData.values*100, yData.values.reshape(1, len(yData))*100, zData.values*100, cmap='Reds')
# xDataForScatter, yDataForScatter = nu.meshgrid(xData.values*100, yData.values*100)
# zDataForScatter = data['varriationRateZ.csv'].values.ravel()*100
# ax.scatter3D(xDataForScatter.ravel(), yDataForScatter.ravel(), zDataForScatter, c=zDataForScatter, cmap='Reds')
pl.show()


# ax = fig.add_subplot(2, 2, 1)
# ax.set_xlabel('l [cm]', fontsize=16, labelpad=16)
# ax.set_ylabel('Variation Rate Z [%]', fontsize=16, labelpad=16)
# ax.tick_params(labelsize=14)
# ax.set_title('h = 5cm', fontsize=18)
# zData = data['varriationRateZ.csv']
# varRateZ = zData.iloc[0, :]
# pl.plot(yData*100, varRateZ*100)
#
# ax = fig.add_subplot(2, 2, 2)
# ax.set_xlabel('l [cm]', fontsize=16, labelpad=16)
# ax.set_ylabel('Mean Bz [mT]', fontsize=16, labelpad=16)
# ax.set_title('h = 5cm', fontsize=18)
# ax.tick_params(labelsize=14)
# zData = data['meanBz.csv']
# meanBz = zData.iloc[0, :]
# pl.plot(yData*100, meanBz*1000)
#
# ax = fig.add_subplot(2, 2, 3)
# ax.set_xlabel('l [cm]', fontsize=16, labelpad=16)
# ax.set_ylabel('Variation Rate Z [%]', fontsize=16, labelpad=16)
# ax.set_title('h = 6cm', fontsize=18)
# ax.tick_params(labelsize=14)
# zData = data['varriationRateZ.csv']
# varRateZ = zData.iloc[11, :]
# pl.plot(yData*100, varRateZ*100)
#
# ax = fig.add_subplot(2, 2, 4)
# ax.set_xlabel('l [cm]', fontsize=16, labelpad=16)
# ax.set_ylabel('Mean Bz [mT]', fontsize=16, labelpad=16)
# ax.set_title('h = 6cm', fontsize=18)
# ax.tick_params(labelsize=14)
# zData = data['meanBz.csv']
# meanBz = zData.iloc[11, :]
# pl.plot(yData*100, meanBz*1000)
# pl.show()
