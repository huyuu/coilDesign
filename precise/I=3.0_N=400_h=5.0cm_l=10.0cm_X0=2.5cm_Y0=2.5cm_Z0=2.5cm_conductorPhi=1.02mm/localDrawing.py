import math as ma
import numpy as nu
from matplotlib import pyplot as pl
from matplotlib import cm
from mpl_toolkits import mplot3d
import csv
import pandas as pd

# Variables

data = {}
sampleFilePaths = ["xSamples.csv", "ySamples.csv"]
resultFilePaths = ['topPlaneBx.csv', 'topPlaneBy.csv', 'topPlaneBz.csv', 'zeroPlaneBx.csv', 'zeroPlaneBy.csv', 'zeroPlaneBz.csv']


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

xData = data['xSamples.csv']
yData = data['ySamples.csv']
zDatas = [ (filePath, data[filePath]) for filePath in resultFilePaths ]

for (name, zData) in zDatas[0:3]:
    fig = pl.figure()
    ax = pl.axes(projection='3d')
    ax.set_xlabel('x [cm]', fontsize=22, labelpad=24)
    ax.set_ylabel('y [cm]', fontsize=22, labelpad=24)
    name = name.split('.')[0]
    ax.set_zlabel(f'{name} [mT]', fontsize=22, labelpad=24)
    ax.tick_params(labelsize=22)
    ax.plot_surface(xData.values*100, yData.values.reshape(1, len(yData))*100, zData.values*1000, cmap='Reds')
    # xDataForScatter, yDataForScatter = nu.meshgrid(xData.values, yData.values)
    # zDataForScatter = data['meanBx.csv'].values.ravel()*1000
    # ax.scatter3D(xDataForScatter, yDataForScatter, zDataForScatter, c=zDataForScatter, cmap='viridis')
    pl.show()

for (name, zData) in zDatas[3:]:
    fig = pl.figure()
    ax = pl.axes(projection='3d')
    ax.set_xlabel('x [cm]', fontsize=22, labelpad=24)
    ax.set_ylabel('y [cm]', fontsize=22, labelpad=24)
    name = name.split('.')[0]
    ax.set_zlabel(f'{name} [mT]', fontsize=22, labelpad=24)
    ax.tick_params(labelsize=22)
    ax.plot_surface(xData.values*100, yData.values.reshape(1, len(yData))*100, zData.values*1000, cmap='Blues')
    # xDataForScatter, yDataForScatter = nu.meshgrid(xData.values, yData.values)
    # zDataForScatter = data['meanBx.csv'].values.ravel()*1000
    # ax.scatter3D(xDataForScatter, yDataForScatter, zDataForScatter, c=zDataForScatter, cmap='viridis')
    pl.show()
