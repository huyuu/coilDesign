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
z = 0.1
y = 0.1  # y=0.1h
x = 0.1  # x=0.1l

lower = 0.56
upper = 0.62
dirName = f'I={I}_N={N}_z={z}h_y={y}h_x={x}l'
fileName = f"variantionRateUnderVariousDFrom{lower}To{upper}.csv"


# Functions

def readVariationRate(fileName):
    return pd.read_csv(f'{dirName}/{fileName}')


# Main

data = readVariationRate(fileName)
distanceArray = data["d(h)"].values
variationRateArray = data["variationRate"].values

pl.plot(distanceArray, variationRateArray)
pl.xlabel('Distance [h]')
pl.ylabel('Variation Rate [%]')
pl.show()
# fig = pl.figure()
# ax = pl.axes(projection='3d')
# ax.set_xlabel('X Distance [cm]', fontsize=22, labelpad=25)
# ax.set_ylabel('Y Distance [cm]', fontsize=22, labelpad=25)
# ax.set_zlabel('Magnetic Field B [mT]', fontsize=22, labelpad=25)
# ax.set_title(f'At z = {preferredZ}d; d = {d}h', fontsize=26)
# ax.tick_params(labelsize=16)
#
# ax.plot_surface(xSamplePoints.values, ySamplePoints.values.reshape(1, len(ySamplePoints)), zElementsOfB.values, rstride=20, cstride=20, cmap='Reds')
# pl.show()
