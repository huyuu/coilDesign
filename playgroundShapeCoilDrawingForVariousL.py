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

lower = 0.1
upper = 2.0
dirName = f'I={I}_N={N}_z={z}h_y={y}h_x={x}l'
fileName = f"variantionRateUnderVariousDFrom{lower}To{upper}"


# Functions

def readVariationRate(fileName):
    return pd.read_csv(f'{dirName}/{fileName}.csv')


# Main

data = readVariationRate(fileName)
distanceArray = data["d(h)"].values
variationRateArray = data["variationRate"].values
meanBArray = data['meanB[mT]'].values

pl.plot(distanceArray, variationRateArray, 'o', markersize=8)
pl.xlabel('Distance [h]', fontsize=40)
pl.ylabel('Variation Rate [%]', fontsize=40)
pl.title('Magnetic Field Variation Rate within (1cm)^2 x length Cube', fontsize=44)
pl.tick_params(labelsize=28)
pl.show()

pl.plot(distanceArray, meanBArray, 'o', markersize=8)
pl.xlabel('Distance [h]', fontsize=40)
pl.ylabel('Average B [mT]', fontsize=40)
pl.title('Average Magnetic Field B within (1cm)^2 x length Cube', fontsize=44)
pl.tick_params(labelsize=28)
pl.show()
