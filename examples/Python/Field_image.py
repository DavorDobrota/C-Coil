import c_coil.coil as c
import c_coil.tensor as ten

import math as m
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import gc
import os

# if using GPU acceleration on Windows, CUDA path needs to be added
# common path is added below, but you should check for differences (e.g. different CUDA version)
# os.add_dll_directory("C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.7/bin")

# image resolution
xDim = 3001
yDim = 3001

# rectangular domain dimensions
xSize = 0.24
ySize = 0.24

# at which point the image is centered
xOffset = 0.0
yOffset = 0.0
zPos = 0.0

# coil used for field computations, orientation set such that it lies in xOy plane
coil = c.Coil(0.03, 0.03, 0.12, 3600, 10)
coil.set_position_and_orientation(ten.Vector3(), m.pi / 2, 0.0)

# setting precision and number of threads used (number of CPU threads is advised)
coil.set_thread_count(8)
precision = c.PrecisionFactor(7.0)

coil.set_default_precision_CPU(precision)
coil.set_default_precision_GPU(precision)

# initialising array for specifying points
positions = ten.Vector3Array()

# initialising image generation parameters, using np arrays to reduce memory use
x = np.zeros((xDim,), dtype=np.float64)
y = np.zeros((yDim,), dtype=np.float64)
pixelValues = np.zeros((xDim, yDim), dtype=np.float64)

# generating point position arguments
for i in range(0, yDim):
    for j in range(0, xDim):
        xPos = xOffset + xSize / (xDim - 1) * (j - (xDim - 1) / 2.0)
        yPos = yOffset + ySize / (yDim - 1) * (i - (yDim - 1) / 2.0)

        positions.append(xPos, yPos, zPos)

# generating positions for image pixels
for i in range(0, xDim):
    x[i] = xSize / (xDim - 1) * (i - (xDim - 1) / 2.0)

for i in range(0, yDim):
    y[i] = ySize / (yDim - 1) * (i - (yDim - 1) / 2.0)

# calling the field compute function
fields = coil.compute_all_B_field_vectors(positions, c.CPU_MT).abs()

# deleting unnecessary data
del positions

# it was noticed that when using the GPU some specific points produce nan or inf
# applying this simple "filter" resolved the issue in many cases
for i in range(0, yDim):
    for j in range(0, xDim):
        temp = fields[i * xDim + j]
        if m.isinf(temp) or m.isnan(temp):
            print(i, j)
            temp = 0.0
            n = 0
            if i - 1 >= 0:
                inter = fields[(i - 1) * xDim + j]
                if not m.isinf(inter) and not m.isnan(inter):
                    temp += inter
                    n += 1
            if i + 1 < yDim:
                inter = fields[(i + 1) * xDim + j]
                if not m.isinf(inter) and not m.isnan(inter):
                    temp += inter
                    n += 1
            if j - 1 >= 0:
                inter = fields[i * xDim + j - 1]
                if not m.isinf(inter) and not m.isnan(inter):
                    temp += inter
                    n += 1
            if j + 1 < xDim:
                inter = fields[i * xDim + j + 1]
                if not m.isinf(inter) and not m.isnan(inter):
                    temp += inter
                    n += 1
            if n > 0:
                temp /= n

        pixelValues[i][j] = temp

# deleting unnecessary data
del fields
gc.collect()

# generating appropriate image to visualise data
fig = plt.figure()
fig.set_size_inches(8, 6.5)
fig.set_dpi(500)

# linear color interpolation to show more than 256 discrete color values (continuous spectrum)
custom = mcolors.LinearSegmentedColormap.from_list("", plt.cm.turbo.colors, N=2**24)

fig.plot = plt.pcolormesh(x, y, pixelValues, cmap=custom)
cbar = plt.colorbar(fig.plot)

# showing and saving the image
fig.show()
fig.savefig('Field_image.png', dpi=500)
