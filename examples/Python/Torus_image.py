import coil_evolution.coil as c
import coil_evolution.coil_group as cg
import coil_evolution.tensor as ten

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import numpy as np
import math as m

import gc
import os

# if using GPU acceleration on Windows, CUDA path needs to be added
# common path is added below, but you should check for differences (e.g. different CUDA version)
# os.add_dll_directory("C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.7/bin")

# image resolution
xDim = 3001
yDim = 3001

# rectangular domain dimensions
xSize = 0.64
ySize = 0.64

# at which point the image is centered
xOffset = 0.0
yOffset = 0.0
zPos = 0.0

# torus parameters
numCoils = 12
torusRadius = 0.2
precision = c.PrecisionFactor(5.0)

# initialising torus group to which member coils are added
torusGroup = cg.CoilGroup()

# adding member coils
for i in range(numCoils):

    torusGroup.add_coil(
        torusRadius / 4,
        torusRadius / 8,
        torusRadius / 4,
        5000,
        100,
        precision,
        16,
        ten.Vector3.get_from_cylindrical_coords(0.0, torusRadius, 2*m.pi * i / numCoils),
        m.pi / 2,
        2*m.pi * i / numCoils + m.pi / 2
    )

# setting precision and number of threads used (number of CPU threads is advised)
torusGroup.set_thread_count(8)
torusGroup.set_default_precision_factor(precision)

# initialising array for specifying points
positions = ten.Vector3Array()

# initialising image generation parameters
x = np.zeros((xDim,), dtype=np.float64)
y = np.zeros((yDim,), dtype=np.float64)
pixelValues = np.zeros((xDim, yDim), dtype=np.float64)

# generating point position arguments
for i in range(0, yDim):
    for j in range(0, xDim):
        xPos = xSize / (xDim - 1) * (j - (xDim - 1) / 2.0)
        yPos = ySize / (yDim - 1) * (i - (yDim - 1) / 2.0)

        positions.append(ten.Vector3(xPos, yPos, zPos))

# generating positions for image pixels
for i in range(0, xDim):
    x[i] = xSize / (xDim - 1) * (i - (xDim - 1) / 2.0)

for i in range(0, yDim):
    y[i] = ySize / (yDim - 1) * (i - (yDim - 1) / 2.0)

# calling the field compute function
fields = torusGroup.compute_all_B_field_vectors(positions, c.CPU_MT).abs()

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
fig.savefig('Torus_image.png', dpi=500)
