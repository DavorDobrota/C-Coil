import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import math as m
import gc

from coil_evolution.coil import *
from coil_evolution.coil_group import *
from coil_evolution.tensor import *

xDim = 6001
yDim = 6001
xSize = 0.64
ySize = 0.64
zPos = 0.0

numCoils = 12

torusGroup = CoilGroup()

torusRadius = 0.2

for i in range(numCoils):

    torusGroup.add_coil(
        torusRadius / 4, torusRadius / 8, torusRadius / 4, 5000, 100, PrecisionFactor(8.0), 16,
        Vector3.get_from_cylindrical_coords(0.0, toTorus_imagerusRadius, 2*m.pi * i / numCoils),
        m.pi / 2, 2*m.pi * i / numCoils + m.pi / 2
    )

torusGroup.set_default_precision_factor(PrecisionFactor(3.0))

positions = Vector3Array()
matrix = np.zeros((xDim, yDim), dtype=np.float64)

x = np.zeros((xDim,), dtype=np.float64)
y = np.zeros((yDim,), dtype=np.float64)

for i in range(0, yDim):
    for j in range(0, xDim):
        xPos = xSize / (xDim - 1) * (j - (xDim - 1) / 2.0)
        yPos = ySize / (yDim - 1) * (i - (yDim - 1) / 2.0)
        positions.append(Vector3(xPos, yPos, zPos))

for i in range(0, xDim):
    x[i] = xSize / (xDim - 1) * (i - (xDim - 1) / 2.0)

for i in range(0, yDim):
    y[i] = ySize / (yDim - 1) * (i - (yDim - 1) / 2.0)

fields = torusGroup.compute_all_B_field_vectors(positions, CPU_MT).abs()

del positions

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

        matrix[i][j] = temp

del fields
gc.collect()

fig = plt.figure()
fig.set_size_inches(8, 6.5)
fig.set_dpi(800)

custom = mcolors.LinearSegmentedColormap.from_list("", plt.cm.turbo.colors, N=2**24)

plot = plt.pcolormesh(x, y, matrix, cmap=custom)
cbar = plt.colorbar(plot)

plt.show()
fig.savefig('torusFieldImg.png', dpi=1000)
