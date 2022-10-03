import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import math as m
import gc

from c_coil.coil import *
from c_coil.tensor import *

xDim = 5001
yDim = 5001
xSize = 0.24
ySize = 0.24
zPos = 0.0

coil = Coil(0.03, 0.03, 0.12, 3600, 10, PrecisionFactor(11.0), 16)
coil.set_position_and_orientation(Vector3(), m.pi / 2, 0.0)

positions = Vector3Array()

matrix = np.zeros((xDim, yDim), dtype=np.float64)

x = np.zeros((xDim,), dtype=np.float64)
y = np.zeros((yDim,), dtype=np.float64)

for i in range(0, yDim):
    for j in range(0, xDim):
        xPos = xSize / (xDim - 1) * (j - (xDim - 1) / 2.0)
        yPos = ySize / (yDim - 1) * (i - (yDim - 1) / 2.0)
        positions.append(xPos, yPos, zPos)

for i in range(0, xDim):
    x[i] = xSize / (xDim - 1) * (i - (xDim - 1) / 2.0)

for i in range(0, yDim):
    y[i] = ySize / (yDim - 1) * (i - (yDim - 1) / 2.0)

fields = coil.compute_all_B_field_vectors(positions, CPU_MT).abs()

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
fig.set_size_inches(8.0, 6.5)
fig.set_dpi(500)

custom = mcolors.LinearSegmentedColormap.from_list("gist_yarg", plt.get_cmap("turbo").colors, N=2**24)

plot = plt.pcolormesh(x, y, matrix, cmap=plt.get_cmap("Greys"))
cbar = plt.colorbar(plot)

plt.show()
fig.savefig('Field_image.png', dpi=500)
