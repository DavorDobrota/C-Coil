import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from c_coil.coil import *
from c_coil.tensor import *
import math as m
import gc


xDim = 301
yDim = 301

xSize = 0.025
ySize = 0.1

xOffset = 0.0025
yOffset = -0.05

zPos = 0.0

prim = Coil(0.01022, 0.011, 0.0022, 20, PrecisionFactor(7.0), 16)
sec = Coil(0.01022, 0.011, 0.0022, 20, PrecisionFactor(7.0), 16)

prim.set_position_and_orientation(Vector3(), m.pi / 2, 0.0)
sec.set_position_and_orientation(Vector3(), m.pi / 2, 0.0)

primPositions = Vector3Array()
secPositions = Vector3Array()

primOrientationY = []
primOrientationZ = []
secOrientationY = []
secOrientationZ = []

matrix = []

x = []
y = []

for i in range(0, yDim):
    for j in range(0, xDim):

        xPos = xOffset + xSize / (xDim - 1) * j
        yPos = yOffset + ySize / (yDim - 1) * i

        primPositions.append(Vector3(0.0, 0.0, 0.0))
        primOrientationY.append(m.pi / 2)
        primOrientationZ.append(0.0)

        secPositions.append(Vector3(xPos, yPos, zPos))
        secOrientationY.append(m.pi / 2)
        secOrientationZ.append(0.0)

for i in range(0, xDim):
    x.append(xOffset + xSize / (xDim - 1) * i)

for i in range(0, yDim):
    y.append(yOffset + ySize / (yDim - 1) * i)

mInductances = Coil.compute_all_mutual_inductance_arrangements(prim, sec,
                                                               primPositions, secPositions,
                                                               primOrientationY, primOrientationZ,
                                                               secOrientationY, secOrientationZ,
                                                               PrecisionFactor(4.0), CPU_MT)


for i in range(0, yDim):
    mat = []
    for j in range(0, xDim):
        mat.append(mInductances[i * xDim + j] * 1000000)
    matrix.append(mat)

del mInductances
gc.collect()

fig = plt.figure()
fig.set_size_inches(8, 6)
fig.set_dpi(800)

custom = mcolors.LinearSegmentedColormap.from_list("", plt.cm.turbo.colors, N=2**24)

plot = plt.pcolormesh(x, y, matrix, cmap=custom)
cbar = plt.colorbar(plot)

plt.show()
fig.savefig('mInductanceImage.png', dpi=1000)