import c_coil.coil as c
import c_coil.tensor as ten

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import math as m
import gc

import os

# if using GPU acceleration on Windows, CUDA path needs to be added
# common path is added below, but you should check for differences (e.g. different CUDA version)
# os.add_dll_directory("C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.7/bin")

# image resolution
xDim = 301
yDim = 301

# rectangular domain dimensions
xSize = 0.025
ySize = 0.1

# at which point the image is centered
xOffset = 0.0022
yOffset = -0.05
zPos = 0.0

# coils used in computations
prim = c.Coil(0.01022, 0.011, 0.0022, 20)
sec = c.Coil(0.01022, 0.011, 0.0022, 20)

# setting precision and number of threads used (number of CPU threads is advised)
prim.set_thread_count(8)
precision = c.PrecisionFactor(4.0)

# initialising arrays for specifying arrangements
primPositions = ten.Vector3Array()
secPositions = ten.Vector3Array()

primOrientationY = []
primOrientationZ = []
secOrientationY = []
secOrientationZ = []

# initialising image generation parameters
x = []
y = []
pixelValues = []

# generating configuration arguments, the secondary coil changes position
for i in range(0, yDim):
    for j in range(0, xDim):

        xPos = xOffset + xSize / (xDim - 1) * j
        yPos = yOffset + ySize / (yDim - 1) * i

        primPositions.append(ten.Vector3(0.0, 0.0, 0.0))
        primOrientationY.append(m.pi / 2)
        primOrientationZ.append(0.0)

        secPositions.append(ten.Vector3(xPos, yPos, zPos))
        secOrientationY.append(m.pi / 2)
        secOrientationZ.append(0.0)

# generating positions for image pixels
for i in range(0, xDim):
    x.append(xOffset + xSize / (xDim - 1) * i)

for i in range(0, yDim):
    y.append(yOffset + ySize / (yDim - 1) * i)

# computing mutual inductance for given arrangements
configurations = c.Coil.compute_all_mutual_inductance_arrangements(
    prim, sec,
    primPositions, secPositions,
    primOrientationY, primOrientationZ,
    secOrientationY, secOrientationZ,
    precision, c.CPU_MT
)
# computing force and torque for given arrangements
# configurations = c.Coil.compute_all_force_torque_arrangements(
#     prim, sec,
#     primPositions, secPositions,
#     primOrientationY, primOrientationZ,
#     secOrientationY, secOrientationZ,
#     precision, c.CPU_MT
# )


for i in range(0, yDim):
    row = []
    for j in range(0, xDim):
        # mutual inductance
        row.append(configurations[i * xDim + j])

        # force magnitude
        # row.append(configurations[i * xDim + j][0].abs())

        # torque magnitude
        # row.append(configurations[i * xDim + j][1].abs())
    pixelValues.append(row)

# deleting unnecessary data and calling gc to reduce memory use
del primPositions
del secPositions

del primOrientationY
del primOrientationZ
del secOrientationY
del secOrientationZ

del configurations

gc.collect()

# generating appropriate image to visualise data
fig = plt.figure()
fig.set_size_inches(8, 6)
fig.set_dpi(500)

# linear color interpolation to show more than 256 discrete color values (continuous spectrum)
custom = mcolors.LinearSegmentedColormap.from_list("", plt.cm.turbo.colors, N=2**24)

fig.plot = plt.pcolormesh(x, y, pixelValues, cmap=custom)
cbar = plt.colorbar(fig.plot)

# showing and saving the image
fig.show()
fig.savefig('Arrangements_image.png', dpi=500)
