
xDim = int32(200);
yDim = int32(200);

xSize = 0.025;
ySize = 0.1;

xOffset = 0.0025;
yOffset = -0.05;

zPos = 0.0;

prim = py.c_coil.coil.Coil(0.01022, 0.011, 0.0022, int32(20));
sec = py.c_coil.coil.Coil(0.01022, 0.011, 0.0022, int32(20));

prim.set_position_and_orientation(py.c_coil.tensor.Vector3(), pi / 2, 0.0);
sec.set_position_and_orientation(py.c_coil.tensor.Vector3(), pi / 2, 0.0);

prim.set_thread_count(int32(16));

primPositions = py.c_coil.tensor.Vector3Array();
secPositions = py.c_coil.tensor.Vector3Array();

primPositions.reserve((xDim + 1) * (yDim + 1));
secPositions.reserve((xDim + 1) * (yDim + 1));

primOrientationY = zeros(1, (xDim + 1) * (yDim + 1));
primOrientationZ = zeros(1, (xDim + 1) * (yDim + 1));
secOrientationY = zeros(1, (xDim + 1) * (yDim + 1));
secOrientationZ = zeros(1, (xDim + 1) * (yDim + 1));

xPos = zeros(1, xDim + 1);
yPos = zeros(1, yDim + 1);

for i = 0:double(yDim)
    for j = 0:double(xDim)
        x = xOffset + xSize / double(xDim) * j;
        y = yOffset + ySize / double(yDim) * i;

        primPositions.append(py.c_coil.tensor.Vector3(0.0, 0.0, 0.0))
        primOrientationY(i * (xDim + 1) + j + 1) = pi / 2.0;
        primOrientationZ(i * (xDim + 1) + j + 1) = 0.0;

        secPositions.append(py.c_coil.tensor.Vector3(x, y, zPos))
        secOrientationY(i * (xDim + 1) + j + 1) = pi / 2.0;
        secOrientationZ(i * (xDim + 1) + j + 1) = 0.0;
    end
end

for i = 0:double(xDim)
    xPos(i + 1) = xOffset + xSize / double(xDim) * i;
end

for i = 0:double(yDim)
    yPos(i + 1) = yOffset + ySize / double(yDim) * i;
end

mInductances = py.c_coil.coil.Coil.compute_all_mutual_inductance_arrangements(...
    prim, sec, primPositions, secPositions,...
    primOrientationY, primOrientationZ,...
    secOrientationY, secOrientationZ,...
    py.c_coil.coil.PrecisionFactor(4.0), py.c_coil.coil.CPU_MT...
);

matInd = zeros(xDim + 1, yDim + 1);

for i = 0:double(yDim)
    for j = 0:double(xDim)
        matInd(i + 1, j + 1) = mInductances{i * (xDim + 1) + j + 1};
    end
end

surf(xPos, yPos, matInd)
