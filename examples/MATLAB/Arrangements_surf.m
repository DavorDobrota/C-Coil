% if using GPU acceleration on Windows, CUDA path needs to be added
% common path is added below, but you should check for differences (e.g. different CUDA version)
% py.os.add_dll_directory("C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.7/bin")

% image resolution
xDim = int32(200);
yDim = int32(200);

% rectangular domain dimensions
xSize = 0.025;
ySize = 0.1;

% at which point the image is centered
xOffset = 0.0025;
yOffset = -0.05;
zPos = 0.0;

% coils used in computations
prim = py.c_coil.coil.Coil(0.01022, 0.011, 0.0022, int32(20));
sec = py.c_coil.coil.Coil(0.01022, 0.011, 0.0022, int32(20));

% setting precision and number of threads used (number of CPU threads is advised)
prim.set_thread_count(int32(8));
precision = py.c_coil.coil.PrecisionFactor(4.0);

% initialising arrays for specifying arrangements
primPositions = py.c_coil.tensor.Vector3Array();
secPositions = py.c_coil.tensor.Vector3Array();

primPositions.reserve((xDim + 1) * (yDim + 1));
secPositions.reserve((xDim + 1) * (yDim + 1));

primOrientationY = zeros(1, (xDim + 1) * (yDim + 1));
primOrientationZ = zeros(1, (xDim + 1) * (yDim + 1));
secOrientationY = zeros(1, (xDim + 1) * (yDim + 1));
secOrientationZ = zeros(1, (xDim + 1) * (yDim + 1));

% initialising image generation parameters
x = zeros(1, xDim + 1);
y = zeros(1, yDim + 1);
matInd = zeros(xDim + 1, yDim + 1);

% generating configuration arguments, the secondary coil changes position
for i = 0:double(yDim)
    for j = 0:double(xDim)
        xPos = xOffset + xSize / double(xDim) * j;
        yPos = yOffset + ySize / double(yDim) * i;

        primPositions.append(py.c_coil.tensor.Vector3(0.0, 0.0, 0.0))
        primOrientationY(i * (xDim + 1) + j + 1) = pi / 2.0;
        primOrientationZ(i * (xDim + 1) + j + 1) = 0.0;

        secPositions.append(py.c_coil.tensor.Vector3(xPos, yPos, zPos))
        secOrientationY(i * (xDim + 1) + j + 1) = pi / 2.0;
        secOrientationZ(i * (xDim + 1) + j + 1) = 0.0;
    end
end

% generating positions for image pixels
for i = 0:double(xDim)
    x(i + 1) = xOffset + xSize / double(xDim) * i;
end

for i = 0:double(yDim)
    y(i + 1) = yOffset + ySize / double(yDim) * i;
end

% computing mutual inductance for given arrangements
configurations = py.c_coil.coil.Coil.compute_all_mutual_inductance_arrangements(...
    prim, sec,...
    primPositions, secPositions,...
    primOrientationY, primOrientationZ,...
    secOrientationY, secOrientationZ,...
    precision, py.c_coil.coil.CPU_MT...
);

% computing force and torque for given arrangements
% configurations = py.c_coil.coil.Coil.compute_all_force_torque_arrangements(...
%     prim, sec,...
%     primPositions, secPositions,...
%     primOrientationY, primOrientationZ,...
%     secOrientationY, secOrientationZ,...
%     precision, py.c_coil.coil.CPU_MT...
% );

for i = 0:double(yDim)
    for j = 0:double(xDim)
%         mutual inductance
        matInd(i + 1, j + 1) = configurations{i * (xDim + 1) + j + 1};

%         force magniture
%         matInd(i + 1, j + 1) = configurations{i * (xDim + 1) + j + 1}{1}.abs();

%         torque magniture
%         matInd(i + 1, j + 1) = configurations{i * (xDim + 1) + j + 1}{2}.abs();
    end
end

surf(xPos, yPos, matInd)
