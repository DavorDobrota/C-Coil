% if using GPU acceleration on Windows, CUDA path needs to be added
% common path is added below, but you should check for differences (e.g. different CUDA version)
% py.os.add_dll_directory("C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.7/bin")

% image resolution
xDim = int32(500);
yDim = int32(500);

% rectangular domain dimensions
xSize = 0.24;
ySize = 0.24;

% at which point the image is centered
xOffset = 0.0;
yOffset = 0.0;
zPos = 0.0;

% coil used for field computations, orientation set such that it lies in xOy plane
coil = py.c_coil.coil.Coil(0.03, 0.03, 0.12, py.int(10000), 10);
coil.set_position_and_orientation(py.c_coil.tensor.Vector3(), pi/2, pi/2);

% setting precision and number of threads used (number of CPU threads is advised)
coil.set_thread_count(py.int(16));
precision = py.c_coil.coil.PrecisionFactor(7.0);

coil.set_default_precision_CPU(py.c_coil.coil.PrecisionFactor(11.0));
coil.set_default_precision_GPU(py.c_coil.coil.PrecisionFactor(7.0));

% initialising array for specifying points
positions = py.c_coil.tensor.Vector3Array();
positions.reserve(int32((xDim + 1) * (yDim + 1)));

% initialising image generation parameters, using np arrays to reduce memory use
x = zeros(1, 501);
y = zeros(1, 501);
fieldValues = zeros(501, 501);

% generating point position arguments
for i = 0:double(yDim)
    for j = 0:double(xDim)
        xPos = xOffset + xSize / double(xDim) * (j - double(xDim) / 2.0);
        yPos = yOffset + ySize / double(yDim) * (i - double(yDim) / 2.0);
        positions.append(py.c_coil.tensor.Vector3(xPos, yPos, zPos));
    end
end

% generating positions for image pixels
for i = 0:double(xDim)
    x(i + 1) = xOffset + xSize / double(xDim) * i;
end

for i = 0:double(yDim)
    y(i + 1) = yOffset + ySize / double(yDim) * i;
end

% calling the field compute function
fields = coil.compute_all_B_field_vectors(positions, py.c_coil.coil.CPU_MT).abs();

for i = 0:500
    for j = 0:500
        fieldValues(i + 1, j + 1) = fields(i * 501 + j + 1);
    end
end

surf(x, y, fieldValues)