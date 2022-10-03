
coil = py.c_coil.coil.Coil(0.1, 0.1, 0.1, py.int(10000), 10);
coil.set_position_and_orientation(py.c_coil.tensor.Vector3(), pi/2, pi/2);
coil.set_thread_count(py.int(16));

coil.set_default_precision_CPU(py.c_coil.coil.PrecisionFactor(11.0));
coil.set_default_precision_GPU(py.c_coil.coil.PrecisionFactor(7.0));

xPos = [];
yPos = [];

positions = py.c_coil.tensor.Vector3Array();

for i = 0:500
    for j = 0:500
        xPos(i + 1) = i * 0.001 - 0.25;
        yPos(j + 1) = j * 0.001 - 0.25;
        positions.append(py.c_coil.tensor.Vector3(xPos(i+1), yPos(j+1), 0.0));
    end
end

fields = coil.compute_all_B_field_vectors(positions, py.c_coil.coil.CPU_MT).abs();

fieldValues = zeros(501, 501);

for i = 0:500
    for j = 0:500
        fieldValues(i + 1, j + 1) = fields(i * 501 + j + 1);
    end
end

surf(xPos, yPos, fieldValues)