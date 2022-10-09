
coil = py.c_coil.coil.Coil(0.1, 0.1, 0.1, py.int(10000), 10);

zPos = [];
positions = py.c_coil.tensor.Vector3Array();

for i = 0:800 
    zPos(i + 1) = i * 0.001 - 0.4;
    positions.append(py.c_coil.tensor.Vector3(0.0, 0.0, zPos(i + 1)));
end

fields = coil.compute_all_B_field_vectors(positions);

plot (zPos, fields.abs());
