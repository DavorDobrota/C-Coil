#include <string>

#include <PyMain.h>
#include <Tensor.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


void initTensor(py::module_ &mainModule)
{
    py::module_ tensorModule = mainModule.def_submodule("tensor");
    
    py::enum_<vec3::CoordinateSystem> coordinateSystem(tensorModule, "CoordinateSystem");
    py::class_<vec3::CoordVector3> coordVector3(tensorModule, "CoordVector3");
    py::class_<vec3::FieldVector3> fieldVector3(tensorModule, "FieldVector3");
    py::class_<vec3::Matrix3> matrix3(tensorModule, "Matrix3");


    // CoordinateSystem

    coordinateSystem.value("CARTESIAN", vec3::CoordinateSystem::CARTESIAN)
        .value("CYLINDRICAL", vec3::CoordinateSystem::CYLINDRICAL)
        .value("SPHERICAL", vec3::CoordinateSystem::SPHERICAL)
        .export_values();


    // CoordVector3

    coordVector3.def_readwrite("elem_1", &vec3::CoordVector3::elem1)
        .def_readwrite("elem_2", &vec3::CoordVector3::elem2)
        .def_readwrite("elem_3", &vec3::CoordVector3::elem3);

    coordVector3.def(py::init<>())
        .def(
            py::init<vec3::CoordinateSystem, double, double, double>(),
            py::arg("system"), py::arg("elem_1"), py::arg("elem_2"), py::arg("elem_3"));

    coordVector3.def("is_cartesian", &vec3::CoordVector3::isCartesian)
        .def("is_cylindrical", &vec3::CoordVector3::isCylindrical)
        .def("is_spherical", &vec3::CoordVector3::isSpherical);

    coordVector3.def("convert_to_cartesian", &vec3::CoordVector3::convertToCartesian)
        .def("convert_to_cylindrical", &vec3::CoordVector3::convertToCylindrical)
        .def("convert_to_spherical", &vec3::CoordVector3::convertToSpherical);

    coordVector3.def_static(
            "convert_all_to_cartesian", &vec3::CoordVector3::convertAllToCartesian, py::arg("vectors"))
        .def_static("convert_all_to_cylindrical", &vec3::CoordVector3::convertAllToCylindrical, py::arg("vectors"))
        .def_static("convert_all_to_spherical", &vec3::CoordVector3::convertAllToSpherical, py::arg("vectors"));

    coordVector3.def_static("convert_to_field_vector", &vec3::CoordVector3::convertToFieldVector, py::arg("vector"))
        .def_static("convert_to_coord_vector", &vec3::CoordVector3::convertToCoordVector, py::arg("vector"));

    coordVector3.def_property_readonly("coordinate_system", &vec3::CoordVector3::getCoordinateSystem);

    coordVector3.def("__repr__", &vec3::CoordVector3::operator std::string);


    // FieldVector3

    fieldVector3.def_readwrite("x", &vec3::FieldVector3::x)
        .def_readwrite("y", &vec3::FieldVector3::y)
        .def_readwrite("z", &vec3::FieldVector3::z);

    fieldVector3.def(py::init<>())
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"));

    fieldVector3.def("__add__", &vec3::FieldVector3::operator+)
        .def("__iadd__", &vec3::FieldVector3::operator+=)
        .def("__sub__", &vec3::FieldVector3::operator-)
        .def("__isub__", &vec3::FieldVector3::operator-=)
        .def("__mul__", &vec3::FieldVector3::operator*)
        .def("__imul__", &vec3::FieldVector3::operator*=);

    fieldVector3.def("magnitude", &vec3::FieldVector3::magnitude);

    fieldVector3.def_static(
        "scalar_product", &vec3::FieldVector3::scalarProduct,
        py::arg("vector1"), py::arg("vector2"));

    fieldVector3.def_static(
        "cross_product", &vec3::FieldVector3::crossProduct,
        py::arg("vector1"), py::arg("vector2"));

    fieldVector3.def("__repr__", &vec3::FieldVector3::operator std::string);


    // Matrix3

    matrix3.def_readwrite("xx", &vec3::Matrix3::xx)
        .def_readwrite("xy", &vec3::Matrix3::xy)
        .def_readwrite("xz", &vec3::Matrix3::xz)
        .def_readwrite("yx", &vec3::Matrix3::yx)
        .def_readwrite("yy", &vec3::Matrix3::yy)
        .def_readwrite("yz", &vec3::Matrix3::yz)
        .def_readwrite("zx", &vec3::Matrix3::zx)
        .def_readwrite("zy", &vec3::Matrix3::zy)
        .def_readwrite("zz", &vec3::Matrix3::zz);

    matrix3.def(py::init<>())
        .def(
            py::init<double, double, double, double, double, double, double, double, double>(),
            py::arg("xx") = 0.0, py::arg("xy") = 0.0, py::arg("xz") = 0.0,
            py::arg("yx") = 0.0, py::arg("yy") = 0.0, py::arg("yz") = 0.0,
            py::arg("zx") = 0.0, py::arg("zy") = 0.0, py::arg("zz") = 0.0);

    matrix3.def("__add__", &vec3::Matrix3::operator+)
        .def("__iadd__", &vec3::Matrix3::operator+=)
        .def(
            "__mul__",
            static_cast<vec3::Matrix3 (vec3::Matrix3::*)(const vec3::Matrix3&) const>(&vec3::Matrix3::operator*)
        )
        .def(
            "__mul__",
            static_cast<vec3::FieldVector3 (vec3::Matrix3::*)(const vec3::FieldVector3&) const>(&vec3::Matrix3::operator*)
        );

    matrix3.def("__repr__", &vec3::Matrix3::operator std::string);
}
