#include <sstream>
#include <string>

#include <PyMain.h>
#include <Tensor.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;


void initTensor(py::module_ &mainModule)
{
    py::module_ tensorModule = mainModule.def_submodule("tensor");

    py::class_<vec3::Vector3> vector3(tensorModule, "Vector3");
    py::class_<vec3::Matrix3> matrix3(tensorModule, "Matrix3");
    py::class_<vec3::Triplet> triplet(tensorModule, "Triplet");

    py::class_<vec3::Vector3Array> vector3Array(tensorModule, "Vector3Array");
    py::class_<vec3::Matrix3Array> matrix3Array(tensorModule, "Matrix3Array");


    // FieldVector3

    vector3.def_readwrite("x", &vec3::Vector3::x)
        .def_readwrite("y", &vec3::Vector3::y)
        .def_readwrite("z", &vec3::Vector3::z);

    vector3.def(py::init<>())
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"));

    vector3.def("__add__", &vec3::Vector3::operator+)
        .def("__iadd__", &vec3::Vector3::operator+=)
        .def("__sub__", &vec3::Vector3::operator-)
        .def("__isub__", &vec3::Vector3::operator-=)
        .def("__mul__", &vec3::Vector3::operator*)
        .def("__imul__", &vec3::Vector3::operator*=);

    vector3.def("abs", &vec3::Vector3::abs);

    vector3.def_static(
        "scalar_product", &vec3::Vector3::scalarProduct,
        py::arg("vector1"), py::arg("vector2"));

    vector3.def_static(
        "cross_product", &vec3::Vector3::crossProduct,
        py::arg("vector1"), py::arg("vector2"));

    vector3.def_static(
        "get_from_cylindrical_coords", &vec3::Vector3::getFromCylindricalCoords,
        py::arg("z"), py::arg("r"), py::arg("phi"));

    vector3.def_static(
        "get_from_spherical_coords", &vec3::Vector3::getFromSphericalCoords,
        py::arg("z"), py::arg("theta"), py::arg("phi"));

    vector3
        .def("get_as_cylindrical_coords", &vec3::Vector3::getAsCylindricalCoords)
        .def("get_as_spherical_coords", &vec3::Vector3::getAsSphericalCoords);

    vector3.def("__repr__", &vec3::Vector3::operator std::string);


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

    matrix3.def("det", &vec3::Matrix3::det);

    matrix3.def("__add__", &vec3::Matrix3::operator+)
        .def("__iadd__", &vec3::Matrix3::operator+=)
        .def("__imul__", &vec3::Matrix3::operator*=)
        .def(
            "__mul__",
            static_cast<vec3::Matrix3 (vec3::Matrix3::*)(double) const>(&vec3::Matrix3::operator*))
        .def(
            "__mul__",
            static_cast<vec3::Matrix3 (vec3::Matrix3::*)(const vec3::Matrix3&) const>(&vec3::Matrix3::operator*)
        )
        .def(
            "__mul__",
            static_cast<vec3::Vector3 (vec3::Matrix3::*)(const vec3::Vector3&) const>(&vec3::Matrix3::operator*)
        );

    matrix3.def("__repr__", &vec3::Matrix3::operator std::string);


    // Triplet

    triplet.def_readwrite("first", &vec3::Triplet::first)
           .def_readwrite("second", &vec3::Triplet::second)
           .def_readwrite("third", &vec3::Triplet::third);

    triplet.def(py::init<>())
           .def(
               py::init<double, double, double>(),
               py::arg("first"), py::arg("second"), py::arg("third")
           );

    triplet.def("__repr__", &vec3::Triplet::operator std::string);


    // Vector3Array

    vector3Array.def(py::init<>())
        .def(py::init<size_t>(), py::arg("init_size"))
        .def(py::init<const std::vector<vec3::Vector3>&>(), py::arg("vector_list"));

    vector3Array
        .def(
            "append",
            static_cast<void (vec3::Vector3Array::*)(const vec3::Vector3&)>(&vec3::Vector3Array::append),
            py::arg("appended_vector"))
        .def(
            "append",
            static_cast<void (vec3::Vector3Array::*)(double, double, double)>(&vec3::Vector3Array::append),
            py::arg("x"), py::arg("y"), py::arg("z"))
        .def("reserve", &vec3::Vector3Array::reserve, py::arg("reserve_size"))
        .def("resize", &vec3::Vector3Array::resize, py::arg("new_size"))
        .def("clear", &vec3::Vector3Array::clear)
        .def("size", &vec3::Vector3Array::size);

    vector3Array.def("items", &vec3::Vector3Array::getItems);

    vector3Array.def("x", &vec3::Vector3Array::x)
        .def("y", &vec3::Vector3Array::y)
        .def("z", &vec3::Vector3Array::z)
        .def("abs", &vec3::Vector3Array::abs);

    vector3Array
        .def(
            "__getitem__",
            [](const vec3::Vector3Array &self, long long index){
                if(index < 0)
                    index += self.size();

                if(index < 0 || index >= self.size())
                    throw std::out_of_range("Array index out of range!");
                return self[index];
            }
        )
        .def(
            "__setitem__",
            [](vec3::Vector3Array &self, long long index, const vec3::Vector3 &item){
                if(index < 0)
                    index += self.size();

                if(index < 0 || index >= self.size())
                    throw std::out_of_range("Array index out of range!");
                self[index] = item;
            }
        )
        .def("__iadd__", &vec3::Vector3Array::operator+=);

    vector3Array.def("__repr__", &vec3::Vector3Array::operator std::string);


    // Matrix3Array

    matrix3Array.def(py::init<>())
            .def(py::init<size_t>(), py::arg("init_size"))
            .def(py::init<const std::vector<vec3::Matrix3>&>(), py::arg("matrix_list"));

    matrix3Array
        .def(
            "append",
            static_cast<void (vec3::Matrix3Array::*)(const vec3::Matrix3&)>(&vec3::Matrix3Array::append),
            py::arg("appended_matrix"))
        .def(
            "append",
            static_cast<
                void (vec3::Matrix3Array::*)
                (double, double, double, double, double, double, double, double, double)
            > (&vec3::Matrix3Array::append),
            py::arg("xx"), py::arg("xy"), py::arg("xz"),
            py::arg("yx"), py::arg("yy"), py::arg("yz"),
            py::arg("zx"), py::arg("zy"), py::arg("zz")
        )
        .def("reserve", &vec3::Matrix3Array::reserve, py::arg("reserve_size"))
        .def("resize", &vec3::Matrix3Array::resize, py::arg("new_size"))
        .def("clear", &vec3::Matrix3Array::clear)
        .def("size", &vec3::Matrix3Array::size);

    matrix3Array.def("items", &vec3::Matrix3Array::getItems);

    matrix3Array
        .def("xx", &vec3::Matrix3Array::xx)
        .def("xy", &vec3::Matrix3Array::xy)
        .def("xz", &vec3::Matrix3Array::xz)
        .def("yx", &vec3::Matrix3Array::yx)
        .def("yy", &vec3::Matrix3Array::yy)
        .def("yz", &vec3::Matrix3Array::yz)
        .def("zx", &vec3::Matrix3Array::zx)
        .def("zy", &vec3::Matrix3Array::zy)
        .def("zz", &vec3::Matrix3Array::zz)
        .def("det", &vec3::Matrix3Array::det);

    matrix3Array
        .def(
            "__getitem__",
            [](const vec3::Matrix3Array &self, long long index){
                if(index < 0)
                    index += self.size();

                if(index < 0 || index >= self.size())
                    throw std::out_of_range("Array index out of range!");
                return self[index];
            }
        )
        .def(
            "__setitem__",
            [](vec3::Matrix3Array &self, long long index, const vec3::Matrix3 &item){
                if(index < 0)
                    index += self.size();

                if(index < 0 || index >= self.size())
                    throw std::out_of_range("Array index out of range!");
                self[index] = item;
            }
        )
        .def("__iadd__", &vec3::Matrix3Array::operator+=);

    matrix3Array.def("__repr__", &vec3::Matrix3Array::operator std::string);
}
