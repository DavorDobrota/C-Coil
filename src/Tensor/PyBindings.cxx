#include <sstream>

#include <Tensor.h>
#include <PyMain.h>
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

    coordVector3.def_readwrite("component1", &vec3::CoordVector3::component1)
        .def_readwrite("component2", &vec3::CoordVector3::component2)
        .def_readwrite("component3", &vec3::CoordVector3::component3);

    coordVector3.def(py::init<>())
        .def(py::init<vec3::CoordinateSystem, double, double, double>());

    coordVector3.def("is_cartesian", &vec3::CoordVector3::isCartesian)
        .def("is_cylindrical", &vec3::CoordVector3::isCylindrical)
        .def("is_spherical", &vec3::CoordVector3::isSpherical);

    coordVector3.def("convert_to_cartesian", &vec3::CoordVector3::convertToCartesian)
        .def("convert_to_cylindrical", &vec3::CoordVector3::convertToCylindrical)
        .def("convert_to_spherical", &vec3::CoordVector3::convertToSpherical);

    coordVector3.def_static("convert_all_to_cartesian", &vec3::CoordVector3::convertAllToCartesian)
            .def_static("convert_all_to_cylindrical", &vec3::CoordVector3::convertAllToCylindrical)
            .def_static("convert_all_to_spherical", &vec3::CoordVector3::convertAllToSpherical);

    coordVector3.def_static("convert_to_field_vector", &vec3::CoordVector3::convertToFieldVector)
        .def_static("convert_to_coord_vector", &vec3::CoordVector3::convertToCoordVector);

    coordVector3.def_property_readonly("coordinateSystem", &vec3::CoordVector3::getCoordinateSystem);

    coordVector3.def("__repr__", [](vec3::CoordVector3 self) -> std::string {
        std::stringstream output;

        switch(self.getCoordinateSystem())
        {
            case vec3::CoordinateSystem::CARTESIAN:
                output << "{x: " << self.component1 << ", y: " << self.component2
                       << ", z: " << self.component3 << "}";
                break;
            case vec3::CoordinateSystem::CYLINDRICAL:
                output << "{z: " << self.component1 << ", r: " << self.component2
                       << ", p: " << self.component3 << "}";
                break;
            case vec3::CoordinateSystem::SPHERICAL:
                output << "{r: " << self.component1 << ", t: " << self.component2
                       << ", p: " << self.component3 << "}";
                break;
        };

        return output.str();
    });


    // FieldVector3

    fieldVector3.def_readwrite("x_component", &vec3::FieldVector3::xComponent)
        .def_readwrite("y_component", &vec3::FieldVector3::yComponent)
        .def_readwrite("z_component", &vec3::FieldVector3::zComponent);

    fieldVector3.def(py::init<>())
        .def(py::init<double, double, double>());

    fieldVector3.def("__add__", &vec3::FieldVector3::operator+)
        .def("__iadd__", &vec3::FieldVector3::operator+=)
        .def("__sub__", &vec3::FieldVector3::operator-)
        .def("__isub__", &vec3::FieldVector3::operator-=)
        .def("__mul__", &vec3::FieldVector3::operator*)
        .def("__imul__", &vec3::FieldVector3::operator*=);

    fieldVector3.def("magnitude", &vec3::FieldVector3::magnitude);

    fieldVector3.def_static("scalar_product", &vec3::FieldVector3::scalarProduct);
    fieldVector3.def_static("cross_product", &vec3::FieldVector3::crossProduct);

    fieldVector3.def("__repr__", [](vec3::FieldVector3 self) -> std::string {
        std::stringstream output;

        output << "{x: " << self.xComponent << ", y: " << self.yComponent
               << ", z: " << self.zComponent << "}";

        return output.str();
    });


    // Matrix3

    matrix3.def_readwrite("xx_element", &vec3::Matrix3::xxElement)
        .def_readwrite("xy_element", &vec3::Matrix3::xyElement)
        .def_readwrite("xz_element", &vec3::Matrix3::xzElement)
        .def_readwrite("yx_element", &vec3::Matrix3::yxElement)
        .def_readwrite("yy_element", &vec3::Matrix3::yyElement)
        .def_readwrite("yz_element", &vec3::Matrix3::yzElement)
        .def_readwrite("zx_element", &vec3::Matrix3::zxElement)
        .def_readwrite("zy_element", &vec3::Matrix3::zyElement)
        .def_readwrite("zz_element", &vec3::Matrix3::zzElement);

    matrix3.def(py::init<>())
        .def(py::init<double, double, double, double, double, double, double, double, double>());

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

    matrix3.def("__repr__", [](vec3::Matrix3 self) -> std::string {
        std::stringstream output;

        output << "[[" << self.xxElement << ", " << self.xyElement << ", " << self.xzElement << "], ["
                       << self.yxElement << ", " << self.yyElement << ", " << self.yzElement << "], ["
                       << self.zxElement << ", " << self.zyElement << ", " << self.zzElement << "]]";

        return output.str();
    });
}
