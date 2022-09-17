#ifndef GENERAL_COIL_PROGRAM_VECTOR3_H
#define GENERAL_COIL_PROGRAM_VECTOR3_H

#include <vector>
#include <iostream>


namespace vec3
{

    class Vector3;
    class Matrix3;

    class Triplet;

    class Vector3Array;
    class Matrix3Array;

    /**
     * @brief Represents a rank 1 tensor of dimension 3, which is a member of the oriented Euclidean vector space.
     * It is commonly referred to as a 3D (Cartesian) Vector.
     * @details Basic operations are supported such as addition (+, +=), subtraction (-, -=), multiplication (*, *=),
     * and magnitude (abs()). The inner (dot) product and vector (cross) product are defined for any two Vector3 objects.
     * For enhanced flexibility, the vector can be defined and obtained in spherical and cylindrical coordinates.
     */
    class Vector3
    {
        public:

            double x;
            double y;
            double z;

            ///@brief Default constructor, creates a null vector (0, 0, 0).
            Vector3();
            ///@brief Creates a vector with coordinates (x, y, z).
            explicit Vector3(double x, double y, double z);

            Vector3 operator+(const Vector3 &otherVec) const;
            Vector3 operator+=(const Vector3 &otherVec);
            Vector3 operator-(const Vector3 &otherVec) const;
            Vector3 operator-=(const Vector3 &otherVec);
            Vector3 operator*(double multiplier) const;
            Vector3 operator*=(double multiplier);

            ///@brief Returns the magnitude (Euclidean norm) of the vector.
            [[nodiscard]] double abs() const;

            ///@brief Returns the inner product of two Vector3 objects.
            static double scalarProduct(Vector3 vector1, Vector3 vector2);
            ///@brief Returns the vector product of two Vector3 objects.
            static Vector3 crossProduct(Vector3 vector1, Vector3 vector2);

            ///@brief Creates a Vector3 from cylindrical coordinates (z, r, phi), phi [0, 2PI].
            static Vector3 getFromCylindricalCoords(double z, double r, double phi);
            ///@brief Creates a Vector3 from spherical coordinates (r,theta, phi), theta [0, PI], phi [0, 2PI].
            static Vector3 getFromSphericalCoords(double r, double theta, double phi);

            ///@brief Returns a Triplet representing the vector in cylindrical coordinates (z, r, phi).
            [[nodiscard]] Triplet getAsCylindricalCoords() const;
            ///@brief Returns a Triplet representing the vector in spherical coordinates (z, r, phi).
            [[nodiscard]] Triplet getAsSphericalCoords() const;

            ///@brief Generates a string object with all components of the 3D Vector.
            explicit operator std::string() const;
    };

    /**
     * @brief Represents a rank 2 tensor of dimension 3, which is a represented as a square matrix.
     * @details Basic operations are supported, such as matrix addition (+, +=), multiplication by a scalar (*, *=),
     * Vector3 transformation (*), matrix multiplication (*), and determinant calculation (det()).
     */
    class Matrix3
    {
        public:

            double xx;
            double xy;
            double xz;

            double yx;
            double yy;
            double yz;

            double zx;
            double zy;
            double zz;

            ///@brief Default constructor, creates a null matrix.
            Matrix3();
            ///@brief Creates matrix with components ((xx, xy, xz), (yx, yy, yz), (zx, zy, zz)).
            explicit Matrix3(double xx, double xy, double xz, double yx, double yy, double yz, double zx, double zy, double zz);

            ///@brief Returns the determinant of the determinant of a square matrix.
            [[nodiscard]] double det() const;

            Matrix3 operator+(const Matrix3 &mat) const;
            Matrix3 operator+=(const Matrix3 &mat);

            Matrix3 operator*(double multiplier) const;
            void operator*=(double multiplier);
            Matrix3 operator*(const Matrix3 &mat) const;
            Vector3 operator*(const Vector3 &vec) const;

            ///@brief Generates a string object with all components of the 3x3 Matrix.
            explicit operator std::string() const;
    };

    /**
     * @brief Represents a general ordered sequence (tuple) with 3 elements.
     * @details Used for representing Vector3 data in different forms, such as cylindrical and spherical coordinates,
     * which are not apt for proper calculations, but are useful for intermediate calculations.
     */
    class Triplet
    {
        public:

            double first;
            double second;
            double third;

            ///@brief Default constructor, returns a tuple of zeros (0, 0, 0).
            Triplet();
            ///@brief Creates a 3-tuple with elements (first, second, third).
            explicit Triplet(double first, double second, double third);

            ///@brief Generates a string object with the values of the triplet.
            explicit operator std::string() const;
    };

    /**
     * @brief Represents std::vector<vec3::Vector3> for easier handling and additional features.
     * Allows only x, y, or z components, as well as abs() values, to be extracted to std::vector<double>.
     * @details A reference to the encapsulated std::vector<vec3::Vector3> can be retrieved for faster C++ calculations.
     * Basic std::vector functionality is implemented (reserve, resize, clear, size)
     * with a python inspired append method for adding elements. Elements can also be added with +=.
     * Reduces memory use in Python.
     */
    class Vector3Array
    {
        private:

            std::vector<Vector3> vectorArray;

        public:

            ///@brief Default constructor, creates an empty encapsulated std::vector<vec3::Vector3>.
            Vector3Array();
            ///@brief Creates an encapsulated std::vector<vec3::Vector3> of given size.
            explicit Vector3Array(size_t initSize);
            ///@brief Creates an encapsulated std::vector<vec3::Vector3> filled with given values.
            explicit Vector3Array(const std::vector<Vector3> &vectorArray);

            ///@brief Applies std::vector push_back with the given Vector3.
            void append(const Vector3 &appendedVector3);
            ///@brief Applies std::vector emplace_back with given 3 values (x, y, z).
            void append(double x, double y, double z);
            ///@brief Applies std::vector reserve with the provided size for faster append operations.
            void reserve(size_t reserveSize);
            ///@brief Applies std::vector resize with the provided size.
            void resize(size_t newSize);
            ///@brief Applies std::vector clear.
            void clear();
            ///@brief Returns the size of encapsulated std::vector.
            [[nodiscard]] size_t size() const;

            ///@brief Returns a reference to encapsulated std::vector<vec3::Vector3>.
            std::vector<Vector3>& getItems();

            ///@brief Returns a std::vector<double> of only x components of Vector3.
            [[nodiscard]] std::vector<double> x() const;
            ///@brief Returns a std::vector<double> of only y components of Vector3.
            [[nodiscard]] std::vector<double> y() const;
            ///@brief Returns a std::vector<double> of only z components of Vector3.
            [[nodiscard]] std::vector<double> z() const;
            ///@brief Returns a std::vector<double> of magnitudes of Vector3.
            [[nodiscard]] std::vector<double> abs() const;

            Vector3& operator[](size_t index);
            const Vector3& operator[](size_t index) const;
            Vector3Array& operator+=(const Vector3 &appendedVector3);

            ///@brief Generates a string object from all elements of encapsulated std::vector.
            explicit operator std::string() const;
    };

    /**
     * @brief Represents std::vector<vec3::Matrix3> for easier handling and additional features.
     * Allows only xx, xy, xz, yx, yy, yz, zx, zy, or zz components, and det() values, to be extracted to a std::vector.
     * @details A reference to the encapsulated std::vector<vec3::Matrix3> can be retrieved for faster C++ calculations.
     * Basic std::vector functionality is implemented (reserve, resize, clear, size)
     * with a python inspired append method for adding elements. Elements can also be added with +=.
     * Reduces memory use in Python.
     */
    class Matrix3Array
    {
        private:

        std::vector<Matrix3> matrixArray;

        public:

            ///@brief Default constructor, creates an empty encapsulated std::vector<vec3::Matrix3>.
            Matrix3Array();
            ///@brief Creates an encapsulated std::vector<vec3::Matrix3> of given size.
            explicit Matrix3Array(size_t initSize);
            ///@brief Creates an encapsulated std::vector<vec3::Vector3> filled with given values.
            explicit Matrix3Array(const std::vector<Matrix3> &matrixArray);

            ///@brief Applies std::vector push_back with the given Matrix3.
            void append(const Matrix3 &appendedMatrix3);
            ///@brief Applies std::vector emplace_back with given 3 values ((xx, xy, xz), (yx, yy, yz), (zx, zy, zz)).
            void append(double xx, double xy, double xz, double yx, double yy, double yz, double zx, double zy, double zz);
            ///@brief Applies std::vector reserve with the provided size for faster append operations.
            void reserve(size_t reserveSize);
            ///@brief Applies std::vector resize with the provided size.
            void resize(size_t newSize);
            ///@brief Returns the size of encapsulated std::vector.
            void clear();
            [[nodiscard]] size_t size() const;

            ///@brief Returns a reference to encapsulated std::vector<vec3::Matrix3>.
            std::vector<Matrix3>& getItems();

            ///@brief Returns a std::vector<double> of only xx components of Matrix3.
            [[nodiscard]] std::vector<double> xx() const;
            ///@brief Returns a std::vector<double> of only xy components of Matrix3.
            [[nodiscard]] std::vector<double> xy() const;
            ///@brief Returns a std::vector<double> of only xz components of Matrix3.
            [[nodiscard]] std::vector<double> xz() const;
            ///@brief Returns a std::vector<double> of only yx components of Matrix3.
            [[nodiscard]] std::vector<double> yx() const;
            ///@brief Returns a std::vector<double> of only yy components of Matrix3.
            [[nodiscard]] std::vector<double> yy() const;
            ///@brief Returns a std::vector<double> of only yz components of Matrix3.
            [[nodiscard]] std::vector<double> yz() const;
            ///@brief Returns a std::vector<double> of only zx components of Matrix3.
            [[nodiscard]] std::vector<double> zx() const;
            ///@brief Returns a std::vector<double> of only zy components of Matrix3.
            [[nodiscard]] std::vector<double> zy() const;
            ///@brief Returns a std::vector<double> of only zz components of Matrix3.
            [[nodiscard]] std::vector<double> zz() const;
            ///@brief Returns a std::vector<double> of determinants of Matrix3.
            [[nodiscard]] std::vector<double> det() const;

            Matrix3& operator[](size_t index);
            const Matrix3& operator[](size_t index) const;
            Matrix3Array& operator+=(const Matrix3 &appendedMatrix3);

            ///@brief Generates a string object from all elements of encapsulated std::vector.
            explicit operator std::string() const;
    };
}


#endif //GENERAL_COIL_PROGRAM_Vector3_H
