#include "Tensor.h"

#include <sstream>


namespace vec3
{
    Matrix3Array::Matrix3Array() { this->matrixArray = {}; }

    Matrix3Array::Matrix3Array(size_t initSize) { this->matrixArray.resize(initSize); }

    Matrix3Array::Matrix3Array(const std::vector<Matrix3> &matrixArray)
    {
        std::copy(matrixArray.begin(), matrixArray.end(), this->matrixArray.begin());
    }


    void Matrix3Array::append(const Matrix3 &appendedMatrix3) { this->matrixArray.push_back(appendedMatrix3); }

    void Matrix3Array::append(double xx, double xy, double xz,
                                    double yx, double yy, double yz,
                                    double zx, double zy, double zz)
    {
        this->matrixArray.emplace_back(xx, xy, xz, yx, yy, yz, zx, zy, zz);
    }

    void Matrix3Array::reserve(size_t reserveSize) { this->matrixArray.reserve(reserveSize); }

    void Matrix3Array::resize(size_t newSize) { this->matrixArray.resize(newSize); }

    size_t Matrix3Array::size() const { return this->matrixArray.size(); }

    void Matrix3Array::clear() { this->matrixArray.clear(); }

    std::vector<Matrix3> & Matrix3Array::getStdVectorRef() { return this->matrixArray; }


    std::vector<double> Matrix3Array::xx() const
    {
        std::vector<double> outputArr;
        outputArr.reserve(this->matrixArray.size());

        for (const auto & i : this->matrixArray)
            outputArr.emplace_back(i.xx);

        return outputArr;
    }

    std::vector<double> Matrix3Array::xy() const
    {
        std::vector<double> outputArr;
        outputArr.reserve(this->matrixArray.size());

        for (const auto & i : this->matrixArray)
            outputArr.emplace_back(i.xy);

        return outputArr;
    }

    std::vector<double> Matrix3Array::xz() const
    {
        std::vector<double> outputArr;
        outputArr.reserve(this->matrixArray.size());

        for (const auto & i : this->matrixArray)
            outputArr.emplace_back(i.xz);

        return outputArr;
    }

    std::vector<double> Matrix3Array::yx() const
    {
        std::vector<double> outputArr;
        outputArr.reserve(this->matrixArray.size());

        for (const auto & i : this->matrixArray)
            outputArr.emplace_back(i.yx);

        return outputArr;
    }

    std::vector<double> Matrix3Array::yy() const
    {
        std::vector<double> outputArr;
        outputArr.reserve(this->matrixArray.size());

        for (const auto & i : this->matrixArray)
            outputArr.emplace_back(i.yy);

        return outputArr;
    }

    std::vector<double> Matrix3Array::yz() const
    {
        std::vector<double> outputArr;
        outputArr.reserve(this->matrixArray.size());

        for (const auto & i : this->matrixArray)
            outputArr.emplace_back(i.yz);

        return outputArr;
    }

    std::vector<double> Matrix3Array::zx() const
    {
        std::vector<double> outputArr;
        outputArr.reserve(this->matrixArray.size());

        for (const auto & i : this->matrixArray)
            outputArr.emplace_back(i.zx);

        return outputArr;
    }

    std::vector<double> Matrix3Array::zy() const
    {
        std::vector<double> outputArr;
        outputArr.reserve(this->matrixArray.size());

        for (const auto & i : this->matrixArray)
            outputArr.emplace_back(i.zy);

        return outputArr;
    }

    std::vector<double> Matrix3Array::zz() const
    {
        std::vector<double> outputArr;
        outputArr.reserve(this->matrixArray.size());

        for (const auto & i : this->matrixArray)
            outputArr.emplace_back(i.zz);

        return outputArr;
    }

    std::vector<double> Matrix3Array::det() const
    {
        std::vector<double> outputArr;
        outputArr.reserve(this->matrixArray.size());

        for (const auto & i : this->matrixArray)
            outputArr.emplace_back(i.det());

        return outputArr;
    }


    Matrix3 & Matrix3Array::operator[](size_t index)
    {
        return static_cast<Matrix3 &>(this->matrixArray[index]);
    }

    const Matrix3 & Matrix3Array::operator[](size_t index) const
    {
        return this->matrixArray[index];
    }

    void Matrix3Array::operator+=(const Matrix3 &appendedMatrix3)
    {
        this->matrixArray.push_back(appendedMatrix3);
    }

    Matrix3Array::operator std::string() const
    {
        std::stringstream output;

        auto stringifyVector = [](auto &ar) -> std::string
        {
            std::stringstream output;

            output << "[";

            for(int i = 0; i < ar.size(); i++)
            {
                if(i != 0)
                    output << ", ";
                output << std::string(ar[i]);
            }

            output << "]";

            return output.str();
        };

        output << "Matrix3Array("
               << "matrix_array=" << stringifyVector(this->matrixArray)
               << ")";

        return output.str();
    }
}
