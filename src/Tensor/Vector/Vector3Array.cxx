#include "Tensor.h"

#include <sstream>


namespace vec3
{
    Vector3Array::Vector3Array() { this->vectorArray = {}; }

    Vector3Array::Vector3Array(size_t initSize) { this->vectorArray.resize(initSize); }

    Vector3Array::Vector3Array(const std::vector<Vector3> &vectorArray)
    {
        std::copy(vectorArray.begin(), vectorArray.end(), this->vectorArray.begin());
    }


    void Vector3Array::append(const Vector3 &appendedVector3) { this->vectorArray.push_back(appendedVector3); }

    void Vector3Array::append(double x, double y, double z) { this->vectorArray.emplace_back(x, y, z); }

    void Vector3Array::reserve(size_t reserveSize) { this->vectorArray.reserve(reserveSize); }

    void Vector3Array::resize(size_t newSize) { this->vectorArray.resize(newSize); }

    size_t Vector3Array::size() const { return this->vectorArray.size(); }

    void Vector3Array::clear() { this->vectorArray.clear(); }

    std::vector<Vector3> & Vector3Array::getStdVectorRef() { return this->vectorArray; }


    std::vector<double> Vector3Array::x() const
    {
        std::vector<double> outputArr;
        outputArr.reserve(this->vectorArray.size());

        for (const auto & i : this->vectorArray)
            outputArr.emplace_back(i.x);

        return outputArr;
    }

    std::vector<double> Vector3Array::y() const
    {
        std::vector<double> outputArr;
        outputArr.reserve(this->vectorArray.size());

        for (const auto & i : this->vectorArray)
            outputArr.emplace_back(i.y);

        return outputArr;
    }

    std::vector<double> Vector3Array::z() const
    {
        std::vector<double> outputArr;
        outputArr.reserve(this->vectorArray.size());

        for (const auto & i : this->vectorArray)
            outputArr.emplace_back(i.z);

        return outputArr;
    }

    std::vector<double> Vector3Array::abs() const
    {
        std::vector<double> outputArr;
        outputArr.reserve(this->vectorArray.size());

        for (const auto & i : this->vectorArray)
            outputArr.emplace_back(i.abs());

        return outputArr;
    }

    Vector3 & Vector3Array::operator[](size_t index)
    {
        return static_cast<Vector3 &>(this->vectorArray[index]);
    }

    const Vector3 &Vector3Array::operator[](size_t index) const
    {
        return this->vectorArray[index];
    }

    void Vector3Array::operator+=(const Vector3 &appendedVector3)
    {
        this->vectorArray.push_back(appendedVector3);
    }


    Vector3Array::operator std::string() const
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

        output << "Vector3Array("
               << "vector_array=" << stringifyVector(this->vectorArray)
               << ")";

        return output.str();
    }
}
