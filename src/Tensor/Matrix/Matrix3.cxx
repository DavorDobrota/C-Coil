#include "Tensor.h"

#include <sstream>


namespace vec3
{
    Matrix3::Matrix3() : Matrix3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) {}

    Matrix3::Matrix3(double xx, double xy, double xz, double yx, double yy, double yz, double zx, double zy, double zz) :
            xx(xx), xy(xy), xz(xz),
            yx(yx), yy(yy), yz(yz),
            zx(zx), zy(zy), zz(zz) {}

    double Matrix3::det() const
    {
        return xx * (yy * zz - yz * zy) - xy * (yx * zz - yz * zx) + xz * (yx * zy - yy * zx);
    }


    Matrix3 Matrix3::operator+(const Matrix3 &mat) const
    {
        return Matrix3(this->xx + mat.xx, this->xy + mat.xy, this->xz + mat.xz,
                             this->yx + mat.yx, this->yy + mat.yy, this->yz + mat.yz,
                             this->zx + mat.zx, this->zy + mat.zy, this->zz + mat.zz);
    }

    Matrix3 Matrix3::operator+=(const Matrix3 &mat)
    {
        this->xx += mat.xx; this->xy += mat.xy; this->xz += mat.xz;
        this->yx += mat.yx; this->yy += mat.yy; this->yz += mat.yz;
        this->zx += mat.zx; this->zy += mat.zy; this->zz += mat.zz;
        return *this;
    }

    Matrix3 Matrix3::operator*(double multiplier) const
    {
        return Matrix3(this->xx * multiplier, this->xy * multiplier, this->xz * multiplier,
                             this->yx * multiplier, this->yy * multiplier, this->yz * multiplier,
                             this->zx * multiplier, this->zy * multiplier, this->zz * multiplier);
    }

    void Matrix3::operator*=(double multiplier)
    {
        this->xx *= multiplier; this->xy *= multiplier; this->xz *= multiplier;
        this->yx *= multiplier; this->yy *= multiplier; this->yz *= multiplier;
        this->zx *= multiplier; this->zy *= multiplier; this->zz *= multiplier;
    }

    Vector3 Matrix3::operator*(const Vector3 &vec) const
    {
        return Vector3(this->xx * vec.x + this->xy * vec.y + this->xz * vec.z,
                             this->yx * vec.x + this->yy * vec.y + this->yz * vec.z,
                             this->zx * vec.x + this->zy * vec.y + this->zz * vec.z);
    }

    Matrix3 Matrix3::operator*(const Matrix3 &mat) const
    {
        return Matrix3(this->xx * mat.xx + this->xy * mat.yx + this->xz * mat.zx,
                             this->xx * mat.xy + this->xy * mat.yy + this->xz * mat.zy,
                             this->xx * mat.xz + this->xy * mat.yz + this->xz * mat.zz,
                             this->yx * mat.xx + this->yy * mat.yx + this->yz * mat.zx,
                             this->yx * mat.xy + this->yy * mat.yy + this->yz * mat.zy,
                             this->yx * mat.xz + this->yy * mat.yz + this->yz * mat.zz,
                             this->zx * mat.xx + this->zy * mat.yx + this->zz * mat.zx,
                             this->zx * mat.xy + this->zy * mat.yy + this->zz * mat.zy,
                             this->zx * mat.xz + this->zy * mat.yz + this->zz * mat.zz);
    }

    Matrix3::operator std::string() const
    {
        std::stringstream output;

        output << "Matrix3([[" << xx << ", " << xy << ", " << xz << "], ["
               << yx << ", " << yy << ", " << yz << "], ["
               << zx << ", " << zy << ", " << zz << "]])";

        return output.str();
    }
}
