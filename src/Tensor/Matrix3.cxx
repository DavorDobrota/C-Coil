#include "Tensor.h"

vec3::Matrix3::Matrix3() : Matrix3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) {}

vec3::Matrix3::Matrix3(double xx, double xy, double xz, double yx, double yy, double yz, double zx, double zy, double zz) :
                        xxElement(xx), xyElement(xy), xzElement(xz),
                        yxElement(yx), yyElement(yy), yzElement(yz),
                        zxElement(zx), zyElement(zy), zzElement(zz) {}


vec3::FieldVector3 vec3::Matrix3::operator*(const FieldVector3 &vec) const
{
    return vec3::FieldVector3(this->xxElement * vec.xComponent + this->xyElement * vec.yComponent + this->xzElement * vec.zComponent,
                              this->yxElement * vec.xComponent + this->yyElement * vec.yComponent + this->yzElement * vec.zComponent,
                              this->zxElement * vec.xComponent + this->zyElement * vec.yComponent + this->zzElement * vec.zComponent);
}