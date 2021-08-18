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

vec3::Matrix3 vec3::Matrix3::operator*(const Matrix3 &mat) const
{
    return vec3::Matrix3(this->xxElement * mat.xxElement + this->xyElement * mat.yxElement + this->xzElement * mat.zxElement,
                         this->xxElement * mat.xyElement + this->xyElement * mat.yyElement + this->xzElement * mat.zyElement,
                         this->xxElement * mat.xzElement + this->xyElement * mat.yzElement + this->xzElement * mat.zzElement,
                         this->yxElement * mat.xxElement + this->yyElement * mat.yxElement + this->yzElement * mat.zxElement,
                         this->yxElement * mat.xyElement + this->yyElement * mat.yyElement + this->yzElement * mat.zyElement,
                         this->yxElement * mat.xzElement + this->yyElement * mat.yzElement + this->yzElement * mat.zzElement,
                         this->zxElement * mat.xxElement + this->zyElement * mat.yxElement + this->zzElement * mat.zxElement,
                         this->zxElement * mat.xyElement + this->zyElement * mat.yyElement + this->zzElement * mat.zyElement,
                         this->zxElement * mat.xzElement + this->zyElement * mat.yzElement + this->zzElement * mat.zzElement);
}
