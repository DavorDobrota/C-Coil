#include "Tensor.h"

vec3::Matrix3::Matrix3() : Matrix3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) {}

vec3::Matrix3::Matrix3(double xx, double xy, double xz, double yx, double yy, double yz, double zx, double zy, double zz) :
                        xxElement(xx), xyElement(xy), xzElement(xz),
                        yxElement(yx), yyElement(yy), yzElement(yz),
                        zxElement(zx), zyElement(zy), zzElement(zz) {}

vec3::FieldVector3 vec3::Matrix3::matrixVectorMultiplication(vec3::Matrix3 matrix, vec3::FieldVector3 vector)
{
    return vec3::FieldVector3(matrix.xxElement * vector.xComponent + matrix.xyElement * vector.yComponent + matrix.xzElement * vector.zComponent,
                              matrix.yxElement * vector.xComponent + matrix.yyElement * vector.yComponent + matrix.yzElement * vector.zComponent,
                              matrix.zxElement * vector.xComponent + matrix.zyElement * vector.yComponent + matrix.zzElement * vector.zComponent);
}