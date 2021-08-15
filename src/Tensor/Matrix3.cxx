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

vec3::Matrix3 vec3::Matrix3::matrixMultiplication(vec3::Matrix3 mat1, vec3::Matrix3 mat2)
{
    return vec3::Matrix3(mat1.xxElement * mat2.xxElement + mat1.xyElement * mat2.yxElement + mat1.xzElement * mat2.zxElement,
                         mat1.xxElement * mat2.xyElement + mat1.xyElement * mat2.yyElement + mat1.xzElement * mat2.zyElement,
                         mat1.xxElement * mat2.xzElement + mat1.xyElement * mat2.yzElement + mat1.xzElement * mat2.zzElement,
                         mat1.yxElement * mat2.xxElement + mat1.yyElement * mat2.yxElement + mat1.yzElement * mat2.zxElement,
                         mat1.yxElement * mat2.xyElement + mat1.yyElement * mat2.yyElement + mat1.yzElement * mat2.zyElement,
                         mat1.yxElement * mat2.xzElement + mat1.yyElement * mat2.yzElement + mat1.yzElement * mat2.zzElement,
                         mat1.zxElement * mat2.xxElement + mat1.zyElement * mat2.yxElement + mat1.zzElement * mat2.zxElement,
                         mat1.zxElement * mat2.xyElement + mat1.zyElement * mat2.yyElement + mat1.zzElement * mat2.zyElement,
                         mat1.zxElement * mat2.xzElement + mat1.zyElement * mat2.yzElement + mat1.zzElement * mat2.zzElement);
}
