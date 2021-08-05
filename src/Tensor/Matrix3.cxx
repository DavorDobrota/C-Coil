#include "Tensor.h"

vec3::Matrix3::Matrix3() : Matrix3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) {}

vec3::Matrix3::Matrix3(double xx, double xy, double xz, double yx, double yy, double yz, double zx, double zy, double zz) :
                        xxElement(xx), xyElement(xy), xzElement(xz),
                        yxElement(yx), yyElement(yy), yzElement(yz),
                        zxElement(zx), zyElement(zy), zzElement(zz) {}