#include <quat.h>
#include <matrix.h>
#include <transform.h>

using namespace ptracer;
using namespace ptracer::math;

Quat Quat::operator+(const Quat& q) const {
    return {x + q.x, y + q.y, z + q.z, w + q.w};
}

Quat& Quat::operator+=(const Quat& quat) {
    x += quat.x;
    y += quat.y;
    z += quat.z;
    w += quat.w;
    return *this;
}

Quat Quat::operator-(const Quat& q) const {
    return {x - q.x, y - q.y, z - q.z, w - q.w};
}

Quat& Quat::operator-=(const Quat& quat) {
    x -= quat.x;
    y -= quat.y;
    z -= quat.z;
    w -= quat.w;
    return *this;
}

Quat Quat::operator*(Float s) const {
    return {x * s, y * s, z * s, w * s};
}

Quat& Quat::operator*=(Float s) {
    x *= s;
    y *= s;
    z *= s;
    w *= s;
    return *this;
}

Quat Quat::operator/(Float s) const {
    assert(s != 0);
    Float recip = 1.0 / s;
    return {x * recip, y * recip, z * recip, w * recip};
}

Quat& Quat::operator/=(Float s) {
    assert(s != 0);
    Float recip = 1.0 / s;
    x *= recip;
    y *= recip;
    z *= recip;
    w *= recip;
    return *this;
}

Transform Quat::transform() const {
    Mat4 mat;

    // Diagonal
    mat.m[0][0] = 1.0 - 2.0 * (y * y + z * z);
    mat.m[1][1] = 1.0 - 2.0 * (x * x + z * z);
    mat.m[2][2] = 1.0 - 2.0 * (x * x + y * y);
    mat.m[3][3] = 1.0;

    mat.m[0][1] = 2.0 * (x * y + z * w);
    mat.m[0][2] = 2.0 * (x * z - y * w);
    mat.m[0][3] = 0;

    mat.m[1][0] = 2.0 * (x * y - z * w);
    mat.m[1][2] = 2.0 * (y * z + x * w);
    mat.m[1][3] = 0;

    mat.m[2][0] = 2.0 * (x * z + y * w);
    mat.m[2][1] = 2.0 * (y * z - x * w);
    mat.m[2][3] = 0;

    mat.m[3][0] = 0;
    mat.m[3][1] = 0;
    mat.m[3][2] = 0;

    return {transpose(mat), mat};
}