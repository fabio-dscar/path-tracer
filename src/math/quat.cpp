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
    DCHECK_NE(s, 0);
    const Float recip = 1.0 / s;
    return {x * recip, y * recip, z * recip, w * recip};
}

Quat& Quat::operator/=(Float s) {
    DCHECK_NE(s, 0);
    const Float recip = 1.0 / s;
    x *= recip;
    y *= recip;
    z *= recip;
    w *= recip;
    return *this;
}

Transform Quat::transform() const {
    Mat4 mat;

    // Diagonal
    mat(0, 0) = 1.0 - 2.0 * (y * y + z * z);
    mat(1, 1) = 1.0 - 2.0 * (x * x + z * z);
    mat(2, 2) = 1.0 - 2.0 * (x * x + y * y);
    mat(3, 3) = 1.0;

    mat(0, 1) = 2.0 * (x * y + z * w);
    mat(0, 2) = 2.0 * (x * z - y * w);
    mat(0, 3) = 0;

    mat(1, 0) = 2.0 * (x * y - z * w);
    mat(1, 2) = 2.0 * (y * z + x * w);
    mat(1, 3) = 0;

    mat(2, 0) = 2.0 * (x * z + y * w);
    mat(2, 1) = 2.0 * (y * z - x * w);
    mat(2, 3) = 0;

    mat(3, 0) = 0;
    mat(3, 1) = 0;
    mat(3, 2) = 0;

    return {Transpose(mat), mat};
}

Quat math::Slerp(Float t, const Quat& q1, const Quat& q2) {
    const Float dotp = Dot(q1, q2);
    if (dotp > 0.999) {
        return Normalize((1.0 - t) * q1 + t * q2);
    } else {
        const Float clampCos = Clamp(dotp, -1, 1);
        const Float theta    = SafeAcos(clampCos);
        const Quat aux       = Normalize(q2 - q1 * dotp);

        return aux * std::sin(theta * t) + q1 * std::cos(theta * t);
    }
}