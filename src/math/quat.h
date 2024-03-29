#ifndef PT_QUAT_H
#define PT_QUAT_H

#include <math.h>
#include <vector.h>

namespace ptracer {

class Transform;

namespace math {

class Quat {
public:
    Float x{0}, y{0}, z{0};
    Float w{1};

    Quat() = default;
    Quat(Float x, Float y, Float z, Float w) : x(x), y(y), z(z), w(w) {}
    Quat(const Vec3& vec, Float w) : x(vec.x), y(vec.y), z(vec.z), w(w) {}
    // Quat(const Transform& transform);

    Quat operator+(const Quat& quat) const;
    Quat& operator+=(const Quat& quat);

    Quat operator-(const Quat& quat) const;
    Quat& operator-=(const Quat& quat);

    Quat operator*(Float scalar) const;
    Quat& operator*=(Float scalar);

    Quat operator/(Float scalar) const;
    Quat& operator/=(Float scalar);

    Transform transform() const;
};

inline Quat operator*(Float scalar, const Quat& quat) {
    return quat * scalar;
}

inline Float Dot(const Quat& q1, const Quat& q2) {
    return q1.x * q2.x + q1.y * q2.y + q1.z * q2.z + q1.w * q2.w;
}

inline Quat Normalize(const Quat& q) {
    const Float len = std::sqrt(Dot(q, q));
    DCHECK_NE(len, 0);
    return q / len;
}

Quat Slerp(Float t, const Quat& q1, const Quat& q2);

} // namespace math
} // namespace ptracer

#endif // PT_QUAT_H