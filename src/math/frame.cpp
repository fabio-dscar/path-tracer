#include <frame.h>

using namespace ptracer;

Frame::Frame(const Normal& normal) {
    _z = normalize(Vec3(normal));
    math::basisFromVector(_z, _x, _y);
    _x.normalize();
    _y.normalize();

    // assert(isConsistent());
}

Frame::Frame(const Vec3& x, const Vec3& y, const Vec3& z) {
    _x = normalize(x);
    _y = normalize(y);
    _z = normalize(z);

    // assert(isConsistent());
}

bool Frame::orthonormal() {
    Float lenX = _x.length();
    Float lenY = _y.length();
    Float lenZ = _z.length();

    if (lenX > (1 + F_EPSILON) || lenY > (1 + F_EPSILON) || lenZ > (1 + F_EPSILON))
        return false;

    Float dot12 = absDot(_x, _y);
    Float dot13 = absDot(_x, _z);
    Float dot23 = absDot(_y, _z);
    return (dot12 < F_EPSILON && dot13 < F_EPSILON && dot23 < F_EPSILON);
}

Vec3 Frame::refract(const Vec3& wi, Float intEta, Float extEta, Float cosT) {
    Float eta = extEta / intEta;
    if (cosT > 0) // If we are leaving the surface, swap IORs
        eta = 1.0 / eta;

    return {-eta * wi.x, -eta * wi.y, cosT};
}

Vec3 Frame::refract(const Vec3& wi, const Normal& n, Float eta, Float cosT) {
    return Vec3(eta * -wi + n * (dot(wi, n) * eta + cosT));
}

Vec3 Frame::refract(const Vec3& wi, const Normal& n, Float eta) {
    Float cosI = dot(n, wi);
    Float sin2I = std::max((Float)0.0, (Float)1.0 - cosI * cosI);
    Float sin2T = eta * eta * sin2I;

    // Handle TIR
    if (sin2T >= 1)
        return Vec3(0);

    Float cosT = std::sqrt(1 - sin2T);

    return eta * -wi + (eta * cosI - cosT) * Vec3(n);
}

bool Frame::sameSide(const Vec3& w1, const Vec3& w2) {
    return cosTheta(w1) * cosTheta(w2) > 0;
}

bool Frame::isPosHemisphere(const Vec3& w) {
    return cosTheta(w) > 0;
}

Float Frame::cosAng(const Vec3& w1, const Vec3& w2) {
    // Project w1 and w2
    Vec2 a{w1.x, w1.y};
    Vec2 b{w2.x, w2.y};

    // Use dot product definition
    return math::clamp(math::dot(a, b) / (a.length() * b.length()), -1, 1);
}