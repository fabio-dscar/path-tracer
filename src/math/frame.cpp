#include <frame.h>

using namespace ptracer;

Frame::Frame(const Normal& normal) {
    z = Normalize(normal);

    std::tie(x, y) = BasisFromVector(z);
    x.normalize();
    y.normalize();
}

Frame::Frame(const Vec3& x, const Vec3& y, const Normal& z)
    : x(Normalize(x)), y(Normalize(y)), z(Normalize(z)) {}

bool Frame::orthonormal() {
    Float lenX = x.length();
    Float lenY = y.length();
    Float lenZ = z.length();

    if (lenX > (1 + Epsilon) || lenY > (1 + Epsilon) || lenZ > (1 + Epsilon))
        return false;

    Float dot12 = AbsDot(x, y);
    Float dot13 = AbsDot(x, z);
    Float dot23 = AbsDot(y, z);
    return (dot12 < Epsilon && dot13 < Epsilon && dot23 < Epsilon);
}

Vec3 Frame::Refract(const Vec3& wi, Float intEta, Float extEta, Float cosT) {
    // Swap IORs if we are leaving surface (cos > 0)
    const Float eta = cosT > 0 ? intEta / extEta : extEta / intEta;
    return {-eta * wi.x, -eta * wi.y, cosT};
}

Vec3 Frame::Refract(const Vec3& wi, const Normal& n, Float eta, Float cosT) {
    return Vec3{eta * -wi + n * (Dot(wi, n) * eta + cosT)};
}

Vec3 Frame::Refract(const Vec3& wi, const Normal& n, Float eta) {
    const Float cosI  = Dot(n, wi);
    const Float sin2I = math::Max(0.0, 1.0 - cosI * cosI);
    const Float sin2T = eta * eta * sin2I;

    // Handle TIR
    if (sin2T >= 1)
        return Vec3{0};

    const Float cosT = std::sqrt(1 - sin2T);
    return eta * -wi + (eta * cosI - cosT) * Vec3{n};
}

bool Frame::SameSide(const Vec3& w1, const Vec3& w2) {
    return CosTheta(w1) * CosTheta(w2) > 0;
}

bool Frame::IsPosHemisphere(const Vec3& w) {
    return CosTheta(w) > 0;
}

Float Frame::CosAng(const Vec3& w1, const Vec3& w2) {
    // Project w1 and w2
    const Vec2 a{w1.x, w1.y};
    const Vec2 b{w2.x, w2.y};

    return math::Clamp(Dot(a, b) / (a.length() * b.length()), -1, 1);
}