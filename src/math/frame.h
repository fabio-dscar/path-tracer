#ifndef __FRAME_H__
#define __FRAME_H__

#include <vector.h>

namespace ptracer {

class Frame {
public:
    Vec3 x{1, 0, 0};
    Vec3 y{0, 1, 0};
    Normal z{0, 0, 1};

    Frame() = default;
    Frame(const Normal& normal);
    Frame(const Vec3& x, const Vec3& y, const Normal& z);

    Vec3 toWorld(const Vec3& vec) const { return x * vec.x + y * vec.y + z * vec.z; }

    Vec3 toLocal(const Vec3& vec) const {
        Float dotU = Dot(x, vec);
        Float dotV = Dot(y, vec);
        Float dotN = Dot(z, vec);

        return {dotU, dotV, dotN};
    }

    const Normal& normal() const { return z; }

    bool orthonormal();

    // Static member functions
    static Vec3 Reflect(const Vec3& wi) { return {-wi.x, -wi.y, wi.z}; }
    static Vec3 Reflect(const Vec3& wi, const Normal& n) {
        return Vec3{2.0 * Dot(wi, n) * n - wi};
    }

    static Vec3 Refract(const Vec3& wi, Float eta, Float cosT) {
        return {-eta * wi.x, -eta * wi.y, cosT};
    }
    static Vec3 Refract(const Vec3& wi, Float intEta, Float extEta, Float cosT);
    static Vec3 Refract(const Vec3& wi, const Normal& n, Float eta, Float cosT);
    static Vec3 Refract(const Vec3& wi, const Normal& n, Float eta);
    // static bool refract(const Vec3& wi, Vec3* wt, const Normal& n, Float eta);

    static bool SameSide(const Vec3& w1, const Vec3& w2);
    static bool IsPosHemisphere(const Vec3& w);

    // Frame spherical functions
    static Float CosTheta(const Vec3& w) { return w.z; }
    static Float CosThetaSqr(const Vec3& w) { return w.z * w.z; }
    static Float AbsCosTheta(const Vec3& w) { return std::abs(w.z); }

    static Float SinTheta(const Vec3& w) { return std::sqrt(SinThetaSqr(w)); }
    static Float SinThetaSqr(const Vec3& w) {
        return math::Max(0.0f, 1.0f - CosThetaSqr(w));
    }

    static Float TanTheta(const Vec3& w) { return SinTheta(w) / CosTheta(w); }
    static Float TanThetaSqr(const Vec3& w) { return SinThetaSqr(w) / CosThetaSqr(w); }

    static Float CosPhi(const Vec3& w) {
        Float sin = SinTheta(w);
        return sin == 0 ? 1 : math::Clamp(w.x / sin, -1.0, 1.0);
    }
    static Float SinPhi(const Vec3& w) {
        Float sin = SinTheta(w);
        return sin == 0 ? 0 : math::Clamp(w.y / sin, -1.0, 1.0);
    }

    static Float CosPhiSqr(const Vec3& w) {
        Float cos = CosPhi(w);
        return cos * cos;
    }
    static Float SinPhiSqr(const Vec3& w) {
        Float sin = SinPhi(w);
        return sin * sin;
    }

    static Float CosAng(const Vec3& w1, const Vec3& w2);
};

};     // namespace ptracer

#endif // __FRAME_H__