#ifndef __FRAME_H__
#define __FRAME_H__

#include <vector.h>

namespace ptracer {

class Frame {
public:
    Frame() = default;
    Frame(const Normal& norm);
    Frame(const Vec3& x, const Vec3& y, const Vec3& z);

    Vec3 toWorld(const Vec3& vec) const { return _x * vec.x + _y * vec.y + _z * vec.z; }
    Vec3 toLocal(const Vec3& vec) const {
        Float dotU = dot(_x, vec);
        Float dotV = dot(_y, vec);
        Float dotN = dot(_z, vec);

        return {dotU, dotV, dotN};
    }

    Normal normal() const { return Normal(_z); }

    const Vec3& x() const { return _x; }
    const Vec3& y() const { return _y; }
    const Vec3& z() const { return _z; }

    bool orthonormal();

    static Vec3 reflect(const Vec3& wi) { return {-wi.x, -wi.y, wi.z}; }
    static Vec3 reflect(const Vec3& wi, const Vec3& n) { return 2 * dot(wi, n) * n - wi; }

    static Vec3 refract(const Vec3& wi, Float eta, Float cosT) {
        return {-eta * wi.x, -eta * wi.y, cosT};
    }
    static Vec3 refract(const Vec3& wi, Float intEta, Float extEta, Float cosT);
    static Vec3 refract(const Vec3& wi, const Normal& n, Float eta, Float cosT);
    static Vec3 refract(const Vec3& wi, const Normal& n, Float eta);
    // static bool refract(const Vec3& wi, Vec3* wt, const Normal& n, Float eta);

    static bool sameSide(const Vec3& w1, const Vec3& w2);
    static bool isPosHemisphere(const Vec3& w);

    // Frame spherical functions
    static Float cosTheta(const Vec3& w) { return w.z; }
    static Float cosThetaSqr(const Vec3& w) { return w.z * w.z; }
    static Float absCosTheta(const Vec3& w) { return std::abs(w.z); }

    static Float sinTheta(const Vec3& w) { return std::sqrt(sinThetaSqr(w)); }
    static Float sinThetaSqr(const Vec3& w) {
        return std::max((Float)0.0, (Float)1.0 - cosThetaSqr(w));
    }

    static Float tanTheta(const Vec3& w) { return sinTheta(w) / cosTheta(w); }
    static Float tanThetaSqr(const Vec3& w) { return sinThetaSqr(w) / cosThetaSqr(w); }

    static Float cosPhi(const Vec3& w) {
        Float sin = sinTheta(w);
        return sin == 0 ? 1 : math::clamp(w.x / sin, -1.0, 1.0);
    }
    static Float sinPhi(const Vec3& w) {
        Float sin = sinTheta(w);
        return sin == 0 ? 0 : math::clamp(w.y / sin, -1.0, 1.0);
    }

    static Float cosPhiSqr(const Vec3& w) {
        Float cos = cosPhi(w);
        return cos * cos;
    }
    static Float sinPhiSqr(const Vec3& w) {
        Float sin = sinPhi(w);
        return sin * sin;
    }

    static Float cosAng(const Vec3& w1, const Vec3& w2);

private:
    Vec3 _x{1, 0, 0};
    Vec3 _y{0, 1, 0};
    Vec3 _z{0, 0, 1};
};

}; // namespace ptracer

#endif // __FRAME_H__