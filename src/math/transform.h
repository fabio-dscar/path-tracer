#ifndef __TRANSFORM_H__
#define __TRANSFORM_H__

#include <vector.h>
#include <matrix.h>

namespace ptracer {

class BBox3;
class Ray;

class Transform {
public:
    Transform() : mat(), invMat() {}
    Transform(const Mat4& matrix) : mat(matrix) {
        auto inv = matrix.inverse();
        assert(inv.has_value());
        invMat = inv.value();
    }
    Transform(const Mat4& matrix, const Mat4& invMatrix)
        : mat(matrix), invMat(invMatrix) {}

    const Mat4& matrix() const { return mat; }
    const Mat4& invMatrix() const { return invMat; }

    Vec3 operator()(const Vec3& vec) const;
    Point3 operator()(const Point3& pt) const;
    Normal operator()(const Normal& norm) const;
    Ray operator()(const Ray& ray) const;
    BBox3 operator()(const BBox3& box) const;

    Transform operator*(const Transform& t) const {
        return {mul(mat, t.mat), mul(t.invMat, invMat)};
    }

    bool flipsOrientation() const { return det3x3(mat) < 0; }

private:
    Mat4 mat;
    Mat4 invMat;
};

inline Transform inverse(const Transform& trans) {
    return Transform(trans.invMatrix(), trans.matrix());
}

Transform translate(const Vec3& trans);
Transform scale(Float sX, Float sY, Float sZ);
Transform rotateX(Float degrees);
Transform rotateY(Float degrees);
Transform rotateZ(Float degrees);
Transform ortho(Float zNear, Float zFar);
Transform lookAt(const Point3& pos, const Point3& at, const Vec3& up);

} // namespace ptracer

#endif // __TRANSFORM_H__