#ifndef PT_TRANSFORM_H
#define PT_TRANSFORM_H

#include <vector.h>
#include <matrix.h>

namespace ptracer {

class BBox3;
class Ray;

class Transform {
public:
    Transform() : mat(), invMat() {}
    Transform(const Mat4& matrix) : mat(matrix) {
        const auto inv = matrix.inverse();
        DCHECK(inv.has_value());
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
        return {mat * t.mat, t.invMat * invMat};
    }

    bool flipsOrientation() const { return Det3x3(mat) < 0; }

private:
    Mat4 mat;
    Mat4 invMat;
};

inline Transform Inverse(const Transform& tform) {
    return {tform.invMatrix(), tform.matrix()};
}

Transform Translate(const Vec3& tl);
Transform Scale(Float sX, Float sY, Float sZ);
Transform RotateX(Float degs);
Transform RotateY(Float degs);
Transform RotateZ(Float degs);
Transform Ortho(Float zNear, Float zFar);
Transform LookAt(const Point3& pos, const Point3& at, const Vec3& up);

} // namespace ptracer

#endif // PT_TRANSFORM_H