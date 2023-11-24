#include <transform.h>

#include <bounds.h>

using namespace ptracer;
using namespace ptracer::math;

/*Ray Transform::operator()(const Ray& ray) const {
    Point3 orig = (*this)(ray.origin());
    Vec3 dir = (*this)(ray.dir());

    Ray r = Ray(orig, dir, ray.minT(), ray.maxT());
    r.setTime(ray.time());

    return r;
}*/

// clang-format off
Bounds3 Transform::operator()(const Bounds3& box) const {
    const Transform& Tr = *this;

    Bounds3 ret({0, 0, 0});
    ret.expand(Tr(Point3(box[0].x, box[0].y, box[0].z)));
    ret.expand(Tr(Point3(box[1].x, box[0].y, box[0].z)));
    ret.expand(Tr(Point3(box[0].x, box[1].y, box[0].z)));
    ret.expand(Tr(Point3(box[0].x, box[0].y, box[1].z)));
    ret.expand(Tr(Point3(box[0].x, box[1].y, box[1].z)));
    ret.expand(Tr(Point3(box[1].x, box[1].y, box[0].z)));
    ret.expand(Tr(Point3(box[1].x, box[0].y, box[1].z)));
    ret.expand(Tr(Point3(box[1].x, box[1].y, box[1].z)));

    return ret;
}

Point3 Transform::operator()(const Point3& pt) const {
    Float x = mat.m[0][0] * pt.x + mat.m[0][1] * pt.y + mat.m[0][2] * pt.z + mat.m[0][3];
    Float y = mat.m[1][0] * pt.x + mat.m[1][1] * pt.y + mat.m[1][2] * pt.z + mat.m[1][3];
    Float z = mat.m[2][0] * pt.x + mat.m[2][1] * pt.y + mat.m[2][2] * pt.z + mat.m[2][3];
    Float w = mat.m[3][0] * pt.x + mat.m[3][1] * pt.y + mat.m[3][2] * pt.z + mat.m[3][3];

    if (w == 1)
        return Point3(x, y, z);
    else
        return Point3(x, y, z) / w;
}

Vec3 Transform::operator()(const Vec3& vec) const {
    Float x = mat.m[0][0] * vec.x + mat.m[0][1] * vec.y + mat.m[0][2] * vec.z;
    Float y = mat.m[1][0] * vec.x + mat.m[1][1] * vec.y + mat.m[1][2] * vec.z;
    Float z = mat.m[2][0] * vec.x + mat.m[2][1] * vec.y + mat.m[2][2] * vec.z;

    return Vec3(x, y, z);
}

Normal Transform::operator()(const Normal& n) const {
    Float x = invMat.m[0][0] * n.x + invMat.m[1][0] * n.y + invMat.m[2][0] * n.z;
    Float y = invMat.m[0][1] * n.x + invMat.m[1][1] * n.y + invMat.m[2][1] * n.z;
    Float z = invMat.m[0][2] * n.x + invMat.m[1][2] * n.y + invMat.m[2][2] * n.z;

    return Normal(x, y, z);
}

Transform ptracer::translate(const Vec3& trans) {
    Mat4 mat(1, 0, 0, trans.x,
             0, 1, 0, trans.y,
             0, 0, 1, trans.z,
             0, 0, 0, 1);

    Mat4 invMat(1, 0, 0, -trans.x,
                0, 1, 0, -trans.y,
                0, 0, 1, -trans.z,
                0, 0, 0, 1);

    return Transform(mat, invMat);
}

Transform ptracer::scale(Float sX, Float sY, Float sZ) {
    Mat4 mat(sX, 0, 0, 0,
             0, sY, 0, 0,
             0, 0, sZ, 0,
             0, 0, 0, 1);

    Mat4 invMat(1 / sX, 0, 0, 0,
                0, 1 / sY, 0, 0,
                0, 0, 1 / sZ, 0,
                0, 0, 0, 1);

    return Transform(mat, invMat);
}

Transform ptracer::rotateX(Float degrees) {
    Float rads = radians(degrees);
    Float sin = std::sin(rads);
    Float cos = std::cos(rads);

    Mat4 mat(1, 0, 0, 0,
             0, cos, -sin, 0,
             0, sin, cos, 0,
             0, 0, 0, 1);

    // The inverse of rotation matrices is their own transpose 
    return Transform(mat, transpose(mat));
}

Transform ptracer::rotateY(Float degrees) {
    Float rads = radians(degrees);
    Float sin = std::sin(rads);
    Float cos = std::cos(rads);

    Mat4 mat(cos, 0, sin, 0,
             0, 1, 0, 0,
             -sin, 0, cos, 0,
             0, 0, 0, 1);

    // The inverse of rotation matrices is their own transpose 
    return Transform(mat, transpose(mat));
}

Transform ptracer::rotateZ(Float degrees) {
    Float rads = radians(degrees);
    Float sin = std::sin(rads);
    Float cos = std::cos(rads);

    Mat4 mat(cos, -sin, 0, 0,
             sin, cos, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 1);

    // The inverse of rotation matrices is their own transpose 
    return Transform(mat, transpose(mat));
}
// clang-format off

Transform ptracer::lookAt(const Point3& pos, const Point3& at, const Vec3& up) {
    Mat4 lookAtMat;

    Vec3 n = -normalize(pos - at);
    Vec3 u = normalize(cross(up + Vec3(0.05, 0, 0.05), n));
    Vec3 v = cross(n, u);

    lookAtMat(0, 0) = u.x;
    lookAtMat(1, 0) = u.y;
    lookAtMat(2, 0) = u.z;
    lookAtMat(3, 0) = 0;

    lookAtMat(0, 1) = v.x;
    lookAtMat(1, 1) = v.y;
    lookAtMat(2, 1) = v.z;
    lookAtMat(3, 1) = 0;

    lookAtMat(0, 2) = n.x;
    lookAtMat(1, 2) = n.y;
    lookAtMat(2, 2) = n.z;
    lookAtMat(3, 2) = 0;

    lookAtMat(0, 3) = pos.x;
    lookAtMat(1, 3) = pos.y;
    lookAtMat(2, 3) = pos.z;
    lookAtMat(3, 3) = 1;

    return Transform(lookAtMat);
}

Transform ptracer::ortho(Float zNear, Float zFar) {
    return scale(1, 1, 1 / (zFar - zNear)) * translate(Vec3(0, 0, -zNear));
}
