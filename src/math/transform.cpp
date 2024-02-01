#include <transform.h>

#include <bbox.h>
#include <ray.h>

using namespace ptracer;
using namespace ptracer::math;

Ray Transform::operator()(const Ray& ray) const {
    const Point3 orig = (*this)(ray.o);
    const Vec3 dir    = (*this)(ray.dir);

    // Ray r  = Ray(orig, dir, ray.minT, ray.maxT);
    Ray r  = Ray{orig, dir};
    r.time = ray.time;

    return r;
}

// clang-format off
BBox3 Transform::operator()(const BBox3& box) const {
    const Transform& Tr = *this;

    BBox3 bbox{{0, 0, 0}};
    bbox.expand(Tr(Point3(box[0].x, box[0].y, box[0].z)));
    bbox.expand(Tr(Point3(box[1].x, box[0].y, box[0].z)));
    bbox.expand(Tr(Point3(box[0].x, box[1].y, box[0].z)));
    bbox.expand(Tr(Point3(box[0].x, box[0].y, box[1].z)));
    bbox.expand(Tr(Point3(box[0].x, box[1].y, box[1].z)));
    bbox.expand(Tr(Point3(box[1].x, box[1].y, box[0].z)));
    bbox.expand(Tr(Point3(box[1].x, box[0].y, box[1].z)));
    bbox.expand(Tr(Point3(box[1].x, box[1].y, box[1].z)));

    return bbox;
}

Point3 Transform::operator()(const Point3& pt) const {
    const Float x = mat(0, 0) * pt.x + mat(0, 1) * pt.y + mat(0, 2) * pt.z + mat(0, 3);
    const Float y = mat(1, 0) * pt.x + mat(1, 1) * pt.y + mat(1, 2) * pt.z + mat(1, 3);
    const Float z = mat(2, 0) * pt.x + mat(2, 1) * pt.y + mat(2, 2) * pt.z + mat(2, 3);
    const Float w = mat(3, 0) * pt.x + mat(3, 1) * pt.y + mat(3, 2) * pt.z + mat(3, 3);

    DCHECK_NE(w, 0);

    if (w == 1)
        return Point3{x, y, z};
    else
        return Point3{x, y, z} / w;
}

Vec3 Transform::operator()(const Vec3& vec) const {
    const Float x = mat(0, 0) * vec.x + mat(0, 1) * vec.y + mat(0, 2) * vec.z;
    const Float y = mat(1, 0) * vec.x + mat(1, 1) * vec.y + mat(1, 2) * vec.z;
    const Float z = mat(2, 0) * vec.x + mat(2, 1) * vec.y + mat(2, 2) * vec.z;

    return Vec3{x, y, z};
}

Normal Transform::operator()(const Normal& n) const {
    const Float x = invMat(0, 0) * n.x + invMat(1, 0) * n.y + invMat(2, 0) * n.z;
    const Float y = invMat(0, 1) * n.x + invMat(1, 1) * n.y + invMat(2, 1) * n.z;
    const Float z = invMat(0, 2) * n.x + invMat(1, 2) * n.y + invMat(2, 2) * n.z;

    return Normal{x, y, z};
}

Transform ptracer::Translate(const Vec3& tform) {
    Mat4 mat{1, 0, 0, tform.x,
             0, 1, 0, tform.y,
             0, 0, 1, tform.z,
             0, 0, 0, 1};

    Mat4 inv{1, 0, 0, -tform.x,
             0, 1, 0, -tform.y,
             0, 0, 1, -tform.z,
             0, 0, 0, 1};

    return {mat, inv};
}

Transform ptracer::Scale(Float sX, Float sY, Float sZ) {
    DCHECK_NE(sX, 0);
    DCHECK_NE(sY, 0);
    DCHECK_NE(sZ, 0);

    Mat4 mat{sX,  0,  0, 0,
              0, sY,  0, 0,
              0,  0, sZ, 0,
              0,  0,  0, 1};

    Mat4 inv{1 / sX, 0, 0, 0,
             0, 1 / sY, 0, 0,
             0, 0, 1 / sZ, 0,
             0, 0, 0, 1};

    return {mat, inv};
}

Transform ptracer::RotateX(Float degrees) {
    const Float rads = Radians(degrees);
    const Float sin  = std::sin(rads);
    const Float cos  = std::cos(rads);

    Mat4 mat{1,   0,   0,  0,
             0, cos, -sin, 0,
             0, sin,  cos, 0,
             0,   0,    0, 1};

    Mat4 inv{1,    0,   0, 0,
             0,  cos, sin, 0,
             0, -sin, cos, 0,
             0,    0,   0, 1};

    return {mat, inv};
}

Transform ptracer::RotateY(Float degrees) {
    const Float rads = Radians(degrees);
    const Float sin  = std::sin(rads);
    const Float cos  = std::cos(rads);

    Mat4 mat{ cos, 0, sin, 0,
                0, 1,   0, 0,
             -sin, 0, cos, 0,
                0, 0,   0, 1};

    Mat4 inv{cos, 0, -sin, 0,
               0, 1,    0, 0,
             sin, 0,  cos, 0,
               0, 0,    0, 1};

    return {mat, inv};
}

Transform ptracer::RotateZ(Float degrees) {
    const Float rads = Radians(degrees);
    const Float sin  = std::sin(rads);
    const Float cos  = std::cos(rads);

    Mat4 mat{cos, -sin, 0, 0,
             sin,  cos, 0, 0,
               0,    0, 1, 0,
               0,    0, 0, 1};

    Mat4 inv{ cos, sin, 0, 0,
             -sin, cos, 0, 0,
                0,   0, 1, 0,
                0,   0, 0, 1};

    return {mat, inv};
}

Transform ptracer::LookAt(const Point3& pos, const Point3& at, const Vec3& up) {
    const Vec3 view = pos - at;
    const Float len = view.lengthSqr();
    
    const Vec3 n = len > 0 ? -Normalize(view) : Vec3{0, 0, 1};
    const Vec3 u = Normalize(Cross(up + Vec3{0.05, 0, 0.05}, n));
    const Vec3 v = Cross(n, u);

    return {{u.x, v.x, n.x, pos.x,
             u.y, v.y, n.y, pos.y,
             u.z, v.z, n.z, pos.z,
               0,   0,   0,     1}};
}
// clang-format on

Transform ptracer::Ortho(Float zNear, Float zFar) {
    CHECK_NE(zFar, zNear);
    return Scale(1, 1, 1 / (zFar - zNear)) * Translate({0, 0, -zNear});
}
