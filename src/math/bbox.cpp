#include <bbox.h>
#include <ray.h>
/*
#include <Sphere.h>*/

#include <optional>
#include <tuple>

using namespace ptracer;
using namespace ptracer::math;

const BBox3 BBox3::Unbounded = BBox3();

bool BBox3::contains(const Point3& pos) const {
    return (pos.x <= bMax.x && pos.x >= bMin.x) && (pos.y <= bMax.y && pos.y >= bMin.y) &&
           (pos.z <= bMax.z && pos.z >= bMin.z);
}

BBox3::IsectResult BBox3::intersectPts(const Ray& ray) const {
    Float tmin = ray.minT;
    Float tmax = ray.maxT;

    for (int i = 0; i < 3; ++i) {
        const Float invDir = 1.0 / ray.dir[i];

        // Compute intersections and sort
        Float tNear = (bMin[i] - ray.o[i]) * invDir;
        Float tFar  = (bMax[i] - ray.o[i]) * invDir;
        if (tNear > tFar)
            std::swap(tNear, tFar);

        // Cut the ray's interval with this slab
        if (tNear > tmin)
            tmin = tNear;

        if (tFar < tmax)
            tmax = tFar;

        // If the interval is not valid, it does not intersect
        if (tmin > tmax)
            return {};
    }

    return std::make_tuple(tmin, tmax);
}

Point3 BBox3::center() const {
    return static_cast<Float>(0.5) * (bMax + bMin);
}

Float BBox3::volume() const {
    const Vec3 sizesVec = sizes();
    return sizesVec.x * sizesVec.y * sizesVec.z;
}

/*
Sphere Bounds3::sphere() const {
    Point3 pos = center();
    Float radius = dist(bMax, pos);
    return Sphere(pos, radius + F_EPSILON);
}*/

bool BBox3::overlaps(const BBox3& box) const {
    return (bMax.x >= box[0].x) && (bMin.x <= box[1].x) && (bMax.y >= box[0].y) &&
           (bMin.y <= box[1].y) && (bMax.z >= box[0].z) && (bMin.z <= box[1].z);
}

bool BBox3::isBounded() const {
    return !(bMin.isInfinity() || bMax.isInfinity());
}

void BBox3::expand(Float size) {
    bMin = bMin + Point3(-size);
    bMax = bMax + Point3(size);
}

void BBox3::expand(const Point3& pt) {
    bMin = Min(bMin, pt);
    bMax = Max(bMax, pt);
}

void BBox3::expand(const BBox3& box) {
    bMin = Min(bMin, box[0]);
    bMax = Max(bMax, box[1]);
}

void BBox3::intersect(const BBox3& box) {
    bMin = Max(bMin, box[0]);
    bMax = Min(bMax, box[1]);
}

const BBox2 BBox2::Unbounded = BBox2();

Point2 BBox2::center() const {
    return static_cast<Float>(0.5) * (bMax + bMin);
}

Float BBox2::area() const {
    const Vec2 sizesVec = sizes();
    return sizesVec.x * sizesVec.y;
}

bool BBox2::contains(const Point2& pos) const {
    return (pos.x <= bMax.x && pos.x >= bMin.x) && (pos.y <= bMax.y && pos.y >= bMin.y);
}

bool BBox2::overlaps(const BBox2& box) const {
    return (bMax.x >= box[0].x) && (bMin.x <= box[1].x) && (bMax.y >= box[0].y) &&
           (bMin.y <= box[1].y);
}

void BBox2::expand(Float size) {
    bMin = bMin + Point2(-size);
    bMax = bMax + Point2(size);
}

void BBox2::expand(const Point2& pt) {
    bMin = Min(bMin, pt);
    bMax = Max(bMax, pt);
}

void BBox2::expand(const BBox2& box) {
    bMin = Min(bMin, box[0]);
    bMax = Max(bMax, box[1]);
}

void BBox2::intersect(const BBox2& box) {
    bMin = Max(bMin, box[0]);
    bMax = Min(bMax, box[1]);
}

bool BBox2::isBounded() const {
    return !(bMin.isInfinity() || bMax.isInfinity());
}