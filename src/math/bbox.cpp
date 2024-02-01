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
    return (pos.x <= pMax.x && pos.x >= pMin.x) && (pos.y <= pMax.y && pos.y >= pMin.y) &&
           (pos.z <= pMax.z && pos.z >= pMin.z);
}

BBox3::IsectResult BBox3::intersectPts(const Ray& ray, Float tMax) const {
    Float tMin = 0;

    for (int i = 0; i < 3; ++i) {
        const Float invDir = 1.0 / ray.dir[i];

        Float tNear = (pMin[i] - ray.o[i]) * invDir;
        Float tFar  = (pMax[i] - ray.o[i]) * invDir;
        if (tNear > tFar)
            std::swap(tNear, tFar);

        // Cut the ray's interval with this slab
        tMin = std::max(tMin, tNear);
        tMax = std::min(tMax, tFar);

        if (tMin > tMax)
            return {};
    }

    return std::make_tuple(tMin, tMax);
}

Point3 BBox3::center() const {
    return 0.5f * (pMax + pMin);
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
    return (pMax.x >= box[0].x) && (pMin.x <= box[1].x) && (pMax.y >= box[0].y) &&
           (pMin.y <= box[1].y) && (pMax.z >= box[0].z) && (pMin.z <= box[1].z);
}

bool BBox3::isBounded() const {
    return !(pMin.isInfinity() || pMax.isInfinity());
}

void BBox3::expand(Float size) {
    pMin = pMin + Point3(-size);
    pMax = pMax + Point3(size);
}

void BBox3::expand(const Point3& pt) {
    pMin = Min(pMin, pt);
    pMax = Max(pMax, pt);
}

void BBox3::expand(const BBox3& box) {
    pMin = Min(pMin, box[0]);
    pMax = Max(pMax, box[1]);
}

void BBox3::intersect(const BBox3& box) {
    pMin = Max(pMin, box[0]);
    pMax = Min(pMax, box[1]);
}

const BBox2 BBox2::Unbounded = BBox2();

Point2 BBox2::center() const {
    return static_cast<Float>(0.5) * (pMax + pMin);
}

Float BBox2::area() const {
    const Vec2 sizesVec = sizes();
    return sizesVec.x * sizesVec.y;
}

bool BBox2::contains(const Point2& pos) const {
    return (pos.x <= pMax.x && pos.x >= pMin.x) && (pos.y <= pMax.y && pos.y >= pMin.y);
}

bool BBox2::overlaps(const BBox2& box) const {
    return (pMax.x >= box[0].x) && (pMin.x <= box[1].x) && (pMax.y >= box[0].y) &&
           (pMin.y <= box[1].y);
}

void BBox2::expand(Float size) {
    pMin = pMin + Point2(-size);
    pMax = pMax + Point2(size);
}

void BBox2::expand(const Point2& pt) {
    pMin = Min(pMin, pt);
    pMax = Max(pMax, pt);
}

void BBox2::expand(const BBox2& box) {
    pMin = Min(pMin, box[0]);
    pMax = Max(pMax, box[1]);
}

void BBox2::intersect(const BBox2& box) {
    pMin = Max(pMin, box[0]);
    pMax = Min(pMax, box[1]);
}

bool BBox2::isBounded() const {
    return !(pMin.isInfinity() || pMax.isInfinity());
}