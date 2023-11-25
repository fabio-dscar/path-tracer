#include <bounds.h>
#include <ray.h>
/*
#include <Sphere.h>*/
#include <tuple>

using namespace ptracer;

const Bounds3 Bounds3::Unbounded = Bounds3();

bool Bounds3::contains(const Point3& pos) const {
    return (pos.x <= bMax.x && pos.x >= bMin.x) &&
           (pos.y <= bMax.y && pos.y >= bMin.y) &&
           (pos.z <= bMax.z && pos.z >= bMin.z);
}

bool Bounds3::intersectPts(const Ray& ray, Float* t0, Float* t1) const {
    Float tmin = ray.minT;
    Float tmax = ray.maxT;

    for (int32 i = 0; i < 3; ++i) {
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
            return false;
    }

    *t0 = tmin;
    *t1 = tmax;
    return true;
}

const Point3& Bounds3::operator[](uint32 i) const {
    if (i == 0)
        return bMin;

    return bMax;
}

Point3 Bounds3::center() const {
    return (Float)0.5 * (bMax + bMin);
}

Float Bounds3::volume() const {
    return sizes().cube();
}

/*
Sphere Bounds3::sphere() const {
    Point3 pos = center();
    Float radius = dist(bMax, pos);
    return Sphere(pos, radius + F_EPSILON);
}*/

bool Bounds3::overlaps(const Bounds3& box) const {
    return (bMax.x >= box[0].x) && (bMin.x <= box[1].x) &&
           (bMax.y >= box[0].y) && (bMin.y <= box[1].y) &&
           (bMax.z >= box[0].z) && (bMin.z <= box[1].z);
}

bool Bounds3::isBounded() const {
    return !(bMin.isInfinity() || bMax.isInfinity());
}

void Bounds3::expand(Float size) {
    bMin = bMin + Point3(-size);
    bMax = bMax + Point3(size);
}

void Bounds3::expand(const Point3& pt) {
    bMin = math::min(bMin, pt);
    bMax = math::max(bMax, pt);
}

void Bounds3::expand(const Bounds3& box) {
    bMin = math::min(bMin, box[0]);
    bMax = math::max(bMax, box[1]);
}

void Bounds3::intersect(const Bounds3& box) {
    bMin = math::max(bMin, box[0]);
    bMax = math::min(bMax, box[1]);
}



const Bounds2 Bounds2::Unbounded = Bounds2();

const Point2& Bounds2::operator[](uint32 i) const {
    if (i == 0)
        return bMin;

    return bMax;
}

Vec2 Bounds2::sizes() const {
    return abs(bMax - bMin);
}

Point2 Bounds2::center() const {
    return (Float)0.5 * (bMax + bMin);
}

Float Bounds2::area() const {
    Vec2 len = sizes();
    return len.x * len.y;
}

bool Bounds2::contains(const Point2& pos) const {
    return (pos.x <= bMax.x && pos.x >= bMin.x) &&
           (pos.y <= bMax.y && pos.y >= bMin.y);
}

bool Bounds2::overlaps(const Bounds2& box) const {
    return (bMax.x >= box[0].x) && (bMin.x <= box[1].x) &&
           (bMax.y >= box[0].y) && (bMin.y <= box[1].y);
}

void Bounds2::expand(Float size) {
    bMin = bMin + Point2(-size);
    bMax = bMax + Point2(size);
}

void Bounds2::expand(const Point2& pt) {
    bMin = math::min(bMin, pt);
    bMax = math::max(bMax, pt);
}

void Bounds2::expand(const Bounds2& box) {
    bMin = math::min(bMin, box[0]);
    bMax = math::max(bMax, box[1]);
}

void Bounds2::intersect(const Bounds2& box) {
    bMin = math::max(bMin, box[0]);
    bMax = math::min(bMax, box[1]);
}

bool Bounds2::isBounded() const {
    return !(bMin.isInfinity() || bMax.isInfinity());
}