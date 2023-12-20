#ifndef __RAY_H__
#define __RAY_H__

#include <math.h>
#include <vector.h>
// #include <Spectral.h>
// #include <Frame.h>

namespace ptracer {

class Geometry;
class SurfaceEvent;
class Shape;

class Ray {
public:
    Point3 o{0};
    Vec3 dir{0};
    Float minT         = FltRayOffset;
    mutable Float maxT = FltInfinity;
    Float time         = 0;
    bool isPrimary     = false;

    Ray() = default;
    Ray(const Point3& origin, const Vec3& dir) : o(origin), dir(dir) {}
    Ray(const Point3& origin, const Vec3& dir, Float minT, Float maxT = FltInfinity)
        : o(origin), dir(dir), minT(minT), maxT(maxT) {}
    Ray(const Point3& origin, const Point3& target);

    Point3 operator()(Float t) const { return o + t * dir; }
    Point3 hitPoint() const { return o + maxT * dir; }

    Float arg(const Point3& point) const { return (point - o).length(); }
    bool inRange(Float t) const { return t > minT && t < maxT; }
};

class RayEvent {
public:
    Point3 pt{0};
    Normal n{0};
    Vec3 wo{0};

    RayEvent() = default;
    RayEvent(const Ray& ray) : pt(ray.hitPoint()), wo(-ray.dir) {}
    RayEvent(const Point3& point, const Normal& normal) : pt(point), n(normal) {}
    RayEvent(const Point3& point, const Normal& normal, const Vec3& wo)
        : pt(point), n(normal), wo(wo) {}

    Ray spawnRay(const Point3& target) const { return {pt, target}; }
    Ray spawnRay(const Vec3& dir) const { return {pt, dir}; }
    Ray spawnRay(const Vec3& dir, Float dist) const {
        return {pt, dir, FltRayOffset, dist - FltEpsilon};
    }
};

/*
    class SurfaceEvent : public RayEvent {
    public:
        const Shape* obj = nullptr;
        Point2 uv;
        Frame  gFrame;
        Frame  sFrame;
        bool   backface;

        SurfaceEvent()
            : RayEvent(), obj(nullptr), backface(false) { }

        SurfaceEvent(const Ray& ray, Shape const* obj)
            : RayEvent(ray), obj(obj), backface(false) { }

        bool hit() const;
        void setEvent(const Ray& ray, Shape const* obj, const Normal& normal);
        Color emission(const Vec3& w) const;

        Vec3 toWorld(const Vec3& w) const;
        Vec3 toLocal(const Vec3& w) const;
    };
*/
} // namespace ptracer

#endif // __RAY_H__