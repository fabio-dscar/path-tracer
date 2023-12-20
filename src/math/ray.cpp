
#include <ray.h>
// #include <Shape.h>
// #include <AreaLight.h>

using namespace ptracer;

Ray::Ray(const Point3& origin, const Point3& target) : o(origin) {
    dir  = Normalize(target - origin);
    maxT = arg(target);

    // Apply offset - avoid self-intersection
    maxT -= maxT * FltRayOffset;
}

/*
void SurfaceEvent::setEvent(const Ray& ray, Shape const* shape, const Normal& n) {
    obj = shape;
    point = ray.hitPoint();
    normal = normalize(n);
    gFrame = Frame(normal);
    wo = gFrame.toLocal(-ray.dir());
    sFrame = gFrame;
    backface = false;
    //    sFrame = gFrame;   // Use shading normal instead
        wo     = sFrame.toLocal(-ray.dir());
}

bool SurfaceEvent::hit() const {
    return obj != nullptr;
}

Color SurfaceEvent::emission(const Vec3& w) const {
    // If surface is emissive, return emission
    if (obj->isLight())
        return obj->areaLight()->evalL(*this, w);

    return Color::BLACK;
}

Vec3 SurfaceEvent::toWorld(const Vec3& w) const {
    return sFrame.toWorld(w);
}

Vec3 SurfaceEvent::toLocal(const Vec3& w) const {
    return sFrame.toLocal(w);
}*/