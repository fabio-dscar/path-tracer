#ifndef __BOUNDS_H__
#define __BOUNDS_H__

#include <math.h>
#include <vector.h>

namespace ptracer {

class Ray;
class Sphere;

class Bounds3 {
public:
    static const Bounds3 Unbounded;

    Bounds3() = default;
    Bounds3(const Point3& pt) : bMin(pt), bMax(pt) {}
    Bounds3(const Point3& min, const Point3& max) : bMin(min), bMax(max) {}

    const Point3& min() const { return bMin; }
    const Point3& max() const { return bMax; }
    const Point3& operator[](uint32 i) const;

    Vec3 sizes() const { return abs(bMax - bMin); }

    Point3 center() const;
    Float volume() const;
    //Sphere sphere() const;

    bool intersectPts(const Ray& ray, Float* t0, Float* t1) const;
    bool contains(const Point3& pos) const;
    bool overlaps(const Bounds3& box) const;
    void expand(Float size);
    void expand(const Point3& pt);
    void expand(const Bounds3& box);
    void intersect(const Bounds3& box);
    bool isBounded() const;

private:
    Point3 bMin{-F_INFINITY};
    Point3 bMax{F_INFINITY};
};

// clang-format off
inline Bounds3 expand(const Bounds3& box, const Point3& pt) {
    return Bounds3(Point3(std::min(box[0].x, pt.x), 
                          std::min(box[0].y, pt.y),
                          std::min(box[0].z, pt.z)),
                   Point3(std::max(box[1].x, pt.x), 
                          std::max(box[1].y, pt.y),
                          std::max(box[1].z, pt.z)));
}

inline Bounds3 expand(const Bounds3& box1, const Bounds3& box2) {
    return Bounds3(Point3(std::min(box1[0].x, box2[0].x), 
                          std::min(box1[0].y, box2[0].y),
                          std::min(box1[0].z, box2[0].z)),
                   Point3(std::max(box1[1].x, box2[1].x), 
                          std::max(box1[1].y, box2[1].y),
                          std::max(box1[1].z, box2[1].z)));
}

inline Bounds3 intersection(const Bounds3& box1, const Bounds3& box2) {
    return Bounds3(Point3(std::max(box1[0].x, box2[0].x), 
                          std::max(box1[0].y, box2[0].y),
                          std::max(box1[0].z, box2[0].z)),
                   Point3(std::min(box1[1].x, box2[1].x), 
                          std::min(box1[1].y, box2[1].y),
                          std::min(box1[1].z, box2[1].z)));
}

inline bool overlaps(const Bounds3& box1, const Bounds3& box2) {
    return (box1[1].x >= box2[0].x) && (box1[0].x <= box2[1].x) &&
           (box1[1].y >= box2[0].y) && (box1[0].y <= box2[1].y) &&
           (box1[1].z >= box2[0].z) && (box1[0].z <= box2[1].z);
}
// clang-format on

class Bounds2 {
public:
    static const Bounds2 Unbounded;

    Bounds2() = default;
    Bounds2(const Point2& pt) : bMin(pt), bMax(pt) {}
    Bounds2(const Point2& min, const Point2& max) : bMin(min), bMax(max) {}
    Bounds2(Float xmin, Float ymin, Float xmax, Float ymax)
        : bMin(xmin, ymin), bMax(xmax, ymax) {}

    const Point2& min() const { return bMin; }
    const Point2& max() const { return bMax; }
    const Point2& operator[](uint32 i) const;

    Vec2 sizes() const;
    Point2 center() const;
    Float area() const;

    bool contains(const Point2& pos) const;
    bool overlaps(const Bounds2& box) const;
    void expand(Float size);
    void expand(const Point2& pt);
    void expand(const Bounds2& box);
    void intersect(const Bounds2& box);
    bool isBounded() const;

private:
    Point2 bMin{-F_INFINITY};
    Point2 bMax{F_INFINITY};
};

inline Bounds2 expand(const Bounds2& box, const Point2& pt) {
    return {std::min(box[0].x, pt.x), std::min(box[0].y, pt.y), 
            std::max(box[1].x, pt.x), std::max(box[1].y, pt.y)};
}

inline Bounds2 expand(const Bounds2& box1, const Bounds2& box2) {
    return {std::min(box1[0].x, box2[0].x), std::min(box1[0].y, box2[0].y),
            std::max(box1[1].x, box2[1].x), std::max(box1[1].y, box2[1].y)};
}

inline Bounds2 intersection(const Bounds2& box1, const Bounds2& box2) {
    return {std::max(box1[0].x, box2[0].x), std::max(box1[0].y, box2[0].y),
            std::min(box1[1].x, box2[1].x), std::min(box1[1].y, box2[1].y)};
}

inline bool overlaps(const Bounds2& box1, const Bounds2& box2) {
    return (box1[1].x >= box2[0].x) && (box1[0].x <= box2[1].x) &&
           (box1[1].y >= box2[0].y) && (box1[0].y <= box2[1].y);
}

} // namespace ptracer

#endif // __BOUNDS_H__