#ifndef PT_BBOX_H
#define PT_BBOX_H

#include <math.h>
#include <vector.h>

#include <optional>

namespace ptracer {

class Ray;
class Sphere;

class BBox3 {
public:
    using IsectResult = std::optional<std::tuple<Float, Float>>;

    static const BBox3 Unbounded;

    BBox3() = default;
    BBox3(const Point3& pt) : pMin(pt), pMax(pt) {}
    BBox3(const Point3& min, const Point3& max) : pMin(min), pMax(max) {}

    const Point3& min() const { return pMin; }
    const Point3& max() const { return pMax; }

    const Point3& operator[](unsigned int idx) const { return (idx == 0) ? pMin : pMax; }

    Vec3 sizes() const { return Abs(pMax - pMin); }

    Point3 center() const;
    Float volume() const;
    // Sphere sphere() const;

    IsectResult intersectPts(const Ray& ray, Float tMax) const;

    bool contains(const Point3& pos) const;
    bool overlaps(const BBox3& box) const;
    void expand(Float size);
    void expand(const Point3& pt);
    void expand(const BBox3& box);
    void intersect(const BBox3& box);
    bool isBounded() const;

private:
    Point3 pMin{-FloatInfinity};
    Point3 pMax{FloatInfinity};
};

inline BBox3 Expand(const BBox3& box, const Point3& pt) {
    return {math::Min(box[0], pt), math::Max(box[1], pt)};
}

inline BBox3 Expand(const BBox3& box1, const BBox3& box2) {
    return {math::Min(box1[0], box2[0]), math::Max(box1[1], box2[1])};
}

inline BBox3 Intersection(const BBox3& box1, const BBox3& box2) {
    return {math::Max(box1[0], box2[0]), math::Min(box1[1], box2[1])};
}

inline bool Overlaps(const BBox3& box1, const BBox3& box2) {
    return (box1[1].x >= box2[0].x) && (box1[0].x <= box2[1].x) &&
           (box1[1].y >= box2[0].y) && (box1[0].y <= box2[1].y) &&
           (box1[1].z >= box2[0].z) && (box1[0].z <= box2[1].z);
}

class BBox2 {
public:
    static const BBox2 Unbounded;

    BBox2() = default;
    BBox2(const Point2& pt) : pMin(pt), pMax(pt) {}
    BBox2(const Point2& min, const Point2& max) : pMin(min), pMax(max) {}
    BBox2(Float xmin, Float ymin, Float xmax, Float ymax)
        : pMin(xmin, ymin), pMax(xmax, ymax) {}

    const Point2& min() const { return pMin; }
    const Point2& max() const { return pMax; }

    const Point2& operator[](unsigned int idx) const { return (idx == 0) ? pMin : pMax; }

    Vec2 sizes() const { return Abs(pMax - pMin); }
    Float area() const;
    Point2 center() const;

    bool contains(const Point2& pos) const;
    bool overlaps(const BBox2& box) const;
    void expand(Float size);
    void expand(const Point2& pt);
    void expand(const BBox2& box);
    void intersect(const BBox2& box);
    bool isBounded() const;

private:
    Point2 pMin{-FloatInfinity};
    Point2 pMax{FloatInfinity};
};

inline BBox2 Expand(const BBox2& box, const Point2& pt) {
    return {math::Min(box[0], pt), math::Max(box[1], pt)};
}

inline BBox2 Expand(const BBox2& box1, const BBox2& box2) {
    return {math::Min(box1[0], box2[0]), math::Max(box1[1], box2[1])};
}

inline BBox2 Intersection(const BBox2& box1, const BBox2& box2) {
    return {math::Max(box1[0], box2[0]), math::Min(box1[1], box2[1])};
}

inline bool Overlaps(const BBox2& box1, const BBox2& box2) {
    return (box1[1].x >= box2[0].x) && (box1[0].x <= box2[1].x) &&
           (box1[1].y >= box2[0].y) && (box1[0].y <= box2[1].y);
}

} // namespace ptracer

#endif // PT_BBOX_H