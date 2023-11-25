#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <ptracer.h>
#include <math.h>

namespace ptracer {
namespace math {

template<template<typename> typename Tup, typename T>
concept HasLength = Tup<T>::HasLength::value;

template<typename T>
struct LenType {
    using type = float;
};

template<>
struct LenType<double> {
    using type = double;
};

template<>
struct LenType<long double> {
    using type = long double;
};

template<template<typename> class Tup, typename T>
class Tuple2 {
public:
    T x{}, y{};

    Tuple2() = default;
    Tuple2(T x, T y) : x(x), y(y) {}
    explicit Tuple2(T scalar) : x(scalar), y(scalar) {}

    template<typename U>
    Tuple2(const Tup<U>& tup) : x(T(tup.x)), y(T(tup.y)) {}

    template<typename U>
    auto operator+(const Tup<U>& tup) const -> Tup<decltype(T{} + U{})> {
        return {x + tup.x, y + tup.y};
    }

    template<typename U>
    Tup<T>& operator+=(const Tup<U>& tup) {
        x += tup.x;
        y += tup.y;
        return static_cast<Tup<T>&>(*this);
    }

    template<typename U>
    auto operator-(const Tup<U>& tup) const -> Tup<decltype(T{} - U{})> {
        return {x - tup.x, y - tup.y};
    }

    template<typename U>
    Tup<T>& operator-=(const Tup<U>& tup) {
        x -= tup.x;
        y -= tup.y;
        return static_cast<Tup<T>&>(*this);
    }

    template<typename U>
    auto operator*(const Tup<U>& tup) const -> Tup<decltype(T{} - U{})> {
        return {x * tup.x, y * tup.y};
    }

    template<typename U>
    Tup<T>& operator*=(const Tup<U>& tup) {
        x *= tup.x;
        y *= tup.y;
        return static_cast<Tup<T>&>(*this);
    }

    template<typename U>
    auto operator*(U scalar) const -> Tup<decltype(T{} * U{})> {
        return {x * scalar, y * scalar};
    }

    template<typename U>
    Tup<T>& operator*=(U scalar) {
        x *= scalar;
        y *= scalar;
        return static_cast<Tup<T>&>(*this);
    }

    template<typename U>
    auto operator/(U scalar) const -> Tup<decltype(T{} / U{})> {
        assert(scalar != 0);
        T recip = 1.0 / scalar;
        return {x * recip, y * recip};
    }

    template<typename U>
    Tup<T>& operator/=(U scalar) {
        assert(scalar != 0);
        T recip = 1.0 / scalar;
        x *= recip;
        y *= recip;
        return static_cast<Tup<T>&>(*this);
    }

    Tup<T> operator-() const { return {-x, -y}; }

    T operator[](uint32 idx) const {
        if (idx == 0)
            return x;

        return y;
    }

    T& operator[](uint32 idx) {
        if (idx == 0)
            return x;

        return y;
    }

    template<typename U>
    bool operator==(Tup<U> tup) const {
        return x == tup.x && y == tup.y;
    }

    Tup<T> recip() const {
        assert(x != 0 && y != 0);
        return {1.0 / x, 1.0 / y};
    }

    auto lengthSqr() const -> typename LenType<T>::type
        requires HasLength<Tup, T>
    {
        return x * x + y * y;
    }
    auto length() const -> typename LenType<T>::type
        requires HasLength<Tup, T>
    {
        return std::sqrt(lengthSqr());
    }

    T cube() const { return x * y; }

    void normalize() { *this /= length(); }

    T min() const { return std::min(x, y); }
    T max() const { return std::max(x, y); }

    uint32 maxDim() const {
        if (x > y)
            return 0;

        return 1;
    }

    uint32 minDim() const {
        if (x < y)
            return 0;

        return 1;
    }

    bool isZero() const { return x == 0 && y == 0; }
    bool isInfinity() const { return (std::isinf(x) || std::isinf(y)); }
};

template<template<typename> class Tup, typename T, typename U>
inline auto operator*(U scalar, const Tuple2<Tup, T>& tup) -> Tup<decltype(T{} * U{})> {
    return tup * scalar;
}

template<template<typename> class Tup, typename T, template<typename> class TupU,
         typename U>
inline auto dot(const Tuple2<Tup, T>& tup1, const Tuple2<TupU, U>& tup2) ->
    typename LenType<decltype(T{} + U{})>::type
    requires HasLength<Tup, T> && HasLength<TupU, U>
{
    return tup1.x * tup2.x + tup1.y * tup2.y;
}

template<template<typename> class Tup, typename T, template<typename> class TupU,
         typename U>
inline auto absDot(const Tuple2<Tup, T>& tup1, const Tuple2<TupU, U>& tup2) ->
    typename LenType<decltype(T{} + U{})>::type
    requires HasLength<Tup, T> && HasLength<TupU, U>
{
    return std::abs(dot(tup1, tup2));
}

template<template<typename> class Tup, typename T>
inline Tup<T> abs(const Tuple2<Tup, T>& tup) {
    return {std::abs(tup.x), std::abs(tup.y)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> normalize(const Tuple2<Tup, T>& tup)
    requires HasLength<Tup, T>
{
    return tup / tup.length();
}

template<template<typename> class Tup, typename T>
inline T min(const Tuple2<Tup, T>& tup) {
    return tup.min();
}

template<template<typename> class Tup, typename T>
inline T max(const Tuple2<Tup, T>& tup) {
    return tup.max();
}

template<template<typename> class Tup, typename T>
inline uint32 maxDim(const Tuple2<Tup, T>& tup) {
    return tup.maxDim();
}

template<template<typename> class Tup, typename T>
inline uint32 minDim(const Tuple2<Tup, T>& tup) {
    return tup.minDim();
}

template<template<typename> class Tup, typename T>
inline Tup<T> sign(const Tuple2<Tup, T>& tup) {
    return {sign(tup.x), sign(tup.y)};
}

template<template<typename> class Tup, typename T, template<typename> class TupU,
         typename U>
inline Tup<decltype(T{} + U{})> min(const Tuple2<Tup, T>& tup1,
                                    const Tuple2<TupU, U>& tup2) {
    return {std::min(tup1.x, tup2.x), std::min(tup1.y, tup2.y)};
}

template<template<typename> class Tup, typename T, template<typename> class TupU,
         typename U>
inline Tup<decltype(T{} + U{})> max(const Tuple2<Tup, T>& tup1,
                                    const Tuple2<TupU, U>& tup2) {
    return {std::max(tup1.x, tup2.x), std::max(tup1.y, tup2.y)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> clamp(const Tuple2<Tup, T>& val, T low, T high) {
    return {clamp(val.x, low, high), clamp(val.y, low, high)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> pow(const Tuple2<Tup, T>& tup, Float exp) {
    return {std::pow(tup.x, exp), std::pow(tup.y, exp)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> floor(const Tuple2<Tup, T>& tup) {
    return {std::floor(tup.x), std::floor(tup.y)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> ceil(const Tuple2<Tup, T>& tup) {
    return {std::ceil(tup.x), std::ceil(tup.y)};
}

template<template<typename> class Tup, typename T>
class Tuple3 {
public:
    T x{}, y{}, z{};

    Tuple3() = default;
    Tuple3(T x, T y, T z) : x(x), y(y), z(z) {}
    explicit Tuple3(T scalar) : x(scalar), y(scalar), z(scalar) {}

    template<typename U>
    Tuple3(const Tup<U>& tup) : x(T(tup.x)), y(T(tup.y)), z(T(tup.z)) {}

    template<typename U>
    auto operator+(const Tup<U>& tup) const -> Tup<decltype(T{} + U{})> {
        return {x + tup.x, y + tup.y, z + tup.z};
    }

    template<typename U>
    Tup<T>& operator+=(const Tup<U>& tup) {
        x += tup.x;
        y += tup.y;
        z += tup.z;
        return static_cast<Tup<T>&>(*this);
    }

    template<typename U>
    auto operator-(const Tup<U>& tup) const -> Tup<decltype(T{} - U{})> {
        return {x - tup.x, y - tup.y, z - tup.z};
    }

    template<typename U>
    Tup<T>& operator-=(const Tup<U>& tup) {
        x -= tup.x;
        y -= tup.y;
        z -= tup.z;
        return static_cast<Tup<T>&>(*this);
    }

    template<typename U>
    auto operator*(const Tup<U>& tup) const -> Tup<decltype(T{} - U{})> {
        return {x * tup.x, y * tup.y, z * tup.z};
    }

    template<typename U>
    Tup<T>& operator*=(const Tup<U>& tup) {
        x *= tup.x;
        y *= tup.y;
        z *= tup.z;
        return static_cast<Tup<T>&>(*this);
    }

    template<typename U>
    auto operator*(U scalar) const -> Tup<decltype(T{} * U{})> {
        return {x * scalar, y * scalar, z * scalar};
    }

    template<typename U>
    Tup<T>& operator*=(U scalar) {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return static_cast<Tup<T>&>(*this);
    }

    template<typename U>
    auto operator/(U scalar) const -> Tup<decltype(T{} / U{})> {
        assert(scalar != 0);
        T recip = 1.0 / scalar;
        return {x * recip, y * recip, z * recip};
    }

    template<typename U>
    Tup<T>& operator/=(U scalar) {
        assert(scalar != 0);
        T recip = 1.0 / scalar;
        x *= recip;
        y *= recip;
        z *= recip;
        return static_cast<Tup<T>&>(*this);
    }

    Tup<T> operator-() const { return {-x, -y, -z}; }

    T operator[](uint32 idx) const {
        if (idx == 0)
            return x;

        if (idx == 1)
            return y;

        return z;
    }

    T& operator[](uint32 idx) {
        if (idx == 0)
            return x;

        if (idx == 1)
            return y;

        return z;
    }

    template<typename U>
    bool operator==(Tup<U> tup) const {
        return x == tup.x && y == tup.y && z = tup.z;
    }

    Tup<T> recip() const {
        assert(x != 0 && y != 0 && z != 0);
        return {1.0 / x, 1.0 / y, 1.0 / z};
    }

    auto lengthSqr() const -> typename LenType<T>::type
        requires HasLength<Tup, T>
    {
        return x * x + y * y + z * z;
    }
    auto length() const -> typename LenType<T>::type
        requires HasLength<Tup, T>
    {
        return std::sqrt(lengthSqr());
    }

    T cube() const { return x * y * z; }

    void normalize() { *this /= length(); }

    T min() const { return std::min(x, std::min(y, z)); }
    T max() const { return std::max(x, std::max(y, z)); }

    uint32 maxDim() const {
        if (x > y) {
            if (x > z)
                return 0;
            else
                return 2;
        } else {
            if (y > z)
                return 1;
            else
                return 2;
        }
    }

    uint32 minDim() const {
        if (x < y) {
            if (x < z)
                return 0;
            else
                return 2;
        } else {
            if (y < z)
                return 1;
            else
                return 2;
        }
    }

    bool isZero() const { return x == 0 && y == 0 && z == 0; }
    bool isInfinity() const { return (std::isinf(x) || std::isinf(y) || std::isinf(z)); }
};

template<template<typename> class Tup, typename T>
inline Tup<T> cross(const Tuple3<Tup, T>& tup1, const Tuple3<Tup, T>& tup2)
    requires HasLength<Tup, T>
{
    T tup1x = tup1.x, tup1y = tup1.y, tup1z = tup1.z;
    T tup2x = tup2.x, tup2y = tup2.y, tup2z = tup2.z;

    return {(tup1y * tup2z) - (tup1z * tup2y), (tup1z * tup2x) - (tup1x * tup2z),
            (tup1x * tup2y) - (tup1y * tup2x)};
}

template<template<typename> class Tup, typename T, typename U>
inline auto operator*(U scalar, const Tuple3<Tup, T>& tup) -> Tup<decltype(T{} * U{})> {
    return tup * scalar;
}

template<template<typename> class Tup, typename T, template<typename> class TupU,
         typename U>
inline auto dot(const Tuple3<Tup, T>& tup1, const Tuple3<TupU, U>& tup2) ->
    typename LenType<decltype(T{} + U{})>::type
    requires HasLength<Tup, T> && HasLength<TupU, U>
{
    return tup1.x * tup2.x + tup1.y * tup2.y + tup1.z * tup2.z;
}

template<template<typename> class Tup, typename T, template<typename> class TupU,
         typename U>
inline auto absDot(const Tuple3<Tup, T>& tup1, const Tuple3<TupU, U>& tup2) ->
    typename LenType<decltype(T{} + U{})>::type
    requires HasLength<Tup, T> && HasLength<TupU, U>
{
    return std::abs(dot(tup1, tup2));
}

template<template<typename> class Tup, typename T>
inline Tup<T> abs(const Tuple3<Tup, T>& tup) {
    return {std::abs(tup.x), std::abs(tup.y), std::abs(tup.z)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> normalize(const Tuple3<Tup, T>& tup)
    requires HasLength<Tup, T>
{
    return tup / tup.length();
}

template<template<typename> class Tup, typename T>
inline T min(const Tuple3<Tup, T>& tup) {
    return tup.min();
}

template<template<typename> class Tup, typename T>
inline T max(const Tuple3<Tup, T>& tup) {
    return tup.max();
}

template<template<typename> class Tup, typename T>
inline uint32 maxDim(const Tuple3<Tup, T>& tup) {
    return tup.maxDim();
}

template<template<typename> class Tup, typename T>
inline uint32 minDim(const Tuple3<Tup, T>& tup) {
    return tup.minDim();
}

template<template<typename> class Tup, typename T>
inline Tup<T> sign(const Tuple3<Tup, T>& tup) {
    return {sign(tup.x), sign(tup.y), sign(tup.z)};
}

template<template<typename> class Tup, typename T, template<typename> class TupU,
         typename U>
inline Tup<decltype(T{} + U{})> min(const Tuple3<Tup, T>& tup1,
                                    const Tuple3<TupU, U>& tup2) {
    return {std::min(tup1.x, tup2.x), std::min(tup1.y, tup2.y), std::min(tup1.z, tup2.z)};
}

template<template<typename> class Tup, typename T, template<typename> class TupU,
         typename U>
inline Tup<decltype(T{} + U{})> max(const Tuple3<Tup, T>& tup1,
                                    const Tuple3<TupU, U>& tup2) {
    return {std::max(tup1.x, tup2.x), std::max(tup1.y, tup2.y), std::max(tup1.z, tup2.z)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> clamp(const Tuple3<Tup, T>& val, T low, T high) {
    return {clamp(val.x, low, high), clamp(val.y, low, high), clamp(val.z, low, high)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> pow(const Tuple3<Tup, T>& tup, Float exp) {
    return {std::pow(tup.x, exp), std::pow(tup.y, exp), std::pow(tup.z, exp)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> floor(const Tuple3<Tup, T>& tup) {
    return {std::floor(tup.x), std::floor(tup.y), std::floor(tup.z)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> ceil(const Tuple3<Tup, T>& tup) {
    return {std::ceil(tup.x), std::ceil(tup.y), std::ceil(tup.z)};
}

template<typename T>
class Normal3;

template<typename T>
class Vector3 : public Tuple3<Vector3, T> {
public:
    using HasLength = std::true_type;

    using Tuple3<Vector3, T>::x;
    using Tuple3<Vector3, T>::y;
    using Tuple3<Vector3, T>::z;

    using Tuple3<Vector3, T>::operator+;
    using Tuple3<Vector3, T>::operator+=;
    using Tuple3<Vector3, T>::operator-;
    using Tuple3<Vector3, T>::operator-=;
    using Tuple3<Vector3, T>::operator*;
    using Tuple3<Vector3, T>::operator*=;

    Vector3() = default;
    Vector3(T x, T y, T z) : Tuple3<Vector3, T>(x, y, z) {}
    explicit Vector3(T s) : Tuple3<Vector3, T>(s) {}

    template<typename U>
    explicit Vector3(Vector3<U> v) : Tuple3<Vector3, T>(T(v.x), T(v.y), T(v.z)) {}

    template<typename U>
    explicit Vector3(Normal3<U> n) : Tuple3<Vector3, T>(T(n.x), T(n.y), T(n.z)) {}

    template<typename U>
    auto operator+(const Normal3<U>& v) const -> Vector3<decltype(T{} + U{})> {
        return {x + v.x, y + v.y, z + v.z};
    }
};

template<typename T>
inline void basisFromVector(const Vector3<T>& vec1, Vector3<T>& vec2, Vector3<T>& vec3) {
    // [Duff et. al, 2017] "Building an Orthonormal Basis, Revisited"
    const Float sign = std::copysign(1.0, vec1.z);
    const Float a = -1.0 / (sign + vec1.z);
    const Float b = vec1.x * vec1.y * a;

    vec2 = Vector3<T>(1.0 + sign * vec1.x * vec1.x * a, sign * b, -sign * vec1.x);
    vec3 = Vector3<T>(b, sign + vec1.y * vec1.y * a, -vec1.y);
}

template<typename T>
class Normal3 : public Tuple3<Normal3, T> {
public:
    using HasLength = std::true_type;

    using Tuple3<Normal3, T>::x;
    using Tuple3<Normal3, T>::y;
    using Tuple3<Normal3, T>::z;

    Normal3() = default;
    Normal3(T x, T y, T z) : Tuple3<Normal3, T>(x, y, z) {}
    explicit Normal3(T s) : Tuple3<Normal3, T>(s) {}

    template<typename U>
    explicit Normal3(const Vector3<U>& v) : Tuple3<Normal3, T>(v.x, v.y, v.z) {}

    template<typename U>
    auto operator+(const Vector3<U>& v) const -> Normal3<decltype(T{} + U{})> {
        return {x + v.x, y + v.y, z + v.z};
    }
};

template<typename T>
class Point3T : public Tuple3<Point3T, T> {
public:
    using HasLength = std::false_type;

    using Tuple3<Point3T, T>::x;
    using Tuple3<Point3T, T>::y;
    using Tuple3<Point3T, T>::z;

    using Tuple3<Point3T, T>::operator+;
    using Tuple3<Point3T, T>::operator+=;

    Point3T() = default;
    Point3T(T x, T y, T z) : Tuple3<Point3T, T>(x, y, z) {}
    explicit Point3T(T s) : Tuple3<Point3T, T>(s) {}

    template<typename U>
    auto operator-(const Point3T<U>& pt) const -> Vector3<decltype(T{} - U{})> {
        return {x - pt.x, y - pt.y, z - pt.z};
    }

    template<typename U>
    auto operator+(const Vector3<U>& v) const -> Point3T<decltype(T{} + U{})> {
        return {x + v.x, y + v.y, z + v.z};
    }
    template<typename U>
    Point3T<T>& operator+=(const Vector3<U>& v) {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    template<typename U>
    auto operator-(const Vector3<U>& v) const -> Point3T<decltype(T{} - U{})> {
        return {x - v.x, y - v.y, z - v.z};
    }
    template<typename U>
    Point3T<T>& operator-=(const Vector3<U>& v) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }
};

template<typename T>
inline Vector3<T> posVec(const Point3T<T>& pt) {
    return {pt.x, pt.y, pt.z};
}

template<typename T>
inline auto dist(const Point3T<T>& pt1, const Point3T<T>& pt2) ->
    typename LenType<T>::type {
    return (pt1 - pt2).length();
}

template<typename T>
inline auto distSqr(const Point3T<T>& pt1, const Point3T<T>& pt2) ->
    typename LenType<T>::type {
    return (pt1 - pt2).lengthSqr();
}

/* ------------------------- */

template<typename T>
class Vector2 : public Tuple2<Vector2, T> {
public:
    using HasLength = std::true_type;

    Vector2() = default;
    Vector2(T x, T y) : Tuple2<Vector2, T>(x, y) {}
    explicit Vector2(T s) : Tuple2<Vector2, T>(s) {}

    template<typename U>
    explicit Vector2(Vector2<U> v) : Tuple2<Vector2, T>(T(v.x), T(v.y)) {}
};

template<typename T>
class Point2T : public Tuple2<Point2T, T> {
public:
    using HasLength = std::false_type;

    using Tuple2<Point2T, T>::x;
    using Tuple2<Point2T, T>::y;

    using Tuple2<Point2T, T>::operator+;
    using Tuple2<Point2T, T>::operator+=;
    using Tuple2<Point2T, T>::operator*;
    using Tuple2<Point2T, T>::operator*=;

    Point2T() = default;
    Point2T(T x, T y) : Tuple2<Point2T, T>(x, y) {}
    explicit Point2T(T s) : Tuple2<Point2T, T>(s) {}

    template<typename U>
    auto operator-(const Point2T<U>& pt) const -> Vector2<decltype(T{} - U{})> {
        return {x - pt.x, y - pt.y};
    }

    template<typename U>
    auto operator+(const Vector2<U>& v) const -> Point2T<decltype(T{} + U{})> {
        return {x + v.x, y + v.y};
    }
    template<typename U>
    Point2T<T>& operator+=(const Vector2<U>& v) const {
        x += v.x;
        y += v.y;
        return *this;
    }

    template<typename U>
    auto operator-(const Vector2<U>& v) const -> Point2T<decltype(T{} - U{})> {
        return {x - v.x, y - v.y};
    }
    template<typename U>
    Point2T<T>& operator-=(const Vector2<U>& v) const {
        x -= v.x;
        y -= v.y;
        return *this;
    }
};

template<typename T>
inline Vector2<T> posVec(const Point2T<T>& pt) {
    return {pt.x, pt.y};
}

template<typename T>
inline auto dist(const Point2T<T>& pt1, const Point2T<T>& pt2) ->
    typename LenType<T>::type {
    return (pt1 - pt2).length();
}

template<typename T>
inline auto distSqr(const Point2T<T>& pt1, const Point2T<T>& pt2) ->
    typename LenType<T>::type {
    return (pt1 - pt2).lengthSqr();
}

} // namespace math

typedef math::Vector3<Float> Vec3;
typedef math::Vector3<int32> Vec3i;
typedef math::Vector3<uint32> Vec3ui;

typedef math::Vector2<Float> Vec2;
typedef math::Vector2<int32> Vec2i;
typedef math::Vector2<uint32> Vec2ui;

typedef math::Point2T<Float> Point2;
typedef math::Point2T<int32> Point2i;
typedef math::Point2T<uint32> Point2ui;

typedef math::Point3T<Float> Point3;
typedef math::Point3T<int32> Point3i;
typedef math::Point3T<uint32> Point3ui;

typedef math::Normal3<Float> Normal;

typedef math::Vector3<Float> Color3;

} // namespace ptracer

#endif