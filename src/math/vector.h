#ifndef PT_VECTOR_H
#define PT_VECTOR_H

#include <ptracer.h>
#include <math.h>

namespace ptracer {
namespace math {

template<template<typename> typename Tup, typename T>
concept HasLength = Tup<T>::HasLength::value;

// Return types for lengths
template<typename T>
struct LenType {
    using type = Float;
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
    auto operator*(const Tup<U>& tup) const -> Tup<decltype(T{} * U{})> {
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
        DCHECK_NE(scalar, 0);
        const T recip = 1.0 / scalar;
        return {x * recip, y * recip};
    }

    template<typename U>
    Tup<T>& operator/=(U scalar) {
        DCHECK_NE(scalar, 0);
        const T recip = 1.0 / scalar;
        x *= recip;
        y *= recip;
        return static_cast<Tup<T>&>(*this);
    }

    Tup<T> operator-() const { return {-x, -y}; }

    T operator[](unsigned int idx) const { return (idx == 0) ? x : y; }
    T& operator[](unsigned int idx) { return (idx == 0) ? x : y; }

    template<typename U>
    bool operator==(Tup<U> tup) const {
        return x == tup.x && y == tup.y;
    }

    Tup<T> recip() const {
        DCHECK_NE(x, 0);
        DCHECK_NE(y, 0);
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

    void normalize() { *this /= length(); }

    T min() const {
        using std::min;
        return min(x, y);
    }

    T max() const {
        using std::max;
        return max(x, y);
    }

    unsigned int maxDim() const { return (x > y) ? 0 : 1; }
    unsigned int minDim() const { return (x < y) ? 0 : 1; }

    bool isZero() const { return x == 0 && y == 0; }
    bool isInfinity() const { return std::isinf(x) || std::isinf(y); }
};

template<template<typename> class Tup, typename T, typename U>
inline auto operator*(U scalar, const Tuple2<Tup, T>& tup) -> Tup<decltype(T{} * U{})> {
    return tup * scalar;
}

template<template<typename> class Tup, typename T>
inline auto Length(const Tuple2<Tup, T>& tup) -> typename LenType<T>::type
requires HasLength<Tup, T>
{
    return tup.length();
}

template<
    template<typename> class Tup, typename T, template<typename> class TupU, typename U>
inline auto Dot(const Tuple2<Tup, T>& tup1, const Tuple2<TupU, U>& tup2) ->
    typename LenType<std::common_type_t<T, U>>::type
requires HasLength<Tup, T> && HasLength<TupU, U>
{
    return tup1.x * tup2.x + tup1.y * tup2.y;
}

template<
    template<typename> class Tup, typename T, template<typename> class TupU, typename U>
inline auto AbsDot(const Tuple2<Tup, T>& tup1, const Tuple2<TupU, U>& tup2) ->
    typename LenType<std::common_type_t<T, U>>::type
requires HasLength<Tup, T> && HasLength<TupU, U>
{
    using std::abs;
    return abs(Dot(tup1, tup2));
}

template<template<typename> class Tup, typename T>
inline Tup<T> Abs(const Tuple2<Tup, T>& tup) {
    using std::abs;
    return {abs(tup.x), abs(tup.y)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> Normalize(const Tuple2<Tup, T>& tup)
requires HasLength<Tup, T>
{
    return tup / tup.length();
}

template<template<typename> class Tup, typename T>
inline T Min(const Tuple2<Tup, T>& tup) {
    return tup.min();
}

template<template<typename> class Tup, typename T>
inline T Max(const Tuple2<Tup, T>& tup) {
    return tup.max();
}

template<template<typename> class Tup, typename T>
inline unsigned int MaxDim(const Tuple2<Tup, T>& tup) {
    return tup.maxDim();
}

template<template<typename> class Tup, typename T>
inline unsigned int MinDim(const Tuple2<Tup, T>& tup) {
    return tup.minDim();
}

template<template<typename> class Tup, typename T>
inline Tup<T> Sign(const Tuple2<Tup, T>& tup) {
    return {Sign(tup.x), Sign(tup.y)};
}

template<
    template<typename> class Tup, typename T, template<typename> class TupU, typename U>
inline auto Min(const Tuple2<Tup, T>& tup1, const Tuple2<TupU, U>& tup2)
    -> Tup<std::common_type_t<T, U>> {
    using std::min;
    return {min(tup1.x, tup2.x), min(tup1.y, tup2.y)};
}

template<
    template<typename> class Tup, typename T, template<typename> class TupU, typename U>
inline auto Max(const Tuple2<Tup, T>& tup1, const Tuple2<TupU, U>& tup2)
    -> Tup<std::common_type_t<T, U>> {
    using std::max;
    return {max(tup1.x, tup2.x), max(tup1.y, tup2.y)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> Clamp(const Tuple2<Tup, T>& val, T low, T high) {
    return {Clamp(val.x, low, high), Clamp(val.y, low, high)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> Pow(const Tuple2<Tup, T>& tup, Float exp) {
    using std::pow;
    return {pow(tup.x, exp), pow(tup.y, exp)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> Floor(const Tuple2<Tup, T>& tup) {
    using std::floor;
    return {floor(tup.x), floor(tup.y)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> Ceil(const Tuple2<Tup, T>& tup) {
    using std::ceil;
    return {ceil(tup.x), ceil(tup.y)};
}

template<template<typename> class Tup, typename T>
inline bool IsNan(const Tuple2<Tup, T>& tup) {
    return IsNaN(tup.x) || IsNan(tup.y);
}

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
    using Tuple2<Point2T, T>::x;
    using Tuple2<Point2T, T>::y;

    using Tuple2<Point2T, T>::operator+;
    using Tuple2<Point2T, T>::operator+=;

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
inline Vector2<T> PosVec(const Point2T<T>& pt) {
    return {pt.x, pt.y};
}

template<typename T>
inline auto Dist(const Point2T<T>& pt1, const Point2T<T>& pt2) ->
    typename LenType<T>::type {
    return (pt1 - pt2).length();
}

template<typename T>
inline auto DistSqr(const Point2T<T>& pt1, const Point2T<T>& pt2) ->
    typename LenType<T>::type {
    return (pt1 - pt2).lengthSqr();
}

/* ---------------------------------------  */

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
    auto operator*(const Tup<U>& tup) const -> Tup<decltype(T{} * U{})> {
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
        DCHECK_NE(scalar, 0);
        const T recip = 1.0 / scalar;
        return {x * recip, y * recip, z * recip};
    }

    template<typename U>
    Tup<T>& operator/=(U scalar) {
        DCHECK_NE(scalar, 0);
        const T recip = 1.0 / scalar;
        x *= recip;
        y *= recip;
        z *= recip;
        return static_cast<Tup<T>&>(*this);
    }

    Tup<T> operator-() const { return {-x, -y, -z}; }

    T operator[](unsigned int idx) const { return (idx == 0) ? x : (idx == 1) ? y : z; }
    T& operator[](unsigned int idx) { return (idx == 0) ? x : (idx == 1) ? y : z; }

    template<typename U>
    bool operator==(Tup<U> tup) const {
        return x == tup.x && y == tup.y && z = tup.z;
    }

    Tup<T> recip() const {
        DCHECK_NE(x, 0);
        DCHECK_NE(y, 0);
        DCHECK_NE(z, 0);
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

    void normalize() { *this /= length(); }

    T min() const {
        using std::min;
        return min(x, min(y, z));
    }
    T max() const {
        using std::max;
        return max(x, max(y, z));
    }

    unsigned int maxDim() const {
        if (x > y)
            return (x > z) ? 0 : 2;

        return (y > z) ? 1 : 2;
    }

    unsigned int minDim() const {
        if (x < y)
            return (x < z) ? 0 : 2;

        return (y < z) ? 1 : 2;
    }

    bool isZero() const { return x == 0 && y == 0 && z == 0; }
    bool isInfinity() const { return std::isinf(x) || std::isinf(y) || std::isinf(z); }
};

template<template<typename> class Tup, typename T, typename U>
inline auto operator*(U scalar, const Tuple3<Tup, T>& tup) -> Tup<decltype(T{} * U{})> {
    return tup * scalar;
}

template<template<typename> class Tup, typename T>
inline auto Length(const Tuple3<Tup, T>& tup) -> typename LenType<T>::type
requires HasLength<Tup, T>
{
    return tup.length();
}

template<template<typename> class Tup, typename T>
inline Tup<T> Cross(const Tuple3<Tup, T>& tup1, const Tuple3<Tup, T>& tup2)
requires HasLength<Tup, T>
{
    return {
        (tup1.y * tup2.z) - (tup1.z * tup2.y), (tup1.z * tup2.x) - (tup1.x * tup2.z),
        (tup1.x * tup2.y) - (tup1.y * tup2.x)};
}

template<
    template<typename> class Tup, typename T, template<typename> class TupU, typename U>
inline auto Dot(const Tuple3<Tup, T>& tup1, const Tuple3<TupU, U>& tup2) ->
    typename LenType<std::common_type_t<T, U>>::type
requires HasLength<Tup, T> && HasLength<TupU, U>
{
    return tup1.x * tup2.x + tup1.y * tup2.y + tup1.z * tup2.z;
}

template<
    template<typename> class Tup, typename T, template<typename> class TupU, typename U>
inline auto AbsDot(const Tuple3<Tup, T>& tup1, const Tuple3<TupU, U>& tup2) ->
    typename LenType<std::common_type_t<T, U>>::type
requires HasLength<Tup, T> && HasLength<TupU, U>
{
    using std::abs;
    return abs(Dot(tup1, tup2));
}

template<template<typename> class Tup, typename T>
inline Tup<T> Abs(const Tuple3<Tup, T>& tup) {
    using std::abs;
    return {abs(tup.x), abs(tup.y), abs(tup.z)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> Normalize(const Tuple3<Tup, T>& tup)
requires HasLength<Tup, T>
{
    return tup / tup.length();
}

template<template<typename> class Tup, typename T>
inline T Min(const Tuple3<Tup, T>& tup) {
    return tup.min();
}

template<template<typename> class Tup, typename T>
inline T Max(const Tuple3<Tup, T>& tup) {
    return tup.max();
}

template<template<typename> class Tup, typename T>
inline unsigned int MaxDim(const Tuple3<Tup, T>& tup) {
    return tup.maxDim();
}

template<template<typename> class Tup, typename T>
inline unsigned int MinDim(const Tuple3<Tup, T>& tup) {
    return tup.minDim();
}

template<template<typename> class Tup, typename T>
inline Tup<T> Sign(const Tuple3<Tup, T>& tup) {
    return {Sign(tup.x), Sign(tup.y), Sign(tup.z)};
}

template<
    template<typename> class Tup, typename T, template<typename> class TupU, typename U>
inline auto Min(const Tuple3<Tup, T>& tup1, const Tuple3<TupU, U>& tup2)
    -> Tup<std::common_type_t<T, U>> {
    using std::min;
    return {min(tup1.x, tup2.x), min(tup1.y, tup2.y), min(tup1.z, tup2.z)};
}

template<
    template<typename> class Tup, typename T, template<typename> class TupU, typename U>
inline auto Max(const Tuple3<Tup, T>& tup1, const Tuple3<TupU, U>& tup2)
    -> Tup<std::common_type_t<T, U>> {
    using std::max;
    return {max(tup1.x, tup2.x), max(tup1.y, tup2.y), max(tup1.z, tup2.z)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> Clamp(const Tuple3<Tup, T>& val, T low, T high) {
    return {Clamp(val.x, low, high), Clamp(val.y, low, high), Clamp(val.z, low, high)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> Pow(const Tuple3<Tup, T>& tup, Float exp) {
    using std::pow;
    return {pow(tup.x, exp), pow(tup.y, exp), pow(tup.z, exp)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> Floor(const Tuple3<Tup, T>& tup) {
    using std::floor;
    return {floor(tup.x), floor(tup.y), floor(tup.z)};
}

template<template<typename> class Tup, typename T>
inline Tup<T> Ceil(const Tuple3<Tup, T>& tup) {
    using std::ceil;
    return {ceil(tup.x), ceil(tup.y), ceil(tup.z)};
}

template<template<typename> class Tup, typename T>
inline bool IsNan(const Tuple3<Tup, T>& tup) {
    return IsNaN(tup.x) || IsNan(tup.y) || IsNan(tup.z);
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

template<template<typename> class Tup, typename T>
inline auto BasisFromVector(const Tuple3<Tup, T>& vec)
    -> std::tuple<Vector3<T>, Vector3<T>>
requires HasLength<Tup, T>
{
    // [Duff et. al, 2017] "Building an Orthonormal Basis, Revisited"
    const Float sign = std::copysign(1.0, vec.z);
    const Float a    = -1.0 / (sign + vec.z);
    const Float b    = vec.x * vec.y * a;

    return {
        Vector3<T>(1.0 + sign * vec.x * vec.x * a, sign * b, -sign * vec.x),
        Vector3<T>(b, sign + vec.y * vec.y * a, -vec.y)};
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

    template<typename U>
    auto operator-(const Vector3<U>& v) const -> Normal3<decltype(T{} + U{})> {
        return {x - v.x, y - v.y, z - v.z};
    }
};

template<typename T>
class Point3T : public Tuple3<Point3T, T> {
public:
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
inline Vector3<T> PosVec(const Point3T<T>& pt) {
    return {pt.x, pt.y, pt.z};
}

template<typename T>
inline auto Dist(const Point3T<T>& pt1, const Point3T<T>& pt2) ->
    typename LenType<T>::type {
    return (pt1 - pt2).length();
}

template<typename T>
inline auto DistSqr(const Point3T<T>& pt1, const Point3T<T>& pt2) ->
    typename LenType<T>::type {
    return (pt1 - pt2).lengthSqr();
}

} // namespace math

using Vec3  = math::Vector3<Float>;
using Vec3i = math::Vector3<int>;
using Vec3u = math::Vector3<unsigned int>;

using Vec2  = math::Vector2<Float>;
using Vec2i = math::Vector2<int>;
using Vec2u = math::Vector2<unsigned int>;

using Point3   = math::Point3T<Float>;
using Point3i  = math::Point3T<int>;
using Point3ui = math::Point3T<unsigned int>;

using Point2   = math::Point2T<Float>;
using Point2i  = math::Point2T<int>;
using Point2ui = math::Point2T<unsigned int>;

using Normal = math::Normal3<Float>;

} // namespace ptracer

#endif // PT_VECTOR_H