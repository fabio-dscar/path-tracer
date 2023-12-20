#ifndef __PTR_MATH_H__
#define __PTR_MATH_H__

#include <cmath>
#include <type_traits>

#define PTR_DOUBLE 0

#include <ptracer.h>
#include <functional>
#include <numbers>

namespace ptracer {

#if PTR_DOUBLE
using Float = double;

static constexpr Float FltEpsilon   = 1e-10;
static constexpr Float FltRayOffset = 1e-6;
#else
using Float = float;

static constexpr Float FltEpsilon   = 1e-6f;
static constexpr Float FltRayOffset = 1e-3f;
#endif

static constexpr Float Pi           = 3.141592653589793238462643383279502884;
static constexpr Float InvPi        = std::numbers::inv_pi;
static constexpr Float Inv2Pi       = 1.0 / (2.0 * Pi);
static constexpr Float Inv4Pi       = 1.0 / (4.0 * Pi);
static constexpr Float PiOver2      = Pi / 2.0;
static constexpr Float PiOver4      = Pi / 4.0;
static constexpr Float PiOver180    = Pi / 180.0;
static constexpr Float InvPiOver180 = 1.0 / PiOver180;
static constexpr Float SqrtInvPi    = std::numbers::inv_sqrtpi;
static constexpr Float Sqrt2        = std::numbers::sqrt2;
static constexpr Float InvLog2      = 1.442695040888963387004650940071;
static constexpr Float FltLowest    = std::numeric_limits<Float>::lowest();
static constexpr Float FltMaximum   = std::numeric_limits<Float>::max();
static constexpr Float FltInfinity  = std::numeric_limits<Float>::infinity();

namespace math {

inline bool IsNaN(Float f) {
    return std::isnan(f);
}

template<typename T, typename U, typename K = decltype(T{} + U{})>
inline constexpr K Max(T v1, U v2)
requires(std::is_arithmetic_v<T> && std::is_arithmetic_v<U>)
{
    return std::max(static_cast<K>(v1), static_cast<K>(v2));
}

template<typename T, typename U, typename V>
inline constexpr T Clamp(T val, U low, V high) {
    if (val < low)
        return low;
    else if (val > high)
        return high;
    else
        return val;
}

template<typename T>
inline constexpr T Clamp(T val, T low, T high) {
    return Clamp(val, low, high);
}

inline constexpr int32 Sign(Float scalar) {
    return std::signbit(scalar) ? -1 : 1;
}

inline constexpr Float SafeAcos(Float x) {
    return std::acos(Clamp(x, -1.0, 1.0));
}

inline constexpr Float SafeSqrt(Float x) {
    return std::sqrt(std::max(static_cast<Float>(0.0), x));
}
inline constexpr Float Radians(Float degrees) {
    return PiOver180 * degrees;
}

inline constexpr Float Degrees(Float radians) {
    return InvPiOver180 * radians;
}

inline constexpr Float Log2(Float x) {
    return std::log(x) * InvLog2;
}

inline constexpr Float Lerp(Float t, Float v1, Float v2) {
    return (1 - t) * v1 + t * v2;
}

} // namespace math
} // namespace ptracer

#endif // __PTR_MATH_H__