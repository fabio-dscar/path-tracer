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

static constexpr Float Epsilon   = 1e-10;
static constexpr Float RayOffset = 1e-6;
#else
using Float = float;

static constexpr Float Epsilon   = 1e-6f;
static constexpr Float RayOffset = 1e-3f;
#endif

static constexpr Float Pi            = 3.141592653589793238462643383279502884;
static constexpr Float InvPi         = std::numbers::inv_pi;
static constexpr Float Inv2Pi        = 1.0 / (2.0 * Pi);
static constexpr Float Inv4Pi        = 1.0 / (4.0 * Pi);
static constexpr Float PiOver2       = Pi / 2.0;
static constexpr Float PiOver4       = Pi / 4.0;
static constexpr Float PiOver180     = Pi / 180.0;
static constexpr Float InvPiOver180  = 1.0 / PiOver180;
static constexpr Float SqrtInvPi     = std::numbers::inv_sqrtpi;
static constexpr Float Sqrt2         = std::numbers::sqrt2;
static constexpr Float InvLog2       = 1.442695040888963387004650940071;
static constexpr Float FloatLowest   = std::numeric_limits<Float>::lowest();
static constexpr Float FloatMaximum  = std::numeric_limits<Float>::max();
static constexpr Float FloatInfinity = std::numeric_limits<Float>::infinity();

namespace math {

inline constexpr bool IsNaN(Float f) {
    return std::isnan(f);
}

template<typename T, typename U, typename R = std::common_type_t<T, U>>
inline constexpr R Max(T v1, U v2)
requires(std::is_arithmetic_v<T> && std::is_arithmetic_v<U>)
{
    return std::max(static_cast<R>(v1), static_cast<R>(v2));
}

template<typename T, typename U, typename R = std::common_type_t<T, U>>
inline constexpr R Min(T v1, U v2)
requires(std::is_arithmetic_v<T> && std::is_arithmetic_v<U>)
{
    return std::min(static_cast<R>(v1), static_cast<R>(v2));
}

template<typename T, typename U, typename V>
inline constexpr T Clamp(T val, U low, V high) {
    return Max(low, Min(val, high));
}

inline constexpr int32 Sign(Float scalar) {
    return std::signbit(scalar) ? -1 : 1;
}

inline constexpr Float SafeAcos(Float x) {
    return std::acos(Clamp(x, -1.0, 1.0));
}

inline constexpr Float SafeSqrt(Float x) {
    return std::sqrt(Max(x, 0));
}

inline constexpr Float Radians(Float degs) {
    return PiOver180 * degs;
}

inline constexpr Float Degrees(Float rads) {
    return InvPiOver180 * rads;
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