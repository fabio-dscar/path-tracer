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
typedef double Float;

static constexpr Float F_EPSILON = 1e-10;
static constexpr Float F_RAY_OFFSET = 1e-6;
#else
typedef float Float;

static const Float F_EPSILON = 1e-6f;
static const Float F_RAY_OFFSET = 1e-3f;
#endif

static constexpr Float PI = std::numbers::pi;
static constexpr Float INVPI = std::numbers::inv_pi;
static constexpr Float INV2PI = 1.0 / (2.0 * PI);
static constexpr Float INV4PI = 1.0 / (4.0 * PI);
static constexpr Float PIOVER2 = PI / 2.0;
static constexpr Float PIOVER4 = PI / 4.0;
static constexpr Float PIOVER180 = PI / 180.0;
static constexpr Float INVPIOVER180 = 1.0 / PIOVER180;
static constexpr Float SQRTINVPI = std::numbers::inv_sqrtpi;
static constexpr Float SQRT2 = std::numbers::sqrt2;
static constexpr Float INVLOG2 = 1.442695040888963387004650940071;
static constexpr Float F_LOWEST = std::numeric_limits<Float>::lowest();
static constexpr Float F_MAXIMUM = std::numeric_limits<Float>::max();
static constexpr Float F_INFINITY = std::numeric_limits<Float>::infinity();

namespace math {

inline bool isNaN(Float f) {
    return std::isnan(f);
}

template<typename T, typename U, typename K = decltype(T{} + U{})>
inline constexpr auto max(T v1, U v2) -> K
    requires(std::is_arithmetic_v<T> && std::is_arithmetic_v<U>)
{
    return std::max(static_cast<K>(v1), static_cast<K>(v2));
}

template<typename T, typename U, typename V>
inline T clamp(T val, U low, V high) {
    if (val < low)
        return low;
    else if (val > high)
        return high;
    else
        return val;
}

template<typename T>
inline T clamp(T val, T low, T high) {
    return clamp(val, low, high);
}

inline constexpr int32 sign(Float scalar) {
    if (std::signbit(scalar))
        return -1;
    return 1;
}

inline constexpr Float acosSafe(Float x) {
    return std::acos(clamp(x, -1.0, 1.0));
}

inline constexpr Float sqrtSafe(Float x) {
    return std::sqrt(std::max(static_cast<Float>(0.0), x));
}
inline constexpr Float radians(Float degrees) {
    return PIOVER180 * degrees;
}

inline constexpr Float degrees(Float radians) {
    return INVPIOVER180 * radians;
}

inline constexpr Float log2(Float x) {
    return std::log(x) * INVLOG2;
}

inline constexpr Float lerp(Float t, Float v1, Float v2) {
    return (1 - t) * v1 + t * v2;
}

} // namespace math
} // namespace ptracer

#endif // __PTR_MATH_H__