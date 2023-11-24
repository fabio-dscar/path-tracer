#ifndef __PTRACER_H__
#define __PTRACER_H__

// #define PTR_DEBUG 0

#if defined(_MSC_VER)
#define NOMINMAX
#endif

#if defined(_WIN32) || defined(_WIN64)
#define PTR_WINDOWS
#elif defined(__linux__)
#define PTR_LINUX
#endif

#include <type_traits>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <cstring>
#include <vector>
#include <cassert>
#include <cstdint>

namespace ptracer {

typedef std::uint8_t uint8;
typedef std::uint16_t uint16;
typedef std::uint32_t uint32;
typedef std::uint64_t uint64;

typedef std::int8_t int8;
typedef std::int16_t int16;
typedef std::int32_t int32;
typedef std::int64_t int64;

} // namespace ptr

#endif