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

using uint8  = std::uint8_t;
using uint16 = std::uint16_t;
using uint32 = std::uint32_t;
using uint64 = std::uint64_t;

using int8  = std::int8_t;
using int16 = std::int16_t;
using int32 = std::int32_t;
using int64 = std::int64_t;

} // namespace ptracer

#endif