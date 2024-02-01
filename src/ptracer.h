#ifndef PTRACER_H
#define PTRACER_H

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

#include <check.h>

namespace ptracer {



} // namespace ptracer

#endif