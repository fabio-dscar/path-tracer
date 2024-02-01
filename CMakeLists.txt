cmake_minimum_required(VERSION 3.20)

project(
  ptracer
  VERSION 0.1
  LANGUAGES CXX C)

option(PTR_BUILD_TESTS "Build tests" ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

# ---------------------------------------------------------------------------------------
#     Third party libs
# ---------------------------------------------------------------------------------------
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ext)

# ---------------------------------------------------------------------------------------
#     OpenGL
# ---------------------------------------------------------------------------------------
find_package(OpenGL REQUIRED)

set(PTR_SOURCES
    src/math/transform.cpp
    src/math/bbox.cpp
    src/math/quat.cpp
    src/math/ray.cpp
    src/math/frame.cpp
    src/util/log.cpp
    ${JSON_SOURCES})

set(PTR_INCLUDE_DIRS
    src
    src/math
    src/util
    ext
    ${SPDLOG_INCLUDE_DIR}
    ${BACKWARD_INCLUDE_DIRS}
    ${TINYOBJ_INCLUDE}
    ${GLAD_INCLUDE}
)

set(DEBUG_FLAGS -Wall -Wextra -Wpedantic)
set(RELEASE_FLAGS -O3 -march=native)

add_library(ptracer_lib STATIC ${PTR_SOURCES})
target_compile_features(ptracer_lib PUBLIC cxx_std_20)
target_compile_definitions(ptracer_lib PUBLIC "$<$<CONFIG:Debug>:DEBUG>")
target_include_directories(ptracer_lib PUBLIC ${PTR_INCLUDE_DIRS})
target_link_libraries(ptracer_lib PRIVATE OpenGL::GL glad ${GLFW_LIBRARIES}
                                      Backward::Interface spdlog)
target_compile_options(
  ptracer_lib
  PRIVATE "$<$<OR:$<CONFIG:Release>,$<CONFIG:MinSizeRel>>:${RELEASE_FLAGS}>"
          "$<$<CONFIG:Debug>:${DEBUG_FLAGS}>")

add_executable(ptracer src/main.cpp)
target_compile_features(ptracer PUBLIC cxx_std_20)
target_compile_definitions(ptracer PUBLIC "$<$<CONFIG:Debug>:DEBUG>")
target_include_directories(ptracer PUBLIC ${PTR_INCLUDE_DIRS})
target_link_libraries(ptracer PRIVATE ptracer_lib)

target_compile_options(
  ptracer
  PRIVATE "$<$<OR:$<CONFIG:Release>,$<CONFIG:MinSizeRel>>:${RELEASE_FLAGS}>"
          "$<$<CONFIG:Debug>:${DEBUG_FLAGS}>")

if(PTR_BUILD_TESTS)
  set(PTR_TEST_SOURCES
    src/math/tests/test_vector.cpp
  )

  add_executable(ptr_tests ${PTR_TEST_SOURCES})
  target_compile_features(ptr_tests PUBLIC cxx_std_20)
  target_link_libraries(ptr_tests PRIVATE ptracer_lib Catch2::Catch2WithMain)
  target_include_directories(ptr_tests PRIVATE src src/math src/util)

  enable_testing()
  list(APPEND CMAKE_MODULE_PATH ${CATCH_EXTRAS_DIR})
  include(CTest)
  include(Catch)
  catch_discover_tests(ptr_tests)
endif()