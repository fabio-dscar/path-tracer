cmake_minimum_required(VERSION 3.20)

project(ptracer VERSION 0.1 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

##################################
## Third party libs
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ext)

##################################
## OpenGL
find_package(OpenGL REQUIRED)

set(RT_PBR_SOURCES
  src/math/transform.cpp
  src/math/bbox.cpp
  src/math/quat.cpp
  src/math/ray.cpp
  src/math/frame.cpp
  ${JSON_SOURCES}
)

add_executable(ptracer src/main.cpp ${RT_PBR_SOURCES})
target_include_directories(ptracer PUBLIC
  src
  src/math
  ext  
  ${TINYOBJ_INCLUDE}
  ${GLAD_INCLUDE}
)
target_link_libraries(ptracer PRIVATE
  OpenGL::GL
  glad
  ${GLFW_LIBRARIES}
) 