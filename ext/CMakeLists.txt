cmake_minimum_required (VERSION 3.20)

#########################################
## freeglut
set(FREEGLUT_BUILD_DEMOS OFF CACHE BOOL " " FORCE)
add_subdirectory(freeglut)

#########################################
## glad
set(GLAD_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/glad/include PARENT_SCOPE)
add_library(glad STATIC ${CMAKE_CURRENT_SOURCE_DIR}/glad/src/glad.c)

target_include_directories(glad PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/glad/include)

set_property(TARGET glad PROPERTY FOLDER "ext")

#########################################
## tinyobj
set(TINYOBJ_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/tinyobj PARENT_SCOPE)

#########################################
## json
set(JSON_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/json/json.hpp PARENT_SCOPE)

