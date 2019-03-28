cmake_minimum_required(VERSION 3.1)

project(pbrt)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

###########################
# Sources
file(GLOB_RECURSE pbrt_SOURCES "src/*.cpp")
file(GLOB_RECURSE pbrt_HEADERS "src/*.h")


SET(pbrt_INCLUDE "src" "src/core" "src/shapes" "src/accelerators" "src/materials")

find_package(OpenEXR)
find_package(TIFF)

if(NOT ${OPENEXR_FOUND})
	message(FATAL_ERROR "OPENEXR_FOUND")
endif()

if(NOT ${TIFF_FOUND})
	message(FATAL_ERROR "TIFF_FOUND")
endif()

add_library(pbrt ${pbrt_SOURCES} ${pbrt_SOURCES})
target_include_directories(pbrt PUBLIC ${pbrt_INCLUDE} ${OPENEXR_INCLUDE_DIR} ${TIFF_INCLUDE_DIR})
target_link_libraries(pbrt PUBLIC ${OPENEXR_LIBRARIES} ${ILMBASE_LIBRARIES} ${TIFF_LIBRARIES})
target_compile_definitions(pbrt PUBLIC "-DPBRT_PROBES_NONE")
