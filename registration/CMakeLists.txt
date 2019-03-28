cmake_minimum_required(VERSION 3.1)

project(registration)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

find_package(Cholmod)


###########################
# Sources
file(GLOB_RECURSE trimesh2_SOURCES "trimesh2/libsrc/*.c" "trimesh2/libsrc/*.cc" "trimesh2/optsrc/*.c" "trimesh2/optsrc/*.cc" "trimesh2/utilsrc/*.c" "trimesh2/utilsrc/*.cc" )
file(GLOB_RECURSE trimesh2_HEADERS "trimesh2/*.h")

SET(trimesh2_INCLUDE "trimesh2/include")


add_library(trimesh2 ${trimesh2_SOURCES} ${trimesh2_HEADERS})
target_include_directories(trimesh2 PUBLIC ${trimesh2_INCLUDE} ${CHOLMOD_INCLUDES})
target_link_libraries(trimesh2 PUBLIC ${CHOLMOD_LIBRARIES})


# Sources
file(GLOB_RECURSE tps_SOURCES "tps_alignment/*.c" "tps_alignment/*.cc")
file(GLOB_RECURSE tps_HEADERS "tps_alignment/*.h")

SET(tps_INCLUDE "tps_alignment/include" "tnt")


add_library(tps ${tps_SOURCES} ${tps_HEADERS})
target_include_directories(tps PUBLIC ${tps_INCLUDE})
target_link_libraries(tps PUBLIC trimesh2)
