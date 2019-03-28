cmake_minimum_required(VERSION 3.1)

project(OpenMesh)



###########################
# Sources
file(GLOB_RECURSE OpenMesh_SOURCES "OpenMesh/*.cc")
file(GLOB_RECURSE OpenMesh_HEADERS "OpenMesh/*.hh")

SET(OpenMesh_INCLUDE "./")


add_library(OpenMesh ${OpenMesh_SOURCES} ${OpenMesh_SOURCES})
target_include_directories(OpenMesh PUBLIC ${OpenMesh_INCLUDE})
