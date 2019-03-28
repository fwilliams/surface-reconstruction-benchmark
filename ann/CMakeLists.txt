cmake_minimum_required(VERSION 3.1)

project(ann)



###########################
# Sources
file(GLOB_RECURSE ANN_SOURCES "src/*.cpp")
file(GLOB_RECURSE ANN_HEADERS "src/*.h")

SET(ANN_MAIN "sample/ann_sample.cpp")
SET(ANN_INCLUDE "include")

add_library(ann ${ANN_SOURCES} ${ANN_HEADERS})
target_include_directories(ann PUBLIC ${ANN_INCLUDE})


add_executable(ann_exec ${ANN_MAIN})
target_link_libraries(ann_exec ann)