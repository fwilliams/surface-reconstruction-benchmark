cmake_minimum_required(VERSION 3.1)

project(ReconBench)

add_subdirectory(ann)
add_subdirectory(OpenMesh)
add_subdirectory(pbrt)
add_subdirectory(registration)

find_package(LAPACK)
find_package(PNG)

if(NOT ${LAPACK_FOUND})
	message(FATAL_ERROR "LAPACK_FOUND")
endif()

if(NOT ${PNG_FOUND})
	message(FATAL_ERROR "PNG_FOUND")
endif()


###########################
# Sources
file(GLOB_RECURSE Modeling_SOURCES "modeling/*.cpp")
file(GLOB_RECURSE Modeling_HEADERS "modeling/*.h")
SET(Modeling_INCLUDE "modeling")

add_library(Modeling ${Modeling_SOURCES} ${Modeling_HEADERS})
target_include_directories(Modeling PUBLIC ${Modeling_INCLUDE})
target_link_libraries(Modeling PUBLIC ann OpenMesh pbrt tps ${LAPACK_LIBRARIES})
target_compile_definitions(Modeling PUBLIC ${LAPACK_LINKER_FLAGS})




# Sources
file(GLOB_RECURSE sampler_SOURCES "sampler/*.cpp")
file(GLOB_RECURSE sampler_HEADERS "sampler/*.h")

SET(sampler_INCLUDE "sampler")


add_library(Sampler ${sampler_SOURCES} ${sampler_HEADERS})
target_include_directories(Sampler PUBLIC ${sampler_INCLUDE} ${PNG_INCLUDE_DIRS})
target_link_libraries(Sampler PUBLIC ann OpenMesh pbrt tps ${PNG_LIBRARIES} ${PNG_LIBRARY} Modeling)
target_compile_definitions(Sampler PUBLIC ${PNG_DEFINITIONS})






# Sources
file(GLOB_RECURSE evaluator_HEADERS "evaluator/*.h")
set(evaluator_SOURCES
	evaluator/GlobalStats.cpp
	evaluator/ShortestDistanceMap.cpp
	evaluator/implicit_sampler.cpp
	evaluator/UniformGrid.cpp
	evaluator/mesh_sampler.cpp
	evaluator/recon_evaluation.cpp
	evaluator/KdOpenMesh.cpp
	evaluator/kd_tree.cpp
	evaluator/shape_distribution.cpp
)
SET(evaluator_INCLUDE "evaluator")



add_executable(implicit_uniform ${evaluator_SOURCES} ${evaluator_HEADERS} evaluator/implicit_uniform.cpp)
target_include_directories(implicit_uniform PUBLIC ${evaluator_INCLUDE})
target_link_libraries(implicit_uniform PUBLIC Modeling Sampler)

add_executable(mesh_fixer ${evaluator_SOURCES} ${evaluator_HEADERS} evaluator/mesh_fixer.cpp)
target_include_directories(mesh_fixer PUBLIC ${evaluator_INCLUDE})
target_link_libraries(mesh_fixer PUBLIC Modeling Sampler)

add_executable(run_evaluation ${evaluator_SOURCES} ${evaluator_HEADERS} evaluator/run_evaluation.cpp)
target_include_directories(run_evaluation PUBLIC ${evaluator_INCLUDE})
target_link_libraries(run_evaluation PUBLIC Modeling Sampler)

add_executable(sample_mesh ${evaluator_SOURCES} ${evaluator_HEADERS} evaluator/sample_mesh.cpp)
target_include_directories(sample_mesh PUBLIC ${evaluator_INCLUDE})
target_link_libraries(sample_mesh PUBLIC Modeling Sampler)

add_executable(aggregate_distribution ${evaluator_SOURCES} ${evaluator_HEADERS} evaluator/plot/aggregate_distribution.cpp)
target_include_directories(aggregate_distribution PUBLIC ${evaluator_INCLUDE})
target_link_libraries(aggregate_distribution PUBLIC Modeling Sampler)

add_executable(aggregation ${evaluator_SOURCES} ${evaluator_HEADERS} evaluator/plot/aggregation.cpp)
target_include_directories(aggregation PUBLIC ${evaluator_INCLUDE})
target_link_libraries(aggregation PUBLIC Modeling Sampler)

add_executable(dump_barplot ${evaluator_SOURCES} ${evaluator_HEADERS} evaluator/plot/dump_barplot.cpp)
target_include_directories(dump_barplot PUBLIC ${evaluator_INCLUDE})
target_link_libraries(dump_barplot PUBLIC Modeling Sampler)

add_executable(dump_distribution ${evaluator_SOURCES} ${evaluator_HEADERS} evaluator/plot/dump_distribution.cpp)
target_include_directories(dump_distribution PUBLIC ${evaluator_INCLUDE})
target_link_libraries(dump_distribution PUBLIC Modeling Sampler)

add_executable(dump_lines ${evaluator_SOURCES} ${evaluator_HEADERS} evaluator/plot/dump_lines.cpp)
target_include_directories(dump_lines PUBLIC ${evaluator_INCLUDE})
target_link_libraries(dump_lines PUBLIC Modeling Sampler)

add_executable(dump_meta ${evaluator_SOURCES} ${evaluator_HEADERS} evaluator/plot/dump_meta.cpp)
target_include_directories(dump_meta PUBLIC ${evaluator_INCLUDE})
target_link_libraries(dump_meta PUBLIC Modeling Sampler)

add_executable(dump_single ${evaluator_SOURCES} ${evaluator_HEADERS} evaluator/plot/dump_single.cpp)
target_include_directories(dump_single PUBLIC ${evaluator_INCLUDE})
target_link_libraries(dump_single PUBLIC Modeling Sampler)

add_executable(single_distribution ${evaluator_SOURCES} ${evaluator_HEADERS} evaluator/plot/single_distribution.cpp)
target_include_directories(single_distribution PUBLIC ${evaluator_INCLUDE})
target_link_libraries(single_distribution PUBLIC Modeling Sampler)


add_executable(sample_implicit ${evaluator_SOURCES} ${evaluator_HEADERS} sample_implicit.cpp)
target_include_directories(sample_implicit PUBLIC ${evaluator_INCLUDE})
target_link_libraries(sample_implicit PUBLIC Modeling Sampler)


add_executable(mesh_to_implicit modeling/mesh_to_implicit.cpp)
target_include_directories(mesh_to_implicit PUBLIC ${modeling_INCLUDE})
target_link_libraries(mesh_to_implicit PUBLIC Modeling Sampler)

add_executable(isosurface modeling/isosurface.cpp)
target_include_directories(isosurface PUBLIC ${modeling_INCLUDE})
target_link_libraries(isosurface PUBLIC Modeling Sampler)

add_executable(pc_generator sampler/pc_generator.cpp)
target_include_directories(pc_generator PUBLIC ${sampler_INCLUDE})
target_link_libraries(pc_generator PUBLIC Sampler)

add_executable(uniform_sampler ${evaluator_SOURCES}  ${evaluator_HEADERS} ${Modeling_HEADERS} ${sampler_HEADERS} sampler/run_uniform.cpp)
target_include_directories(uniform_sampler PUBLIC ${sampler_INCLUDE} ${Modeling_INCLUDE} ${evaluator_INCLUDE})
target_link_libraries(uniform_sampler PUBLIC Modeling Sampler)
