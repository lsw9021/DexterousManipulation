cmake_minimum_required(VERSION 2.8.6)
project(fem)
add_compile_options(-std=gnu++11)
add_compile_options(-Wdeprecated-declarations)
SET(CMAKE_BUILD_TYPE Release CACHE STRING
	"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
#	FORCE
	)

include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR})
find_package(Eigen3 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${OPENGL_INCLUDE_DIR})
include_directories(${GLUT_INCLUDE_DIR})

file(GLOB srcs "*.h" "*.cpp")
file(GLOB cst_srcs "Constraint/*.h" "Constraint/*.cpp")
file(GLOB mesh_srcs "Mesh/*.h" "Mesh/*.cpp")

add_library(fem ${srcs} ${cst_srcs} ${mesh_srcs})
target_link_libraries(fem ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})