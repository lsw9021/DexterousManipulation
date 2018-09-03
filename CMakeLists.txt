cmake_minimum_required(VERSION 2.8.6)
project(volumetric)
add_compile_options(-std=gnu++11)
add_compile_options(-Wdeprecated-declarations)

SET(CMAKE_BUILD_TYPE Release CACHE STRING
	"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	# FORCE
	)
SET( ENABLE_VMCON ON CACHE BOOL "Enable vmcon3d" )
SET( ENABLE_RENDER ON CACHE BOOL "Enable vmcon3d render" )
# SET( USE_JOINT_TORQUE ON CACHE BOOL "Enable joint torque.")
SET( USE_MUSCLE ON CACHE BOOL "Use Muscle.")

add_subdirectory( fem )

include(FindOpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(Eigen3 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(Boost 1.45.0 COMPONENTS filesystem)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${OPENGL_INCLUDE_DIR})
include_directories(${GLUT_INCLUDE_DIR})

include_directories(${CMAKE_SOURCE_DIR})
link_directories(./fem/)

file(GLOB vmcon_srcs "vmcon/*.h" "vmcon/*.cpp" "vmcon/iLQR/*.h" "vmcon/iLQR/*.cpp" "vmcon/FSM/*.h" "vmcon/FSM/*.cpp")

file(GLOB render_srcs "render/*.h" "render/*.cpp" ${vmcon_srcs})
# SET(render_src_without_main ${render_srcs})
list(REMOVE_ITEM render_srcs ${CMAKE_CURRENT_SOURCE_DIR}/vmcon/main.cpp)
# SET(render_srcs ${render_src_without_main})
if (ENABLE_VMCON)
	project(vmcon)
	find_package(DART REQUIRED COMPONENTS utils-urdf gui optimizer-ipopt CONFIG)
	add_compile_options(-DHAVE_CSTDDEF)
if(USE_MUSCLE)
	add_compile_options(-DUSE_MUSCLE)
endif()
	include_directories(${Boost_INCLUDE_DIRS})
	include_directories(${DART_INCLUDE_DIRS})

	add_executable(vmcon ${vmcon_srcs})
	target_link_libraries(vmcon ${DART_LIBRARIES} ${Boost_LIBRARIES} fem)
endif()

if (ENABLE_RENDER)
	project(render)
	find_package(DART REQUIRED COMPONENTS utils-urdf gui optimizer-ipopt CONFIG)
	add_compile_options(-DHAVE_CSTDDEF)
if(USE_MUSCLE)
	add_compile_options(-DUSE_MUSCLE)
endif()
	include_directories(${Boost_INCLUDE_DIRS})
	include_directories(${DART_INCLUDE_DIRS})

	add_executable(render ${render_srcs})
	target_link_libraries(render ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${DART_LIBRARIES} ${Boost_LIBRARIES} fem)
endif()
