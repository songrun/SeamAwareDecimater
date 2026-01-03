cmake_minimum_required(VERSION 3.10)
project(SeamAwareDecimater)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

## We need Eigen
find_package(Eigen REQUIRED)
include_directories( "${Eigen_INCLUDE_DIR}" )

## We need libigl
find_package(LIBIGL REQUIRED)
include_directories( "${LIBIGL_INCLUDE_DIR}" )

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    # Try to find Homebrew libomp on macOS
    if(APPLE)
        execute_process(COMMAND brew --prefix libomp OUTPUT_VARIABLE LIBOMP_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(LIBOMP_PREFIX)
            set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
            set(OpenMP_CXX_LIBRARIES "-L${LIBOMP_PREFIX}/lib -lomp")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
            set(LIBIGL_LIBRARIES ${LIBIGL_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
            message(STATUS "Found Homebrew libomp at ${LIBOMP_PREFIX}")
        endif()
    endif()
endif()

## We don't have/want MOSEK
add_definitions(-DIGL_NO_MOSEK)

link_directories(
	/usr/local/lib
) 

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_TIMING "Enable performance timing" OFF)
if(ENABLE_TIMING)
  add_definitions(-DENABLE_TIMING)
endif()

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug, Release, RelWithDebInfo, MinSizeRel." FORCE)
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  add_compile_options(-stdlib=libc++)
endif()

add_library(DEC_LIBS
    OBJECT
    quadric_error_metric.cpp
    cost_and_placement.cpp
    half_edge.cpp
    neighbor_faces_and_boundary.cpp
    detect_foldover.cpp
    decimate.cpp
    collapse_edge_seam.cpp
    )
    
add_executable(decimater
	decimater.cpp
	$<TARGET_OBJECTS:DEC_LIBS>
	)
target_link_libraries ( decimater
	${LIBIGL_LIBRARIES}
)
