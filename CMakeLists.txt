cmake_minimum_required(VERSION 2.8)
project(chrono_models)
include(SBELUtils.cmake)

#set path for local find modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

SET(ENABLE_UNIT_GPU FALSE   CACHE   BOOL    "Turn on for GPU models")
SET(ENABLE_UNIT_MPI FALSE   CACHE   BOOL    "Turn on for GPU models")

# Additional packages
enable_openmp_support()
#enable_mpi_support()
#FIND_PACKAGE(hdf5)
FIND_PACKAGE(GLM)
# User specified path to ChronoEngine
SET(CHRONO_INC       ""   CACHE PATH "Where is ChronoEngine Source located")
SET(CHRONO_PARALLEL_INC       ""   CACHE PATH "Where is ChronoEngine Source located")
SET(CHRONO_OPENGL_INC  ""   CACHE PATH "Where is CHRONOOPENGL located")

# Include directories based on user path to ChronoEngine
SET(INC_DIRS ${CHRONO_INC} "${CHRONO_INC}/collision" "${CHRONO_INC}/collision/bullet" ${CHRONO_PARALLEL_INC} ${CHRONO_OPENGL_INC} ${HDF5_INCLUDE_DIRS})

# User specified path to ChronoEngine Libraries
SET(CHRONO_LIB_PATH  ""   CACHE PATH "Where is CHRONOLIB located")
SET(CHRONO_PARALLEL_LIB_PATH  ""   CACHE PATH "Where is CHRONOPARALLEL located")
SET(CHRONO_OPENGL_LIB_PATH  ""   CACHE PATH "Where is CHRONOOPENGL located")
# Find required Libraries
FIND_LIBRARY(CHRONO_LIB NAMES ChronoEngine PATHS ${CHRONO_LIB_PATH} REQUIRED)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	SET(GCC_OPT "sse" CACHE STRING "gcc optimizations")
	SET_PROPERTY(CACHE GCC_OPT PROPERTY STRINGS "mmx" "sse" "sse2" "sse3" "ssse3" "sse4.1" "sse4.2" "avx")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m${GCC_OPT}")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m${GCC_OPT}")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	SET(INTEL_OPT "SSE2" CACHE STRING "intel optimizations")
	SET_PROPERTY(CACHE INTEL_OPT PROPERTY STRINGS "SSE2" "SSE3" "SSSE3" "SSE4.1" "SSE4.2" "AVX" "CORE-AVX-I" "CORE-AVX2")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -x${INTEL_OPT}")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -x${INTEL_OPT}")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()

# Set unit-specific options
IF(ENABLE_UNIT_GPU)
    enable_cuda_support()
    LIST(APPEND INC_DIRS "${CUDA_TOOLKIT_ROOT_DIR}/include" "${CUDA_SDK_ROOT_DIR}/common/inc" )
    FIND_LIBRARY(CHRONO_LIB_PARALLEL NAMES ChronoEngine_Parallel PATHS ${CHRONO_PARALLEL_LIB_PATH} REQUIRED)
    #FIND_LIBRARY(CHRONO_LIB_POSTPROCESS NAMES ChronoEngine_POSTPROCESS PATHS ${CHRONO_LIB_PATH} REQUIRED)
    FIND_LIBRARY(CHRONO_LIB_OPENGL NAMES ChronoEngine_OPENGL PATHS ${CHRONO_OPENGL_LIB_PATH} REQUIRED)

    SET(GPU_LIBS
        #${CHRONO_LIB_POSTPROCESS}
        ${CHRONO_LIB_PARALLEL}
        ${CHRONO_LIB_OPENGL}
        ${HDF5_LIBRARIES}
        ${MPI_LIBRARIES}
    )
ENDIF()

IF(ENABLE_UNIT_MPI)
    enable_mpi_support()

    FIND_LIBRARY(CHRONO_LIB_MPI NAMES ChronoEngine_MPI PATHS ${CHRONO_LIB_PATH} REQUIRED)

    SET(MPI_LIBS
        ${CHRONO_LIB_MPI}
    )
ENDIF()

include_directories(${INC_DIRS} ${GLM_INCLUDE_DIRS})

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ")

ADD_EXECUTABLE(demo_m113 demo_m113.cpp)      
TARGET_LINK_LIBRARIES(demo_m113 ${CHRONO_LIB} ${GPU_LIBS})

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
