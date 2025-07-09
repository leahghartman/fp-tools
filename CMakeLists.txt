# CMake version requirement
cmake_minimum_required(VERSION 3.16)

# Project information
project(fp-tools LANGUAGES CXX)
set(FPTOOLS_VERSION_MAJOR 1)
set(FPTOOLS_VERSION_MINOR 0)

add_compile_definitions(SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ========= Configurable build options =========
option(USE_MPI     "Enable MPI support"     OFF)
option(USE_OMP     "Enable OpenMP support"  OFF)
option(USE_CUDA    "Enable CUDA support"    OFF)

# ========= Include directories =========
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/ext)

# ========= Source files and Executable Directory =========
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
file(GLOB_RECURSE SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp)
add_executable(fptools ${SOURCES})

# ========= MPI =========
if(USE_MPI)
    find_package(MPI REQUIRED)
    if(MPI_FOUND)
        message(STATUS "MPI found. Enabling MPI support.")
        add_definitions(-DUSE_MPI)
        include_directories(${MPI_INCLUDE_PATH})
        set(EXTRA_LIBS ${EXTRA_LIBS} MPI::MPI_CXX)
    else()
        message(FATAL_ERROR "USE_MPI is ON but MPI was not found.")
    endif()
endif()

# ========= OpenMP =========
if(USE_OMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_FOUND)
        message(STATUS "OpenMP found. Enabling OpenMP support.")
        add_definitions(-DUSE_OMP)
        set(EXTRA_LIBS ${EXTRA_LIBS} OpenMP::OpenMP_CXX)
    else()
        message(FATAL_ERROR "USE_OMP is ON but OpenMP was not found.")
    endif()
endif()

# ========= CUDA =========
if(USE_CUDA)
    find_package(CUDA REQUIRED)
    if(CUDA_FOUND)
        message(STATUS "CUDA found. Enabling CUDA support.")
        add_definitions(-DUSE_CUDA)
        include_directories(${CUDA_INCLUDE_DIRS})
        set(EXTRA_LIBS ${EXTRA_LIBS} ${CUDA_LIBRARIES})
    else()
        message(FATAL_ERROR "USE_CUDA is ON but CUDA was not found.")
    endif()
endif()

# ========= FFTW =========

#find_package(FFTW REQUIRED COMPONENTS DOUBLE)
#target_link_libraries(fptools PRIVATE FFTW::Double)

# ========= Optimization flags (optional) =========
target_compile_options(fptools PRIVATE -O3 -Wall -Wno-comment)
