# =============================================================================
# TRexSelector_Simulations — Root CMakeLists.txt
# =============================================================================
# Standalone CMake project for C++ demo and simulation executables.
# Consumes the installed TRexSelector library via find_package().
#
# Prerequisites:
#   cmake --install <TRexSelector-build-dir> --prefix <prefix>
#   cmake -S . -B build/debug -DCMAKE_PREFIX_PATH=<prefix>
#
# The CMakePresets.json in this directory sets CMAKE_PREFIX_PATH to
# ../TRexSelector/install by default.
# =============================================================================
cmake_minimum_required(VERSION 3.23...4.3)


# ============================================================================
# Compiler Selection (must match TRexSelector — Homebrew LLVM on macOS)
# ============================================================================
if(APPLE AND NOT DEFINED CMAKE_CXX_COMPILER)
    set(HOMEBREW_LLVM_PREFIX "/opt/homebrew/opt/llvm")
    if(EXISTS "${HOMEBREW_LLVM_PREFIX}/bin/clang++")
        set(CMAKE_C_COMPILER   "${HOMEBREW_LLVM_PREFIX}/bin/clang"   CACHE FILEPATH "C compiler")
        set(CMAKE_CXX_COMPILER "${HOMEBREW_LLVM_PREFIX}/bin/clang++" CACHE FILEPATH "C++ compiler")
        message(STATUS "*** Using Homebrew LLVM: ${CMAKE_CXX_COMPILER}")
    else()
        message(WARNING "*** Homebrew LLVM not found — falling back to system compiler")
    endif()
endif()


# ============================================================================
# Output Directories
# ============================================================================
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)


# ============================================================================
# Project Definition
# ============================================================================
project(TRexSelector_Simulations
    VERSION 1.0.0
    DESCRIPTION "Demo and simulation executables for the TRexSelector library"
    LANGUAGES CXX
)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# ============================================================================
# Find TRexSelector (installed library)
# ============================================================================
find_package(TRexSelector REQUIRED)


# ============================================================================
# Find remaining dependencies
# ============================================================================

# OpenMP
if(APPLE)
    if(NOT DEFINED HOMEBREW_LLVM_PREFIX)
        set(HOMEBREW_LLVM_PREFIX "/opt/homebrew/opt/llvm")
    endif()
    set(OpenMP_CXX_FLAGS       "-fopenmp=libomp -I${HOMEBREW_LLVM_PREFIX}/include")
    set(OpenMP_CXX_LIB_NAMES   "libomp")
    set(OpenMP_libomp_LIBRARY  "${HOMEBREW_LLVM_PREFIX}/lib/libomp.dylib"
        CACHE FILEPATH "libomp path" FORCE)
endif()
find_package(OpenMP REQUIRED)

# Eigen3
find_package(Eigen3 REQUIRED NO_MODULE)

# Boost
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()
find_package(Boost 1.80 REQUIRED COMPONENTS headers)

# BLAS/LAPACK
if(APPLE)
    set(BLA_VENDOR Apple)
    find_package(BLAS)
    find_package(LAPACK)
    if(NOT BLAS_FOUND)
        set(BLA_VENDOR OpenBLAS)
        find_package(BLAS REQUIRED)
        find_package(LAPACK REQUIRED)
    endif()
else()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif()

# Matplot++ (only needed here, no longer in the core library)
find_package(Matplot++ REQUIRED)


# ============================================================================
# Configuration Summary
# ============================================================================
message(STATUS "==============================================")
message(STATUS "TRexSelector_Simulations Build Configuration")
message(STATUS "==============================================")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler:   ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "")
message(STATUS "Libraries:")
if(TRexSelector_FOUND)
    message(STATUS "  - TRexSelector: ${TRexSelector_VERSION} (${TRexSelector_DIR})")
else()
    message(STATUS "  - TRexSelector: NOT FOUND")
endif()
if(Matplot++_FOUND)
    message(STATUS "  - Matplot++: ${Matplot++_VERSION}")
else()
    message(STATUS "  - Matplot++: NOT FOUND")
endif()
if(Eigen3_FOUND)
    message(STATUS "  - Eigen3: ${Eigen3_VERSION}")
endif()
if(BLAS_FOUND)
    message(STATUS "  - BLAS Vendor: ${BLA_VENDOR}")
endif()
message(STATUS "==============================================")


# ============================================================================
# C++ demo and simulation executables
# ============================================================================
add_subdirectory(cpp)
