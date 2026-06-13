# TRexSelector_Simulations/cpp/CMakeLists.txt
# ===================================================================================
# CMake configurations for T-Algorithm and TRex demo/simulation executables.
# Links against the installed TRexSelector library via namespaced targets.
#
# Organization:
# -------------
# |-- cpp/
# |-- applications/
# |-- ml_methods/
# |-- trex_selector_methods/
# |-- tsolvers/
# ===================================================================================


# ===================================================================================
# 1. Library Sets
# ===================================================================================
# Base libraries for all demo executables.
set(BASE_LIBS
    Boost::headers
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
)

# Libraries for Data Preprocessing demos
set(DATA_PREPROCESSING_LIBS
    ${BASE_LIBS}
    TRexSelector::trex_ml_methods
)

# Libraries for T-Algorithm based demos
set(T_ALGORITHM_LIBS
    ${BASE_LIBS}
    TRexSelector::trex_ml_methods
    TRexSelector::trex_tsolvers
)

# Libraries for TRex Selector based demos
# (trex_selector_methods transitively pulls in tsolvers, ml_methods, utils)
set(T_REX_SELECTOR_LIBS
    ${BASE_LIBS}
    TRexSelector::trex_selector_methods
)


# ==================================================================================
# Macro: add_demo_target
# ==================================================================================
# Creates a demo executable with proper include paths, linking, and compile flags.
#
# Usage:
#   add_demo_target(<target_name> <demo_file> <library_set> [additional_libraries...])
#
# Arguments:
#   target_name          - Name of the demo executable (e.g., demo_ts_01_tlars)
#   demo_file            - Source file path relative to current CMakeLists.txt
#   library_set          - Library variable (DATA_PREPROCESSING_LIBS, T_ALGORITHM_LIBS,
#                          T_REX_SELECTOR_LIBS)
#   additional_libraries - (Optional) Additional libraries to link against
#
# ==================================================================================
macro(add_demo_target target_name demo_file library_set)

    add_executable(${target_name} ${demo_file})

    # Include directories: installed headers are carried via INTERFACE_INCLUDE_DIRECTORIES
    # on the TRexSelector:: targets. Additional paths can be injected via ADDITIONAL_INCLUDES.
    target_include_directories(${target_name} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${ADDITIONAL_INCLUDES}
    )

    target_link_libraries(${target_name} PRIVATE
        Eigen3::Eigen
        OpenMP::OpenMP_CXX
        ${${library_set}}
        ${ARGN}
    )

    target_compile_definitions(${target_name} PRIVATE
        EIGEN_USE_BLAS
    )

    # NOTE: -ffast-math implies -ffinite-math-only, which breaks std::isnan/isfinite.
    target_compile_options(${target_name} PRIVATE
        -Wno-deprecated-declarations
        $<$<CONFIG:Debug>:-g -O0 -Wall -Wextra>
        $<$<CONFIG:Release>:-O3 -DNDEBUG -march=native -ffast-math -funroll-loops>
    )

    set(ADDITIONAL_INCLUDES "")
endmacro()
# ==================================================================================


# ===================================================================================
# Demo/Simulation Targets
# ===================================================================================
add_subdirectory(ml_methods)
add_subdirectory(tsolvers)
add_subdirectory(trex_selector_methods)
add_subdirectory(applications)
#
