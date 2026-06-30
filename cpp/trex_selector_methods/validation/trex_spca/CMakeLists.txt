# src/demos/trex_selector_methods/validation/trex_spca/CMakeLists.txt
# ==================================================================================
# Validation simulations for the Sparse PCA T-Rex selector.
# ==================================================================================

# trex_spca_sim_utils.hpp lives in the sibling trex_spca/ folder.
# NOTE: add_demo_target resets ADDITIONAL_INCLUDES at the end, so it must be set
#       before every call.

# TODO: validation_trex_spca_01_lambda2_probe depends on removed API:
#   ridge_cv_glmnet, ridge_cv_direct, elastic_net_cv_gaussian
#   (all superseded by enet_cv_ccd / ridge_cv_svd). Disabled until rewritten.
# set(ADDITIONAL_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/../../trex_spca")
# add_demo_target(validation_trex_spca_01_lambda2_probe
#     validation_trex_spca_01_lambda2_probe.cpp
#     T_REX_SELECTOR_LIBS
# )

set(ADDITIONAL_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/../../trex_spca")
add_demo_target(validation_trex_spca_02_solver_comparison
    validation_trex_spca_02_solver_comparison.cpp
    T_REX_SELECTOR_LIBS
)

set(ADDITIONAL_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/../../trex_spca")
add_demo_target(validation_trex_spca_03_scaling_comparison
    validation_trex_spca_03_scaling_comparison.cpp
    T_REX_SELECTOR_LIBS
)

set(ADDITIONAL_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/../../trex_spca")
add_demo_target(validation_trex_spca_04_rdump_pipeline
    validation_trex_spca_04_rdump_pipeline.cpp
    T_REX_SELECTOR_LIBS
)

set(ADDITIONAL_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/../../trex_spca")
add_demo_target(validation_trex_spca_05_lambda2_foldmatch
    validation_trex_spca_05_lambda2_foldmatch.cpp
    T_REX_SELECTOR_LIBS
)

# All targets: redirect DEMO_OUTPUT_DIR to validation_results/.
foreach(_t
    validation_trex_spca_02_solver_comparison
    validation_trex_spca_03_scaling_comparison
    validation_trex_spca_04_rdump_pipeline
    validation_trex_spca_05_lambda2_foldmatch
)
    get_target_property(_defs ${_t} COMPILE_DEFINITIONS)
    string(REPLACE "simulation_results/" "validation_results/" _defs "${_defs}")
    set_target_properties(${_t} PROPERTIES COMPILE_DEFINITIONS "${_defs}")
endforeach()
