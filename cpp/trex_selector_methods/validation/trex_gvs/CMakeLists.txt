# src/demos/trex_selector_methods/validation/trex_gvs/CMakeLists.txt
# ==================================================================================
# Validation simulations for the GVS T-Rex selector.
# ==================================================================================

# trex_gvs_dgps.hpp and trex_gvs_simulation_utils.hpp live in the sibling trex_gvs/ folder.
set(ADDITIONAL_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/../../trex_gvs")

add_demo_target(validation_trex_gvs_01_scaling_solver_comparison
    validation_trex_gvs_01_scaling_solver_comparison.cpp
    T_REX_SELECTOR_LIBS
)

# Output files (if any) go to validation_results/ instead of simulation_results/.
get_target_property(_defs validation_trex_gvs_01_scaling_solver_comparison COMPILE_DEFINITIONS)
string(REPLACE "simulation_results/" "validation_results/" _defs "${_defs}")
set_target_properties(validation_trex_gvs_01_scaling_solver_comparison PROPERTIES COMPILE_DEFINITIONS "${_defs}")
