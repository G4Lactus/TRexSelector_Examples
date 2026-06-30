# src/demos/trex_selector_methods/validation/trex/CMakeLists.txt
# ==================================================================================
# Validation simulations for the classical T-Rex selector.
# ==================================================================================

# trex_sim_utils.hpp lives in the sibling trex/ folder.
set(ADDITIONAL_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/../../trex")

add_demo_target(validation_trex_01_scaling_comparison
    validation_trex_01_scaling_comparison.cpp
    T_REX_SELECTOR_LIBS
)

# Output files (if any) go to validation_results/ instead of simulation_results/.
get_target_property(_defs validation_trex_01_scaling_comparison COMPILE_DEFINITIONS)
string(REPLACE "simulation_results/" "validation_results/" _defs "${_defs}")
set_target_properties(validation_trex_01_scaling_comparison PROPERTIES COMPILE_DEFINITIONS "${_defs}")
