# src/demos/trex_selector_methods/validation/trex_da/CMakeLists.txt
# ==================================================================================
# Validation simulations for the Dependency-Aware T-Rex selector.
# ==================================================================================

# trex_da_sim_common.hpp lives in the sibling trex_da/ folder.
set(ADDITIONAL_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/../../trex_da")

add_demo_target(validation_trex_da_01_bt_dendro_diag
    validation_trex_da_01_bt_dendro_diag/validation_trex_da_01_bt_dendro_diag.cpp
    T_REX_SELECTOR_LIBS
)

# Output files go to validation_results/ instead of simulation_results/.
get_target_property(_defs validation_trex_da_01_bt_dendro_diag COMPILE_DEFINITIONS)
string(REPLACE "simulation_results/" "validation_results/" _defs "${_defs}")
set_target_properties(validation_trex_da_01_bt_dendro_diag PROPERTIES COMPILE_DEFINITIONS "${_defs}")
