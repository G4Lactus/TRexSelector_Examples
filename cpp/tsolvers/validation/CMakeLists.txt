# cpp/tsolvers/validation/CMakeLists.txt
# ==================================================================================
# Validation simulations: correctness and invariance checks for T-Algorithm solvers.
# ==================================================================================

add_demo_target(validation_ts_01_tenet_aug_comparison
    validation_ts_01_tenet_aug_comparison.cpp
    T_ALGORITHM_LIBS
)

# Output files (if any) go to validation_results/ instead of simulation_results/.
get_target_property(_defs validation_ts_01_tenet_aug_comparison COMPILE_DEFINITIONS)
string(REPLACE "simulation_results/" "validation_results/" _defs "${_defs}")
set_target_properties(validation_ts_01_tenet_aug_comparison PROPERTIES COMPILE_DEFINITIONS "${_defs}")

add_demo_target(validation_ts_02_tlars_tlasso_rcompare
    validation_ts_02_tlars_tlasso_rcompare.cpp
    T_ALGORITHM_LIBS
)

# Output files (if any) go to validation_results/ instead of simulation_results/.
get_target_property(_defs validation_ts_02_tlars_tlasso_rcompare COMPILE_DEFINITIONS)
string(REPLACE "simulation_results/" "validation_results/" _defs "${_defs}")
set_target_properties(validation_ts_02_tlars_tlasso_rcompare PROPERTIES COMPILE_DEFINITIONS "${_defs}")
