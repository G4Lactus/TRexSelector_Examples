# cpp/ml_methods/validation/hac_clustering/CMakeLists.txt
# ==================================================================================
# Validation simulations for HAC clustering.
# ==================================================================================

add_demo_target(validation_mlm_hac_01_rcompare
    validation_mlm_hac_01_rcompare.cpp
    DATA_PREPROCESSING_LIBS
)

# Output files (if any) go to validation_results/ instead of simulation_results/.
get_target_property(_defs validation_mlm_hac_01_rcompare COMPILE_DEFINITIONS)
string(REPLACE "simulation_results/" "validation_results/" _defs "${_defs}")
set_target_properties(validation_mlm_hac_01_rcompare PROPERTIES COMPILE_DEFINITIONS "${_defs}")

add_demo_target(validation_mlm_hac_02_gvs_dummy_rcompare
    validation_mlm_hac_02_gvs_dummy_rcompare.cpp
    DATA_PREPROCESSING_LIBS
)

# Output files (if any) go to validation_results/ instead of simulation_results/.
get_target_property(_defs validation_mlm_hac_02_gvs_dummy_rcompare COMPILE_DEFINITIONS)
string(REPLACE "simulation_results/" "validation_results/" _defs "${_defs}")
set_target_properties(validation_mlm_hac_02_gvs_dummy_rcompare PROPERTIES COMPILE_DEFINITIONS "${_defs}")
