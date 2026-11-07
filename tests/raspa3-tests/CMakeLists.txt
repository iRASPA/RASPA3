add_executable(unit_tests_raspa3
  mc_component_moves.cpp
  mc_fractional_molecule_order.cpp
  mc_gibbs_drift.cpp
  mc_identity_switch_drift.cpp
  mc_thermodynamic_integration_drift.cpp
  mc_muvt_drift.cpp
  mc_multi_trial_cbmc.cpp
  mc_muvt_grids_drift.cpp
  mc_npt_drift.cpp
  mc_nvt_drift.cpp
  mc_nvt_polarization_drift.cpp
  mc_nvt_grids_drift.cpp
  mc_reaction_drift.cpp
  mc_semi_flexible_cbmc.cpp
  mc_cyclic_cbmc.cpp
  mc_u0_sanity.cpp
  mc_concerted_rotation.cpp
  mc_pivot_cbmc.cpp
  mc_double_bridging.cpp
  mc_driver_bookkeeping.cpp
  mc_acceptance_log_space.cpp
  mc_bead_moves.cpp
  mc_cross_links.cpp
  md_parallel_tempering.cpp
  cbmc_api.cpp
  cbmc_base_sampler.cpp
  cbmc_fixed_endpoint.cpp
  cbmc_recoil_growth.cpp
  end_to_end_histogram.cpp
  molecule_shape.cpp
  molecule_backbone.cpp
  transition_matrix.cpp
  input_reader_tmmc.cpp
  reweighted_histogram_wham.cpp
  parallel_tmmc_windows.cpp
  equation_of_state_critical_constants.cpp
  macrostate_phase_split.cpp
  tmmc_logpi_reconstruction.cpp
  dctst_mirror_reinsertion.cpp
  dctst_polarization_drift.cpp
  tmmc_component_tracking.cpp
  tmmc_restart_state.cpp
  parallel_tmmc_wang_landau_switch.cpp
  transition_matrix_bookkeeping.cpp
  system_mole_fractions.cpp
  mc_omit_inter_interactions_drift.cpp
  group_fractional_slots.cpp
  serial_reaction_fractional_layout.cpp
  swap_cfcmc_blocked_pockets.cpp
  enthalpy_of_adsorption.cpp
  lambda_histogram_free_energy.cpp
  cbmc_first_bead_log_space.cpp
  mc_cfcmc_flexible_charged_drift.cpp
  parallel_tempering_unsupported.cpp
  cbmc_thread_pool_grids.cpp
  gibbs_empty_box_statistics.cpp
  mc_move_probabilities_normalized.cpp
  binary_restart_round_trip.cpp
  external_field_grid.cpp
  simulationbox_block_average.cpp
  uninitialized_members_and_table_bounds.cpp
  spatial_decomposition_bond_bend_split.cpp
  polarization_minimization_derivatives.cpp
  mixed_framework_rigid_group_exclusions.cpp
  intramolecular_potential_fixes.cpp
  dctst_transmission_polarization_drift.cpp
  tmmc_property_output.cpp
  main.cpp)

raspa_target_link_libraries(unit_tests_raspa3)
target_link_libraries(unit_tests_raspa3 PUBLIC
                      raspakit
                      GTest::gtest
                      GTest::gtest_main
                      GTest::gmock
                      GTest::gmock_main)

gtest_add_tests(TARGET unit_tests_raspa3 TEST_SUFFIX .noArgs TEST_LIST noArgsTests)

if(BUILD_MAC_PACKAGE)
  install(CODE "execute_process(COMMAND
      codesign --force --options runtime --timestamp --sign \"Developer ID Application: David Dubbeldam (24U2ZRZ6SC)\"  
      \"tests/raspa3-tests/unit_tests_raspa3\")" COMPONENT unit_tests)
endif()

install(TARGETS unit_tests_raspa3 DESTINATION share/raspa3/tests COMPONENT unit_tests)
