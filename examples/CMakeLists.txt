configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/1_mc_methane_in_box/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/2_mc_co2_n2_in_two_independent_boxes/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/3_mc_binarymixture_in_two_boxes/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/4_mc_binary_mixture_propane_butane_in_box/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/5_md_methane_in_box_msd/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/6_mc_enthalpy_of_adsorption_in_mfi_at_zero_loading/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/7_mc_henry_coefficient_of_methane_in_mfi/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/8_mc_adsorption_of_methane_in_mfi/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/9_mc_adsorption_of_butane_in_mfi/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/10_mc_adorption_co2_in_mfi/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/11_mc_adsorption_of_co2_in_cu-btc/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/12_mc_henry_coefficient_of_co2_n2_methane_in_mfi/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/13_md_rdf_water_box/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/14_mc_rdf_water_box/run)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/non-basic/1_mc_adsorption_binary_mixture_co2_ch4_in_irmof-1/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/non-basic/2_mc_npt_co2/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/non-basic/3_mc_npt_propane/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/non-basic/4_mc_gibbs_co2/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/non-basic/5_md_benzene_in_irmof-1/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/non-basic/6_mc_adsorption_co2_in_lta4a_sodium/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/non-basic/7_md_diffusion_co2_in_lta4a_sodium/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/non-basic/8_mc_adsorption_c6_isomers_in_mfi/run)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/advanced/1_mc_gibbs_co2/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/advanced/2_mc_cfcmc_co2_in_mfi/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/advanced/4_hmc_cfcmc_co2_in_mfi/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/advanced/5_tmmc_methane_in_tobacco-667/0/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/advanced/5_tmmc_methane_in_tobacco-667/1/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/advanced/5_tmmc_methane_in_tobacco-667/2/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/advanced/5_tmmc_methane_in_tobacco-667/3/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/advanced/5_tmmc_methane_in_tobacco-667/4/run)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/auxiliary/1_mc_ideal_gas_rosenbluth_weight_butane/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/auxiliary/2_mc_ideal_gas_rosenbluth_weight_c5-c9-573k/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/auxiliary/3_mc_ideal_gas_rosenbluth_weights_c6_isomers-433k/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/auxiliary/4_charge_equilibration/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/auxiliary/5_make_grids/run)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/reduced_units/1_particles_in_box_nvt/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/reduced_units/2_particles_in_box_nvt_cfcmc/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/reduced_units/3_gibbs/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/reduced_units/4_gibbs_cfcmc/run)


install(DIRECTORY basic DESTINATION share/raspa3/examples 
        COMPONENT examples
        PATTERN "run" EXCLUDE
        PERMISSIONS OWNER_READ GROUP_READ)

install(DIRECTORY basic DESTINATION share/raspa3/examples 
        COMPONENT examples
        PATTERN "run" 
        PERMISSIONS  OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)

install(DIRECTORY non-basic DESTINATION share/raspa3/examples 
        COMPONENT examples
        PATTERN "run" EXCLUDE
        PERMISSIONS OWNER_READ GROUP_READ)

install(DIRECTORY non-basic DESTINATION share/raspa3/examples 
        COMPONENT examples
        PATTERN "run" 
        PERMISSIONS  OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)

install(DIRECTORY advanced DESTINATION share/raspa3/examples 
        COMPONENT examples
        PATTERN "run" EXCLUDE
        PERMISSIONS OWNER_READ GROUP_READ)

install(DIRECTORY advanced DESTINATION share/raspa3/examples 
        COMPONENT examples
        PATTERN "run" 
        PERMISSIONS  OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)

install(DIRECTORY auxiliary DESTINATION share/raspa3/examples 
        COMPONENT examples
        PATTERN "run" EXCLUDE
        PERMISSIONS OWNER_READ GROUP_READ)

install(DIRECTORY auxiliary DESTINATION share/raspa3/examples 
        COMPONENT examples
        PATTERN "run" 
        PERMISSIONS  OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)

install(DIRECTORY reduced_units DESTINATION share/raspa3/examples 
        COMPONENT examples
        PATTERN "run" EXCLUDE
        PERMISSIONS OWNER_READ GROUP_READ)

install(DIRECTORY reduced_units DESTINATION share/raspa3/examples 
        COMPONENT examples
        PATTERN "run" 
        PERMISSIONS  OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
