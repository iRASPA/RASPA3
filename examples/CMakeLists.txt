configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/1_mc_methane_in_box/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/2_mc_co2_n2_in_two_independent_boxes/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/4_mc_binarymixture_in_two_boxes/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/6_enthalpy_of_adsorption_in_MFI_at_zero_loading/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/7_adsorption_of_methane_in_MFI/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/basic/8_adsorption_of_co2_in_cu-btc/run)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/non-basic/2_npt_mc_methane-btc/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/non-basic/7_gibbs_co2-btc/run)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/advanced/7_gibbs_co2/run)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/run.cmake.in ${CMAKE_CURRENT_SOURCE_DIR}/advanced/9_cfcmc_co2_in_MFI/run)


install(DIRECTORY basic DESTINATION share/raspa3/examples 
        COMPONENT examples-basic
        PATTERN "run" EXCLUDE
        PERMISSIONS OWNER_READ GROUP_READ)

install(DIRECTORY basic DESTINATION share/raspa3/examples 
        COMPONENT examples-basic
        PATTERN "run" 
        PERMISSIONS  OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
