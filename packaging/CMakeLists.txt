
if(BUILD_PACKAGE)

# rpm-settings
set(CPACK_PACKAGE_VERSION ${VERSION})
set(CPACK_GENERATOR RPM)
set(CPACK_RPM_PACKAGE_NAME "raspa3")
set(CPACK_RPM_PACKAGE_LICENSE "MIT")
set(CPACK_RPM_PACKAGE_DESCRIPTION "Molecular simulation code for adsorption and diffusion in nanoporours materials.")
set(CPACK_RPM_FILE_NAME "RPM-DEFAULT")
set(CPACK_RPM_PACKAGE_RELEASE_DIST ON)
set(CPACK_RPM_PACKAGE_RELEASE 1)
set(CPACK_RPM_PACKAGE_CONTACT "David Dubbeldam (d.dubbeldam@uva.nl)")
set(CPACK_RPM_PACKAGE_VENDOR "David Dubbeldam, Sofia Calero, Thijs J.H. Vlugt, R.Q. Snurr")
set(CPACK_RPM_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})


include(CPack)

endif()
