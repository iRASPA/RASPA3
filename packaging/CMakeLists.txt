
if(BUILD_MAC_PACKAGE)

set(CPACK_PACKAGE_VERSION ${CMAKE_PROJECT_VERSION})
set(CPACK_GENERATOR "productbuild")
set(CPACK_PACKAGE_NAME "raspa")
set(CPACK_PRODUCTBUILD_IDENTIFIER nl.darkwing.raspa)
set(CPACK_PACKAGING_INSTALL_PREFIX /usr/local)
set(CPACK_PACKAGE_DESCRIPTION "Molecular simulation code for adsorption and diffusion in nanoporours materials.")
set(CPACK_PACKAGE_CONTACT "David Dubbeldam (d.dubbeldam@uva.nl)")
set(CPACK_PACKAGE_VENDOR "David Dubbeldam, Sofia Calero, Thijs J.H. Vlugt, R.Q. Snurr")
set(CPACK_PACKAGE_FILE_NAME "raspa-${CMAKE_PROJECT_VERSION}-${CMAKE_OSX_ARCHITECTURES}")
set(CPACK_PKGBUILD_IDENTITY_NAME "3rd Party Mac Developer Installer: David Dubbeldam (24U2ZRZ6SC)")
set(CPACK_RESOURCE_FILE_WELCOME ${CMAKE_SOURCE_DIR}/packaging/welcome.txt)
set(CPACK_RESOURCE_FILE_README ${CMAKE_SOURCE_DIR}/packaging/readme.txt)
set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_SOURCE_DIR}/packaging/license.txt)
set(CPACK_STRIP_FILES ON)

include(CPack)

endif()


if(BUILD_CENTOS6_PACKAGE)

set(CPACK_PACKAGE_VERSION ${VERSION})
set(CPACK_GENERATOR "RPM")
set(CPACK_RPM_PACKAGE_NAME "raspa")
set(CPACK_RPM_PACKAGE_LICENSE "MIT")
set(CPACK_RPM_PACKAGE_DESCRIPTION "Molecular simulation code for adsorption and diffusion in nanoporours materials.")
set(CPACK_RPM_FILE_NAME "RPM-DEFAULT")
set(CPACK_RPM_PACKAGE_RELEASE_DIST ON)
set(CPACK_RPM_PACKAGE_RELEASE 1)
set(CPACK_RPM_PACKAGE_CONTACT "David Dubbeldam (d.dubbeldam@uva.nl)")
set(CPACK_RPM_PACKAGE_VENDOR "David Dubbeldam, Sofia Calero, Thijs J.H. Vlugt, R.Q. Snurr")
set(CPACK_RPM_PACKAGE_DEPENDS "blas, lapack, libgfortran, libquadmath, libgomp, glibc, libgcc")
set(CPACK_STRIP_FILES ON)

include(CPack)

endif()


if(BUILD_RPM_PACKAGE)

set(CPACK_PACKAGE_VERSION ${VERSION})
set(CPACK_GENERATOR "RPM")
set(CPACK_RPM_PACKAGE_NAME "raspa")
set(CPACK_RPM_PACKAGE_LICENSE "MIT")
set(CPACK_RPM_PACKAGE_DESCRIPTION "Molecular simulation code for adsorption and diffusion in nanoporours materials.")
set(CPACK_RPM_FILE_NAME "RPM-DEFAULT")
set(CPACK_RPM_PACKAGE_RELEASE_DIST ON)
set(CPACK_RPM_PACKAGE_RELEASE 1)
set(CPACK_RPM_PACKAGE_CONTACT "David Dubbeldam (d.dubbeldam@uva.nl)")
set(CPACK_RPM_PACKAGE_VENDOR "David Dubbeldam, Sofia Calero, Thijs J.H. Vlugt, R.Q. Snurr")
set(CPACK_RPM_PACKAGE_DEPENDS "blas, lapack, libgfortran, libquadmath, libomp, glibc, libgcc")
set(CPACK_STRIP_FILES ON)

include(CPack)

endif()


if(BUILD_DEB_PACKAGE)

set(CPACK_GENERATOR "DEB")
set(CPACK_DEBIAN_PACKAGE_NAME "raspa")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "d.dubbeldam@uva.nl")
set(CPACK_PACKAGE_NAME ${CMAKE_PROJECT_NAME})
set(CPACK_set_DESTDIR TRUE)
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE.md")
set(CPACK_DEBIAN_PACKAGE_VERSION ${PACKAGE_VERSION})
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libblas64-3, liblapack64-3, libquadmath0, libgfortran5, libgcc-s1, libomp5, libc6")
set(CPACK_PACKAGE_VENDOR "David Dubbeldam, Sofia Calero, Thijs J.H. Vlugt, R.Q. Snurr")
set(CPACK_STRIP_FILES ON)

include(CPack)

endif()
