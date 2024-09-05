# Install script for directory: /Users/dialpuri/Development/nucleofind/package/clipper/core

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/dialpuri/Development/nucleofind/package/clipper/core/libclipper-core.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libclipper-core.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libclipper-core.a")
    execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libclipper-core.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/clipper/core" TYPE FILE FILES
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/container_hkl.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/fftmap.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/hkl_compute.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/test_data.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/clipper_test.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/coords.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/clipper_sysdep.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/nxmap_operator.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/derivs.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/spacegroup_data.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/clipper_precision.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/container_map.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/container_types.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/hkl_lookup.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/container.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/clipper_message.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/cell.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/resol_basisfn.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/hkl_data.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/test_core.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/rotation.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/hkl_operators.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/fftmap_sparse.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/clipper_types.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/clipper_memory.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/symop.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/resol_targetfn.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/clipper_stats.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/clipper_thread.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/spacegroup.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/clipper_instance.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/ramachandran.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/hkl_info.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/hkl_datatypes.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/xmap.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/nxmap.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/map_interp.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/map_utils.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/resol_fn.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/atomsf.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/core/clipper_util.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/clipper" TYPE FILE FILES
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/clipper-cctbx.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/cctbx/clipper_cctbx.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/clipper-mmdb.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/clipper-cif.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/clipper-mmdbold.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/clipper-ccp4.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/minimal-clipper-hkl.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/clipper-phs.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/minimal-clipper-map.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/clipper-minimol.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/clipper.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/clipper-cns.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/clipper-contrib.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/cif/cif_data_io.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper//clipper/phs/phs_io.h"
    )
endif()

