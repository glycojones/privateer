# Install script for directory: /Users/dialpuri/Development/nucleofind/package/clipper/minimol

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/dialpuri/Development/nucleofind/package/clipper/minimol/libclipper-minimol.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libclipper-minimol.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libclipper-minimol.a")
    execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libclipper-minimol.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/clipper/minimol" TYPE FILE FILES
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper/clipper/minimol/minimol_data.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper/clipper/minimol/minimol.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper/clipper/minimol/container_minimol.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper/clipper/minimol/minimol_io_gemmi.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper/clipper/minimol/minimol_io_mmdb.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper/clipper/minimol/minimol_io_seq.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper/clipper/minimol/minimol_seq.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper/clipper/minimol/minimol_utils.h"
    "/Users/dialpuri/Development/nucleofind/package/checkout/clipper/clipper/clipper-minimol.h"
    )
endif()

