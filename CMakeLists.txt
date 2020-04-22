## Project settings
#
cmake_minimum_required(VERSION 3.10)
project(privateer C CXX)
option(BUILD_SHARED_LIBS "Build as a shared library" ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT DEFINED LIB_INSTALL_DIR)
    set(LIB_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)
endif()

if (DEVKIT)
  message(STATUS "Building with system-wide Python3")
endif()

## Required CCP4 CCP4_LIBRARIES
#
find_library(CLIPPERCORELIB clipper-core)
if(NOT CLIPPERCORELIB)
  message(FATAL_ERROR "clipper-core library not found")
endif()
if(CLIPPERCORELIB)
  message(FATAL_ERROR "clipper-core library found")
endif()