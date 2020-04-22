## Project settings
#
cmake_minimum_required(VERSION 3.12)
project(privateer VERSION 0.4 LANGUAGES C CXX)
option(BUILD_SHARED_LIBS "Build as a shared library" ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT DEFINED LIB_INSTALL_DIR)
    set(LIB_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)
endif()


## Find Privateer's dependencies that have been compiled locally. 
#
find_library(MMDB2DEP mmdb2)
find_library(CCP4CDEP ccp4c)
find_library(CCP4SRSDEP ccp4srs)
find_library(CLIPPERCOREDEP clipper-core)
find_library(CLIPPERMMDBDEP clipper-mmdb)
find_library(CLIPPERMINIMOLDEP clipper-minimol)
find_library(CLIPPERCONTRIBDEP clipper-contrib)
find_library(CLIPPERCCP4DEP clipper-ccp4)
find_library(CLIPPERCIFDEP clipper-cif)

find_path(DEPENDENCYHEADERS ${CMAKE_SOURCE_DIR}/dependencies/include)

add_subdirectory(${CMAKE_SOURCE_DIR}/dependencies/pybind11)




