## Project settings
#
cmake_minimum_required(VERSION 3.1)
project(privateer C CXX)
option(BUILD_SHARED_LIBS "Build as a shared library" ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT DEFINED LIB_INSTALL_DIR)
    set(LIB_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)
endif()

## Required CCP4 CCP4_LIBRARIES
#
find_package(CCP4 REQUIRED mmdb2 ccp4c ccp4srs clipper-core clipper-mmdb clipper-minimol
                           clipper-contrib clipper-ccp4 clipper-cif )

## Required Python libraries
#
if(APPLE)
    find_program(PYTHON NAMES $ENV{PYTHON} python3.6
                        HINTS ${CMAKE_INSTALL_PREFIX}/bin)
else()
    find_program(PYTHON NAMES $ENV{PYTHON} python3.6
                        HINTS ${CMAKE_INSTALL_PREFIX}/bin)
endif()

## Configure Python
#
execute_process(COMMAND ${PYTHON} -c "import sys; sys.stdout.write(sys.prefix)"
                OUTPUT_VARIABLE _PY_PREFIX
                RESULT_VARIABLE failure)
if (failure)
    message(FATAL_ERROR "Failed to run ${PYTHON}: ${failure}")
else()
    message(STATUS "Python interpreter: ${PYTHON}")
    message(STATUS "Python sys.prefix: ${_PY_PREFIX}")
endif()
# Get Python libraries
find_library(PYTHON_LIBRARY NAMES python3.6m
             HINTS ${_PY_PREFIX}/lib ${_PY_PREFIX}/lib64 ${_PY_PREFIX}/libs
             NO_DEFAULT_PATH)
if (NOT PYTHON_LIBRARY)
    message(FATAL_ERROR "Python library not found.")
else ()
    message(STATUS "Python library: ${PYTHON_LIBRARY}")
endif()
# Get Python includes
find_path(PYTHON_INCLUDE_DIR NAMES Python.h
          HINTS ${_PY_PREFIX}/include
                ${_PY_PREFIX}/include/python3.6m
          NO_DEFAULT_PATH)
if (NOT PYTHON_INCLUDE_DIR)
    message(FATAL_ERROR "Python header not found.")
else ()
    message(STATUS "Python.h is in: ${PYTHON_INCLUDE_DIR}")
    set(PYTHON_LIBRARIES ${PYTHON_LIBRARY})
endif()

## Get Pybind11 through cmake
#
find_package(pybind11 CONFIG REQUIRED)
MESSAGE( STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
MESSAGE( STATUS "pybind11_INCLUDE_DIRS = ${pybind11_INCLUDE_DIRS}" )
MESSAGE( STATUS "pybind11_LIBRARIES    = ${pybind11_LIBRARIES}" )

## Add all the include directories
#
include_directories(${CCP4_INCLUDE_DIRS}
                    ${CMAKE_CURRENT_SOURCE_DIR}
                    ${pybind11_INCLUDE_DIRS}
                    ${PYTHON_INCLUDE_DIR})

## Create libprivateer
#
add_library(privateer_lib
            clipper-glyco.cpp
            clipper-glyco_data.cpp
            privateer-restraints.cpp
            privateer-lib.cpp)
target_link_libraries ( privateer_lib ${CCP4_LIBRARIES} )
set_target_properties ( privateer_lib    PROPERTIES OUTPUT_NAME privateer )

## Create Privateer executable
#
add_executable(privateer_exec
               privateer.cpp)
target_link_libraries ( privateer_exec privateer_lib ${CCP4_LIBRARIES} )
set_target_properties ( privateer_exec   PROPERTIES OUTPUT_NAME privateer )

## Make hybrid Python module
#
pybind11_add_module(privateer_module privateer-pybind11.cpp)
set_target_properties ( privateer_module PROPERTIES OUTPUT_NAME privateer_core )
target_link_libraries ( privateer_module PRIVATE privateer_lib ${CCP4_LIBRARIES} ${pybind11_LIBRARIES} )

## And finally, install everything in their reet places
#
set     ( PYTHON_MODULE_PATH ${CMAKE_INSTALL_PREFIX}/lib/py3 )
message ( STATUS "Installing Python module to: ${PYTHON_MODULE_PATH}")
file    ( COPY "python/privateer" DESTINATION ${PYTHON_MODULE_PATH})
install ( FILES clipper-glyco.h clipper-glyco_data.h privateer-lib.h
          DESTINATION include/privateer)
install ( TARGETS privateer_lib
          LIBRARY DESTINATION "${LIB_INSTALL_DIR}")
install ( TARGETS privateer_exec DESTINATION bin)
install ( TARGETS privateer_module DESTINATION ${PYTHON_MODULE_PATH}/privateer )
