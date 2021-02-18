## Project settings
#
cmake_minimum_required(VERSION 3.12)
project(privateer VERSION 0.4 LANGUAGES C CXX)

SET(PRIVATEER_CC_FLAGS "-g -O2 -w")
SET(PRIVATEER_CXX_FLAGS "-g -O2 -w")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PRIVATEER_CC_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PRIVATEER_CXX_FLAGS}")

option(BUILD_SHARED_LIBS "Build as a shared library" ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)


set(PRIVATEER_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

#Add nlohmann::json and gemmi that need to be integrated into privateer's shared library
# set(JSON_BuildTests OFF CACHE INTERNAL "")
# add_subdirectory(${CMAKE_SOURCE_DIR}/dependencies/json)

if(NOT DEFINED LIB_INSTALL_DIR)
    set(LIB_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)
endif()

## Required CCP4 CCP4_LIBRARIES
#
find_package(CCP4 REQUIRED mmdb2 ccp4c ccp4srs clipper-core clipper-mmdb clipper-minimol
                           clipper-contrib clipper-ccp4 clipper-cif)

# Required for std::thread in C++11 standard. 
find_package(Threads REQUIRED)                          

## Required Python stuff - link against system Python3 for testing on YSBL machines
#
    # Link against CCP4-distributed Python
if(APPLE)
    find_program(PYTHON NAMES $ENV{PYTHON} ccp4-python
                        HINTS ${CMAKE_INSTALL_PREFIX}/bin )
else()
    find_program(PYTHON NAMES $ENV{PYTHON} ccp4-python python.exe
                        HINTS ${CMAKE_INSTALL_PREFIX}/bin
                        PATHS [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.7\\InstallPath])
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

set(PYTHON_EXECUTABLE ${PYTHON})

# Get Python libraries
find_library(PYTHON_LIBRARY NAMES $ENV{PYTHON} python2.7m python2.7 python
             HINTS ${_PY_PREFIX}/lib64 ${_PY_PREFIX}/lib ${_PY_PREFIX}/libs
             NO_DEFAULT_PATH)
if (NOT PYTHON_LIBRARY)
    message(FATAL_ERROR "Python library not found.")
else ()
    message(STATUS "Python library: ${PYTHON_LIBRARY}")
    set(PYTHON_LIBRARIES ${PYTHON_LIBRARY})
endif()

# Get Python includes
find_path(PYTHON_INCLUDE_DIR NAMES Python.h
          HINTS ${_PY_PREFIX}/include
                ${_PY_PREFIX}/include/python2.7
                ${_PY_PREFIX}/include/python2.7m
          NO_DEFAULT_PATH)


if (NOT PYTHON_INCLUDE_DIR)
    message(FATAL_ERROR "Python header not found.")
else ()
    message(STATUS "Python.h is in: ${PYTHON_INCLUDE_DIR}")
endif()

## Run pytest through 'make test'
#
enable_testing()
add_test(NAME privateer_tests COMMAND ${PYTHON} -m pytest )

## Get pybind11 through cmake
#
set(PYBIND11_PYTHON_VERSION 2.7 CACHE STRING "")
set(LIB_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)
find_package(pybind11 CONFIG REQUIRED)
MESSAGE( STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
MESSAGE( STATUS "pybind11_INCLUDE_DIRS = ${pybind11_INCLUDE_DIRS}" )
MESSAGE( STATUS "pybind11_LIBRARIES    = ${pybind11_LIBRARIES}" )

set(PRIVATEER_SOURCE ${PRIVATEER_SOURCE_DIR}/privateer.cpp)

## Add all the include directories
#
include_directories(${CCP4_INCLUDE_DIRS}
                    ${pybind11_INCLUDE_DIRS}
                    ${PYTHON_INCLUDE_DIR}
                    ${PRIVATEER_SOURCE_DIR})

## Create libprivateer
#
add_library(privateer_lib
            ${PRIVATEER_SOURCE_DIR}/clipper-glyco.cpp
            ${PRIVATEER_SOURCE_DIR}/clipper-glyco_data.cpp
            ${PRIVATEER_SOURCE_DIR}/privateer-blobs.cpp
            ${PRIVATEER_SOURCE_DIR}/privateer-composition.cpp
            ${PRIVATEER_SOURCE_DIR}/privateer-dbquery.cpp
            ${PRIVATEER_SOURCE_DIR}/privateer-cryo_em.cpp
            ${PRIVATEER_SOURCE_DIR}/privateer-parallelism.h
            ${PRIVATEER_SOURCE_DIR}/privateer-xray.cpp
            ${PRIVATEER_SOURCE_DIR}/privateer-restraints.cpp
            ${PRIVATEER_SOURCE_DIR}/privateer-lib.cpp)
  
target_link_libraries ( privateer_lib 
                        ${CCP4_LIBRARIES} 
                        ${pybind11_LIBRARIES} 
                        ${PYTHON_LIBRARY})
                        
set_target_properties ( privateer_lib    PROPERTIES OUTPUT_NAME privateer )

## Create Privateer executable
#
add_executable(privateer_exec
              ${PRIVATEER_SOURCE})

target_link_libraries ( privateer_exec 
                        privateer_lib 
                        ${CCP4_LIBRARIES} 
                        ${pybind11_LIBRARIES}
                        Threads::Threads)
set_target_properties ( privateer_exec   PROPERTIES OUTPUT_NAME privateer )

## Make hybrid Python module
#
pybind11_add_module(privateer_core ${PRIVATEER_SOURCE} ${PRIVATEER_SOURCE_DIR}/privateer-pybind11.cpp ${PRIVATEER_SOURCE_DIR}/privateer-restraints.cpp)
target_link_libraries ( privateer_core PRIVATE 
                        privateer_lib
                        ${CCP4_LIBRARIES} 
                        ${pybind11_LIBRARIES})


## And finally, install everything in their reet places
#


set     ( PYTHON_MODULE_PATH ${CMAKE_INSTALL_PREFIX}/lib/py2 )
message ( STATUS "Installing Python module to: ${PYTHON_MODULE_PATH}")
install ( DIRECTORY ${PRIVATEER_SOURCE_DIR}/privateer DESTINATION ${PYTHON_MODULE_PATH})
install ( FILES ${PRIVATEER_SOURCE_DIR}/clipper-glyco.h
                ${PRIVATEER_SOURCE_DIR}/clipper-glyco_data.h
                ${PRIVATEER_SOURCE_DIR}/privateer-blobs.h
                ${PRIVATEER_SOURCE_DIR}/privateer-composition.h
                ${PRIVATEER_SOURCE_DIR}/privateer-dbquery.h
                ${PRIVATEER_SOURCE_DIR}/privateer-cryo_em.h
                ${PRIVATEER_SOURCE_DIR}/privateer-parallelism.h
                ${PRIVATEER_SOURCE_DIR}/privateer-xray.h
                ${PRIVATEER_SOURCE_DIR}/privateer-restraints.h 
                ${PRIVATEER_SOURCE_DIR}/privateer-lib.h
                ${PRIVATEER_SOURCE_DIR}/database.json
          DESTINATION ${CMAKE_INSTALL_PREFIX}/include/privateer)
install ( TARGETS privateer_lib
          LIBRARY DESTINATION ${LIB_INSTALL_DIR})
install ( TARGETS privateer_exec DESTINATION bin)
install ( TARGETS privateer_core DESTINATION ${PYTHON_MODULE_PATH} )

