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

set(PRIVATEER_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/privateer)



## Find Privateer's dependencies that have been compiled locally. 
#
find_library(MMDB2DEP NAMES mmdb2
            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib)
find_library(CCP4CDEP NAMES ccp4c
            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib)
find_library(CCP4SRSDEP NAMES ccp4srs
            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib)
find_library(CLIPPERCOREDEP NAMES clipper-core
            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib)
find_library(CLIPPERMMDBDEP NAMES clipper-mmdb
            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib)
find_library(CLIPPERMINIMOLDEP NAMES clipper-minimol
            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib)
find_library(CLIPPERCONTRIBDEP NAMES clipper-contrib
            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib)
find_library(CLIPPERCCP4DEP NAMES clipper-ccp4
            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib)
find_library(CLIPPERCIFDEP NAMES clipper-cif
            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib)

#Add nlohmann::json and gemmi that need to be integrated into privateer's shared library
set(JSON_BuildTests OFF CACHE INTERNAL "")
add_subdirectory(${CMAKE_SOURCE_DIR}/dependencies/json)
add_subdirectory(${CMAKE_SOURCE_DIR}/dependencies/pybind11)

find_package(Threads REQUIRED)


#Find Privateer's main.cpp file
set(PRIVATEER_SOURCE ${PRIVATEER_SOURCE_DIR}/cpp/privateer.cpp)

#Location of include files to be linked
include_directories(${CMAKE_SOURCE_DIR}/dependencies/include
                    ${CMAKE_SOURCE_DIR}/dependencies/gemmi/include
                    ${PYBIND11_INCLUDE_DIR}
                    ${PYTHON_INCLUDE_DIRS}
                    ${PRIVATEER_SOURCE_DIR}/cpp
                    ${PRIVATEER_SOURCE_DIR}/cpp/pybind11)

add_library(privateer_lib
            ${PRIVATEER_SOURCE_DIR}/cpp/clipper-glyco.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/clipper-glyco_data.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-blobs.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-composition.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-dbquery.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-cryo_em.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-parallelism.h
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-xray.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pybind11.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-ccp4mg.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-restraints.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-lib.cpp)

target_link_libraries ( privateer_lib 
                        ${MMDB2DEP} 
                        ${CCP4CDEP} 
                        ${CCP4SRSDEP}
                        ${CLIPPERCOREDEP}
                        ${CLIPPERMMDBDEP}
                        ${CLIPPERMINIMOLDEP}
                        ${CLIPPERCONTRIBDEP}
                        ${CLIPPERCCP4DEP}
                        ${CLIPPERCIFDEP}
                        ${PYTHON_LIBRARY}
                        nlohmann_json::nlohmann_json
                        Threads::Threads)

set_target_properties ( privateer_lib    PROPERTIES OUTPUT_NAME privateer )

add_executable(privateer_exec
               ${PRIVATEER_SOURCE})

target_link_libraries ( privateer_exec 
                        privateer_lib 
                        ${MMDB2DEP} 
                        ${CCP4CDEP} 
                        ${CCP4SRSDEP}
                        ${CLIPPERCOREDEP}
                        ${CLIPPERMMDBDEP}
                        ${CLIPPERMINIMOLDEP}
                        ${CLIPPERCONTRIBDEP}
                        ${CLIPPERCCP4DEP}
                        ${CLIPPERCIFDEP}
                        ${PYTHON_LIBRARY}
                        nlohmann_json::nlohmann_json)

set_target_properties ( privateer_exec   PROPERTIES OUTPUT_NAME privateer )

pybind11_add_module(libprivateer     ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pybind11.cpp
                                    ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-ccp4mg.cpp
                                    ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-restraints.cpp)
target_link_libraries ( libprivateer PRIVATE 
                        privateer_lib 
                        ${MMDB2DEP} 
                        ${CCP4CDEP} 
                        ${CCP4SRSDEP}
                        ${CLIPPERCOREDEP}
                        ${CLIPPERMMDBDEP}
                        ${CLIPPERMINIMOLDEP}
                        ${CLIPPERCONTRIBDEP}
                        ${CLIPPERCCP4DEP}
                        ${CLIPPERCIFDEP}
                        ${PYTHON_LIBRARY}
                        nlohmann_json::nlohmann_json)


install ( TARGETS privateer_lib LIBRARY DESTINATION ${PROJECT_SOURCE_DIR})
install ( TARGETS privateer_exec DESTINATION ${PROJECT_SOURCE_DIR})
install ( TARGETS libprivateer DESTINATION ${PROJECT_SOURCE_DIR} )