## Project settings
#
cmake_minimum_required(VERSION 3.12)
project(privateer VERSION 0.5 LANGUAGES C CXX)

SET(PRIVATEER_CC_FLAGS "-g -O2 -w")
SET(PRIVATEER_CXX_FLAGS "-g -O2 -w")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PRIVATEER_CC_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PRIVATEER_CXX_FLAGS}")

option(BUILD_SHARED_LIBS "Build as a shared library" ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)

set(PRIVATEER_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/privateer)
set(PRIVATEERDATA ${CMAKE_SOURCE_DIR}/data)

set(CLIBDENV $ENV{CLIBD})

## Find Privateer's dependencies that have been compiled locally. 
if (MODE STREQUAL "TESTING")
    set (MMDB2DEP ${CMAKE_SOURCE_DIR}/privateer_dependencies/lib/libmmdb2.so)
    set (CCP4CDEP ${CMAKE_SOURCE_DIR}/privateer_dependencies/lib/libccp4srs.so)
    set (CCP4SRSDEP ${CMAKE_SOURCE_DIR}/privateer_dependencies/lib/libccp4c.so)
    set (CLIPPERCOREDEP ${CMAKE_SOURCE_DIR}/privateer_dependencies/lib/libclipper-core.so)
    set (CLIPPERMMDBDEP ${CMAKE_SOURCE_DIR}/privateer_dependencies/lib/libclipper-mmdb.so)
    set (CLIPPERMINIMOLDEP ${CMAKE_SOURCE_DIR}/privateer_dependencies/lib/libclipper-minimol.so)
    set (CLIPPERCONTRIBDEP ${CMAKE_SOURCE_DIR}/privateer_dependencies/lib/libclipper-contrib.so)
    set (CLIPPERCCP4DEP ${CMAKE_SOURCE_DIR}/privateer_dependencies/lib/libclipper-ccp4.so)
    set (CLIPPERCIFDEP ${CMAKE_SOURCE_DIR}/privateer_dependencies/lib/libclipper-cif.so)
else()
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
endif()

if (MODE STREQUAL "TESTING")
    add_subdirectory(${CMAKE_SOURCE_DIR}/privateer_dependencies/pybind11)
else()
    add_subdirectory(${CMAKE_SOURCE_DIR}/dependencies/pybind11)
endif()
find_package(Threads REQUIRED)

#Find Privateer's main.cpp file
set(PRIVATEER_SOURCE ${PRIVATEER_SOURCE_DIR}/cpp/privateer.cpp)

#Location of include files to be linked
if (MODE STREQUAL "TESTING")
    include_directories(${CMAKE_SOURCE_DIR}/privateer_dependencies/include
                    ${CMAKE_SOURCE_DIR}/privateer_dependencies/gemmi/include
                    ${PYBIND11_INCLUDE_DIR}
                    ${PYTHON_INCLUDE_DIRS}
                    ${PRIVATEER_SOURCE_DIR}/cpp
                    ${PRIVATEER_SOURCE_DIR}/cpp/pybind11)
else()
    include_directories(${CMAKE_SOURCE_DIR}/dependencies/include
                    ${CMAKE_SOURCE_DIR}/dependencies/gemmi/include
                    ${PYBIND11_INCLUDE_DIR}
                    ${PYTHON_INCLUDE_DIRS}
                    ${PRIVATEER_SOURCE_DIR}/cpp
                    ${PRIVATEER_SOURCE_DIR}/cpp/pybind11)
endif()

file(COPY ${PRIVATEERDATA}/glycomics/privateer_glycomics_database.json 
    DESTINATION ${CLIBDENV})

file(COPY ${PRIVATEERDATA}/linkage_torsions/privateer_torsion_database.json 
    DESTINATION ${CLIBDENV})

add_library(privateer_lib
            ${PRIVATEER_SOURCE_DIR}/cpp/third-party/gemmi_compilation_unit.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-error.h
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-json.h
            ${PRIVATEER_SOURCE_DIR}/cpp/clipper-glyco.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-db.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/clipper-glyco_data.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-blobs.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-composition.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-dbquery.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-cryo_em.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-xray.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-modelling.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-interactions.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-lib.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pybind11.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-ccp4mg.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-restraints.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pyanalysis.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pymodelling.cpp)

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
                        Threads::Threads)

set_target_properties ( privateer_exec   PROPERTIES OUTPUT_NAME privateer )

pybind11_add_module(privateer_core          ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pybind11.cpp
                                            ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-ccp4mg.cpp
                                            ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-restraints.cpp
                                            ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pyanalysis.cpp
                                            ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pymodelling.cpp)

pybind11_add_module(privateer_modelling     ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pybind11.cpp
                                            ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-ccp4mg.cpp
                                            ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-restraints.cpp
                                            ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pyanalysis.cpp
                                            ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pymodelling.cpp)

target_link_libraries ( privateer_core PRIVATE 
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
                        ${PYTHON_LIBRARY})

target_link_libraries ( privateer_modelling PRIVATE
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
                        ${PYTHON_LIBRARY})


install ( TARGETS privateer_lib LIBRARY DESTINATION ${PROJECT_SOURCE_DIR})
install ( TARGETS privateer_exec DESTINATION ${PROJECT_SOURCE_DIR})
install ( TARGETS privateer_core DESTINATION ${PROJECT_SOURCE_DIR} )
install ( TARGETS privateer_modelling DESTINATION ${PROJECT_SOURCE_DIR} )