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

set(JSON_BuildTests OFF CACHE INTERNAL "")
add_subdirectory(${CMAKE_SOURCE_DIR}/dependencies/json)
add_subdirectory(${CMAKE_SOURCE_DIR}/dependencies/pybind11)
enable_testing()
add_test(NAME privateer_tests COMMAND ${PYTHON} -m pytest )

include_directories(${CMAKE_SOURCE_DIR}/dependencies/include
                    ${CMAKE_SOURCE_DIR}/include}
                    ${PYBIND11_INCLUDE_DIR}
                    ${PYTHON_INCLUDE_DIRS}
                    ${CMAKE_SOURCE_DIR}/dependencies/gemmi/include)
add_library(privateer_lib
            ${CMAKE_SOURCE_DIR}/src/privateer-blobs.cpp
            ${CMAKE_SOURCE_DIR}/src/clipper-glyco.cpp
            ${CMAKE_SOURCE_DIR}/src/clipper-glyco_data.cpp
            ${CMAKE_SOURCE_DIR}/src/privateer-restraints.cpp
            ${CMAKE_SOURCE_DIR}/src/privateer-lib.cpp)
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
                        nlohmann_json::nlohmann_json)
set_target_properties ( privateer_lib    PROPERTIES OUTPUT_NAME privateer )

add_executable(privateer_exec
               privateer.cpp)
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

pybind11_add_module(privateer_module ${CMAKE_SOURCE_DIR}/src/privateer-pybind11.cpp ${CMAKE_SOURCE_DIR}/src/privateer-restraints.cpp)
set_target_properties ( privateer_module PROPERTIES OUTPUT_NAME privateer_core )
target_compile_options( privateer_module PRIVATE "-Wno-shadow")
target_link_libraries ( privateer_module PRIVATE 
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

set     ( PYTHON_MODULE_PATH ${CMAKE_INSTALL_PREFIX}/lib/py3 )
install ( DIRECTORY "python/privateer" DESTINATION ${PYTHON_MODULE_PATH})
# install ( FILES include/clipper-glyco.h include/clipper-glyco_data.h include/privateer-restraints.h include/privateer-lib.h include/privateer-blobs.h
#             DESTINATION ${CMAKE_INSTALL_PREFIX}/include/privateer)
install ( TARGETS privateer_lib
            LIBRARY DESTINATION ${LIB_INSTALL_DIR})
install ( TARGETS privateer_exec DESTINATION ${CMAKE_SOURCE_DIR})
install ( TARGETS privateer_module DESTINATION ${PYTHON_MODULE_PATH}/privateer )
