cmake_minimum_required(VERSION 3.4)
project(ccp4srs CXX)

include(GNUInstallDirs)

option(BUILD_SHARED_LIBS "build shared libs (true)" TRUE)

include(CheckLibraryExists)

CHECK_LIBRARY_EXISTS(m cos "" HAVE_LIB_M)                                                                                                
                                                                                                                                         
if (HAVE_LIB_M)                                                                                                                          
    set(EXTRA_LIBS ${EXTRA_LIBS} m)                                                                                                      
endif (HAVE_LIB_M)

include( CheckIncludeFileCXX )
include( CheckFunctionExists )
include( CheckTypeSize)

check_include_file_cxx( "dirent.h" HAVE_DIRENT_H)

set(PACKAGE "ccp4srs")
set(VERSION "1.0")

set(ssm_src ${CMAKE_CURRENT_SOURCE_DIR}/../checkout/ccp4srs)

configure_file(${ssm_src}/config.h.cmakein
               ${CMAKE_CURRENT_BINARY_DIR}/config.h @ONLY)

set(prefix      ${CMAKE_INSTALL_PREFIX})
set(exec_prefix ${CMAKE_INSTALL_PREFIX})
set(libdir      ${CMAKE_INSTALL_PREFIX}/lib)
set(includedir  ${CMAKE_INSTALL_PREFIX}/include)
configure_file(${ssm_src}/ccp4srs.pc.cmakein
   ${CMAKE_CURRENT_BINARY_DIR}/ccp4srs.pc.in @ONLY)
 file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/ccp4srs.pc
   INPUT ${CMAKE_CURRENT_BINARY_DIR}/ccp4srs.pc.in)

set(ccp4srs_sources
	${ssm_src}/ccp4srs/ccp4srs_angle.cpp
        ${ssm_src}/ccp4srs/ccp4srs_manager.cpp
	${ssm_src}/ccp4srs/ccp4srs_atom.cpp
        ${ssm_src}/ccp4srs/ccp4srs_monomer.cpp
	${ssm_src}/ccp4srs/ccp4srs_base.cpp
        ${ssm_src}/ccp4srs/ccp4srs_plane.cpp
	${ssm_src}/ccp4srs/ccp4srs_bond.cpp
        ${ssm_src}/ccp4srs/ccp4srs_torsion.cpp
	${ssm_src}/ccp4srs/ccp4srs_chem.cpp
        ${ssm_src}/ccp4srs/ccp4srs_tree.cpp
	${ssm_src}/ccp4srs/ccp4srs_chicenter.cpp
        ${ssm_src}/ccp4srs/ccp4srs_types.cpp
	${ssm_src}/ccp4srs/ccp4srs_container.cpp
        ${ssm_src}/ccp4srs/memio_.cpp 
	${ssm_src}/ccp4srs/ccp4srs_entry.cpp
        ${ssm_src}/ccp4srs/mpfile_.cpp
	${ssm_src}/ccp4srs/ccp4srs_index.cpp
        ${ssm_src}/ccp4srs/ccp4srs_smiles.cpp
)

set(ccp4srs_headers
	${ssm_src}/ccp4srs/ccp4srs_angle.h
        ${ssm_src}/ccp4srs/ccp4srs_index.h
	${ssm_src}/ccp4srs/ccp4srs_atom.h
        ${ssm_src}/ccp4srs/ccp4srs_manager.h
	${ssm_src}/ccp4srs/ccp4srs_base.h
        ${ssm_src}/ccp4srs/ccp4srs_monomer.h
	${ssm_src}/ccp4srs/ccp4srs_bond.h
        ${ssm_src}/ccp4srs/ccp4srs_plane.h 
	${ssm_src}/ccp4srs/ccp4srs_chem.h
        ${ssm_src}/ccp4srs/ccp4srs_torsion.h
	${ssm_src}/ccp4srs/ccp4srs_chicenter.h
        ${ssm_src}/ccp4srs/ccp4srs_tree.h
	${ssm_src}/ccp4srs/ccp4srs_container.h
        ${ssm_src}/ccp4srs/ccp4srs_types.h
	${ssm_src}/ccp4srs/ccp4srs_defs.h
        ${ssm_src}/ccp4srs/memio_.h
	${ssm_src}/ccp4srs/ccp4srs_entry.h
        ${ssm_src}/ccp4srs/mpfile_.h
        ${ssm_src}/ccp4srs/ccp4srs_smiles.h
)

add_library(ccp4srs STATIC ${ccp4srs_sources} ${ccp4srs_headers})
add_library(ccp4srs::ccp4srs ALIAS ccp4srs)

target_link_directories(ccp4srs PUBLIC ../mmdb2 ../ccp4)
target_link_libraries(ccp4srs PRIVATE 
                             z
                             mmdb2 ccp4
                             ${EXTRA_LIBS}
)

target_include_directories(ccp4srs PRIVATE .
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/ccp4srs
    .. ../checkout/mmdb2 ../checkout/libccp/ccp4 ../checkout/zlib
    ${CMAKE_BINARY_DIR}
)

set_target_properties(ccp4srs PROPERTIES PUBLIC_HEADER  "${ccp4srs_headers}" )
set_target_properties(ccp4srs PROPERTIES SOVERSION ${VERSION})
if(MSVC)
    set_target_properties(ccp4srs PROPERTIES DEFINE_SYMBOL "CCP4SRS_DLL_EXPORT")
endif()

set(srsgen_sources
    ${ssm_src}/srsgen.cpp
)
if(MSVC)
    set(srsgen_sources ${srsgen_sources} ${ssm_src}/windows/dirent.cpp )
endif(MSVC)

add_executable(srsgen ${srsgen_sources})
target_link_directories(srsgen PUBLIC ../mmdb2 ../ccp4)
target_link_libraries(srsgen PRIVATE 
                             ccp4srs::ccp4srs
                             mmdb2 ccp4 nodefs.js
                             ${EXTRA_LIBS}
)
target_include_directories(srsgen PRIVATE .
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/ccp4srs
    .. ../checkout/mmdb2 ../checkout/libccp4/ccp4 ../checkout/zlib
    ${CMAKE_BINARY_DIR}
)
                                 
add_executable(srsview ${ssm_src}/srsview.cpp)
target_link_directories(srsview PUBLIC ../mmdb2 ../ccp4)
target_link_libraries(srsview PRIVATE 
                             ccp4srs::ccp4srs
                             mmdb2 ccp4 nodefs.js
                             ${EXTRA_LIBS}
)
target_include_directories(srsview PRIVATE .
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/ccp4srs
    .. ../checkout/mmdb2 ../checkout/libccp4/ccp4 ../checkout/zlib
    ${CMAKE_BINARY_DIR}
)

add_executable(pdb2to3 ${ssm_src}/srsrem.cpp)

target_include_directories(pdb2to3 PRIVATE .
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/ccp4srs
    .. ../checkout/mmdb2 ../checkout/libccp4/ccp4 ../checkout/zlib
    ${CMAKE_BINARY_DIR}
)

target_link_directories(pdb2to3 PUBLIC ../mmdb2 ../ccp4 ../checkout/zlib)
target_link_libraries(pdb2to3 PRIVATE 
                             ccp4srs::ccp4srs
                             mmdb2 ccp4 nodefs.js
                             
                             ${EXTRA_LIBS}
)
target_include_directories(srsgen PRIVATE .
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/ccp4srs
    .. ../checkout/mmdb2 ../checkout/libccp4/ccp4
    ${CMAKE_BINARY_DIR}
)

install(TARGETS ccp4srs
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/ccp4srs
)
