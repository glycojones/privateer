# project(ccp4)

# add_library(ccp4 STATIC 
# ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_array.c
# ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_general.c
# ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_parser.c
# ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_program.c
# ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_unitcell.c
# ${WRK_DIR}/checkout/libccp4/ccp4/cmap_accessor.c
# ${WRK_DIR}/checkout/libccp4/ccp4/cmap_close.c
# ${WRK_DIR}/checkout/libccp4/ccp4/cmap_data.c
# ${WRK_DIR}/checkout/libccp4/ccp4/cmap_header.c
# ${WRK_DIR}/checkout/libccp4/ccp4/cmap_labels.c
# ${WRK_DIR}/checkout/libccp4/ccp4/cmap_open.c
# ${WRK_DIR}/checkout/libccp4/ccp4/cmap_skew.c
# ${WRK_DIR}/checkout/libccp4/ccp4/cmap_stats.c
# ${WRK_DIR}/checkout/libccp4/ccp4/cmap_symop.c
# ${WRK_DIR}/checkout/libccp4/ccp4/cmtzlib.c
# ${WRK_DIR}/checkout/libccp4/ccp4/csymlib.c
# ${WRK_DIR}/checkout/libccp4/ccp4/cvecmat.c
# ${WRK_DIR}/checkout/libccp4/ccp4/library_err.c
# ${WRK_DIR}/checkout/libccp4/ccp4/library_file.c
# ${WRK_DIR}/checkout/libccp4/ccp4/library_utils.c
# ${WRK_DIR}/checkout/libccp4/ccp4/pack_c.c
# )

# set(ccp4c_HEADERS
#  ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_file_err.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_program.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_unitcell.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/cmap_errno.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/cmap_stats.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/csymlib.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/library_file.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/w32mvs.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_fortran.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_spg.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_utils.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/cmap_header.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/cmaplib.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/cvecmat.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/mtzdata.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_array.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_general.h   
#  ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_vars.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/cmap_labels.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/cmaplib_f.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/overview.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_errno.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_parser.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_types.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/cmap_data.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/cmap_skew.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/cmtzlib.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/pack_c.h
#  ${WRK_DIR}/checkout/libccp4/ccp4/ccp4_sysdep.h
# )

# set_target_properties(ccp4 PROPERTIES PUBLIC_HEADER "${ccp4c_HEADERS}")

# install(TARGETS ccp4
# LIBRARY DESTINATION lib
# PUBLIC_HEADER DESTINATION include/ccp4
# )


cmake_minimum_required(VERSION 3.15)
project(ccp4)
option(BUILD_SHARED_LIBS "Build shared library" ON)
option(INSTALL_DEV_FILES "Install headers, cmake configs" ON)
file(GLOB ccp4_src "${ccp4-dependencies_SOURCE_DIR}/ccp4/*.c")
#message(STATUS ${mmdb_src})
file(GLOB ccp4_headers "${ccp4-dependencies_SOURCE_DIR}/ccp4/*.h")
list(REMOVE_ITEM ccp4_src "${ccp4-dependencies_SOURCE_DIR}/ccp4/vmslibrary.c")
add_library(ccp4 STATIC ${ccp4_src} ${ccp4_headers})
#set(ccp4c_HEADERS ${ccp4_headers})

set_target_properties(ccp4 PROPERTIES PUBLIC_HEADER "${ccp4c_HEADERS}")
if (INSTALL_DEV_FILES)
    install(TARGETS ccp4
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include/ccp4
    #LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    #PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include/ccp4
    )
endif()