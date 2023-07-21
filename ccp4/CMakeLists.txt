project(ccp4)

add_library(ccp4 STATIC 
../checkout/libccp4/ccp4/ccp4_array.c
../checkout/libccp4/ccp4/ccp4_general.c
../checkout/libccp4/ccp4/ccp4_parser.c
../checkout/libccp4/ccp4/ccp4_program.c
../checkout/libccp4/ccp4/ccp4_unitcell.c
../checkout/libccp4/ccp4/cmap_accessor.c
../checkout/libccp4/ccp4/cmap_close.c
../checkout/libccp4/ccp4/cmap_data.c
../checkout/libccp4/ccp4/cmap_header.c
../checkout/libccp4/ccp4/cmap_labels.c
../checkout/libccp4/ccp4/cmap_open.c
../checkout/libccp4/ccp4/cmap_skew.c
../checkout/libccp4/ccp4/cmap_stats.c
../checkout/libccp4/ccp4/cmap_symop.c
../checkout/libccp4/ccp4/cmtzlib.c
../checkout/libccp4/ccp4/csymlib.c
../checkout/libccp4/ccp4/cvecmat.c
../checkout/libccp4/ccp4/library_err.c
../checkout/libccp4/ccp4/library_file.c
../checkout/libccp4/ccp4/library_utils.c
../checkout/libccp4/ccp4/pack_c.c
)

set(ccp4c_HEADERS
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/ccp4_file_err.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/ccp4_program.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/ccp4_unitcell.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/cmap_errno.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/cmap_stats.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/csymlib.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/library_file.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/w32mvs.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/ccp4_fortran.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/ccp4_spg.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/ccp4_utils.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/cmap_header.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/cmaplib.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/cvecmat.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/mtzdata.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/ccp4_array.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/ccp4_general.h   
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/ccp4_vars.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/cmap_labels.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/cmaplib_f.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/overview.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/ccp4_errno.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/ccp4_parser.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/ccp4_types.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/cmap_data.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/cmap_skew.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/cmtzlib.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/pack_c.h
 ${CMAKE_SOURCE_DIR}/checkout/libccp4/ccp4/ccp4_sysdep.h
)

set_target_properties(ccp4 PROPERTIES PUBLIC_HEADER "${ccp4c_HEADERS}")

install(TARGETS ccp4
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/ccp4
)
