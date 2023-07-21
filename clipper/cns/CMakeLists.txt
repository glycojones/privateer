project(clipper-cns)

add_library(clipper-cns STATIC 
../../checkout/clipper/clipper/cns/cns_hkl_io.cpp
../../checkout/clipper/clipper/cns/cns_map_io.cpp
)

target_include_directories(clipper-cns PRIVATE ../.. ../../checkout/clipper ../../fftw ../../rfftw)

set(clipper-cns_HEADERS
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/cns/cns_map_io.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/cns/cns_hkl_io.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/clipper-cns.h
)

target_compile_options(clipper-cns PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-cns PROPERTIES PUBLIC_HEADER "${clipper-cns_HEADERS}")

install(TARGETS clipper-cns
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/clipper/cns
)
