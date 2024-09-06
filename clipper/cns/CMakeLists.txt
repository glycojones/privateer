project(clipper-cns)

add_library(clipper-cns STATIC 
${clipper-dependencies_SOURCE_DIR}/clipper/cns/cns_hkl_io.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/cns/cns_map_io.cpp
)

#target_include_directories(clipper-cns PRIVATE ../.. ${WRK_DIR}/checkout/clipper ${WRK_DIR}/fftw ${WRK_DIR}/rfftw)
target_include_directories(clipper-cns PRIVATE
        ${clipper-dependencies_SOURCE_DIR}/clipper
        ${fftw-dependencies_SOURCE_DIR}
)

set(clipper-cns_HEADERS
${clipper-dependencies_SOURCE_DIR}/clipper/cns/cns_map_io.h
${clipper-dependencies_SOURCE_DIR}/clipper/cns/cns_hkl_io.h
${clipper-dependencies_SOURCE_DIR}/clipper/clipper-cns.h
)

target_compile_options(clipper-cns PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-cns PROPERTIES PUBLIC_HEADER "${clipper-cns_HEADERS}")

install(TARGETS clipper-cns
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/clipper/cns
)
