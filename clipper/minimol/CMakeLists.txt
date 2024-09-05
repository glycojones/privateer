project(clipper-minimol)

add_library(clipper-minimol STATIC
        ${clipper-dependencies_SOURCE_DIR}/clipper/minimol/container_minimol.cpp
        ${clipper-dependencies_SOURCE_DIR}/clipper/minimol/minimol.cpp
        ${clipper-dependencies_SOURCE_DIR}/clipper/minimol/minimol_data.cpp
        ${clipper-dependencies_SOURCE_DIR}/clipper/minimol/minimol_io.cpp
        ${clipper-dependencies_SOURCE_DIR}/clipper/minimol/minimol_seq.cpp
        ${clipper-dependencies_SOURCE_DIR}/clipper/minimol/minimol_utils.cpp
)

#target_include_directories(clipper-minimol PRIVATE ../../checkout/mmdb2/ ../../checkout/clipper ../../checkout/fftw-2.1.5/fftw ../../checkout/fftw-2.1.5/rfftw)

target_include_directories(clipper-minimol PRIVATE
        ${mmdb2-dependencies_SOURCE_DIR}
        ${clipper-dependencies_SOURCE_DIR}/clipper
        ${fftw-dependencies_SOURCE_DIR}
)

set(clipper-minimol_HEADERS
        ${clipper-dependencies_SOURCE_DIR}/clipper/minimol/minimol_data.h
        ${clipper-dependencies_SOURCE_DIR}/clipper/minimol/minimol.h
        ${clipper-dependencies_SOURCE_DIR}/clipper/minimol/container_minimol.h
        ${clipper-dependencies_SOURCE_DIR}/clipper/minimol/minimol_io.h
        ${clipper-dependencies_SOURCE_DIR}/clipper/minimol/minimol_seq.h
        ${clipper-dependencies_SOURCE_DIR}/clipper/minimol/minimol_utils.h
        ${clipper-dependencies_SOURCE_DIR}/clipper/clipper-minimol.h
)

target_compile_options(clipper-minimol PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-minimol PROPERTIES PUBLIC_HEADER "${clipper-minimol_HEADERS}")

install(TARGETS clipper-minimol
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/clipper/minimol
)