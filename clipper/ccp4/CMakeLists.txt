project(clipper-ccp4)

add_library(clipper-ccp4 STATIC 
${clipper-dependencies_SOURCE_DIR}/clipper/ccp4/ccp4_map_io.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/ccp4/ccp4_mtz_io.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/ccp4/ccp4_mtz_types.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/ccp4/ccp4_utils.cpp
)

target_include_directories(clipper-ccp4 PRIVATE
        ${ccp4-dependencies_SOURCE_DIR}
        ${clipper-dependencies_SOURCE_DIR}/clipper
        ${fftw-dependencies_SOURCE_DIR}
)