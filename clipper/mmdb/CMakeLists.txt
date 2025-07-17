project(clipper-mmdb)

add_library(clipper-mmdb STATIC 
${clipper-dependencies_SOURCE_DIR}/clipper/mmdb/clipper_mmdb.cpp
)

#target_include_directories(clipper-mmdb PRIVATE ${WRK_DIR}/checkout/mmdb2/ ${WRK_DIR}/checkout/clipper ${WRK_DIR}/fftw ${WRK_DIR}/rfftw)

target_include_directories(clipper-mmdb PRIVATE
        ${mmdb2-dependencies_SOURCE_DIR}
        ${clipper-dependencies_SOURCE_DIR}/clipper
        ${fftw-dependencies_SOURCE_DIR}
)
set(clipper-mmdb_HEADERS
        ${clipper-dependencies_SOURCE_DIR}/clipper/mmdb/clipper_mmdb.h
)

target_compile_options(clipper-mmdb PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-mmdb PROPERTIES PUBLIC_HEADER "${clipper-mmdb_HEADERS}")

install(TARGETS clipper-mmdb
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/clipper/mmdb
)
