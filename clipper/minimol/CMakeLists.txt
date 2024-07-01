project(clipper-minimol)

add_library(clipper-minimol STATIC 
../../checkout/clipper/clipper/minimol/container_minimol.cpp
../../checkout/clipper/clipper/minimol/minimol.cpp
../../checkout/clipper/clipper/minimol/minimol_data.cpp
../../checkout/clipper/clipper/minimol/minimol_io_gemmi.cpp
../../checkout/clipper/clipper/minimol/minimol_io_mmdb.cpp
../../checkout/clipper/clipper/minimol/minimol_seq.cpp
../../checkout/clipper/clipper/minimol/minimol_utils.cpp
)

target_include_directories(clipper-minimol PRIVATE
 ${CMAKE_SOURCE_DIR}/install/include/
 ${CMAKE_SOURCE_DIR}/install/include/gemmi
 ../../checkout/mmdb2/ .
 ./../checkout/clipper
  ../../checkout/fftw-2.1.5/fftw
  ../../checkout/fftw-2.1.5/rfftw
    ../../checkout/gemmi/include

  )

set(clipper-minimol_HEADERS
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/minimol/minimol_data.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/minimol/minimol.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/minimol/container_minimol.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/minimol/minimol_io_gemmi.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/minimol/minimol_io_mmdb.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/minimol/minimol_seq.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/minimol/minimol_utils.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/clipper-minimol.h
)

target_compile_options(clipper-minimol PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-minimol PROPERTIES PUBLIC_HEADER "${clipper-minimol_HEADERS}")

install(TARGETS clipper-minimol
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/clipper/minimol
)