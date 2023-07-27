project(clipper-gemmi)

add_library(clipper-gemmi STATIC 
../../checkout/clipper/clipper/gemmi/clipper_gemmi.cpp)

target_include_directories(clipper-gemmi PRIVATE
  ../../checkout/clipper/clipper
  ../../checkout/clipper/
  ../../checkout/gemmi/include)

set(clipper-gemmi_HEADERS
${CMAKE_SOURCE_DIR}/checkout/clipper/clipper/gemmi/clipper_gemmi.h
)

target_compile_options(clipper-gemmi PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-gemmi PROPERTIES PUBLIC_HEADER "${clipper-gemmi_HEADERS}")

install(TARGETS clipper-gemmi
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/clipper/gemmi
)
