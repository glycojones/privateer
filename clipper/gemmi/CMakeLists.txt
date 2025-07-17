project(clipper-gemmi)

add_library(clipper-gemmi STATIC 
${WRK_DIR}/checkout/clipper/clipper/gemmi/clipper_gemmi.cpp)

target_include_directories(clipper-gemmi PRIVATE
  ${WRK_DIR}/checkout/clipper/clipper
  ${WRK_DIR}/checkout/clipper/
  ${WRK_DIR}/checkout/gemmi/include)

set(clipper-gemmi_HEADERS
${WRK_DIR}/checkout/clipper/clipper/gemmi/clipper_gemmi.h
)

target_compile_options(clipper-gemmi PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-gemmi PROPERTIES PUBLIC_HEADER "${clipper-gemmi_HEADERS}")

install(TARGETS clipper-gemmi
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/clipper/gemmi
)
