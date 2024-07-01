project(clipper-gemmi)

add_library(clipper-gemmi STATIC
../../checkout/clipper/clipper/gemmi/clipper_gemmi.cpp
../../checkout/clipper/clipper/gemmi/clipper_gemmi_model.cpp
)
#../../checkout/clipper/clipper/gemmi/clipper_gemmi_map.cpp


target_include_directories(clipper-gemmi PRIVATE
  ../../checkout/clipper/clipper
  ../../checkout/clipper/
  ../../checkout/gemmi/include)

target_link_directories(clipper-gemmi PRIVATE
  ../../gemmi
 )

 target_link_libraries(clipper-gemmi PRIVATE gemmi_cpp)

set(clipper-gemmi_HEADERS
${CMAKE_SOURCE_DIR}/checkout/clipper/clipper/gemmi/clipper_gemmi.h
# ${CMAKE_SOURCE_DIR}/checkout/clipper/clipper/gemmi/clipper_gemmi_map.h
${CMAKE_SOURCE_DIR}/checkout/clipper/clipper/gemmi/clipper_gemmi_model.h
)

target_compile_options(clipper-gemmi PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-gemmi PROPERTIES PUBLIC_HEADER "${clipper-gemmi_HEADERS}")

install(TARGETS clipper-gemmi
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/clipper/gemmi
)