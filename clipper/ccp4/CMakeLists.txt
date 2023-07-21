project(clipper-ccp4)

add_library(clipper-ccp4 STATIC 
../../checkout/clipper/clipper/ccp4/ccp4_map_io.cpp
../../checkout/clipper/clipper/ccp4/ccp4_mtz_io.cpp
../../checkout/clipper/clipper/ccp4/ccp4_mtz_types.cpp
../../checkout/clipper/clipper/ccp4/ccp4_utils.cpp
)

target_include_directories(clipper-ccp4 PRIVATE ../../checkout/libccp4 ../../checkout/clipper ../../fftw ../../rfftw)
