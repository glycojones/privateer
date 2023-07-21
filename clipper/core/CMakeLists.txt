project(clipper-core)

add_library(clipper-core STATIC 
../../checkout/clipper/clipper/core/atomsf.cpp
../../checkout/clipper/clipper/core/cell.cpp
../../checkout/clipper/clipper/core/clipper_instance.cpp
../../checkout/clipper/clipper/core/clipper_memory.cpp
../../checkout/clipper/clipper/core/clipper_message.cpp
../../checkout/clipper/clipper/core/clipper_stats.cpp
../../checkout/clipper/clipper/core/clipper_test.cpp
../../checkout/clipper/clipper/core/clipper_thread.cpp
../../checkout/clipper/clipper/core/clipper_types.cpp
../../checkout/clipper/clipper/core/clipper_util.cpp
../../checkout/clipper/clipper/core/container.cpp
../../checkout/clipper/clipper/core/container_hkl.cpp
../../checkout/clipper/clipper/core/container_map.cpp
../../checkout/clipper/clipper/core/container_types.cpp
../../checkout/clipper/clipper/core/coords.cpp
../../checkout/clipper/clipper/core/derivs.cpp
../../checkout/clipper/clipper/core/fftmap.cpp
../../checkout/clipper/clipper/core/fftmap_sparse.cpp
../../checkout/clipper/clipper/core/hkl_compute.cpp
../../checkout/clipper/clipper/core/hkl_data.cpp
../../checkout/clipper/clipper/core/hkl_datatypes.cpp
../../checkout/clipper/clipper/core/hkl_info.cpp
../../checkout/clipper/clipper/core/hkl_lookup.cpp
../../checkout/clipper/clipper/core/hkl_operators.cpp
../../checkout/clipper/clipper/core/map_interp.cpp
../../checkout/clipper/clipper/core/map_utils.cpp
../../checkout/clipper/clipper/core/nxmap.cpp
../../checkout/clipper/clipper/core/nxmap_operator.cpp
../../checkout/clipper/clipper/core/ramachandran.cpp
../../checkout/clipper/clipper/core/resol_basisfn.cpp
../../checkout/clipper/clipper/core/resol_fn.cpp
../../checkout/clipper/clipper/core/resol_targetfn.cpp
../../checkout/clipper/clipper/core/rotation.cpp
../../checkout/clipper/clipper/core/spacegroup.cpp
../../checkout/clipper/clipper/core/spacegroup_data.cpp
../../checkout/clipper/clipper/core/symop.cpp
../../checkout/clipper/clipper/core/test_core.cpp
../../checkout/clipper/clipper/core/test_data.cpp
../../checkout/clipper/clipper/core/xmap.cpp
)

target_include_directories(clipper-core PRIVATE ../../checkout/fftw-2.1.5/fftw ../../checkout/fftw-2.1.5/rfftw)
target_compile_options(clipper-core PUBLIC "-DFFTW_ENABLE_FLOAT")

set(clipper-core_HEADERS
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/container_hkl.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/fftmap.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/hkl_compute.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/test_data.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/clipper_test.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/coords.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/clipper_sysdep.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/nxmap_operator.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/derivs.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/spacegroup_data.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/clipper_precision.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/container_map.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/container_types.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/hkl_lookup.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/container.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/clipper_message.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/cell.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/resol_basisfn.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/hkl_data.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/test_core.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/rotation.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/hkl_operators.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/fftmap_sparse.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/clipper_types.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/clipper_memory.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/symop.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/resol_targetfn.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/clipper_stats.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/clipper_thread.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/spacegroup.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/clipper_instance.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/ramachandran.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/hkl_info.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/hkl_datatypes.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/xmap.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/nxmap.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/map_interp.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/map_utils.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/resol_fn.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/atomsf.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/core/clipper_util.h
)

target_compile_options(clipper-core PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-core PROPERTIES PUBLIC_HEADER "${clipper-core_HEADERS}")

install(TARGETS clipper-core
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/clipper/core
)

set(clipper-headers_HEADERS
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/clipper-cctbx.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/cctbx/clipper_cctbx.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/clipper-mmdb.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/clipper-cif.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/clipper-mmdbold.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/clipper-ccp4.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/minimal-clipper-hkl.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/clipper-phs.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/minimal-clipper-map.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/clipper-minimol.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/clipper.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/clipper-cns.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/clipper-contrib.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/cif/cif_data_io.h
${CMAKE_SOURCE_DIR}/checkout/clipper//clipper/phs/phs_io.h
)

install(FILES 
${clipper-headers_HEADERS}
DESTINATION include/clipper
)
