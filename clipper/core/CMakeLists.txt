project(clipper-core)

add_library(clipper-core STATIC 
${clipper-dependencies_SOURCE_DIR}/clipper/core/atomsf.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/cell.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_instance.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_memory.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_message.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_stats.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_test.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_thread.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_types.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_util.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/container.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/container_hkl.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/container_map.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/container_types.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/coords.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/derivs.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/fftmap.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/fftmap_sparse.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/hkl_compute.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/hkl_data.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/hkl_datatypes.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/hkl_info.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/hkl_lookup.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/hkl_operators.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/map_interp.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/map_utils.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/nxmap.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/nxmap_operator.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/ramachandran.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/resol_basisfn.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/resol_fn.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/resol_targetfn.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/rotation.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/spacegroup.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/spacegroup_data.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/symop.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/test_core.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/test_data.cpp
${clipper-dependencies_SOURCE_DIR}/clipper/core/xmap.cpp
)

target_include_directories(clipper-core PRIVATE
        ${fftw-dependencies_SOURCE_DIR}/fftw
        ${fftw-dependencies_SOURCE_DIR}/rfftw)
target_compile_options(clipper-core PUBLIC "-DFFTW_ENABLE_FLOAT")

set(clipper-core_HEADERS
${clipper-dependencies_SOURCE_DIR}/clipper/core/container_hkl.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/fftmap.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/hkl_compute.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/test_data.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_test.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/coords.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_sysdep.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/nxmap_operator.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/derivs.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/spacegroup_data.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_precision.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/container_map.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/container_types.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/hkl_lookup.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/container.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_message.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/cell.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/resol_basisfn.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/hkl_data.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/test_core.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/rotation.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/hkl_operators.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/fftmap_sparse.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_types.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_memory.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/symop.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/resol_targetfn.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_stats.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_thread.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/spacegroup.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_instance.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/ramachandran.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/hkl_info.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/hkl_datatypes.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/xmap.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/nxmap.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/map_interp.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/map_utils.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/resol_fn.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/atomsf.h
${clipper-dependencies_SOURCE_DIR}/clipper/core/clipper_util.h
)

target_compile_options(clipper-core PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-core PROPERTIES PUBLIC_HEADER "${clipper-core_HEADERS}")

install(TARGETS clipper-core
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/clipper/core
)

set(clipper-headers_HEADERS
${clipper-dependencies_SOURCE_DIR}//clipper/clipper-cctbx.h
${clipper-dependencies_SOURCE_DIR}//clipper/cctbx/clipper_cctbx.h
${clipper-dependencies_SOURCE_DIR}//clipper/clipper-mmdb.h
${clipper-dependencies_SOURCE_DIR}//clipper/clipper-cif.h
${clipper-dependencies_SOURCE_DIR}//clipper/clipper-mmdbold.h
${clipper-dependencies_SOURCE_DIR}//clipper/clipper-ccp4.h
${clipper-dependencies_SOURCE_DIR}//clipper/minimal-clipper-hkl.h
${clipper-dependencies_SOURCE_DIR}//clipper/clipper-phs.h
${clipper-dependencies_SOURCE_DIR}//clipper/minimal-clipper-map.h
${clipper-dependencies_SOURCE_DIR}//clipper/clipper-minimol.h
${clipper-dependencies_SOURCE_DIR}//clipper/clipper.h
${clipper-dependencies_SOURCE_DIR}//clipper/clipper-cns.h
${clipper-dependencies_SOURCE_DIR}//clipper/clipper-contrib.h
${clipper-dependencies_SOURCE_DIR}//clipper/cif/cif_data_io.h
${clipper-dependencies_SOURCE_DIR}//clipper/phs/phs_io.h
)

install(FILES 
${clipper-headers_HEADERS}
DESTINATION include/clipper
)
