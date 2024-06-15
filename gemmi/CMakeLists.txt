cmake_minimum_required(VERSION 3.13)

project(gemmi LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
message(STATUS "Compiling with C++ standard: ${CMAKE_CXX_STANDARD}")
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(gemmi_src ../checkout/gemmi)

find_package(ZLIB)
if (ZLIB_FOUND)
    include_directories("${ZLIB_INCLUDE_DIR}")
    message(STATUS "The build will use zlib code from emscripten")
endif()

include_directories("${gemmi_src}/include"
        "${gemmi_src}/third_party")

add_library(gemmi_cpp STATIC
        ${gemmi_src}/src/assembly.cpp
        ${gemmi_src}/src/calculate.cpp
        ${gemmi_src}/src/crd.cpp
        ${gemmi_src}/src/eig3.cpp
        ${gemmi_src}/src/mmcif.cpp
        ${gemmi_src}/src/mmread_gz.cpp
        ${gemmi_src}/src/mtz.cpp
        ${gemmi_src}/src/mtz2cif.cpp
        ${gemmi_src}/src/polyheur.cpp
        ${gemmi_src}/src/read_cif.cpp
        ${gemmi_src}/src/resinfo.cpp
        ${gemmi_src}/src/riding_h.cpp
        ${gemmi_src}/src/sprintf.cpp
        ${gemmi_src}/src/to_mmcif.cpp
        ${gemmi_src}/src/to_pdb.cpp
        ${gemmi_src}/src/monlib.cpp
        ${gemmi_src}/src/topo.cpp
        ${gemmi_src}/src/xds_ascii.cpp)

#set_target_properties(gemmi_cpp PROPERTIES COMPILE_FLAGS "-fexceptions -pthread -s USE_ZLIB=1 -s USE_PTHREADS=1" )

set(gemmi_HEADERS
        ${gemmi_src}/include/gemmi/addends.hpp
        ${gemmi_src}/include/gemmi/align.hpp
        ${gemmi_src}/include/gemmi/assembly.hpp
        ${gemmi_src}/include/gemmi/asudata.hpp
        ${gemmi_src}/include/gemmi/asumask.hpp
        ${gemmi_src}/include/gemmi/atof.hpp
        ${gemmi_src}/include/gemmi/atox.hpp
        ${gemmi_src}/include/gemmi/bessel.hpp
        ${gemmi_src}/include/gemmi/binner.hpp
        ${gemmi_src}/include/gemmi/blob.hpp
        ${gemmi_src}/include/gemmi/bond_idx.hpp
        ${gemmi_src}/include/gemmi/c4322.hpp
        ${gemmi_src}/include/gemmi/calculate.hpp
        ${gemmi_src}/include/gemmi/ccp4.hpp
        ${gemmi_src}/include/gemmi/cellred.hpp
        ${gemmi_src}/include/gemmi/chemcomp.hpp
        ${gemmi_src}/include/gemmi/chemcomp_xyz.hpp
        ${gemmi_src}/include/gemmi/cif.hpp
        ${gemmi_src}/include/gemmi/cif2mtz.hpp
        ${gemmi_src}/include/gemmi/cifdoc.hpp
        ${gemmi_src}/include/gemmi/contact.hpp
        ${gemmi_src}/include/gemmi/crd.hpp
        ${gemmi_src}/include/gemmi/ddl.hpp
        ${gemmi_src}/include/gemmi/dencalc.hpp
        ${gemmi_src}/include/gemmi/dirwalk.hpp
        ${gemmi_src}/include/gemmi/ecalc.hpp
        ${gemmi_src}/include/gemmi/eig3.hpp
        ${gemmi_src}/include/gemmi/elem.hpp
        ${gemmi_src}/include/gemmi/enumstr.hpp
        ${gemmi_src}/include/gemmi/fail.hpp
        ${gemmi_src}/include/gemmi/fileutil.hpp
        ${gemmi_src}/include/gemmi/floodfill.hpp
        ${gemmi_src}/include/gemmi/formfact.hpp
        ${gemmi_src}/include/gemmi/fourier.hpp
        ${gemmi_src}/include/gemmi/fprime.hpp
        ${gemmi_src}/include/gemmi/fstream.hpp
        ${gemmi_src}/include/gemmi/grid.hpp
        ${gemmi_src}/include/gemmi/gz.hpp
        ${gemmi_src}/include/gemmi/input.hpp
        ${gemmi_src}/include/gemmi/interop.hpp
        ${gemmi_src}/include/gemmi/it92.hpp
        ${gemmi_src}/include/gemmi/iterator.hpp
        ${gemmi_src}/include/gemmi/json.hpp
        ${gemmi_src}/include/gemmi/levmar.hpp
        ${gemmi_src}/include/gemmi/linkhunt.hpp
        ${gemmi_src}/include/gemmi/math.hpp
        ${gemmi_src}/include/gemmi/merge.hpp
        ${gemmi_src}/include/gemmi/metadata.hpp
        ${gemmi_src}/include/gemmi/mmcif.hpp
        ${gemmi_src}/include/gemmi/mmcif_impl.hpp
        ${gemmi_src}/include/gemmi/mmdb.hpp
        ${gemmi_src}/include/gemmi/mmread.hpp
        ${gemmi_src}/include/gemmi/mmread_gz.hpp
        ${gemmi_src}/include/gemmi/model.hpp
        ${gemmi_src}/include/gemmi/modify.hpp
        ${gemmi_src}/include/gemmi/monlib.hpp
        ${gemmi_src}/include/gemmi/mtz.hpp
        ${gemmi_src}/include/gemmi/mtz2cif.hpp
        ${gemmi_src}/include/gemmi/neighbor.hpp
        ${gemmi_src}/include/gemmi/neutron92.hpp
        ${gemmi_src}/include/gemmi/numb.hpp
        ${gemmi_src}/include/gemmi/pdb.hpp
        ${gemmi_src}/include/gemmi/pdb_id.hpp
        ${gemmi_src}/include/gemmi/pirfasta.hpp
        ${gemmi_src}/include/gemmi/polyheur.hpp
        ${gemmi_src}/include/gemmi/qcp.hpp
        ${gemmi_src}/include/gemmi/read_cif.hpp
        ${gemmi_src}/include/gemmi/read_map.hpp
        ${gemmi_src}/include/gemmi/recgrid.hpp
        ${gemmi_src}/include/gemmi/reciproc.hpp
        ${gemmi_src}/include/gemmi/refln.hpp
        ${gemmi_src}/include/gemmi/remarks.hpp
        ${gemmi_src}/include/gemmi/resinfo.hpp
        ${gemmi_src}/include/gemmi/riding_h.hpp
        ${gemmi_src}/include/gemmi/scaling.hpp
        ${gemmi_src}/include/gemmi/select.hpp
        ${gemmi_src}/include/gemmi/seqalign.hpp
        ${gemmi_src}/include/gemmi/seqid.hpp
        ${gemmi_src}/include/gemmi/seqtools.hpp
        ${gemmi_src}/include/gemmi/sfcalc.hpp
        ${gemmi_src}/include/gemmi/small.hpp
        ${gemmi_src}/include/gemmi/smcif.hpp
        ${gemmi_src}/include/gemmi/solmask.hpp
        ${gemmi_src}/include/gemmi/span.hpp
        ${gemmi_src}/include/gemmi/sprintf.hpp
        ${gemmi_src}/include/gemmi/stats.hpp
        ${gemmi_src}/include/gemmi/symmetry.hpp
        ${gemmi_src}/include/gemmi/to_chemcomp.hpp
        ${gemmi_src}/include/gemmi/to_cif.hpp
        ${gemmi_src}/include/gemmi/to_json.hpp
        ${gemmi_src}/include/gemmi/to_mmcif.hpp
        ${gemmi_src}/include/gemmi/to_pdb.hpp
        ${gemmi_src}/include/gemmi/topo.hpp
        ${gemmi_src}/include/gemmi/twin.hpp
        ${gemmi_src}/include/gemmi/unitcell.hpp
        ${gemmi_src}/include/gemmi/utf.hpp
        ${gemmi_src}/include/gemmi/util.hpp
        ${gemmi_src}/include/gemmi/version.hpp
        ${gemmi_src}/include/gemmi/xds_ascii.hpp


)

set(gemmi_third_party-headers_HEADERS
        ${gemmi_src}/include/gemmi/third_party/fast_float.h
        ${gemmi_src}/include/gemmi/third_party/pocketfft_hdronly.h
        ${gemmi_src}/include/gemmi/third_party/sajson.h
        ${gemmi_src}/include/gemmi/third_party/tinydir.h
)

set(gemmi_third_party_tao-headers_HEADERS
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl.hpp
)

set(gemmi_third_party_tao_pegtl-headers_HEADERS
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/analyze.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/apply_mode.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/argv_input.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/ascii.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/buffer_input.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/config.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/cstream_input.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/eol.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/eol_pair.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/file_input.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/input_error.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/istream_input.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/memory_input.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/mmap_input.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/normal.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/nothing.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/parse.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/parse_error.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/position.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/read_input.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/rewind_mode.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/rules.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/string_input.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/tracking_mode.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/utf16.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/utf32.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/utf8.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/version.hpp
)

set(gemmi_third_party_tao_pegtl_analysis-headers_HEADERS
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/analysis/analyze_cycles.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/analysis/counted.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/analysis/generic.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/analysis/grammar_info.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/analysis/insert_guard.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/analysis/insert_rules.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/analysis/rule_info.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/analysis/rule_type.hpp
)

set(gemmi_third_party_tao_pegtl_internal-headers_HEADERS
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/action.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/action_input.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/alnum.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/alpha.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/any.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/apply.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/apply0.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/apply0_single.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/apply_single.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/at.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/bof.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/bol.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/bump_help.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/bump_impl.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/bytes.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/control.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/cr_crlf_eol.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/cr_eol.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/crlf_eol.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/cstream_reader.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/cstring_reader.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/demangle.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/demangle_cxxabi.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/demangle_nop.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/demangle_sanitise.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/disable.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/discard.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/dusel_mode.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/duseltronik.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/enable.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/endian.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/endian_gcc.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/endian_win.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/eof.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/eol.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/eolf.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/file_mapper.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/file_opener.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/file_reader.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/has_apply.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/has_apply0.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/identifier.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/if_apply.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/if_must.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/if_must_else.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/if_then_else.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/input_pair.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/integer_sequence.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/istream_reader.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/istring.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/iterator.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/lf_crlf_eol.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/lf_eol.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/list.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/list_must.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/list_tail.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/list_tail_pad.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/marker.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/minus.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/must.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/not_at.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/one.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/opt.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/pad.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/pad_opt.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/peek_char.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/peek_utf16.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/peek_utf32.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/peek_utf8.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/pegtl_string.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/plus.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/raise.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/range.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/ranges.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/rep.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/rep_min.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/rep_min_max.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/rep_opt.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/require.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/result_on_found.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/rule_conjunction.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/rules.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/seq.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/skip_control.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/sor.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/star.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/star_must.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/state.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/string.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/trivial.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/try_catch_type.hpp
        ${gemmi_src}/include/gemmi/third_party/tao/pegtl/internal/until.hpp
)

set_target_properties(gemmi_cpp PROPERTIES PUBLIC_HEADER "${gemmi_HEADERS}")

install(TARGETS gemmi_cpp
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/gemmi
)

install(FILES
        ${gemmi_third_party-headers_HEADERS}
        DESTINATION include/gemmi/third_party
)

install(FILES
        ${gemmi_third_party_tao-headers_HEADERS}
        DESTINATION include/gemmi/third_party/tao
)

install(FILES
        ${gemmi_third_party_tao_pegtl-headers_HEADERS}
        DESTINATION include/gemmi/third_party/tao/pegtl
)

install(FILES
        ${gemmi_third_party_tao_pegtl_analysis-headers_HEADERS}
        DESTINATION include/gemmi/third_party/tao/pegtl/analysis
)

install(FILES
        ${gemmi_third_party_tao_pegtl_internal-headers_HEADERS}
        DESTINATION include/gemmi/third_party/tao/pegtl/internal
)