## Project settings
#
cmake_minimum_required(VERSION 3.12)
project(privateer VERSION 0.5 LANGUAGES C CXX)

SET(PRIVATEER_CC_FLAGS "-O2 -w")
SET(PRIVATEER_CXX_FLAGS "-O2 -w -fwasm-exceptions -I${CMAKE_INSTALL_PREFIX}/include")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PRIVATEER_CC_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PRIVATEER_CXX_FLAGS} ")

option(BUILD_SHARED_LIBS "Build as a shared library" ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)

set(PRIVATEER_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/privateer)
set(PRIVATEERDATA ${CMAKE_SOURCE_DIR}/data)

add_link_options(-sDEFAULT_PTHREAD_STACK_SIZE=2MB -sSTACK_SIZE=5MB)
add_compile_options(--preload-file data/linkage_torsions@/)

set(CLIBDENV $ENV{CLIBD})

  add_subdirectory(gemmi)
  add_subdirectory(ccp4)
  add_subdirectory(rfftw)
  add_subdirectory(fftw)
  add_subdirectory(mmdb2)
  add_subdirectory(clipper/minimol)
  add_subdirectory(clipper/core)
  add_subdirectory(clipper/ccp4)
  add_subdirectory(clipper/contrib)
  add_subdirectory(clipper/cns)
  add_subdirectory(clipper/mmdb)
  add_subdirectory(clipper/cif)
  add_subdirectory(clipper/phs)
  add_subdirectory(clipper/gemmi)

set(PRIVATEER_SOURCE ${PRIVATEER_SOURCE_DIR}/cpp/privateer.cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}
 checkout/mmdb2 
 checkout/ssm/ssm 
 checkout/libccp4/ccp4 
 checkout/clipper 
 checkout/gemmi/include
 dependencies/gemmi
 checkout/zlib
 ${PRIVATEER_SOURCE_DIR}/cpp
 ${PRIVATEER_SOURCE_DIR}/cpp/pybind11
  ${PRIVATEER_SOURCE_DIR}/cpp/third-party/simdjson/

#  checkout/ccp4srs   
 )

 link_directories( checkout/mmdb2 
 checkout/ssm/ssm 
 checkout/libccp4/ccp4 
 checkout/clipper 
#  checkout/ccp4srs 
 checkout/gemmi/include
 ${PRIVATEER_SOURCE_DIR}/cpp
 ${PRIVATEER_SOURCE_DIR}/cpp/pybind11)


# add_subdirectory(${CMAKE_SOURCE_DIR}/dependencies/pybind11)

find_package(Threads REQUIRED)

#Find Privateer's main.cpp file

# #Location of include files to be linked
# include_directories(
#   # ${CMAKE_SOURCE_DIR}/dependencies/include
#                     # ${CMAKE_SOURCE_DIR}/dependencies/gemmi/include
#                     # ${PYBIND11_INCLUDE_DIR}
#                     # ${PYTHON_INCLUDE_DIRS}
#                     # ${PRIVATEER_SOURCE_DIR}/cpp
#                     # ${PRIVATEER_SOURCE_DIR}/cpp/pybind11
#                     )

# file(COPY ${PRIVATEERDATA}/glycomics/privateer_glycomics_database.json 
#     DESTINATION ${CLIBDENV})

# file(COPY ${PRIVATEERDATA}/linkage_torsions/privateer_torsion_database.json 
#     DESTINATION ${CLIBDENV})

add_library(privateer_lib
            ${PRIVATEER_SOURCE_DIR}/cpp/third-party/gemmi_compilation_unit.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-error.h
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-json.h
            ${PRIVATEER_SOURCE_DIR}/cpp/clipper-glyco.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-db.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-blobs.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-composition.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-dbquery.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-cryo_em.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-xray.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-modelling.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-interactions.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/privateer-lib.cpp
            # ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pybind11.cpp
            # ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-ccp4mg.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-restraints.cpp
            # ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pyanalysis.cpp
            # ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pymodelling.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/clipper-glyco_data.cpp
           # ${PRIVATEER_SOURCE_DIR}/cpp/privateer-bind.cpp
            ${PRIVATEER_SOURCE_DIR}/cpp/third-party/simdjson/simdjson.cpp
            )


add_executable(privateer_exec ${PRIVATEER_SOURCE_DIR}/cpp/privateer-bind.cpp)
               

set_target_properties (privateer_exec  PROPERTIES OUTPUT_NAME privateer )

target_link_options(privateer_exec PRIVATE --preload-file data/linkage_torsions@/ -sDEMANGLE_SUPPORT=0 -sEXIT_RUNTIME=1 -sASSERTIONS=2 -fexceptions -sALLOW_MEMORY_GROWTH=1 --bind -sFORCE_FILESYSTEM=1 -sMODULARIZE=1 -sEXPORT_NAME=privateer_module -sEXPORTED_RUNTIME_METHODS=['FS'] -sMODULARIZE=1 -sEXPORT_ES6=1 )

target_link_directories(privateer_exec
 PUBLIC 
 mmdb2 
 clipper/mmdb 
 clipper/cif 
 clipper/core 
 clipper/contrib 
 clipper/ccp4 
 clipper/gemmi
 clipper/minimol
 fftw 
 rfftw 
 ccp4 
 ../checkout/gemmi
 )


target_link_libraries(
  privateer_exec
  privateer_lib
 mmdb2
 clipper-mmdb
 clipper-cif
 clipper-core 
 clipper-contrib 
 clipper-ccp4 
 clipper-gemmi
 clipper-minimol
 fftw 
 rfftw 
 ccp4 
 nodefs.js 
 Threads::Threads
 )
                               

# set_target_properties ( privateer_exec  PROPERTIES OUTPUT_NAME privateer )
# SET_SOURCE_FILES_PROPERTIES(privateer-frontend_scripting.i PROPERTIES CPLUSPLUS ON)

set(privateer_HEADERS
${PRIVATEER_SOURCE_DIR}/cpp/privateer-restraints.h
${PRIVATEER_SOURCE_DIR}/cpp/privateer-xray.h
${PRIVATEER_SOURCE_DIR}/cpp/privateer-blobs.h
${PRIVATEER_SOURCE_DIR}/cpp/privateer-dbquery.h
${PRIVATEER_SOURCE_DIR}/cpp/privateer-error.h
${PRIVATEER_SOURCE_DIR}/cpp/clipper-glyco_data.h
${PRIVATEER_SOURCE_DIR}/cpp/privateer-composition.h
${PRIVATEER_SOURCE_DIR}/cpp/privateer-lib.h
${PRIVATEER_SOURCE_DIR}/cpp/privateer-json.h
${PRIVATEER_SOURCE_DIR}/cpp/clipper-glyco.h
${PRIVATEER_SOURCE_DIR}/cpp/privateer-parallelism.h
${PRIVATEER_SOURCE_DIR}/cpp/privateer-cryo_em.h
${PRIVATEER_SOURCE_DIR}/cpp/third-party/simdjson/simdjson.h
)

set_target_properties(privateer_exec PROPERTIES PUBLIC_HEADER "${privateer_HEADERS}")
# target_compile_options(privateer_exec PUBLIC --preload-file /home/jordan/dev/privateer_wasm/privateer/data/linkage_torsions@/)

set_target_properties(privateer_exec PROPERTIES COMPILE_FLAGS "-s ENVIRONMENT='web'")
install(TARGETS privateer_lib
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/privateer
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/privateer.wasm
${CMAKE_CURRENT_BINARY_DIR}/privateer.js
DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/web/privateer_wasm/wasm)

set(PRIVATEER_SOURCE ${PRIVATEER_SOURCE_DIR}/cpp/privateer.cpp)

# install ( TARGETS privateer_lib LIBRARY DESTINATION ${PROJECT_SOURCE_DIR})
# install ( TARGETS privateer_exec DESTINATION ${PROJECT_SOURCE_DIR})
# install ( TARGETS privateer_core DESTINATION ${PROJECT_SOURCE_DIR} )
# install ( TARGETS privateer_modelling DESTINATION ${PROJECT_SOURCE_DIR} )