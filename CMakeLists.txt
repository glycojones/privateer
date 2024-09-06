# New Build System CMake File - Jordan Dialpuri - 05/09/24

cmake_minimum_required(VERSION 3.15)
project(
  privateer
  VERSION 0.5
  LANGUAGES C CXX)

set(PRIVATEER_CC_FLAGS "-O2 -w")
set(PRIVATEER_CXX_FLAGS "-O2 -w -I${CMAKE_INSTALL_PREFIX}/include")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PRIVATEER_CC_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PRIVATEER_CXX_FLAGS} ")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(BUILD_SHARED_LIBS "Build as a shared library" ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)

set(PRIVATEER_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/privateer)
set(PRIVATEERDATA ${CMAKE_SOURCE_DIR}/data)
set(CLIBDENV $ENV{CLIBD})

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

set(WRK_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

include(FetchContent)
# Download Clipper
FetchContent_Declare(
  clipper-dependencies
  URL http://www.ysbl.york.ac.uk/jsd523/clipper-2.1.20201109.tar.gz)

# Download MMDB
FetchContent_Declare(mmdb2-dependencies
                     URL http://www.ysbl.york.ac.uk/jsd523/mmdb2-2.0.22.tar.gz)

# Download FFTW
FetchContent_Declare(fftw-dependencies
                     URL http://www.fftw.org/fftw-2.1.5.tar.gz)

# Download CCP4
FetchContent_Declare(ccp4-dependencies
                     URL http://www.ysbl.york.ac.uk/jsd523/libccp4-8.0.0.tar.gz)

# Download Gemmi
FetchContent_Declare(gemmi-dependencies
                     URL http://www.ysbl.york.ac.uk/jsd523/gemmi-0.6.5.tar.gz)

FetchContent_MakeAvailable(
  clipper-dependencies mmdb2-dependencies fftw-dependencies ccp4-dependencies
  gemmi-dependencies)

add_subdirectory(${WRK_DIR}/clipper/minimol)
add_subdirectory(${WRK_DIR}/clipper/core)
add_subdirectory(${WRK_DIR}/clipper/contrib)
add_subdirectory(${WRK_DIR}/clipper/cns)
add_subdirectory(${WRK_DIR}/clipper/mmdb)
add_subdirectory(${WRK_DIR}/clipper/cif)
add_subdirectory(${WRK_DIR}/clipper/phs)
add_subdirectory(${WRK_DIR}/clipper/ccp4)
add_subdirectory(${WRK_DIR}/fftw)
add_subdirectory(${WRK_DIR}/rfftw)
add_subdirectory(${WRK_DIR}/ccp4)
add_subdirectory(${WRK_DIR}/mmdb2)
add_subdirectory(${WRK_DIR}/gemmi)
# add_subdirectory(${WRK_DIR}/dependencies/pybind11)

set(PRIVATEER_SOURCE ${PRIVATEER_SOURCE_DIR}/cpp/privateer.cpp)

# include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ${WRK_DIR}/checkout/mmdb2
# ${WRK_DIR}/checkout/clipper ${WRK_DIR}/checkout/libccp4/ccp4
# ${WRK_DIR}/checkout/gemmi/include ${WRK_DIR}/checkout/ssm/ssm
# ${WRK_DIR}/checkout/zlib ${WRK_DIR}/checkout/pybind11/include
# ${PRIVATEER_SOURCE_DIR}/cpp ${PRIVATEER_SOURCE_DIR}/cpp/pybind11
# ${Python_INCLUDE_DIRS} ${pybind11_INCLUDE_DIRS} )

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${clipper-dependencies_SOURCE_DIR}
  ${gemmi-dependencies_SOURCE_DIR}/include
  ${mmdb2-dependencies_SOURCE_DIR}
  ${ccp4-dependencies_SOURCE_DIR}
  ${PRIVATEER_SOURCE_DIR}/cpp
  ${PRIVATEER_SOURCE_DIR}/cpp/pybind11)

link_directories(
  ${CMAKE_BINARY_DIR}/mmdb2
  ${CMAKE_BINARY_DIR}/clipper
  ${CMAKE_BINARY_DIR}/ccp4
  ${WRK_DIR}/src
  ${CMAKE_BINARY_DIR}
  ${PRIVATEER_SOURCE_DIR}/cpp
  ${PRIVATEER_SOURCE_DIR}/cpp/pybind11
  ${Python_LIBRARY_DIRS})

# link_directories( ${WRK_DIR}/mmdb2 ${WRK_DIR}/clipper ${WRK_DIR}/ccp4
# ${WRK_DIR}/gemmi ${WRK_DIR}/src ${CMAKE_BINARY_DIR}
# ${PRIVATEER_SOURCE_DIR}/cpp ${PRIVATEER_SOURCE_DIR}/cpp/pybind11
# ${Python_LIBRARY_DIRS} )

# add_custom_target(copy_resources ALL COMMAND ${CMAKE_COMMAND} -E
# copy_directory ${CMAKE_SOURCE_DIR}/data ${CLIBDENV}/privateer_data COMMENT
# "Copying Privateer Data files to ${CLIBDENV}/privateer_data")

find_package(Threads REQUIRED)

# file(COPY ${PRIVATEERDATA}/glycomics/privateer_glycomics_database.json
# DESTINATION ${CLIBDENV})

# file(COPY ${PRIVATEERDATA}/linkage_torsions/privateer_torsion_database.json
# DESTINATION ${CLIBDENV})

add_library(
  privateer_lib STATIC
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
  ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pybind11.cpp
  ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-ccp4mg.cpp
  ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-restraints.cpp
  ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pyanalysis.cpp
  ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pymodelling.cpp
  ${PRIVATEER_SOURCE_DIR}/cpp/clipper-glyco_data.cpp)

target_link_directories(
  privateer_lib
  PUBLIC
  ${CMAKE_BINARY_DIR}/clipper/mmdb
  ${CMAKE_BINARY_DIR}/clipper/cif
  ${CMAKE_BINARY_DIR}/clipper/ccp4
  ${CMAKE_BINARY_DIR}/clipper/minimol
  ${CMAKE_BINARY_DIR}/clipper/core
  ${CMAKE_BINARY_DIR}/clipper/contrib
  ${CMAKE_BINARY_DIR}/fftw
  ${CMAKE_BINARY_DIR}/rfftw
  ${CMAKE_BINARY_DIR}/ccp4
  ${CMAKE_BINARY_DIR}/mmdb2
  ${CMAKE_BINARY_DIR}/gemmi)

target_link_libraries(
  privateer_lib
  PUBLIC clipper-contrib
         clipper-mmdb
         clipper-cif
         clipper-ccp4
         clipper-minimol
         clipper-core
         rfftw
         fftw
         ccp4
         mmdb2
         gemmi_cpp
         Threads::Threads
         pybind11::module)

add_executable(privateer_exec ${PRIVATEER_SOURCE_DIR}/cpp/privateer.cpp)

set_target_properties(privateer_exec PROPERTIES OUTPUT_NAME privateer)

pybind11_add_module(
  privateer_core
  ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pybind11.cpp
  ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-ccp4mg.cpp
  ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-restraints.cpp
  ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pyanalysis.cpp
  ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pymodelling.cpp)

pybind11_add_module(
  privateer_modelling
  ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pybind11.cpp
  ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-ccp4mg.cpp
  ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-restraints.cpp
  ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pyanalysis.cpp
  ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-pymodelling.cpp)

target_link_directories(
  privateer_exec
  PUBLIC
  ${CMAKE_BINARY_DIR}/clipper/mmdb
  ${CMAKE_BINARY_DIR}/clipper/cif
  ${CMAKE_BINARY_DIR}/clipper/ccp4
  ${CMAKE_BINARY_DIR}/clipper/minimol
  ${CMAKE_BINARY_DIR}/clipper/core
  ${CMAKE_BINARY_DIR}/clipper/contrib
  ${CMAKE_BINARY_DIR}/fftw
  ${CMAKE_BINARY_DIR}/rfftw
  ${CMAKE_BINARY_DIR}/ccp4
  ${CMAKE_BINARY_DIR}/mmdb2
  ${CMAKE_BINARY_DIR}/gemmi
  ${Python_LIBRARY_DIRS})

target_link_libraries(
  privateer_exec
  privateer_lib
  clipper-contrib
  clipper-mmdb
  clipper-cif
  clipper-ccp4
  clipper-minimol
  clipper-core
  rfftw
  fftw
  ccp4
  mmdb2
  gemmi_cpp
  Threads::Threads
  pybind11::module
  ${Python_LIBRARIES})

target_link_libraries(
  privateer_core
  PRIVATE privateer_lib
          clipper-contrib
          clipper-mmdb
          clipper-cif
          clipper-ccp4
          clipper-minimol
          clipper-core
          rfftw
          fftw
          ccp4
          mmdb2
          gemmi_cpp
          Threads::Threads
          ${Python_LIBRARIES}
          pybind11::module)

target_link_libraries(
  privateer_modelling
  PRIVATE privateer_lib
          clipper-contrib
          clipper-mmdb
          clipper-cif
          clipper-ccp4
          clipper-minimol
          clipper-core
          clipper-minimol
          rfftw
          fftw
          ccp4
          mmdb2
          gemmi_cpp
          Threads::Threads
          ${Python_LIBRARIES}
          pybind11::module)

set(privateer_HEADERS
    ${PRIVATEER_SOURCE_DIR}/cpp/pybind11/privateer-restraints.h
    ${PRIVATEER_SOURCE_DIR}/cpp/privateer-xray.h
    ${PRIVATEER_SOURCE_DIR}/cpp/privateer-blobs.h
    ${PRIVATEER_SOURCE_DIR}/cpp/privateer-dbquery.h
    ${PRIVATEER_SOURCE_DIR}/cpp/privateer-error.h
    ${PRIVATEER_SOURCE_DIR}/cpp/clipper-glyco_data.h
    ${PRIVATEER_SOURCE_DIR}/cpp/privateer-composition.h
    ${PRIVATEER_SOURCE_DIR}/cpp/privateer-lib.h
    ${PRIVATEER_SOURCE_DIR}/cpp/privateer-json.h
    ${PRIVATEER_SOURCE_DIR}/cpp/clipper-glyco.h
    # ${PRIVATEER_SOURCE_DIR}/cpp/privateer-parallelism.h
    ${PRIVATEER_SOURCE_DIR}/cpp/privateer-cryo_em.h)

# set_target_properties(privateer_core PROPERTIES LINK_FLAGS "-undefined
# dynamic_lookup") set_target_properties(privateer_modelling PROPERTIES
# LINK_FLAGS "-undefined dynamic_lookup")

# set_target_properties(privateer_exec PROPERTIES PUBLIC_HEADER
# "${privateer_HEADERS}")
set(PRIVATEER_SOURCE ${PRIVATEER_SOURCE_DIR}/cpp/privateer.cpp)

install(TARGETS privateer_lib LIBRARY DESTINATION privateer)
install(TARGETS privateer_exec DESTINATION ${CMAKE_SOURCE_DIR})
install(TARGETS privateer_core LIBRARY DESTINATION privateer)
install(TARGETS privateer_modelling LIBRARY DESTINATION privateer)
