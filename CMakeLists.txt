cmake_minimum_required(VERSION 2.8)
project(privateer C CXX)

find_package(CCP4 REQUIRED mmdb ccp4c clipper-core clipper-mmdb clipper-minimol
                           clipper-contrib clipper-ccp4 clipper-cif)

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(${CCP4_INCLUDE_DIRS})

add_executable(privateer-validate
               clipper-glyco.cpp clipper-glyco_data.cpp
               privateer-lib.cpp privateer-validate.cpp)

target_link_libraries(privateer-validate ${CCP4_LIBRARIES})

install(TARGETS privateer-validate DESTINATION bin)

