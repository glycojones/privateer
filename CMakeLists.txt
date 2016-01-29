cmake_minimum_required(VERSION 2.8)
project(privateer C CXX)

option(BUILD_SHARED_LIBS "Build as a shared library" ON)

if(NOT DEFINED LIB_INSTALL_DIR)
    set(LIB_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)
endif()

find_package(CCP4 REQUIRED mmdb2 ccp4c ccp4srs clipper-core clipper-mmdb clipper-minimol
                           clipper-contrib clipper-ccp4 clipper-cif)

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(${CCP4_INCLUDE_DIRS})

add_library(privateer
            clipper-glyco.cpp
            clipper-glyco_data.cpp
            privateer-lib.cpp)
target_link_libraries(privateer ${CCP4_LIBRARIES})

add_executable(privateer-validate
               privateer.cpp)

target_link_libraries(privateer-validate privateer)

set_target_properties(privateer-validate PROPERTIES OUTPUT_NAME privateer)

install(FILES clipper-glyco.h clipper-glyco_data.h privateer-lib.h
        DESTINATION include/privateer)
install(TARGETS privateer
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION "${LIB_INSTALL_DIR}"
        LIBRARY DESTINATION "${LIB_INSTALL_DIR}")
install(TARGETS privateer-validate DESTINATION bin)
if (UNIX)
  install(CODE "
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink privateer
            \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/privateer-validate\")
        ")
endif()
