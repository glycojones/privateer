cmake_minimum_required(VERSION 3.1)
project(privateer C CXX)

option(BUILD_SHARED_LIBS "Build as a shared library" ON)

set(CMAKE_CXX_STANDARD 11)

if(NOT DEFINED LIB_INSTALL_DIR)
    set(LIB_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib)
endif()

find_package(CCP4 REQUIRED mmdb2 ccp4c ccp4srs clipper-core clipper-mmdb clipper-minimol
                           clipper-contrib clipper-ccp4 clipper-cif)

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

if(APPLE)
    find_program(PYTHON NAMES $ENV{PYTHON} ccp4-python python2.7 python.exe
                        HINTS ${CMAKE_INSTALL_PREFIX}/bin
                        PATHS [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.7\\InstallPath])
else()
    find_program(PYTHON NAMES $ENV{PYTHON} python2.7 python.exe
                        HINTS ${CMAKE_INSTALL_PREFIX}/bin
                        PATHS [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.7\\InstallPath])
endif()


execute_process(COMMAND ${PYTHON} -c "import sys; sys.stdout.write(sys.prefix)"
                OUTPUT_VARIABLE _PY_PREFIX
                RESULT_VARIABLE failure)

if (failure)
    message(FATAL_ERROR "Failed to run ${PYTHON}: ${failure}")
else()
    message(STATUS "Python interpreter: ${PYTHON}")
endif()

message(STATUS "Python sys.prefix: ${_PY_PREFIX}")
find_library(PYTHON_LIBRARY NAMES python2.7 python27
             HINTS ${_PY_PREFIX}/lib ${_PY_PREFIX}/lib64 ${_PY_PREFIX}/libs
             NO_DEFAULT_PATH)
                                                          
if (NOT PYTHON_LIBRARY)
    message(FATAL_ERROR "Python library not found.")
endif()

message(STATUS "Python library: ${PYTHON_LIBRARY}")

find_path(PYTHON_INCLUDE_DIR NAMES Python.h
          HINTS ${_PY_PREFIX}/include
                ${_PY_PREFIX}/include/python2.7
          NO_DEFAULT_PATH)

message(STATUS "Python.h is in: ${PYTHON_INCLUDE_DIR}")

if (NOT PYTHON_INCLUDE_DIR)
    message(FATAL_ERROR "Python header not found.")
endif()

set(PYTHON_LIBRARIES ${PYTHON_LIBRARY})

SET(CMAKE_SWIG_FLAGS "")

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(${CCP4_INCLUDE_DIRS}
                    ${CMAKE_CURRENT_SOURCE_DIR}
                    "${PYTHON_INCLUDE_DIR}")

add_library(privateer
            clipper-glyco.cpp
            clipper-glyco_data.cpp
            privateer-lib.cpp)
target_link_libraries(privateer ${CCP4_LIBRARIES})

add_executable(privateer-validate
               privateer.cpp)

target_link_libraries(privateer-validate privateer)

set_target_properties(privateer-validate PROPERTIES OUTPUT_NAME privateer)

SET_SOURCE_FILES_PROPERTIES(privateer-frontend_scripting.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(privateer python privateer-frontend_scripting.i )

if(WIN32)
    set(SWIG_CLIPPER_MODULE_SUFFIX "pyd")
else()
    set(SWIG_CLIPPER_MODULE_SUFFIX "so")
endif()

SWIG_LINK_LIBRARIES(privateer privateer ${PYTHON_LIBRARIES} ${CLIPPER-MMDB_LIBRARY} ${CLIPPER-CIF_LIBRARY} ${CLIPPER-CCP4_LIBRARY}
        ${CLIPPER-CONTRIB_LIBRARY} ${CLIPPER-CORE_LIBRARIES} ${CLIPPER-MINIMOL_LIBRARY} ${MMDB2_LIBRARY} ${CCP4C_LIBRARY} m)

set ( PYTHON_MODULE_PATH ${CMAKE_INSTALL_PREFIX}/lib/py2)

message(STATUS "Python module path: ${PYTHON_MODULE_PATH}")


install(FILES clipper-glyco.h clipper-glyco_data.h privateer-lib.h
        DESTINATION include/privateer)
install(TARGETS privateer
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION "${LIB_INSTALL_DIR}"
        LIBRARY DESTINATION "${LIB_INSTALL_DIR}")
install(TARGETS privateer-validate DESTINATION bin)

install ( FILES ${CMAKE_CURRENT_BINARY_DIR}/_privateer.${SWIG_CLIPPER_MODULE_SUFFIX} DESTINATION ${PYTHON_MODULE_PATH} )
install ( FILES ${CMAKE_CURRENT_BINARY_DIR}/privateer.py DESTINATION ${PYTHON_MODULE_PATH} )


if (UNIX)
  install(CODE "
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink privateer
            \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/privateer-validate\")
        ")
endif()
