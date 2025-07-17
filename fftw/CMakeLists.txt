project(fftw)

add_library(fftw STATIC 
${fftw-dependencies_SOURCE_DIR}//fftw/config.c
${fftw-dependencies_SOURCE_DIR}//fftw/executor.c
${fftw-dependencies_SOURCE_DIR}//fftw/fftwf77.c
${fftw-dependencies_SOURCE_DIR}//fftw/fftwnd.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_1.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_10.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_11.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_12.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_13.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_14.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_15.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_16.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_2.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_3.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_32.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_4.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_5.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_6.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_64.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_7.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_8.c
${fftw-dependencies_SOURCE_DIR}//fftw/fn_9.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_1.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_10.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_11.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_12.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_13.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_14.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_15.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_16.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_2.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_3.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_32.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_4.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_5.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_6.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_64.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_7.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_8.c
${fftw-dependencies_SOURCE_DIR}//fftw/fni_9.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftw_10.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftw_16.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftw_2.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftw_3.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftw_32.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftw_4.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftw_5.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftw_6.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftw_64.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftw_7.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftw_8.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftw_9.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftwi_10.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftwi_16.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftwi_2.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftwi_3.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftwi_32.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftwi_4.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftwi_5.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftwi_6.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftwi_64.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftwi_7.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftwi_8.c
${fftw-dependencies_SOURCE_DIR}//fftw/ftwi_9.c
${fftw-dependencies_SOURCE_DIR}//fftw/generic.c
${fftw-dependencies_SOURCE_DIR}//fftw/malloc.c
${fftw-dependencies_SOURCE_DIR}//fftw/planner.c
${fftw-dependencies_SOURCE_DIR}//fftw/putils.c
${fftw-dependencies_SOURCE_DIR}//fftw/rader.c
${fftw-dependencies_SOURCE_DIR}//fftw/timer.c
${fftw-dependencies_SOURCE_DIR}//fftw/twiddle.c
${fftw-dependencies_SOURCE_DIR}//fftw/wisdom.c
${fftw-dependencies_SOURCE_DIR}//fftw/wisdomio.c
)

set(fftw_HEADERS
${fftw-dependencies_SOURCE_DIR}///fftw/fftw.h
${fftw-dependencies_SOURCE_DIR}///fftw/fftw-int.h
${fftw-dependencies_SOURCE_DIR}///fftw/config.h
)

set(fftw_HEADERS
${fftw-dependencies_SOURCE_DIR}//fftw/fftw.h
${fftw-dependencies_SOURCE_DIR}//fftw/fftw-int.h
${fftw-dependencies_SOURCE_DIR}//fftw/config.h
)

target_compile_options(fftw PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(fftw PROPERTIES PUBLIC_HEADER "${fftw_HEADERS}")

install(TARGETS fftw
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/fftw
)
