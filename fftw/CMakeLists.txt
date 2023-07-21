project(fftw)

add_library(fftw STATIC 
../checkout/fftw-2.1.5/fftw/config.c
../checkout/fftw-2.1.5/fftw/executor.c
../checkout/fftw-2.1.5/fftw/fftwf77.c
../checkout/fftw-2.1.5/fftw/fftwnd.c
../checkout/fftw-2.1.5/fftw/fn_1.c
../checkout/fftw-2.1.5/fftw/fn_10.c
../checkout/fftw-2.1.5/fftw/fn_11.c
../checkout/fftw-2.1.5/fftw/fn_12.c
../checkout/fftw-2.1.5/fftw/fn_13.c
../checkout/fftw-2.1.5/fftw/fn_14.c
../checkout/fftw-2.1.5/fftw/fn_15.c
../checkout/fftw-2.1.5/fftw/fn_16.c
../checkout/fftw-2.1.5/fftw/fn_2.c
../checkout/fftw-2.1.5/fftw/fn_3.c
../checkout/fftw-2.1.5/fftw/fn_32.c
../checkout/fftw-2.1.5/fftw/fn_4.c
../checkout/fftw-2.1.5/fftw/fn_5.c
../checkout/fftw-2.1.5/fftw/fn_6.c
../checkout/fftw-2.1.5/fftw/fn_64.c
../checkout/fftw-2.1.5/fftw/fn_7.c
../checkout/fftw-2.1.5/fftw/fn_8.c
../checkout/fftw-2.1.5/fftw/fn_9.c
../checkout/fftw-2.1.5/fftw/fni_1.c
../checkout/fftw-2.1.5/fftw/fni_10.c
../checkout/fftw-2.1.5/fftw/fni_11.c
../checkout/fftw-2.1.5/fftw/fni_12.c
../checkout/fftw-2.1.5/fftw/fni_13.c
../checkout/fftw-2.1.5/fftw/fni_14.c
../checkout/fftw-2.1.5/fftw/fni_15.c
../checkout/fftw-2.1.5/fftw/fni_16.c
../checkout/fftw-2.1.5/fftw/fni_2.c
../checkout/fftw-2.1.5/fftw/fni_3.c
../checkout/fftw-2.1.5/fftw/fni_32.c
../checkout/fftw-2.1.5/fftw/fni_4.c
../checkout/fftw-2.1.5/fftw/fni_5.c
../checkout/fftw-2.1.5/fftw/fni_6.c
../checkout/fftw-2.1.5/fftw/fni_64.c
../checkout/fftw-2.1.5/fftw/fni_7.c
../checkout/fftw-2.1.5/fftw/fni_8.c
../checkout/fftw-2.1.5/fftw/fni_9.c
../checkout/fftw-2.1.5/fftw/ftw_10.c
../checkout/fftw-2.1.5/fftw/ftw_16.c
../checkout/fftw-2.1.5/fftw/ftw_2.c
../checkout/fftw-2.1.5/fftw/ftw_3.c
../checkout/fftw-2.1.5/fftw/ftw_32.c
../checkout/fftw-2.1.5/fftw/ftw_4.c
../checkout/fftw-2.1.5/fftw/ftw_5.c
../checkout/fftw-2.1.5/fftw/ftw_6.c
../checkout/fftw-2.1.5/fftw/ftw_64.c
../checkout/fftw-2.1.5/fftw/ftw_7.c
../checkout/fftw-2.1.5/fftw/ftw_8.c
../checkout/fftw-2.1.5/fftw/ftw_9.c
../checkout/fftw-2.1.5/fftw/ftwi_10.c
../checkout/fftw-2.1.5/fftw/ftwi_16.c
../checkout/fftw-2.1.5/fftw/ftwi_2.c
../checkout/fftw-2.1.5/fftw/ftwi_3.c
../checkout/fftw-2.1.5/fftw/ftwi_32.c
../checkout/fftw-2.1.5/fftw/ftwi_4.c
../checkout/fftw-2.1.5/fftw/ftwi_5.c
../checkout/fftw-2.1.5/fftw/ftwi_6.c
../checkout/fftw-2.1.5/fftw/ftwi_64.c
../checkout/fftw-2.1.5/fftw/ftwi_7.c
../checkout/fftw-2.1.5/fftw/ftwi_8.c
../checkout/fftw-2.1.5/fftw/ftwi_9.c
../checkout/fftw-2.1.5/fftw/generic.c
../checkout/fftw-2.1.5/fftw/malloc.c
../checkout/fftw-2.1.5/fftw/planner.c
../checkout/fftw-2.1.5/fftw/putils.c
../checkout/fftw-2.1.5/fftw/rader.c
../checkout/fftw-2.1.5/fftw/timer.c
../checkout/fftw-2.1.5/fftw/twiddle.c
../checkout/fftw-2.1.5/fftw/wisdom.c
../checkout/fftw-2.1.5/fftw/wisdomio.c
)

set(fftw_HEADERS
${CMAKE_SOURCE_DIR}/checkout/fftw-2.1.5//fftw/fftw.h
${CMAKE_SOURCE_DIR}/checkout/fftw-2.1.5//fftw/fftw-int.h
${CMAKE_SOURCE_DIR}/checkout/fftw-2.1.5//fftw/config.h
)

set(fftw_HEADERS
${CMAKE_SOURCE_DIR}/checkout/fftw-2.1.5/fftw/fftw.h
${CMAKE_SOURCE_DIR}/checkout/fftw-2.1.5/fftw/fftw-int.h
${CMAKE_SOURCE_DIR}/checkout/fftw-2.1.5/fftw/config.h
)

target_compile_options(fftw PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(fftw PROPERTIES PUBLIC_HEADER "${fftw_HEADERS}")

install(TARGETS fftw
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/fftw
)
