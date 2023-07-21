project(rfftw)

add_library(rfftw STATIC 
../checkout/fftw-2.1.5/rfftw/fcr_1.c
../checkout/fftw-2.1.5/rfftw/fcr_10.c
../checkout/fftw-2.1.5/rfftw/fcr_11.c
../checkout/fftw-2.1.5/rfftw/fcr_12.c
../checkout/fftw-2.1.5/rfftw/fcr_128.c
../checkout/fftw-2.1.5/rfftw/fcr_13.c
../checkout/fftw-2.1.5/rfftw/fcr_14.c
../checkout/fftw-2.1.5/rfftw/fcr_15.c
../checkout/fftw-2.1.5/rfftw/fcr_16.c
../checkout/fftw-2.1.5/rfftw/fcr_2.c
../checkout/fftw-2.1.5/rfftw/fcr_3.c
../checkout/fftw-2.1.5/rfftw/fcr_32.c
../checkout/fftw-2.1.5/rfftw/fcr_4.c
../checkout/fftw-2.1.5/rfftw/fcr_5.c
../checkout/fftw-2.1.5/rfftw/fcr_6.c
../checkout/fftw-2.1.5/rfftw/fcr_64.c
../checkout/fftw-2.1.5/rfftw/fcr_7.c
../checkout/fftw-2.1.5/rfftw/fcr_8.c
../checkout/fftw-2.1.5/rfftw/fcr_9.c
../checkout/fftw-2.1.5/rfftw/fhb_10.c
../checkout/fftw-2.1.5/rfftw/fhb_16.c
../checkout/fftw-2.1.5/rfftw/fhb_2.c
../checkout/fftw-2.1.5/rfftw/fhb_3.c
../checkout/fftw-2.1.5/rfftw/fhb_32.c
../checkout/fftw-2.1.5/rfftw/fhb_4.c
../checkout/fftw-2.1.5/rfftw/fhb_5.c
../checkout/fftw-2.1.5/rfftw/fhb_6.c
../checkout/fftw-2.1.5/rfftw/fhb_7.c
../checkout/fftw-2.1.5/rfftw/fhb_8.c
../checkout/fftw-2.1.5/rfftw/fhb_9.c
../checkout/fftw-2.1.5/rfftw/fhf_10.c
../checkout/fftw-2.1.5/rfftw/fhf_16.c
../checkout/fftw-2.1.5/rfftw/fhf_2.c
../checkout/fftw-2.1.5/rfftw/fhf_3.c
../checkout/fftw-2.1.5/rfftw/fhf_32.c
../checkout/fftw-2.1.5/rfftw/fhf_4.c
../checkout/fftw-2.1.5/rfftw/fhf_5.c
../checkout/fftw-2.1.5/rfftw/fhf_6.c
../checkout/fftw-2.1.5/rfftw/fhf_7.c
../checkout/fftw-2.1.5/rfftw/fhf_8.c
../checkout/fftw-2.1.5/rfftw/fhf_9.c
../checkout/fftw-2.1.5/rfftw/frc_1.c
../checkout/fftw-2.1.5/rfftw/frc_10.c
../checkout/fftw-2.1.5/rfftw/frc_11.c
../checkout/fftw-2.1.5/rfftw/frc_12.c
../checkout/fftw-2.1.5/rfftw/frc_128.c
../checkout/fftw-2.1.5/rfftw/frc_13.c
../checkout/fftw-2.1.5/rfftw/frc_14.c
../checkout/fftw-2.1.5/rfftw/frc_15.c
../checkout/fftw-2.1.5/rfftw/frc_16.c
../checkout/fftw-2.1.5/rfftw/frc_2.c
../checkout/fftw-2.1.5/rfftw/frc_3.c
../checkout/fftw-2.1.5/rfftw/frc_32.c
../checkout/fftw-2.1.5/rfftw/frc_4.c
../checkout/fftw-2.1.5/rfftw/frc_5.c
../checkout/fftw-2.1.5/rfftw/frc_6.c
../checkout/fftw-2.1.5/rfftw/frc_64.c
../checkout/fftw-2.1.5/rfftw/frc_7.c
../checkout/fftw-2.1.5/rfftw/frc_8.c
../checkout/fftw-2.1.5/rfftw/frc_9.c
../checkout/fftw-2.1.5/rfftw/rconfig.c
../checkout/fftw-2.1.5/rfftw/rexec.c
../checkout/fftw-2.1.5/rfftw/rexec2.c
../checkout/fftw-2.1.5/rfftw/rfftwf77.c
../checkout/fftw-2.1.5/rfftw/rfftwnd.c
../checkout/fftw-2.1.5/rfftw/rgeneric.c
../checkout/fftw-2.1.5/rfftw/rplanner.c
)

target_include_directories(rfftw PRIVATE ../checkout/fftw-2.1.5/fftw)

set(rfftw_HEADERS
${CMAKE_SOURCE_DIR}/checkout/fftw-2.1.5/rfftw/rfftw.h
)
target_compile_options(rfftw PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(rfftw PROPERTIES PUBLIC_HEADER "${rfftw_HEADERS}")

install(TARGETS rfftw
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/rfftw
)
