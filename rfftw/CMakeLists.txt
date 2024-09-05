project(rfftw)

add_library(rfftw STATIC 
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_1.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_10.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_11.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_12.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_128.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_13.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_14.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_15.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_16.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_2.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_3.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_32.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_4.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_5.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_6.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_64.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_7.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_8.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fcr_9.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhb_10.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhb_16.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhb_2.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhb_3.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhb_32.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhb_4.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhb_5.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhb_6.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhb_7.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhb_8.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhb_9.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhf_10.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhf_16.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhf_2.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhf_3.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhf_32.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhf_4.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhf_5.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhf_6.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhf_7.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhf_8.c
${fftw-dependencies_SOURCE_DIR}//rfftw/fhf_9.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_1.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_10.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_11.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_12.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_128.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_13.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_14.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_15.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_16.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_2.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_3.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_32.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_4.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_5.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_6.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_64.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_7.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_8.c
${fftw-dependencies_SOURCE_DIR}//rfftw/frc_9.c
${fftw-dependencies_SOURCE_DIR}//rfftw/rconfig.c
${fftw-dependencies_SOURCE_DIR}//rfftw/rexec.c
${fftw-dependencies_SOURCE_DIR}//rfftw/rexec2.c
${fftw-dependencies_SOURCE_DIR}//rfftw/rfftwf77.c
${fftw-dependencies_SOURCE_DIR}//rfftw/rfftwnd.c
${fftw-dependencies_SOURCE_DIR}//rfftw/rgeneric.c
${fftw-dependencies_SOURCE_DIR}//rfftw/rplanner.c
)

target_include_directories(rfftw PRIVATE ${fftw-dependencies_SOURCE_DIR}//fftw)

set(rfftw_HEADERS
${fftw-dependencies_SOURCE_DIR}//rfftw/rfftw.h
)
target_compile_options(rfftw PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(rfftw PROPERTIES PUBLIC_HEADER "${rfftw_HEADERS}")

install(TARGETS rfftw
LIBRARY DESTINATION lib
PUBLIC_HEADER DESTINATION include/rfftw
)
