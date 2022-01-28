AC_INIT(ccp4srs, 1.0.0)

AC_CONFIG_AUX_DIR(build-aux)
AC_CONFIG_MACRO_DIR(m4)
AC_CONFIG_SRCDIR(ccp4srs/ccp4srs_angle.cpp)
AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE([1.11 foreign subdir-objects silent-rules -Wall])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
LT_INIT([disable-shared win32-dll])
AC_PROG_CXX
AC_PROG_INSTALL

AM_MAINTAINER_MODE dnl disable (by default) maintainer mode

AC_SEARCH_LIBS([cos], [m], [], [AC_MSG_ERROR([need math library])])
AC_CHECK_LIB(z, crc32, , AC_MSG_ERROR(zlib library not found))
AC_CHECK_HEADER(zlib.h, :, AC_MSG_ERROR(zlib.h header not found))

AC_LANG_PUSH([C++])
AC_CHECK_HEADER(mmdb2/mmdb_manager.h, :,
                AC_MSG_ERROR([mmdb2 headers not found.]))
dnl to make ccp4/configure script work, missing library at this stage
dnl is not an error
AC_CHECK_LIB(mmdb2, MMDB_CCIF_Init, :,
             AC_MSG_WARN([mmdb2 library not found.]))
AC_LANG_POP([C++])

AC_OUTPUT([Makefile])
