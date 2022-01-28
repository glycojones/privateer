/*
     ccp4_sysdep.h: System-dependent definitions
     Copyright (C) 2001  CCLRC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
*/

/** @file ccp4_sysdep.h
 *
 *  @brief System-dependent definitions.
 *
 *  @author Charles Ballard, based in part on earlier versions
 */

#ifndef __CCP4_BITS
#define __CCP4_BITS

#ifdef __sgi   /* in ANSI mode */
#  ifndef sgi
#    define sgi
#  endif
#endif

#ifndef VMS
#  if defined (vms) || defined (__vms) || defined (__VMS)
#    define VMS
#  endif
#endif

#if defined (sun) || defined (__sun)
#  if !defined(__STDC__) || defined(__GNUC__)
#    if !defined(G77)
      extern char *sys_errlist [];
#     define strerror(i) sys_errlist[i] /* k&r compiler doesn't have it */
#    endif
#  endif
#endif

#if defined (_AIX) || defined(___AIX) || defined (__hpux)
#  define CALL_LIKE_HPUX 1
#elif defined (VMS)
#  define CALL_LIKE_VMS 1
#elif defined (_MSC_VER) && (_MSC_VER >= 800)
#  define CALL_LIKE_MVS 2
#elif defined (_MSC_VER) || (defined (_WIN32) && !defined(__MINGW32__))
#  define CALL_LIKE_MVS 1
#else
#  define CALL_LIKE_SUN 1
#endif

#ifndef _POSIX_SOURCE
#define _POSIX_SOURCE
#endif

#include <stdio.h>

#if defined (VMS)
#  include <descrip.h>          /* non-POSIX */
#  define NOUNISTD
#else
#  include <sys/types.h>
#  include <sys/stat.h>
#  if !defined (_WIN32) && !defined (_MSC_VER)
#    include <sys/times.h>
#  endif
#  ifdef _MSC_VER
#    define NOUNISTD
#  endif
#endif

#include <stddef.h>
#include <string.h>

#ifndef NOUNISTD
#  include <unistd.h>
#else
#  ifndef VMS 
#    ifndef _MSC_VER
#      include <sys/file.h>     /* ESV, old Concentrix */ /* non-POSIX */
#    endif
#  endif
#endif
#ifndef NOSTDLIB                /* for TitanOS 4.2, at least? */
#  include <stdlib.h>
#endif

#include <errno.h>
#include <ctype.h>

#if defined(_AIX) || defined (__hpux) || defined(F2C) ||\
    defined(G77) || defined(_WIN32) || defined (sun) /* would do no harm on others, though */
#  include <time.h>
#endif

#include <limits.h>
#include <float.h>

#if defined (F2C)
#  define Skip_f2c_Undefs
#  include "f2c.h"
#endif
#if defined (G77)
#  define Skip_f2c_Undefs       /* g2c.h infelicity... */
#  if defined (HAVE_G2C_H)
#    include "g2c.h"
#  endif
#endif

/* Using MSVC need __declspec */
#if defined(__WIN32__) || defined(_WIN32)
#  if defined(_MSC_VER) && defined(DLL_EXPORT)
#    define CCP4_DL_IMPORT(type) __declspec(dllexport) type
#    define CCP4_DL_EXPORT __declspec(dllexport)
#  elif defined(_MSC_VER)
#    define CCP4_DL_IMPORT(type) __declspec(dllimport) type
#    define CCP4_DL_EXPORT
#  else
#    define CCP4_DL_IMPORT(type) type
#    define CCP4_DL_EXPORT
#  endif
#else
#  define CCP4_DL_IMPORT(type) type
#  define CCP4_DL_EXPORT
#endif

/* defined in library_utils.c */
#if defined (_MSC_VER) &&  _MSC_VER < 1800
  double rint(double x);
#endif

#ifdef _MSC_VER
#define  M_PI            3.14159265358979323846
#endif

#ifdef _WIN32
#  define PATH_SEPARATOR '\\'
#  define EXT_SEPARATOR '.'
#else
#  define PATH_SEPARATOR '/'
#  define EXT_SEPARATOR '.'
#endif

#define MAXFLEN       512    /**< the maximum length of a filename in CCP4 */
#define MAXFILES       16    /**< maximum number of files open symultaneously */
#define DEFMODE         2    /**< default mode access for random access files */

#define IRRELEVANT_OP   0
#define READ_OP         1
#define WRITE_OP        2

#include<fcntl.h>
#ifndef SEEK_SET
#  define SEEK_SET 0
#  define SEEK_CUR 1
#  define SEEK_END 2
#endif /* ! SEEK_SET */
#ifndef O_WRONLY
#define O_RDONLY 0x0000       /**< i/o mode: read-only */
#define O_WRONLY 0x0001       /**< i/o mode: write-only  */
#define O_RDWR   0x0002       /**< i/o mode: read and write  */
#define O_APPEND 0x0008       /**< i/o mode: append to existing file  */
#define O_CREAT  0x0200       /**< i/o mode: create file  */
#define O_TRUNC  0x0400       /**< i/o mode: truncate existing file  */
#endif
#define O_TMP    0x0010       /**< i/o mode: scratch file */

/* Before version 6.3 we defined BYTE, INT16 and INT32 (without the CCP4_
 * prefix). The prefix has been added to avoid name conflicts.
 */
#define CCP4_BYTE  0
#define CCP4_INT16 1
#define CCP4_INT32 6
#define FLOAT32 2
#define COMP32  3
#define COMP64  4

#define DFNTI_MBO       1       /**< Motorola byte order 2's compl */
#define DFNTI_IBO       4       /**< Intel byte order 2's compl */

#define DFNTF_BEIEEE    1       /**< big endian IEEE (canonical) */
#define DFNTF_VAX       2       /**< Vax format */
#define DFNTF_CONVEXNATIVE 5    /**< Convex native floats */
#define DFNTF_LEIEEE    4       /**< little-endian IEEE format */

#if defined (VAX) || defined (vax) /* gcc seems to use vax */
#  define NATIVEFT DFNTF_VAX
#  define NATIVEIT DFNTI_IBO
#endif

#if defined(MIPSEL) || defined(i386) || defined(i860) || defined(__ia64__) || defined(__amd64__) || defined(__x86_64__) || defined(_M_AMD64)
#  define NATIVEIT DFNTI_IBO
#  define NATIVEFT DFNTF_LEIEEE
#endif

#if defined (powerpc) || defined (__powerpc__) || defined (__ppc__) || \
      defined __PPC || defined (__s390__) || defined (__s390x__) || \
      defined (__hppa__)
#  define NATIVEIT DFNTI_MBO
#  define NATIVEFT DFNTF_BEIEEE
#endif

#ifdef __alpha
#  ifdef VMS
#    if __IEEE_FLOAT == 1
#      define NATIVEFT DFNTF_LEIEEE
#    else
#      define NATIVEFT DFNTF_VAX
#    endif
#  else                       /* assume OSF/1 */
#    define NATIVEFT DFNTF_LEIEEE
#  endif
#  define NATIVEIT DFNTI_IBO
#endif

#if defined(MIPSEB) || defined(__hpux) || defined(_AIX) || defined(m68k) || defined(mc68000) || defined(sparc) || defined (__sparc__)
#  define NATIVEIT DFNTI_MBO
#  define NATIVEFT DFNTF_BEIEEE
#endif

#if defined(__ARM__) || defined(__arm__) || defined(__aarch64__)
# if defined(__ARMEB__) || defined (__AARCH64EB__)
#  define NATIVEIT DFNTI_MBO
#  define NATIVEFT DFNTF_BEIEEE
# elif defined(__ARMEL__) || defined (__AARCH64EL__)
#  define NATIVEIT DFNTI_IBO
#  define NATIVEFT DFNTF_LEIEEE
# endif
#endif

/* From time to time new architectures are added here, often because Linux
 * packagers want to build it on all platforms supported by their distro. 
 * Here we try to catch machines not listed explicitely above, under
 * assumption that endianness is the same for floating point numbers
 * as for integers. Which is safe assumption on modern standard computers
 * (not embedded systems), according to
 * http://en.wikipedia.org/wiki/Endianness#Floating-point_and_endianness
 */
#if !defined(NATIVEIT) && !defined(NATIVEFT) && defined(__BYTE_ORDER)
# if __BYTE_ORDER == __LITTLE_ENDIAN
#  define NATIVEIT DFNTI_IBO
#  define NATIVEFT DFNTF_LEIEEE
# elif __BYTE_ORDER == __BIG_ENDIAN
#  define NATIVEIT DFNTI_MBO
#  define NATIVEFT DFNTF_BEIEEE
# endif
#endif

#ifndef NATIVEFT
#  error "Can't determine machine number format"
#endif

#define DFNT_UINT       0       /**< unsigned int */
#define DFNT_SINT       1       /**< short int */
#define DFNT_INT        2       /**< int */
#define DFNT_UCHAR      3       /**< unsigned char */
#define DFNT_CHAR       4       /**< char */
#define DFNT_FLOAT      5       /**< float */
#define DFNT_DOUBLE     6       /**< double */

#endif /* __CCP4_BITS */
