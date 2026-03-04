/*
     ccp4_utils.h: headers for utility functions.
     Copyright (C) 2001  CCLRC, Charles Ballard

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

/** @file ccp4_utils.h
 *  @brief   Utility functions.
 *  @author  Charles Ballard
 */

 /* Modifications by Tristan Croll, 2016-2019:
  *
  * - Native Windows compatibility
  */

#ifndef __CCP4_UTILS
#define __CCP4_UTILS

#include <string.h>
#include "ccp4_types.h"
#include "library_file.h"
/* rcsidh[] = "$Id$" */

#include "imex.h"

#ifdef __cplusplus
namespace CCP4 {
extern "C" {
#endif

/****************************************************************************
 * Function prototypes                                                      *
 ****************************************************************************/

CCP4_IMEX size_t ccp4_utils_flength (char *, int);

CCP4_IMEX int ccp4_utils_translate_mode_float(float *, const void *, int, int);

CCP4_IMEX void ccp4_utils_fatal (const char *);

CCP4_IMEX void ccp4_utils_print (const char *message);

CCP4_IMEX int ccp4_utils_setenv (char *);

/* turn on line buffering for stdout */
CCP4_IMEX int ccp4_utils_outbuf (void);

/* turn off any buffering on stdin */
CCP4_IMEX int ccp4_utils_noinpbuf (void);

CCP4_IMEX union float_uint_uchar ccp4_nan ();

CCP4_IMEX int ccp4_utils_isnan (const union float_uint_uchar *);

CCP4_IMEX void ccp4_utils_bml (int, union float_uint_uchar *);

CCP4_IMEX void ccp4_utils_wrg (int, union float_uint_uchar *, float *);

CCP4_IMEX void ccp4_utils_hgetlimits (int *, float *);

CCP4_IMEX int ccp4_utils_mkdir (const char *, const char *);

CCP4_IMEX int ccp4_utils_chmod (const char *, const char *);

CCP4_IMEX void *ccp4_utils_malloc(size_t);

CCP4_IMEX void *ccp4_utils_realloc(void *, size_t);

CCP4_IMEX void *ccp4_utils_calloc(size_t, size_t);

CCP4_IMEX int ccp4_file_size(const char *);

CCP4_IMEX char *ccp4_utils_username(void);

CCP4_IMEX char *ccp4_utils_basename(const char *filename);

CCP4_IMEX char *ccp4_utils_pathname(const char *filename);

CCP4_IMEX char *ccp4_utils_extension(const char *filename);

CCP4_IMEX char *ccp4_utils_joinfilenames(const char *dir, const char *file);

CCP4_IMEX void ccp4_utils_idate (int *);

CCP4_IMEX char *ccp4_utils_date(char *);

CCP4_IMEX void ccp4_utils_itime (int *);

CCP4_IMEX char *ccp4_utils_time(char *);

CCP4_IMEX float ccp4_utils_etime (float *);

#if defined (_MSC_VER)
CCP4_IMEX double ccp4_erfc( double x );
#endif

/****************************************************************************
*  End of prototypes                                                        *
*****************************************************************************/
#ifdef __cplusplus
}
}
#endif

#endif  /* __CCP4_UTILS */
