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

#ifndef __CCP4_UTILS
#define __CCP4_UTILS

#include <string.h>
#include "ccp4_types.h"
#include "library_file.h"
/* rcsidh[] = "$Id$" */

#ifdef __cplusplus
namespace CCP4 {
extern "C" {
#endif

/****************************************************************************
 * Function prototypes                                                      *
 ****************************************************************************/

size_t ccp4_utils_flength (char *, int);

int ccp4_utils_translate_mode_float(float *, const void *, int, int);

void ccp4_utils_fatal (const char *);

void ccp4_utils_print (const char *message);

int ccp4_utils_setenv (char *);

/* turn on line buffering for stdout */
int ccp4_utils_outbuf (void);

/* turn off any buffering on stdin */
int ccp4_utils_noinpbuf (void);

union float_uint_uchar ccp4_nan ();

int ccp4_utils_isnan (const union float_uint_uchar *);

void ccp4_utils_bml (int, union float_uint_uchar *);

void ccp4_utils_wrg (int, union float_uint_uchar *, float *);

void ccp4_utils_hgetlimits (int *, float *);

int ccp4_utils_mkdir (const char *, const char *);

int ccp4_utils_chmod (const char *, const char *);

void *ccp4_utils_malloc(size_t);

void *ccp4_utils_realloc(void *, size_t);

void *ccp4_utils_calloc(size_t, size_t);

int ccp4_file_size(const char *);

char *ccp4_utils_username(void);

char *ccp4_utils_basename(const char *filename);

char *ccp4_utils_pathname(const char *filename);

char *ccp4_utils_extension(const char *filename);

char *ccp4_utils_joinfilenames(const char *dir, const char *file);

void ccp4_utils_idate (int *);

char *ccp4_utils_date(char *);

void ccp4_utils_itime (int *);

char *ccp4_utils_time(char *);

float ccp4_utils_etime (float *);

#if defined (_MSC_VER)
double ccp4_erfc( double x );
#endif

/****************************************************************************
*  End of prototypes                                                        *
*****************************************************************************/
#ifdef __cplusplus
}
}
#endif

#endif  /* __CCP4_UTILS */
