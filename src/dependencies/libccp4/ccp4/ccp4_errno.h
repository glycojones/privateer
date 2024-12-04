/*
     ccp4_errno.h: Header file for error handling routines
     Copyright (C) 2001  CCLRC, Charles Ballard and Martyn Winn

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
/** @file ccp4_errno.h
 *  Header file for error handling routines
 *  base error codes on system errors. 
 */

#ifndef __CCP4_ERROR_GUARD
#define __CCP4_ERROR_GUARD

#include <errno.h>

/* rcsidhe[] = "$Id$" */

#ifndef CCP4_ERRSYSTEM
#define CCP4_ERRSYSTEM(x) (((x)&0xfff)<<24)
#endif
#ifndef CCP4_ERRLEVEL
#define CCP4_ERRLEVEL(x)  (((x)&0xf)<<16)
#endif
#ifndef CCP4_ERRSETLEVEL
#define CCP4_ERRSETLEVEL(y,x) ((y) & (~CCP4_ERRLEVEL(0xf)) | CCP4_ERRLEVEL(x)))
#endif
#ifndef CCP4_ERRGETSYS   
#define CCP4_ERRGETSYS(x)   (((x)>>24)&0xfff)
#endif
#ifndef CCP4_ERRGETLEVEL
#define CCP4_ERRGETLEVEL(x) (((x)>>16)&0xf)
#endif
#ifndef CCP4_ERRGETCODE
#define CCP4_ERRGETCODE(x)  ((x)&0xffff)
#endif

#define CCP4_ERR_SYS CCP4_ERRSYSTEM(0x0)
#define CCP4_ERR_FILE CCP4_ERRSYSTEM(0x1)
#define CCP4_ERR_COORD CCP4_ERRSYSTEM(0x2)
#define CCP4_ERR_MTZ CCP4_ERRSYSTEM(0x3)
#define CCP4_ERR_MAP CCP4_ERRSYSTEM(0x4)
#define CCP4_ERR_UTILS CCP4_ERRSYSTEM(0x5)
#define CCP4_ERR_PARS CCP4_ERRSYSTEM(0x6)
#define CCP4_ERR_SYM CCP4_ERRSYSTEM(0x7)
#define CCP4_ERR_GEN CCP4_ERRSYSTEM(0x8)

#define CCP4_COUNT(x) sizeof(x)/sizeof(x[0])

/** @global ccp4_errno: global variable that stores the error last error
 *           code from the ccp4 libraries 
 * | 12 bits - library | 4 bits - level | 16 bits - code |
 *
 *  associated macros
 *    CCP4_ERR_SYS   0     OS error
 *    CCP4_ERR_FILE  1     io library
 *    CCP4_ERR_COORD 2     mmdb
 *    CCP4_ERR_MTZ   3     cmtz
 *    CCP4_ERR_MAP   4     map io
 *    CCP4_ERR_UTILS 5     utility routines
 *    CCP4_ERR_PARS  6     parser routines
 *    CCP4_ERR_SYM   7     csymlib 
 *
 * and bit manipulation
 *    CCP4_ERRSYSTEM   system mask for setting
 *    CCP4_ERRLEVEL    error level mask
 *    CCP4_ERRSETLEVEL error level mask for setting error level
 *    CCP4_ERRGETSYS   mask for returning system
 *    CCP4_ERRGETLEVEL mask for returning level
 *    CCP4_ERRGETCODE  mask for returning the code
 *
 * error levels
 *    0  Success
 *    1  Informational
 *    2  Warning 
 *    3  Error
 *    4  Fatal
 */ 
#ifdef __cplusplus
extern "C" {
#endif
extern int ccp4_errno;
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
namespace CCP4 {
extern "C" {
#endif

/** Print out passed message and internal message based upon
 *  ccp4_errno
 *            "message : error message "
 * @param message (const char *)
 * @return void
 */
void ccp4_error( const char *);

/** Obtain character string based upon error code.
 *  Typical use ccp4_strerror(ccp4_errno)
 *  The returned string is statically allocated in the
 *  library_err.c file and should not be freed.
 * @param error code (int)
 * @return const pointer to error message (const char *) 
 */
const char *ccp4_strerror( int);

/** Wrapper for ccp4_error which also calls exit(1)
 * @param message (const char *)
 * @return void
 */
void ccp4_fatal(const char *);

/** Function to set verbosity level for messages from
 *  ccp4_signal. Currently just off (0) and on (1).
 *  It should be generalised to be able to switch
 *  individual components on and off, i.e. replace 1 by
 *  a mask.
 *  cf. ccp4VerbosityLevel which sets the verbosity level
 *  for ccp4printf  These are separate as they may be used
 *  differently.
 * @param iverb If >= 0 then set the verbosity level to the
 *  value of iverb. If < 0 (by convention -1) then report
 *  current level. 
 * @return current verbosity level
 */
int ccp4_liberr_verbosity(int iverb);

/** Routine to set ccp4_errno and print out message for
 *  error tracing. This should be the only way in
 *  which ccp4_errno is set.
 *  See error codes above for levels and systems.
 *  A callback with prototype void function(void)
 *  may also be passed to the routine.
 *  Note: FATAL calls exit(1). 
 *  If ccp4_liberr_verbosity returns 0, then ccp4_signal sets
 *  ccp4_errno and returns without doing anything else.
 * @param error code (int)
 * @param message (const char * const) 
 * @param callback (point to routine void (*)(void) )
 * @return void
 */
void ccp4_signal(const int, const char *const, void (*)());

int cfile_perror(const char *);

#ifdef __cplusplus
}
}
#endif

#endif  /*!CCP4_ERROR_GUARD */
