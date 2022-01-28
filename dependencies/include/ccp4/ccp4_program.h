/*
     ccp4_program.h: Headers to utilies to set and fetch program information.
     Copyright (C) 2001  CCLRC, Peter Briggs

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


/** @file ccp4_program.h
 *  Utilies to set and fetch program information.
 *  Peter Briggs CCP4 May 2001
 */

/*------------------------------------------------------------------*/

/* Macro definitions */

/*------------------------------------------------------------------*/

#ifndef __CCP4Program__
#define __CCP4Program__

/* rcsidhp[] = "$Id$" */

#ifdef  __cplusplus
namespace CCP4 {
extern "C" {
#endif

#define CCP4_VERSION_NO "7.1"
#define CCP4_PATCH_LEVEL "7.1.0"

/* Maximum lengths of strings holding program names and versions */
#define MAXLEN_PROGNAME    80
#define MAXLEN_PROGVERSION 80
#define MAXLEN_RCSDATE     80

/*------------------------------------------------------------------*/

/* Type Definitions */

/*------------------------------------------------------------------*/

/* Define a type which is a pointer to a function taking an integer
   and a pointer to character, and returning an integer */
typedef int (*CCP4INTFUNCPTR)(int, const char *);

/*------------------------------------------------------------------*/

/* Function Prototypes */

/*------------------------------------------------------------------*/

/** Register or query program version.
 * @param progvers Program version string, or NULL to query existing value.
 * @return Program version string.
 */
char *ccp4_prog_vers(const char *progvers);

/** Query ccp4 version.
 * @return CCP4 version string.
 */
char *ccp4_vers_no(void);

/** Set or return program name.
 * @param progname Program name, or NULL to query existing value.
 * @return Program name
 */
char *ccp4ProgramName(const char *progname);

/** Print program info for -i option.
 */ 
void ccp4_prog_info(void);

/** Set or return program RCS date
 * @param rcs_string Date string, or NULL to query existing value.
 * @return Date string
 */
char *ccp4RCSDate(const char *rcs_string);

/** Set or print program time information
 * @param init
 */
void ccp4ProgramTime(int init);

/** Set or return the reference verbosity level
 * Always return the verbosity level - if verboselevel is
 * between 0 and 9 then reset the verbosity level to
 * verboselevel
 * @param level Verbosity level, or -1 to query existing value.
 * @return Verbosity level
 */
int ccp4VerbosityLevel(int level);

/** Set or invoke a user-defined callback function
 * The callback must be of the form "function(const int, const char *)"
 * This is essentially an internal function which operates in one of two
 * modes - in "set" mode the named function is stored and the remaining
 * arguments are discarded; in "invoke" mode the stored function is
 * executed with the supplied values (the supplied name is discarded).
 * @param mycallback Callback function (discarded in "invoke" mode)
 * @param mode Either "set" or "invoke"
 * @param ierr An error level equivalent to that used in ccperror
 * @param message A message string equivalent to that used in ccperror
 * @return Result of the executed function (invoke mode)
 */
int ccp4Callback(CCP4INTFUNCPTR mycallback, char *mode, int ierr,
                 const char *message);

/** Set a user-defined callback function
 * This is a wrapper to ccp4Callback - it stores a user-defined
 * callback function which must be of the form 
 * "function(const int, const char *)"
 * @param mycallback Callback function
 * @return 1 (if the function is stored), 0 (if it is not) 
 */
int ccp4SetCallback(CCP4INTFUNCPTR mycallback);

/** Invoke the user-defined callback function
 * This is a wrapper to ccp4Callback - it executes the user-defined
 * callback function previously stored.
 * @param ierr An error level equivalent to that used in ccperror
 * @param message A message string equivalent to that used in ccperror
 * @return Result of the executed function
 */
int ccp4InvokeCallback(int ierr, const char *message);

/** A dummy callback function used by default in ccp4CallOnExit
 * Internal function. This function does nothing.
 * @param level Severity level supplied from ccperror
 * @param message Message text supplied from ccperror
 * @return Always returns 1
*/      
int ccp4NullCallback(int level, const char *message);

/** Check existence of licence agreement
 * @param name Name of licence, e.g. "CCP4".
 * @return 1 for licence exists, else 0.
 */
int ccp4_licence_exists(const char *name);

/** Register or query html output level.
 * @param ihtml_in 0 = turn off html output, 1 = turn on html output, -1 = query existing value
 * @return 0 = no html output, 1 = html output
 */
int html_log_output(int ihtml_in);

/** Register or query summary output level.
 * @param isumm_in 0 = turn off summary output, 1 = turn on summary output, -1 = query existing value
 * @return 0 = no summary output, 1 = summary output
 */
int summary_output(int isumm_in);

#ifdef __cplusplus
} 
} 
#endif

#endif   /* __CCP4Program__ */
