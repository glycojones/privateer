/*
     library_err.c: Error handling library
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

/** @file library_err.c
 *  Error handling library.
 *  Charles Ballard
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ccp4_errno.h"
/* rcsid[] = "$Id$" */

/** @global ccp4_errno: global to store data 
*/
int ccp4_errno = 0;

/* error_levels: error level descriptions */
static const char * const error_levels[] =
  {
    "Success",                                   /* 0 */
    "Informational",                             /* 1 */
    "Warning",                                   /* 2 */
    "Error",                                     /* 3 */
    "FATAL ERROR"                                /* 4 */
};

/* file io errors */
static const char *const cfile_errlist[] =
  {
    "Error 0",                                   /* 0 = CIO_Ok */
    "Bad mode",                                  /* 1 = CIO_BadMode */
    "Cannot open file",                          /* 2 = CIO_CantOpenFile */
    "Too many open files",                       /* 3 = CIO_MaxFile */
    "Read failed",                               /* 4 = CIO_ReadFail */
    "Write failed",                              /* 5 = CIO_WriteFail */
    "Close fail",                                /* 6 = CIO_CloseFail */
    "Seek fail",                                 /* 7 = CIO_SeekFail */
    "Null pointer passed",                       /* 8 = CIO_NullPtr */
    "End of File",                               /* 9 = CIO_EOF */
    "No file"                                    /* 10 = CIO_NoFile */
    "File not open",                             /* 11 = CIO_NotOpen */
    "Unlink failed"                              /* 12 = CIO_UnlinkFail */
  };

/* map library errors */
static const char *const cmap_errlist[] =
  {
    "Error 0",                                   /* 0 = CMERR_Ok */
    "Unassigned unit",                           /* 1 = CMERR_NoChannel */
    "Unassigned unit or disposed file",          /* 2 = CMERR_NoFile */
    "Logical name does not exist",               /* 3 = CMERR_NoLogicalName */
    "Cannot open file",                          /* 4 = CMERR_CantOpenFile */
    "No associated header",                      /* 5 = CMERR_NoHeader */
    "Read failed",                               /* 6 = CMERR_ReadFail */
    "Write failed",                              /* 7 = CMERR_WriteFail */
    "Parameter or dimension is incorrect    ",   /* 8 = CMERR_ParamError */
    "Unrecognised keyword",                      /* 9 = CMERR_UnrecognK */
    "File stamp error",                          /* 10 = CMERR_FileStamp */
    "Symmetry",                                  /* 11 = CMERR_SymErr */
    "Cannot allocate memory",                    /* 12 = CMERR_AllocFail */
    "Too many open files",                       /* 13 = CMERR_MaxFile */
  };

/* mtz library errrors */
static const char *const cmtz_errlist[] =
  {
    "Error 0",                                   /* 0 = CMTZERR_Ok */
    "Unassigned unit",                           /* 1 = CMTZERR_NoChannel */
    "Null file handle, file not opened",         /* 2 = CMTZERR_NoFile */
    "Logical name does not exist",               /* 3 = CMTZERR_NoLogicalName */
    "Cannot open file",                          /* 4 = CMTZERR_CantOpenFile */
    "No associated header",                      /* 5 = CMTZERR_NoHeader */
    "Read failed",                               /* 6 = CMTZERR_ReadFail */
    "Write failed",                              /* 7 = CMTZERR_WriteFail */
    "Function parameter is incorrect",           /* 8 = CMTZERR_ParamError */
    "Invalid cell dimensions",                   /* 9 = CMTZERR_Cellerr */
    "File stamp error",                          /* 10 = CMTZERR_FileStamp */
    "Symmetry",                                  /* 11 = CMTZERR_SymErr */
    "Cannot allocate memory",                    /* 12 = CMTZERR_AllocFail */
    "Too many open files",                       /* 13 = CMTZERR_MaxFile */
    "Failed to initialise parser",               /* 14 = CMTZERR_ParserFail */
    "File not identified as MTZ",                /* 15 = CMTZERR_NotMTZ */
    "Missing or incomplete dataset information in input file.", /* 16 = CMTZERR_DatasetIncomplete */
    "No architecture information in file.",      /* 17 = CMTZERR_NoArch */
    "Attempt to access unallocated dataset",     /* 18 = CMTZERR_NullDataset */
    "Input MTZ file has incorrect major version for current library",      /* 19 = CMTZERR_BadVersion */
    "MTZ header is corrupted: missing tokens in SYMINF record", /* 20 = CMTZERR_SYMINFIncomplete */
    "MTZ header is corrupted: missing tokens in COLUMN record", /* 21 = CMTZERR_COLUMNIncomplete */
    "Batch headers corrupted",                   /* 22 = CMTZERR_BadBatchHeader */
    "Input MTZ file has different minor version to that supported by current library",      /* 23 = CMTZERR_DifferentVersion */
    "File column type different from type expected by program",      /* 24 = CMTZERR_ColTypeMismatch */
    "MTZ header: error in column group specification",   /* 25 = CMTZERR_ColGroupError */
    "MTZ header: error in column source specification",   /* 26 = CMTZERR_ColSourceError */
  };

/* parser library errors */
static const char *const cpars_errlist[] =
  {
    "Error 0",                                   /* 0 = CPARSERR_Ok */
    "Maximum number of tokens exceeded",         /* 1 = CPARSERR_MaxTokExceeded */
    "Cannot allocate memory",                    /* 2 = CPARSERR_AllocFail */
    "Null pointer",                              /* 3 = CPARSERR_NullPointer */
    "Line is longer than allocated length, so truncated", /* 4 = CPARSERR_LongLine */
    "Failed to open external command file",      /* 5 = CPARSERR_CantOpenFile */
    "Failed to get name for external file",      /* 6 = CPARSERR_NoName */
    "Overflow - exponent is too big to be evaluated",  /* 7 = CPARSERR_ExpOverflow */
    "Underflow - exponent is too small to be evaluated", /* 8 = CPARSERR_ExpUnderflow */
    "Problem in mat4_to_symop",                  /* 9 = CPARSERR_MatToSymop */
    "Failed to interpret symop string",          /* 10 = CPARSERR_SymopToMat */
  };

/* symmetry library errors */
static const char *const csym_errlist[] =
  {
    "Error 0",                                   /* 0 = CSYMERR_Ok */
    "Failed to initialise parser",               /* 1 = CSYMERR_ParserFail */
    "Cannot find SYMINFO file - no symmetry information", /* 2 = CSYMERR_NoSyminfoFile */
    "Pointer to spacegroup structure is NULL",   /* 3 = CSYMERR_NullSpacegroup */
    "ASU definition not found for this spacegroup", /* 4 = CSYMERR_NoAsuDefined */
    "Undefined Laue code for this spacegroup",   /* 5 = CSYMERR_NoLaueCodeDefined */
    "Not enough tokens on SYMINFO line",         /* 6 = CSYMERR_SyminfoTokensMissing */
  };

static const char *const cgen_errlist[] =
  {
    "Error 0",                                   /* 0 = CGENERR_Ok */
    "Cannot allocate memory",                    /* 1 = CGENERR_AllocFail */
    "Cannot set environment variable",           /* 2 = CGENERR_CantSetEnvironment */
    "Maximum number of logical names exceeded",  /* 3 = CGENERR_MaxNamesExceeded */
    "Use: -e filename",                          /* 4 = CGENERR_EOptionUseError */
    "Use: -d filename",                          /* 5 = CGENERR_DOptionUseError */
    "Use: <logical name> <file name>",           /* 6 = CGENERR_LogicalNameUseError */
    "Cannot open environ.def",                   /* 7 = CGENERR_CantOpenEnvFile */
    "Cannot open default.def",                   /* 8 = CGENERR_CantOpenDefFile */
    "Cannot parse environ.def file",             /* 9 = CGENERR_ParseEnvFail */
    "Cannot parse default.def file",             /* 10= CGENERR_ParseDefFail */
    "Cannot find input file",                    /* 11= CGENERR_CantFindInFile */
    "Failed to set path for environ.def file",   /* 12= CGENERR_EnvPathFail */
    "Failed to set path for default.def file",   /* 13= CGENERR_DefPathFail */
    "Cannot get CLIBD from environment",         /* 14= CGENERR_CantGetClibd */
    "Cannot get CCP4_SCR from environment",      /* 15= CGENERR_CantGetCcp4Scr */
  };

struct error_system {
  char system[32];
  int system_nerr;
  const char * const *error_list;
};

/* construct error list */
static const struct error_system ccp4_errlist[] = {
    {"system", 0, 0, },
    {"library_file", CCP4_COUNT(cfile_errlist), cfile_errlist,},
    {"mmdb", 0, 0,},
    {"mtz", CCP4_COUNT(cmtz_errlist), cmtz_errlist,},
    {"ccp4_map", CCP4_COUNT(cmap_errlist), cmap_errlist,},
    {"utils", 0, 0},
    {"ccp4_parser", CCP4_COUNT(cpars_errlist), cpars_errlist,},
    {"csym", CCP4_COUNT(csym_errlist), csym_errlist,},
    {"ccp4_general", CCP4_COUNT(cgen_errlist), cgen_errlist,}
};

static const int ccp4_system_nerr = CCP4_COUNT(ccp4_errlist);

/* Obtain character string based upon error code.
    Typical use ccp4_strerror(ccp4_errno)
    The returned string is statically allocated in the
    library_err.c file and should not be freed.
    param error code (int)
    returns const pointer to error message.
*/
const char *ccp4_strerror(int error)
{
  int system = CCP4_ERRGETSYS(error);
  /* int level = CCP4_ERRGETLEVEL(error); */
  int code = CCP4_ERRGETCODE(error);

  if (error == -1 || system == 0)
    return strerror(errno);

  if (system >= ccp4_system_nerr)
    return ("bad system error");
  if (code >= ccp4_errlist[system].system_nerr)
    return ("bad error code");
  return (ccp4_errlist[system].error_list[code]);
}

/* Print out passed message and internal message based upon
    ccp4_errno
    "message : error message "
    param message (const char *)
    return void      
*/
void ccp4_error (const char *msg)
{
  const char *colon;

  if (msg == 0 || *msg == '\0')
    colon = "";
  else
    colon = ": ";
 
  fprintf (stderr, "%s%s%s\n",
           msg, colon, ccp4_strerror(ccp4_errno));
  if(ccp4_errno != -1 && CCP4_ERRGETSYS(ccp4_errno)) {
    fprintf (stderr, "System: %s\nLevel: %d\n", 
             ccp4_errlist[CCP4_ERRGETSYS(ccp4_errno)].system,
             CCP4_ERRGETLEVEL(ccp4_errno));
    if (errno)
      fprintf (stderr, "%s%s\n",
               "Last system message: ",strerror(errno)); }
}

/* Wrapper for ccp4_error which also calls exit(1)
   param message (const char *)
*/
void ccp4_fatal (const char *message)
{
  ccp4_error(message);
  exit(1);
}

int CFile_Perror(const char *msg)
{
  const char * colon;
  int error = CCP4_ERRGETCODE(ccp4_errno);
  int cfile_nerr = ccp4_errlist[4].system_nerr;

  if (msg == NULL || msg == '\0') colon = "";
  else colon = ": ";

  if (error > 0 && error <= cfile_nerr) {
    fprintf(stderr,"%s%s%s \n",
            msg,colon,cfile_errlist[error]);
    return error; }
  fprintf(stderr,"Unknown error code");
  return -1;
}

int ccp4_liberr_verbosity(int iverb) {
  static int verbosity_level=1;

  if (iverb >= 0)
    verbosity_level = iverb;

  return verbosity_level;
}

/* Routine to set ccp4_errno and print out message for
    error tracing. This should be the only way in
    which ccp4_errno is set.
    See error codes above for levels and systems.
    A callback with prototype void function(void)
    may also be passed to the routine.
    Note: FATAL calls exit(1).
    param error code (int)
    param message (const char * const)
    param callback (point to routine)
*/
void ccp4_signal(const int code, const char * const msg, 
		 void (*callback) ())
{
  int severity = CCP4_ERRGETLEVEL(code),
    system     = CCP4_ERRGETSYS(code),
    msg_no     = CCP4_ERRGETCODE(code),
    fatal_err  = (severity == 4);
  static const char msg_fmt[] = ">>>>>> CCP4 library signal %s:%s (%s)\n\t raised in %s <<<<<<\n";
   static const char sys_fmt[] = ">>>>>> System signal %d:%s (%s)\n\t raised in %s <<<<<<\n";
  ccp4_errno = code;

  /* use this function to control whether error messages are printed */
  if (!ccp4_liberr_verbosity(-1)) return;

  if (system == 0) {
    if (msg) 
      printf(sys_fmt,
	     errno,
	     strerror(errno),
	     error_levels[severity],
	     msg);
    else 
      printf(">>>>>> System signal %d:%s (%s) <<<<<<", 
	     errno, 
	     strerror(errno), 
	     error_levels[severity]);
    ccp4_errno = errno; }
  else 
    if (msg) 
      printf(msg_fmt,
	     ccp4_errlist[system].system,
	     ccp4_errlist[system].error_list[msg_no],
	     error_levels[severity],
	     msg);
    else
      printf(">>>>>> CCP4 library signal %s:%s (%s) <<<<<<\n",
	     ccp4_errlist[system].system,
	     ccp4_errlist[system].error_list[msg_no],
	     error_levels[severity]);
  
  if (callback)
    (*callback)();

  if (fatal_err) exit(1);
}


