/*
     ccp4_program.c: Utilies to set and fetch program information.
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

/** @file ccp4_program.c
 *  Utilies to set and fetch program information.
 *  Peter Briggs CCP4 May 2001
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ccp4_program.h"
#include "ccp4_parser.h"
#include "ccp4_utils.h"
#include "ccp4_general.h"
/* rcsid[] = "$Id$" */

static char ccp4version[MAXLEN_PROGVERSION];

char *ccp4_prog_vers(const char *progvers) 
{
  static char programversion[MAXLEN_PROGVERSION]="";

  if (progvers) {
    strncpy(programversion, progvers, MAXLEN_PROGVERSION);
    programversion[MAXLEN_PROGVERSION-1] = '\0';
  }
  return programversion;
}

char *ccp4_vers_no(void)
{
  static int init=0;

  char *filepath=NULL, *filename=NULL;
  char *vfile="/lib/ccp4/MAJOR_MINOR";
  FILE *cfile;
  int i;

  if (!init) {
  strcpy(ccp4version,CCP4_VERSION_NO);

  filepath = (char *) getenv("CCP4");
  if (filepath) {
      filename = (char *) ccp4_utils_malloc(sizeof(char)*(strlen(filepath)+strlen(vfile))+1);
      strcpy(filename,filepath);
      strcat(filename,vfile);
      if (ccpexists(filename)) {
         cfile=fopen(filename,"r");
         if (cfile) {
           fgets(ccp4version,MAXLEN_PROGVERSION,cfile);
           i = strlen(ccp4version)-1;
           while (isspace(ccp4version[i]) ) {
             ccp4version[i--]='\0';
           }
         }
      }
      /* Make sure that we clean up */
      if (filename) free(filename);
    }
    init=1;
  }
  return ccp4version;
}
/*------------------------------------------------------------------*/

/* ccp4ProgramName

   Set or return program name

   Always returns a pointer to the program name
   If progname is not NULL then set the program name to
   progname.

   NB Default program name will be returned as "CCP4",
   until reset by the calling subprogram.
*/
char *ccp4ProgramName(const char *progname)
{
  static char programname[MAXLEN_PROGNAME]="CCP4";
  int         i;

  if (progname) {
    i = 0;
    while (progname[i] != '\0' && i < MAXLEN_PROGNAME) {
      programname[i] = progname[i];
      ++i;
    }
    if (i == MAXLEN_PROGNAME) {
      programname[MAXLEN_PROGNAME-1] = '\0';
    } else {
      programname[i] = '\0';
    }
  }
  return programname;
}

/* ccp4_prog_info

   Print program info for -i option.
 */
void ccp4_prog_info(void)
{
    printf("CCP4 software suite: library version %s\n",ccp4_vers_no());
    printf("CCP4 software suite: patch level     %s\n",ccp4_vers_no());
    printf("Program:             %s",ccp4ProgramName(NULL));
    if (ccp4_prog_vers(NULL) && strlen(ccp4_prog_vers(NULL))) 
      printf("; version %s",ccp4_prog_vers(NULL));
    printf("\n");
}

/* ccp4RCSDate

   Set or return program RCS date

   If the input string is not a NULL pointer then
   it is assumed to be an RCS string
   This is processed to extract a date string in
   the form "DD/MM/YY" (day/month/year), which is
   then stored.

   ccp4RCSDate always returns the currently
   stored date string.
*/
char *ccp4RCSDate(const char *rcs_string)
{
  static char RCSDate[MAXLEN_RCSDATE]="";
  char        tmpstr1[8],tmpstr2[3];

  /* Deconstruct the RCS string passed to this
     function */
  if (rcs_string) {
    /* Extract useful data from RCS string for examination */
    strncpy(tmpstr1,rcs_string,7);
    tmpstr1[7] = '\0';
    strncpy(tmpstr2,rcs_string,2);
    tmpstr2[2] = '\0';
    if (strncmp(tmpstr1,"$Date: ",7) == 0) {
      /* Raw form of RCS string (not exported) i.e.:
	 "$Date$"
      */
      /* Build the date string in the form DD/MM/YY */
      strncpy(RCSDate,rcs_string+15,2);
      strncat(RCSDate,"/",1);
      strncat(RCSDate,rcs_string+12,2);
      strncat(RCSDate,"/",1);
      strncat(RCSDate,rcs_string+9,2);
    } else if (strlen(rcs_string) > 10 &&
	       (strncmp(tmpstr2,"19",2) == 0 || strncmp(tmpstr2,"20",2)) ) {
      /* RCS string after export i.e.:
	 "2003/05/14 11:45:13 ..." */
      /* Build the date string in the form DD/MM/YY */
      strncpy(RCSDate,rcs_string+8,2);
      strncat(RCSDate,"/",1);
      strncat(RCSDate,rcs_string+5,2);
      strncat(RCSDate,"/",1);
      strncat(RCSDate,rcs_string+2,2);
    } else {
      /* Fallback */
      strncpy(RCSDate,"",1);
    }
  }
  /* Always return the stored date */
  return RCSDate;
}

/* ccp4ProgramTime

   Set or print program time information
*/
void ccp4ProgramTime(int init)
{
  static int elaps0=0;
  static float tarray0[2];
  int elaps;
  float tarray[2];

  if (init || !elaps0 ) {
    elaps0 = time(NULL);
    ccp4_utils_etime(tarray0);
  } else {
    elaps = time(NULL) - elaps0;
    ccp4_utils_etime(tarray);

    printf("Times: User: %9.1fs System: %6.1fs Elapsed: %5d:%2.2d  \n",
	   tarray[0]-tarray0[0],tarray[1]-tarray0[1],elaps/60,elaps%60);
  }

}

/* ccp4VerbosityLevel

   Set or return the reference verbosity level

   Always return the verbosity level - if verboselevel is
   between 0 and 9 then reset the verbosity level to
   verboselevel
*/
int ccp4VerbosityLevel(int level)
{
  /* The default level is 1 */
  static int verbositylevel=1;

  if (level > -1 && level < 10)
    verbositylevel = level;
  return verbositylevel;
}

/*------------------------------------------------------------------*/

/*
  Callback functionality
 */

/* ccp4Callback

   Set or invoke a user-defined callback function.

   Internal function: applications should use the API functions
   ccp4SetCallback and ccp4InvokeCallback
 */
int ccp4Callback(CCP4INTFUNCPTR mycallback, char *mode, int ierr,
                 const char *message)
{
  static CCP4INTFUNCPTR callback=ccp4NullCallback;

  if (strncmp(mode,"set",3) == 0) {
    /* Set callback
       Store the pointer to the callback function */
    callback=mycallback;
    return 1;
  } else if (strncmp(mode,"invoke",3) == 0) {
    /* Invoke callback
       Execute the callback function */
    return callback(ierr,message);
  }
  /* Unrecognised mode */
  return 0;
}

/* ccp4SetCallback

   Store a pointer to a user-defined callback function of
   the form "int func(int, char *)"

   This is a wrapper to ccp4Callback in "set" mode.
 */
int ccp4SetCallback(CCP4INTFUNCPTR mycallback)
{
  return ccp4Callback(mycallback,"set",-1,"No message");
}

/* ccp4InvokeCallback

   Execute the user-defined callback function (previously
   set up using ccp4SetCallback) with the supplied
   arguments.

   This is a wrapper to ccp4Callback in "invoke" mode.
 */
int ccp4InvokeCallback(int ierr, const char *message)
{
  return ccp4Callback(ccp4NullCallback,"invoke",ierr,message);
}

/* Default null callback function

   Internal function: this is the default callback function
   used by ccp4Callback if no user-defined function has been
   specified.
 */
int ccp4NullCallback(int level, const char *message)
{
  /* This is the default callback function which takes no
     action */
  return 1;
}

/*------------------------------------------------------------------*/

/* check existence of licence agreement */

int ccp4_licence_exists(const char *name)
{
  int sue=1,lpath;
  char *filepath=NULL,*filename=NULL,tmp_string[20];

  strtoupper(tmp_string,name);
  if (strmatch(tmp_string,"CCP4")) {
    filepath = (char *) getenv("CCP4");
    if (filepath) {
      lpath = strlen(filepath);
      filename = (char *) ccp4_utils_malloc(sizeof(char)*(lpath+12));
      strcpy(filename,filepath);
      strcpy(filename+lpath,"/.agree2ccp4");
      if (ccpexists(filename)) sue = 0;
      /* Make sure that we clean up */
      if (filename) free(filename);
    }
    if (sue == 1) {
      filepath = (char *) getenv("HOME");
      if (filepath) {
        lpath = strlen(filepath);
        filename = (char *) ccp4_utils_malloc(sizeof(char)*(lpath+12));
        strcpy(filename,filepath);
        strcpy(filename+lpath,"/.agree2ccp4");
        if (ccpexists(filename)) sue = 0;
	/* Make sure that we clean up */
	if (filename) free(filename);
      }
    }
    if (sue == 1) {
      ccperror(1,"Cannot find required license agreements!");
      return 0;
    }
  }
  return 1;
}

/* html_log_output and summary_output currently only used by ccperror to
   tidy up Fortran program output. Defaults are 0 for C programs. */
int html_log_output(int ihtml_in) {
  static int ihtml=0;

  if (ihtml_in >= 0)
    ihtml = ihtml_in;
  return ihtml;
}

int summary_output(int isumm_in) {
  static int isumm=0;

  if (isumm_in >= 0)
    isumm = isumm_in;
  return isumm;
}
