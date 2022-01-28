/*
     ccp4_general.c: General library functions and utilities.
     Copyright (C) 2001  CCLRC, Peter Briggs et al

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

/** @file ccp4_general.c
 *  General library functions and utilities.
 *  Peter Briggs et al
 */

/*   ccp4_general.c
     Peter Briggs CCP4 May 2001/Feb 2003

     General library functions and utilities

     ccperror(level,message)
     Error reporting and program termination

     ccp4printf(level,format,printargs)
     A wrapper for vprintf - acts like printf with an additional
     argument which is checked against the reference verbosity
     level.

     ccp4fyp(argc,argv)
     Initialise environment for CCP4 programs and parse the
     command line arguments

     ccp4setenv(logical_name,value,envname,envtype,envext,ienv,
     no_overwrt)
     Set up file names and associate with logical_name in
     environment

     ccpexists(filename)
     Check if file exists and can be opened for reading

     ccpputenv(logical_name,file_name)
     Set environment variable logical_name to have value file_name
*/

/*------------------------------------------------------------------*/

/* CCP4 library clones */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

#ifdef HAVE_CONFIG_H
#include <config.h> /* PACKAGE_ROOT */
#endif
/* Library header files */
#include "ccp4_fortran.h"
#include "ccp4_utils.h"
#include "ccp4_parser.h"
#include "ccp4_general.h"
#include "ccp4_program.h"
#include "cmtzlib.h"
#include "csymlib.h"
#include "ccp4_errno.h"
/* rcsid[] = "$Id$" */

/*------------------------------------------------------------------*/

/* ccperror

   Error reporting and program termination

   ccperror is called with a message level ierr and a message string
   The exact behaviour is determined by the value of ierr:

   0 : normal termination and stop
   1 : fatal error and stop
   2 : report severe warning
   3 : report information
   4 : report from library

   ierr=-1 also reports last system error and terminates.

   The text in message is sent to standard output, and to standard
   error for ierr=1 (or -1).

*/
int ccperror(int ierr, const char *message)
{
  /* Execute user-defined function callback */
  ccp4InvokeCallback(ierr,message);
  /* Do messaging */
  ccperror_noexit(ierr, message);
  /* Exit if necessary */
  if (ierr==0) {
    exit(0);
  } else if (ierr==1 || ierr==-1) {
    exit(1);
  }
  return 0;
}

/* ccperror_noexit

   As above, but doesn't call exit()

*/

int ccperror_noexit(int ierr, const char *message)
{
  char *prog_name=NULL;

  /* Get the program name */
  prog_name = ccp4ProgramName(NULL);
  if (!prog_name) prog_name = strdup("CCP4");

  if (ierr==0) {
    /* Level 0 : normal termination */
    ccp4printf(0," %s:  %s\n",prog_name,message);
    /* Get the amount of time elapsed since start of
       program. Initialised by ccp4fyp */
    ccp4ProgramTime(0);

    /* closing tags - these are simplified w.r.t. Fortranic 
       versions for this special case */
    if (html_log_output(-1)) {
      printf("</pre>\n");
      printf("</html>\n");
    }
    if (summary_output(-1)) {
      if (html_log_output(-1)) {
        printf("<!--SUMMARY_END--></FONT></B>\n");
      } else {
        printf("<!--SUMMARY_END-->\n");
      }
    }

  } else if (ierr==1 || ierr==-1) {
    /* Level 1 (-1) : fatal error */
    /* If ierr=-1 then print last system error
       N.B. Use of perror in this context is untested by me */
    if (ierr < 0) perror("Last system error message");
    /* Send the message to the standard error */
    fprintf(stderr," %s:  %s\n",prog_name,message);
    /* Also to the standard out */
    ccp4printf(0," %s:  %s\n",prog_name,message);
    /* Get the amount of time elapsed since start of
       program. Initialised by ccp4fyp */
    ccp4ProgramTime(0);

    /* closing tags - these are simplified w.r.t. Fortranic 
       versions for this special case */
    if (html_log_output(-1)) {
      printf("</pre>\n");
      printf("</html>\n");
    }
    if (summary_output(-1)) {
      if (html_log_output(-1)) {
        printf("<!--SUMMARY_END--></FONT></B>\n");
      } else {
        printf("<!--SUMMARY_END-->\n");
      }
    }

  } else if (ierr==2) {
    /* Level 2 : severe warning */
    ccp4printf(0," \n $TEXT:Warning: $$ comment $$ \n WARNING: %s\n $$\n",
	   message);

  } else if (ierr>2)  {
    /* Levels higher than 2 : report information */
    ccp4printf(0,"%s\n",message);
  }
  return 0;
}

/*------------------------------------------------------------------*/

/* ccp4printf

   This is a wrapper for vprintf

   If the supplied message is less than or equal to the reference
   verbosity level then the format string and remaining arguments
   are passed to vprintf to be printed on stdout.
   Otherwise nothing is printed.

   The format string has the same format as the format strings
   passed to printf, with the remaining arguments being the values
   (if any) to be inserted into placeholders in the format string.

   Returns the number of bytes printed, or zero if no printing was
   performed, or a negative number for an error from vprintf.
*/
int ccp4printf(int level, char *format, ...)
{
  int     nbytes=0;
  va_list args;

  /* Check the level against the refence verbosity level */
  if (level <= ccp4VerbosityLevel(-1)) { 
    /* Use the vprint function to print */
    va_start(args,format);
    nbytes = vprintf(format,args);
    va_end(args);
  }
  /* Return to calling function */
  return nbytes;
}

/*------------------------------------------------------------------*/

/* return parent of directory that contains arg, e.g.
 * /foo/bin/prog        -> /foo    (actually /foo/bin/..)
 * C:\CCP4\bin\prog.exe -> C:\CCP4 (actually C:\CCP4\bin\..)
 * foo/prog             -> .       (actually foo/..)
 * prog                 -> ..
 * ../prog              -> ../..
 */
static const char* get_parent_directory(const char* arg)
{
    const char *last_sep = strrchr(arg, PATH_SEPARATOR);
    const char *basename = last_sep != NULL ? last_sep + 1 : arg;
    int dir_len = basename - arg;
    char* parent_dir = (char*) ccp4_utils_malloc(dir_len + 3);
    strncpy(parent_dir, arg, dir_len);
    strcpy(parent_dir+dir_len, "..");
    return parent_dir;
}

/*------------------------------------------------------------------*/

/* ccp4fyp

   Initialise environment for CCP4 programs and parse the
   command line arguments

   Background:
   This function processes the command line arguments supplied via
   the argv array, and performs initialisations of the CCP4
   environment within the program based on these arguments.

   CCP4 programs which use ccp4fyp should be called with the following
   syntax:

   <progname> [switches] <logical_name1> <value1> <logical_name2> <value2>...

   The command switches all start with a hyphen (-) and immediately
   follow the program name. The remainer of the arguments are logical
   name-value pairs.

   On successful completion ccp4fyp returns 0 (zero).
   On encountering an error, ccp4fyp registers the error condition via
   ccp4_signal and returns a non-zero value.
*/
int ccp4fyp(int argc, char **argv)
{
  /* Diagnostics/debug output */
  int diag=0;

  int  i,iarg=1,arg_end=0,ienv=0;
  char *envname[CCP4_MAXNAMES],*envtype[CCP4_MAXNAMES],*envext[CCP4_MAXNAMES];

  /* Flags for processing command line switches */
  int  info=0,ihelp=1,idefault=0,ienviron=0,nohtml=0,nosummary=0;
  char *testarg=NULL;

  /* Filenames, directories etc */
#ifdef PACKAGE_ROOT
  const char *pkgroot=PACKAGE_ROOT;
#else
  const char *pkgroot=get_parent_directory(argv[0]);
#endif
  char *basename=NULL,*cinclude=NULL,*home=NULL;
  char *dir=NULL,*std_dir=NULL,*tmpstr=NULL;

  /* Decoding environ/defaults files */
  char line[CCP4_MAXLINE];
  char *logical_name=NULL,*file_name=NULL,*file_type=NULL,*file_ext=NULL;
  CCP4PARSERARRAY *parser=NULL;

  /* Environ.def file */
  int  env_init=1;
  char *env_file=NULL;
  char *env_logical_name=NULL,*env_file_type=NULL,*env_file_ext=NULL;
  FILE *envfp;

  /* Default.def file */
  int  def_init=1;
  char *def_file=NULL;
  char *def_logical_name=NULL,*def_file_name=NULL;
  FILE *deffp;

  /* Used for error messages from ccp4setenv */
  int  ierr;

  /* Begin */
  if (diag) printf("CCP4FYP: starting\n");

  if (diag) 
   for (i = 0; i < argc; ++i) 
      printf("ccp4fyp: comand line argument %d %s\n",i,argv[i]);

  /* ------------------------------------------------------ */
  /* Initialise program name and timing information */
  /* ------------------------------------------------------ */
  ccp4ProgramTime(1);

  /*ccp4ProgramName(ccp4_utils_basename(argv[0])); */
  basename = ccp4_utils_basename(argv[0]);
  ccp4ProgramName(basename);
  free(basename);
  basename = NULL;

  /* ------------------------------------------------------ */
  /* Process command line option switches */
  /* ------------------------------------------------------ */
  /* Nb ignore the first argument (iarg=0),
     because this is just the executable name */
  while (iarg < argc && !arg_end) {
    if (diag) printf("CCP4FYP: command line argument %d = \"%s\"\n",iarg,argv[iarg]);

    if (argv[iarg][0] == '-') {
      /* An argument of the form -option */
      if (diag) printf("--> This is an option switch\n");
      
      /* Remove the leading hyphen */
      testarg = (char *) ccp4_utils_realloc(testarg,(strlen(argv[iarg])+1)*sizeof(char));
      strtoupper(testarg,&(argv[iarg][1]));
      if (diag) printf("--> Truncated and uppercased it is now \"%s\"\n",testarg);

      /* Test for each possible option:
	 -v(erbose) 0-9  verbose output level
	 -h(elp) 0-9     (alias for -v)
	 -n              don't read environ.def and default.def
	 -d <filename>   use <filename> instead of default.def
	 -e <filename>   use <filename> instead of environ.def
	 -i              print CCP4 library version, program name and
	                 program version to standard output, and exit. 
	 -nohtml         suppress printing of html tags
	 -nosummary      suppress printing of summary tags
      */
      if (testarg[0] == 'V' || testarg[0] == 'H') {
	/* -v|h <0-9>: Verbose option */
	if (diag) printf("--> Identified -v option:");
	/* Set ihelp to point to the argument which should
	   hold the verbosity level, and process later */
	ihelp = ++iarg;

      } else if (testarg[0] == 'N') {
	/* -n, -nohtml, -nosummary */
	if (strlen(testarg) == 1) {
	  /* -n: Don't read environ.def and default.def */
	  if (diag) printf("--> Identified -n option\n");
	  def_init = 0;
	  env_init = 0;
	} else if (strncmp("NOHTML",testarg,3)==0) {
	  /* -nohtml: Don't write html tags into logfile */
	  nohtml = -1;
	} else if (strncmp("NOSUMMARY",testarg,3)==0) {
	  /* -nosummary: Don't write summary tags into logfile */
	  nosummary = -1;
	}
      } else if (testarg[0] == 'D') {
	/* -d <filename>: Use non-default default.def file */
	if (diag) printf("--> Identified -d option:");
	/* Set idefault to point to the argument which should
	   hold the default.def file, and process later */
	idefault = ++iarg;
      } else if (testarg[0] == 'E') {
	/* -e <filename>: Use non-default environ.def file */
	if (diag) printf("--> Identified -e option:");
	/* Set ienviron to point to the argument which should
	   hold the environ.def file, and process later */
	ienviron = ++iarg;
      } else if (testarg[0] == 'I') {
	/* -i(nfo): Info option */
	if (diag) printf("--> Identified -i(nfo) option\n");
	info = 1;

	/* Unrecognised switch */
      } else {
	ccp4printf(1,"Ignoring unrecognised switch \"%s\"\n",argv[iarg]);
      }

      /* Next argument */
      iarg++;

    } else {
      /* Found an argument which doesn't start with a "-"
	 This is the end of the command line switches */
      if (diag) printf("CCP4FYP: end of command line switches\n");
      arg_end = 1;
    }
  }
  /* Release the memory associated with the temporary argument pointer */
  if (testarg) {
    free(testarg);
    testarg = NULL;
  }

  /* ------------------------------------------------------ */
  /* Finished processing command line options */
  /* ------------------------------------------------------ */

  /* At this point we may need to perform some actions based
     on the switches supplied by the user */

  /* Program information requested by -i(nfo) option */
  if (info) {
    /* Print program information and stop */
    ccp4_prog_info();
    exit(0);
  }

  /* Initialise debug (verbosity) level
     Level 0 switches off all output (silent mode)
     Level 9 prints everything
     Level 1 is normal I guess
     It is not clear from documentation precisely what is/is not
     output for the levels between 0 and 9 */
  if (ihelp > 0) {
    /* Extract the level from the argument list - ihelp
       points to which argument should hold it */
    if (ihelp < argc) {
      if (strlen(argv[ihelp]) == 1 && isdigit(argv[ihelp][0])) {
	ihelp = atoi(argv[ihelp]);
      } else {
	ihelp = 1;
      }
      if (diag) printf("Verbose level %d\n",ihelp);
    } else {
      ihelp = 1;
      if (diag) puts("No verbose level - defaulting to 1.");
    }
  }
  ccp4VerbosityLevel(ihelp);

  /* ------------------------------------------------------ */
  /* Initialise HTML and SUMMARY tags */
  /* ------------------------------------------------------ */

  /* Initialise html, summary tags by setting environment variables
     for Fortran html library to pick up. No direct call, as
     C programs may not link to Fortran library. */
  if (nohtml) ccpputenv("CCP_SUPPRESS_HTML","1");
  if (nosummary) ccpputenv("CCP_SUPPRESS_SUMMARY","1");

  /* ------------------------------------------------------ */
  /* Get useful directories (CINCL and HOME) */
  /* ------------------------------------------------------ */

  /* Get value of CINCL variable -
     Note that (1) getenv returns a null pointer if no name is
     found, (2) we should not free the resulting address */
  cinclude = (char *) getenv("CINCL");
  if (!cinclude) {
    if (diag) printf("--> CINCL env var has no value assigned\n");
  } else {
    if (diag) printf("--> CINCL is \"%s\"\n",cinclude);
  }

  /* Get value of HOME variable (same notes apply as for CINCL) */
  home = (char *) getenv("HOME");
  if (!home) {
    if (diag) printf("--> HOME env var has no value assigned\n");
  } else {
    if (diag) printf("--> HOME is \"%s\"\n",home);
  }    

  /* ------------------------------------------------------ */
  /* Environ.def file */
  /* ------------------------------------------------------ */

  /* ------------------------------------------------------ */
  /* Use non-standard environ.def file */
  /* ------------------------------------------------------ */
  /* The search logic is:
     1. If a directory is specified as part of the filename, then
	use the filename as is, otherwise
     2. if the HOME variable is defined, then use $HOME/filename,
        otherwise
     3. use the filename as is (in current directory).
  */
  if (ienviron > 0) {
    /* Non-standard environ.def was specified */
    if (ienviron < argc) {
      if (env_file) free(env_file);
      env_file = (char *) ccp4_utils_malloc(sizeof(char)*(strlen(argv[ienviron])+1));
      if (!env_file) {
	/* Couldn't allocate memory to store filename 
	   Do clean up and exit */
	ccp4fyp_cleanup(ienv,envname,envtype,envext,logical_name,file_name,
			file_type,file_ext,env_file,def_file,dir,parser);
	ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_AllocFail),"ccp4fyp",NULL);
	return 1;
      }
      strcpy(env_file,argv[ienviron]);
      if (diag) printf("CCP4FYP: env file is \"%s\"\n",env_file);
      env_init = 1;
      /* Check whether this includes the path */
      if (dir) free(dir);
      dir = ccp4_utils_pathname(env_file);
      if (dir && dir[0] == '\0') {
	/* No path name - try and use $HOME/file_name */
	if (diag) puts("CCP4FYP: env file has no path.");
	if (home) {
	  tmpstr = ccp4_utils_joinfilenames(home,env_file);
	  if (diag) printf("CCP4FYP: HOME exists, joined filename \"%s\"\n",tmpstr);
	  if (tmpstr) {
	    if (env_file) free(env_file);
	    env_file = tmpstr;
	  } else {
	    /* Failed to make complete filename
	       Do clean up and exit */
	    ccp4fyp_cleanup(ienv,envname,envtype,envext,logical_name,file_name,
			    file_type,file_ext,env_file,def_file,dir,parser);
	    ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_EnvPathFail),"ccp4fyp",NULL);
	    return 1;
	  }
	}
      }
      if (diag) printf(" environ.def file is \"%s\"\n",env_file);
    } else {
      /* Not enough arguments in the arg list
	 Do clean up and exit */
      if (diag) printf(" no filename found\n");
      ccp4fyp_cleanup(ienv,envname,envtype,envext,logical_name,file_name,
		      file_type,file_ext,env_file,def_file,dir,parser);
      ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_EOptionUseError),"ccp4fyp",NULL);
      return 1;
    }
  }

  /* ------------------------------------------------------ */
  /* Use standard environ.def file */
  /* ------------------------------------------------------ */
  /* The search logic is:
     1. If the CINCL variable is defined, then use $CINCL/filename,
	otherwise
     2. if the HOME variable is defined, then use $HOME/filename,
        otherwise
     3. use the filename as is (in current directory).
  */
  if (!env_file) {
    /* Use the standard environ.def file in CINCL */
    if (diag) printf("--> use standard environ.def file\n");
    std_dir = NULL;
    if (cinclude) {
      std_dir = cinclude;
    } else if (home) {
      std_dir = home;
    }
    /* Set the full path for the environ.def file */
    if (env_file) free(env_file);
    if (std_dir) {
      if (diag) printf("--> leading directory is \"%s\"\n",std_dir);
      env_file = ccp4_utils_joinfilenames(std_dir,"environ.def");
    } else {
      env_file = ccp4_utils_joinfilenames(".","environ.def");
    }
    if (!env_file) {
      /* Failed to make full path name for environ.def
	 Do clean up and exit */
      ccp4fyp_cleanup(ienv,envname,envtype,envext,logical_name,file_name,
		      file_type,file_ext,env_file,def_file,dir,parser);
      ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_EnvPathFail),"ccp4fyp",NULL);
      return 1;
    }
    if (diag) printf("--> Full path for environ.def is \"%s\"\n",env_file);
  }

  /* ------------------------------------------------------ */
  /* Read in environ.def */
  /* ------------------------------------------------------ */
  /* environ.def contains lines of the form
     LOGICALNAME=type.ext # comments
     where type is "in", "out" or "inout"
     and ext is the default extension, e.g. "mtz" or "scr"
  */
  if (env_init && env_file) {
    if (diag) printf("CCP4FYP: reading environ.def file\n");
    if (diag) printf("--> Full path for environ.def is \"%s\"\n",env_file);
    /* Open the environ.def file as read-only*/
    ccp4printf(2,"Opening file \"%s\"\n",env_file);
    envfp = fopen(env_file,"r");
    /* did not find the environ.def file in std_dir
       so try "PACKAGE_ROOT" */
    if (!envfp) {
      free(env_file);
      env_file = ccp4_utils_joinfilenames(ccp4_utils_joinfilenames(
                  ccp4_utils_joinfilenames(pkgroot,"share"),"ccp4"),"environ.def");
      if (diag) printf("CCP4FYP: reading environ.def file\n");
      if (diag) printf("--> Full path for environ.def is \"%s\"\n",env_file);
      /* Open the environ.def file as read-only*/
      ccp4printf(2,"Opening file \"%s\"\n",env_file);
      
      envfp = fopen(env_file,"r");
    } 
    /* multiple failures */ 
    if (!envfp) {
      /* Failed to open the file
       Do clean up and exit */
      ccp4fyp_cleanup(ienv,envname,envtype,envext,logical_name,file_name,
		      file_type,file_ext,env_file,def_file,dir,parser);
      ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_CantOpenEnvFile),"ccp4fyp",NULL);
      return 1;
    } else {
      /* Set up a ccp4_parser array to deal with the contents of
	 environ.def */
      parser = (CCP4PARSERARRAY *) ccp4_parse_start(CCP4_MAXTOKS);
      /* Set the delimiters to whitespace, = and .
	 This should split lines into the three components */
      ccp4_parse_delimiters(parser," \t=.",NULL);
      
      /* Read from the file until EOF*/
      while (fgets(line,CCP4_MAXLINE,envfp)) {
	/* Remove the trailing newline from fgets */
	line[strlen(line)-1] = '\0';
	/* Use ccp4_parse to get the tokens on each line */
	ccp4_parse_reset(parser);
	if (ccp4_parse(line,parser) == 3) {
	  env_logical_name = parser->token[0].fullstring;
	  env_file_type    = parser->token[1].fullstring;
	  env_file_ext     = parser->token[2].fullstring;
	  /* Check that we have values for all three components */
	  if (!env_logical_name || !env_file_type || !env_file_ext) {
	    /* Error parsing the line 
	       Do clean up and exit */
	    ccp4fyp_cleanup(ienv,envname,envtype,envext,logical_name,file_name,
			    file_type,file_ext,env_file,def_file,dir,parser);
	    if (envfp) fclose(envfp);
	    ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_ParseEnvFail),"ccp4fyp",NULL);
	    return -1;
	  } else {
	    /* Store in arrays for use when decoding default.def
	       and logical names on command line */
	    if (ienv+1 == CCP4_MAXNAMES) {
	      /* Exceeded the allowed number of logical names
		 Do clean up and exit */
	      ccp4fyp_cleanup(ienv,envname,envtype,envext,logical_name,file_name,
			      file_type,file_ext,env_file,def_file,dir,parser);
	      if (envfp) fclose(envfp);
	      ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_MaxNamesExceeded),"ccp4fyp",NULL);
	      return 1;
	    } else {
	      /* Store logical name in envname */
	      envname[ienv] = (char *)
		ccp4_utils_malloc(sizeof(char)*(strlen(env_logical_name)+1));
	      strcpy(envname[ienv],env_logical_name);
	      /* Store file type in envtype */
	      envtype[ienv] = (char *)
		ccp4_utils_malloc(sizeof(char)*(strlen(env_file_type)+1));
	      strcpy(envtype[ienv],env_file_type);
	      /* File extension in envext */
	      envext[ienv] = (char *)
		ccp4_utils_malloc(sizeof(char)*(strlen(env_file_ext)+1));
	      strcpy(envext[ienv],env_file_ext);
	
	      if (diag) printf("Decoded line: %s = %s.%s\n",envname[ienv],
			       envtype[ienv],envext[ienv]);

	      /* Increment ienv counter for number of name-pairs
		 read in */
	      ienv++;
	    }
	  }
	  /* Reset number of tokens before reading next line */
	  ccp4_parse_reset(parser);
	}
	/* End of loop over lines in file */
      }
      /* Close the environ.def file and free memory storing
	 the filename*/
      fclose(envfp);
      if (env_file) {
	free(env_file);
	env_file = NULL;
      }

      /* Finished with the parser array */
      ccp4_parse_end(parser);
      parser = NULL;
    }
  }

  /* ------------------------------------------------------ */
  /* Default.def file */
  /* ------------------------------------------------------ */

  /* ------------------------------------------------------ */
  /* Non-standard default.def file */
  /* ------------------------------------------------------ */
  /* The search logic is:
     1. If a directory is specified as part of the filename, then
	use the filename as is, otherwise
     2. if the HOME variable is defined, then use $HOME/filename,
        otherwise
     3. use the filename as is (in current directory).
  */
  if (idefault > 0) {
    /* Extract the filename from the argument list - idefault
       points to which argument should hold it */
    if (idefault < argc) {
      if (def_file) free(def_file);
      def_file = (char *) ccp4_utils_malloc(sizeof(char)*(strlen(argv[idefault])+1));
      if (!def_file) {
	/* Couldn't allocate memory to store filename
	   Do clean up and exit */
	ccp4fyp_cleanup(ienv,envname,envtype,envext,logical_name,file_name,
			file_type,file_ext,env_file,def_file,dir,parser);
	ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_AllocFail),"ccp4fyp",NULL);
	return 1;
      }
      strcpy(def_file,argv[idefault]);
      def_init = 1;
      if (diag) printf("CCP4FYP: def file is \"%s\"\n",def_file);
      /* Check whether this includes the path */
      if (dir) free(dir);
      dir = ccp4_utils_pathname(def_file);
      if (dir && dir[0] == '\0') {
	/* No path name - try and use $HOME/file_name */
	if (diag) puts("CCP4FYP: def file has no path.");
	if (home) {
	  tmpstr = ccp4_utils_joinfilenames(home,def_file);
	  if (diag) printf("CCP4FYP: HOME exists, joined filename \"%s\"\n",tmpstr);
	  if (tmpstr) {
	    if (def_file) free(def_file);
	    def_file = tmpstr;
	  } else {
	    /* Failed to make complete filename
	       Do clean up and exit */
	    ccp4fyp_cleanup(ienv,envname,envtype,envext,logical_name,file_name,
			    file_type,file_ext,env_file,def_file,dir,parser);
	    ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_DefPathFail),"ccp4fyp",NULL);
	    return 1;
	  }
	}
      }
      if (diag) printf(" default.def file is \"%s\"\n",def_file);
    } else {
      /* Not enough arguments in the arg list
	 Do clean up and exit */
      if (diag) printf(" no filename found\n");
      ccp4fyp_cleanup(ienv,envname,envtype,envext,logical_name,file_name,
		      file_type,file_ext,env_file,def_file,dir,parser);
      ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_DOptionUseError),"ccp4fyp",NULL);
      return 1;
    }
  }

  /* ------------------------------------------------------ */
  /* Standard default.def file */
  /* ------------------------------------------------------ */
  /* The search logic is:
     1. If the CINCL variable is defined, then use $CINCL/filename,
	otherwise
     2. if the HOME variable is defined, then use $HOME/filename,
        otherwise
     3. use the filename as is (in current directory).
  */
  if (!def_file) {
    /* Use the standard default.def */
    if (diag) printf("--> use standard default.def file\n");
    std_dir = NULL;
    if (cinclude) {
      std_dir = cinclude;
    } else if (home) {
      std_dir = home;
    }
    /* Set the full path for the default.def file */
    if (def_file) free(def_file);
    if (std_dir) {
      if (diag) printf("--> leading directory is \"%s\"\n",std_dir);
      def_file = ccp4_utils_joinfilenames(std_dir,"default.def");
    } else {
      def_file = ccp4_utils_joinfilenames(".","default.def");
    }
    if (!def_file) {
      /* Unable to set the filename
	 Do clean up and exit */
      ccp4fyp_cleanup(ienv,envname,envtype,envext,logical_name,file_name,
		      file_type,file_ext,env_file,def_file,dir,parser);
      ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_DefPathFail),"ccp4fyp",NULL);
      return 1;
    }
    if (diag) printf("--> Full path for default.def is \"%s\"\n",def_file);
  }

  /* ------------------------------------------------------ */
  /* Read in default.def */
  /* ------------------------------------------------------ */
  /* default.def contains lines of the form
	   LOGICALNAME=FILENAME # comments
  */
  if (def_init && def_file) {
    if (diag) printf("CCP4FYP: reading default.def file\n");
    if (diag) printf("--> Full path for default.def is \"%s\"\n",def_file);
    /* Open the default.def file as read-only*/
    ccp4printf(2,"Opening file \"%s\"\n",def_file);
    deffp = fopen(def_file,"r");
    /* did not find the default.def file in std_dir
       so try "PACKAGE_ROOT" */
    if (!deffp) {
      free(def_file);
      def_file = ccp4_utils_joinfilenames(ccp4_utils_joinfilenames(
                  ccp4_utils_joinfilenames(pkgroot,"share"),"ccp4"),"default.def");
      if (diag) printf("CCP4FYP: reading default.def file\n");
      if (diag) printf("--> Full path for default.def is \"%s\"\n",def_file);
      /* Open the default.def file as read-only*/
      ccp4printf(2,"Opening file \"%s\"\n",def_file);
        
      deffp = fopen(def_file,"r"); 
    }
    /* multiple failures */
    if (!deffp) {
      /* Failed to open the file
	 Do clean up and exit */
      ccp4fyp_cleanup(ienv,envname,envtype,envext,logical_name,file_name,
		      file_type,file_ext,env_file,def_file,dir,parser);
      ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_CantOpenDefFile),"ccp4fyp",NULL);
      return 1;
    }
    /* Set a ccp4_parser array to deal with the contents of default.def */
    parser = (CCP4PARSERARRAY *) ccp4_parse_start(CCP4_MAXTOKS);
    /* Set the delimiters to whitespace and =
       This should split lines into the two components */
    ccp4_parse_delimiters(parser," \t=",NULL);

    /* Read from the file until EOF*/
    while (fgets(line,CCP4_MAXLINE,deffp)) {

      /* Remove the trailing newline from fgets */
      line[strlen(line)-1] = '\0';

      /* Use ccp4_parse to get the tokens on each line */
      ccp4_parse_reset(parser);
      if (ccp4_parse(line,parser) == 2) {
	def_logical_name = parser->token[0].fullstring;
	def_file_name    = parser->token[1].fullstring;
	if (!def_logical_name || !def_file_name) {
	  /* Failed to parse the line - 
	     do clean up and exit */
	  ccp4fyp_cleanup(ienv,envname,envtype,envext,logical_name,file_name,
			  file_type,file_ext,env_file,def_file,dir,parser);
	  if (deffp) fclose(deffp);
	  ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_ParseDefFail),"ccp4fyp",NULL);
	  return -1;
	}
	if (diag) printf("Decoded line: %s = %s\n",def_logical_name,def_file_name);

	/* Set up the environment for this pair
	   Don't overwrite any existing logical name */
	ierr = ccp4setenv(def_logical_name,def_file_name,
			  envname,envtype,envext,&ienv,1);
	if (ierr) {
	  /* An error from ccp4setenv
	     Clean up and exit */
	  if (deffp) fclose(deffp);
	  ccp4fyp_cleanup(ienv,envname,envtype,envext,logical_name,file_name,
			  file_type,file_ext,env_file,def_file,dir,parser);
	  return ierr;
	}
	/* Reset number of tokens before reading next line */
	ccp4_parse_reset(parser);
      }  
    }
    /* Close the default.def file */
    fclose(deffp);
    if (def_file) {
      free(def_file);
      def_file = NULL;
    }

    /* Finished with the parser array */
    ccp4_parse_end(parser);
    parser = NULL;
  }

  /* ------------------------------------------------------ */
  /* Process remaining command line arguments */
  /* ------------------------------------------------------ */

  /* Read in the rest of the command line arguments
     These should consist of pairs of arguments i.e.
     <logical name> <file name>
  */

  ccp4printf(2,"Processing Command Line Arguments\n");
  while (iarg < argc) {
    /* Get logical name and uppercase it */
    if (logical_name) free(logical_name); 
    logical_name = (char *) ccp4_utils_malloc((strlen(argv[iarg])+1)*sizeof(char));
    if (diag) printf("--> Raw logical name: \"%s\"\n",argv[iarg]);
    strtoupper(logical_name,argv[iarg]);
    logical_name[strlen(argv[iarg])] = '\0';
    if (diag) printf("--> Logical name: \"%s\"",logical_name);
    iarg++;
    /* Get associated filename */
    if (iarg < argc) {
      if (file_name) free(file_name);
      file_name = (char *) ccp4_utils_malloc((strlen(argv[iarg])+1)*sizeof(char));
      strcpy(file_name,argv[iarg]);
      if (diag) printf("   file name: \"%s\"\n",file_name);
      /* Set up the environment for this pair
	 Do overwrite any existing logical name */
      ierr = ccp4setenv(logical_name,file_name,envname,envtype,envext,&ienv,0);
      if (diag) printf("CCP4FYP: returned from ccp4setenv, ierr = %d\n",ierr);
      if (ierr) {
	/* An error from ccp4setenv
	   Clean up and exit */
	ccp4fyp_cleanup(ienv,envname,envtype,envext,logical_name,file_name,
			file_type,file_ext,env_file,def_file,dir,parser);
	return ierr;
      }
      iarg++;
    } else {
      /* No associated filename
	 Do clean up and exit with error */
      if (diag) printf("  no associated file name\n");
      ccp4fyp_cleanup(ienv,envname,envtype,envext,logical_name,file_name,
		      file_type,file_ext,env_file,def_file,dir,parser);
      ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_LogicalNameUseError),"ccp4fyp",NULL);
      return 1;
    }
    /* Next pair of arguments */
  }
  ccp4printf(2,"End of pre-processing stage\n");

  /* ------------------------------------------------------ */
  /* Finished processing - do local clean up */
  /* ------------------------------------------------------ */

  ccp4fyp_cleanup(ienv,envname,envtype,envext,logical_name,file_name,
		  file_type,file_ext,env_file,def_file,dir,parser);
  if (diag) printf("CCP4FYP: ending\n");
  return 0;
}

/*------------------------------------------------------------------*/

/* ccp4fyp_cleanup

   Internal function.
   Called on exit from ccp4fyp - free memory associated with
   pointers used inside ccp4fyp.
*/
int ccp4fyp_cleanup(int ienv, char **envname, char **envtype, char **envext,
		    char *logical_name, char *file_name, char *file_type,
		    char *file_ext, char *env_file, char *def_file,
		    char *dir, CCP4PARSERARRAY *parser)
{
  int i;
  /* Parser */
  if (parser) ccp4_parse_end(parser);
  /* Free single valued pointers, if set */
  if (file_type) free(file_type);
  if (file_ext) free(file_ext);
  if (env_file) free(env_file);
  if (def_file) free(def_file);
  if (logical_name) free(logical_name);
  if (file_name) free(file_name);
  if (dir) free(dir);
  /* Free arrays of pointers */
  if (ienv > 0) {
    for (i=0; i<ienv; ++i) {
      if (envname[i]) free(envname[i]);
      if (envtype[i]) free(envtype[i]);
      if (envext[i]) free(envext[i]);
    }
  }
  return 0;
}

/*------------------------------------------------------------------*/

/* ccp4setenv

   Set environment variables

   Associates logical_name with value. It is passed arrays of name,
   type and extension for ienv number of name lines read from
   environ.def, which may be used to modify the value before
   assignment.

   If no_overwrt is true then any existing logical_names are not
   redefined.

   ccp4setenv returns 0 on sucess, non-zero values on failure.
   In the case of failure an error code is registered with ccp4_signal
   and a non-zero value is returned to the calling subprogram.
*/
int ccp4setenv(char *logical_name, char* value, char **envname,
		char **envtype, char **envext, int *ienv, int no_overwrt)
{
  int  diag=0;
  int  icount,lext=0,lroot=0,lpath=0,lname,lprognam,procid;
  char *clibd,*cscr;
  char *file_ext=NULL,*file_root=NULL,*file_path=NULL,*file_name=NULL;
  char *tmpstr1=NULL;

  /* Begin */
  if (diag) printf("CCP4SETENV: started, ienv = %d\n",*ienv);

  /* Exit if the logical name already exists and we are in
     no-overwrite mode */
  if (getenv(logical_name) && no_overwrt) {
    /* Do nothing and return without an error */
    ccp4setenv_cleanup(file_ext,file_root,file_path,file_name);
    return 0;
  }

  /* Look for a match between logical_name and the names from
     environ.def */
  if (diag) printf("CCP4SETENV: looking for a match to logical name \"%s\"\n",
		   logical_name);
  /* compare on strlen(envname[icount]) characters, so that e.g. HKLIN1 will match */
  icount = 0;
  while (icount<*ienv && strncmp(logical_name,envname[icount],strlen(envname[icount]))) {
    icount++;
  }
  if (icount == *ienv) {
    /* Not in the list - non-standard logical name */
    if (diag) printf("%s is non-standard logical name\n",logical_name);
    /* Add it to the list
       Name is logical name
       Type is "undef"
       Extension is the extension of the filename */
    if ((*ienv)+1 == CCP4_MAXNAMES) {
      /* Exceeded the allowed number of logical names
	 Clean up, set message and return an error */
      ccp4setenv_cleanup(file_ext,file_root,file_path,file_name);
      ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_MaxNamesExceeded),"ccp4setenv",NULL);
      return 1;
    } else {
      /* Store logical name in envname */
      if (diag) printf("Storing non-standard logical name ... \n");
      envname[*ienv] = (char *)
	ccp4_utils_malloc(sizeof(char)*(strlen(logical_name)+1));
      strcpy(envname[*ienv],logical_name);
      if (diag) printf("... done logical name:\"%s\" ... \n",envname[*ienv]);
      /* Store file type "undef" in envtype */
      envtype[*ienv] = (char *)
	ccp4_utils_malloc(sizeof(char)*(strlen("undef")+1));
      strcpy(envtype[*ienv],"undef");
      if (diag) printf("... done type:\"%s\" ... \n",envtype[*ienv]);
      /* File extension in envext */
      if (file_ext) free(file_ext);
      file_ext = ccp4_utils_extension(value);
      envext[*ienv] = (char *)
	ccp4_utils_malloc(sizeof(char)*(strlen(file_ext)+1));
      strcpy(envext[*ienv],file_ext);
      if (diag) printf("... done extension:\"%s\" ... \n",envext[*ienv]);
      /* Increment the number of values stored */
      ++(*ienv);
      if (diag) printf("CCP4SETENV: now storing %d entries\n",*ienv);
    }
  } else {
    if (diag) printf("CCP4SETENV: Matched logical name for number %d\n",icount);
  }

  /* Split the supplied value up into components to extract get the
     file extension, the file root (i.e. name less extension) and
     the file path
     It is possible there is no path and/or no extension */
  if (diag) printf("CCP4SETENV: supplied file = \"%s\"\n",value);

  /* Get file path */
  if (file_path) free(file_path);
  file_path = ccp4_utils_pathname(value);
  lpath = strlen(file_path)-1;
  if (diag) printf("CCP4SETENV: path = \"%s\"\n",file_path); 

  /* Get file extension */
  if (file_ext) free(file_ext);
  file_ext = ccp4_utils_extension(value);
  lext = strlen(file_ext);
  if (diag) printf("CCP4SETENV: extension = \"%s\"\n",file_ext); 

  /* Get file root */
  if (file_root) free(file_root);
  file_root = ccp4_utils_basename(value);
  lroot = strlen(file_root);
  if (diag) printf("CCP4SETENV: root = \"%s\"\n",file_root);

  /* Add the appropriate file extension if none was supplied
     The exception is for /dev/null or NL: */
  if (!strmatch(value,"/dev/null") && !strmatch(value,"NL:")) {
    if (lext <= 0) {
      /* Add extension */
      if (icount < *ienv) {
	lext = strlen(envext[icount]); 
	if (file_ext) free(file_ext);
	file_ext = (char *) ccp4_utils_malloc(sizeof(char)*(lext+1));
	strncpy(file_ext,envext[icount],(lext+1));
	if (diag) printf("CCP4SETENV: added extension \"%s\"\n",file_ext);
      } else {
	if (diag) printf("CCP4SETENV: not adding extension for non-standard file name\n");
      }
    }
    /* If no path supplied then get appropriate path depending
       on the extension
       .dic, .lib, .bes, .prt, .cif = $CLIBD
       .scr = $CCP4_SCR
    */
    
    /* Fetch the appropriate path name for the file except when its HKLIN or HKLOUT (could be in CWD) */
    if ((lpath < 0) && ((strcmp(logical_name, "HKLIN") != 0) && (strcmp(logical_name, "HKLOUT") != 0)
         && (strcmp(logical_name, "LIBIN") != 0) && (strcmp(logical_name, "LIB_IN") != 0)
         && (strcmp(logical_name, "LIBOUT") != 0) && (strcmp(logical_name, "LIB_OUT") != 0))) {
      /* Fetch the appropriate path name from the environment */
      
      if (strmatch(file_ext,"lib") || strmatch(file_ext,"dic")
	  || strmatch(file_ext,"bes") || strmatch(file_ext,"prt") || strmatch(file_ext, "cif")) {
	/* Fetch CLIBD */
	clibd = (char *) getenv("CLIBD");
	if (clibd) {
	  if (diag) printf("CCP4SETENV: CLIBD = \"%s\"\n",clibd);
	  /* Store in file_path */
	  lpath = strlen(clibd);
	  if (file_path) free(file_path);
	  file_path = (char *) ccp4_utils_malloc(sizeof(char)*(lpath+1));
	  strncpy(file_path,clibd,(lpath+1));
	  if (diag) printf("CCP4SETENV: set file path to CLIBD = \"%s\"\n",file_path);
	} else {
	  /* Couldn't get CLIBD from the environment
	     Clean up, set message and return an error */
	  ccp4setenv_cleanup(file_ext,file_root,file_path,file_name);
	  ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_CantGetClibd),"ccp4setenv",NULL);
	  printf("Couldn't CLIBD from the environment - check setup\n");
	  return 1;
	}

      } else if (strmatch(file_ext,"scr")) {
        /* Scratch files in a special place
	   actually create <ccp4_scr>/<prognm>_.<pid>
	*/
	/* Fetch CCP4_SCR */
	cscr = (char *) getenv("CCP4_SCR");
	if (cscr) {
	  if (diag) printf("CCP4SETENV: CCP4_SCR = \"%s\"\n",cscr);
	  /* Store in file_path */
	  lpath = strlen(cscr);
	  if (file_path) free(file_path);
	  file_path = (char *) ccp4_utils_malloc(sizeof(char)*(lpath+1));
	  strncpy(file_path,cscr,(lpath+1));
	  if (diag) printf("CCP4SETENV: set file path to CCP4_SCR = \"%s\"\n",file_path);
	} else {
	  /* Couldn't get CCP4_SCR
	     Clean up, set message and return an error */
	  ccp4setenv_cleanup(file_ext,file_root,file_path,file_name);
	  ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_CantGetCcp4Scr),"ccp4setenv",NULL);
	  return 1;
	}
        /* Replace <file_root> with <prognam>_<file_root> */
	lprognam = strlen(ccp4ProgramName(NULL));
        tmpstr1 = ccp4_utils_malloc(sizeof(char)*(lprognam + lroot + 2));
        strncpy(tmpstr1,ccp4ProgramName(NULL),lprognam);
	tmpstr1[lprognam] = '\0';
        strncat(tmpstr1,"_",1);
	strncat(tmpstr1,file_root,lroot);
        if (file_root) free(file_root);
	file_root = tmpstr1;
	lroot = strlen(file_root);
	if (diag) printf("CCP4SETENV: updated file_root = \"%s\"\n",file_root);
	/* Replace scr extension with the process id
	 In fact to guarantee that it is always 5 characters,
	 take the id number modulo 100,000 */
	procid = (int) getpid();
	if (diag) printf("CCP4SETENV: initial procid = %d\n",procid);
        procid = procid % CCP4_MODULO;
	if (diag) printf("CCP4SETENV: procid = %d",procid);
	if (file_ext) free(file_ext);
	file_ext = (char*) ccp4_utils_malloc(sizeof(char)*6);
	sprintf(file_ext,"%05d",procid);
        lext = 5;
	if (diag) printf(" giving file extension \"%s\"\n",file_ext);
      }
      /* No special path for this particular extension */
    }
  } else {
    if (diag) printf("CCP4SETENV: detected dev-null\n");
  }

  /* Build the filename */
  lname = lpath + 1;
  file_name = (char *) ccp4_utils_realloc(file_name,sizeof(char)*(lname + 1));
  if (lpath < 0) {
    file_name[0] = '\0';
  } else if (lpath == 0) {
    file_name[0] = PATH_SEPARATOR;
    file_name[1] = '\0';
  } else {
    strncpy(file_name,file_path,lname);
    file_name[lpath] = PATH_SEPARATOR;
    file_name[lpath+1] = '\0';
  }
  if (diag) printf("CCP4SETENV: building filename = \"%s\"\n",file_name);
  lname = lname + lroot;
  file_name = (char *) ccp4_utils_realloc(file_name,sizeof(char)*(lname + 1));
  if (lroot) {
    strcat(file_name,file_root);
  }
  if (diag) printf("CCP4SETENV: building filename = \"%s\"\n",file_name);
  if (lext > 0) {
    lname = lname + lext + 1;
    file_name = (char *) ccp4_utils_realloc(file_name,sizeof(char)*(lname + 1));
    strcat(file_name,".");
    if (lext) {
      strcat(file_name,file_ext);
    }
    file_name[lname] = '\0';
  }
  if (diag) printf("CCP4SETENV: building filename = \"%s\"\n",file_name);

  /* Test that (non-default) input files exist */
  if (icount < *ienv) {
    if (strmatch(envtype[icount],"in") && !no_overwrt) {
      /* Does the file exist? */
      if (diag) printf("CCP4SETENV: checking for existence of input file\n");
      if (ccpexists(file_name)) {
	if (diag) printf("CCP4SETENV: \"%s\" can be opened for reading\n",file_name);
      } else {
	/* File doesn't exist/cannot be opened for reading
	   Clean up, set message and return an error */
	if (diag) printf("CCP4SETENV: \"%s\" cannot be opened for reading\n",file_name);
	printf("File: \"%s\"\nCannot be opened for reading\n",file_name);
	ccp4setenv_cleanup(file_ext,file_root,file_path,file_name);
	ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_CantFindInFile),"ccp4setenv",NULL);
	return -1;
      }
    }
  } else {
    if (diag) printf("CCP4SETENV: cannot determine file type (input or output)\n");
  }

  /* Set the environment variable */
  if (ccpputenv(logical_name,file_name)) {
    if (diag) printf("CCP4SETENV: ccpputenv returned okay\n");
  } else {
    /* Unable to set environment variable
       Clean up and exit the program */
    ccp4setenv_cleanup(file_ext,file_root,file_path,file_name);
    ccp4_signal(CCP4_ERRLEVEL(3) | CGEN_ERRNO(CGENERR_CantSetEnvironment),"ccp4setenv",NULL);
    printf("Cannot create environment variable: \"%s\"\n",logical_name);
    return -1;
  }

  if (diag) {
    tmpstr1 = (char *) getenv(logical_name);
    if (tmpstr1) {
      printf("CCP4SETENV: logical name %s has value \"%s\"\n",logical_name,tmpstr1);
    } else {
      printf("CCP4SETENV: failed to fetch value for logical_name \"%s\"\n",logical_name);
    }
  }

  /* Free dynamically allocated memory before returning */
  ccp4setenv_cleanup(file_ext,file_root,file_path,file_name);
  return 0;
}

/*------------------------------------------------------------------*/

/* ccp4setenv_cleanup

   Internal function.
   Called on exit from ccp4setenv - free memory associated with
   pointers used inside ccp4setenv.
*/
int ccp4setenv_cleanup(char *file_ext, char *file_root, char *file_path,
		       char *file_name)
{
  if (file_ext) free(file_ext);
  if (file_root) free(file_root);
  if (file_path) free(file_path);
  if (file_name) free(file_name);
  return 1;
}

/*------------------------------------------------------------------*/

/* ccpexists

   Check if named file exists
   The check is performed by attempting to fopen the file for
   read only, then immediately closing the file. If this method
   proves to be unsatisfactory then it may be necessary to investigate
   using access or stat instead.

   Returns 1 if the file can be opened for read (=exists), 0 otherwise.
*/
int ccpexists(char *filename)
{
  FILE *fp;

  if (filename) {
    fp = fopen(filename,"r");
    if (fp) {
      fclose(fp);
      return 1;
    }
  }
  return 0;
}

/*------------------------------------------------------------------*/

/* ccpputenv

   This is a wrapper for the ccp4_utils_setenv command.

   It must be supplied with a logical name (the name of a variable
   which will be set in the environment) and a file name (the value which
   will be assigned to that variable).

   Returns 1 if successful and 0 otherwise.

   Notes:
   1. Platform-dependency is encoded in ccp4_utils_setenv.
   2. Dynamically allocated strings passed to ccpputenv should be
   freed by the calling subprogram to avoid memory leaks.
*/
int ccpputenv(char *logical_name, char *file_name)
{
  int ltmpstr,diag=0;
  char *tmpstr;

  if (logical_name && file_name) {
    /* Allocate memory for temporary string */
    ltmpstr = strlen(logical_name) + strlen(file_name) + 1;
    tmpstr = (char *) ccp4_utils_malloc(sizeof(char)*(ltmpstr+1));
    /* putenv requires a string of the form "logical_name=file_name" */
    if (tmpstr) {
      strcpy(tmpstr,logical_name);
      strcat(tmpstr,"=");
      strcat(tmpstr,file_name);
      tmpstr[ltmpstr] = '\0';
      if (diag) printf("CCPPUTENV: string going into ccp4_utils_setenv is \"%s\"\n",tmpstr);
      /* ccp4_utils_setenv returns 0 on success */
      if (ccp4_utils_setenv(tmpstr) == 0) {
        /* free tmpstr here as ccp4_utils_setenv does separate malloc */
        free (tmpstr);
        return 1;
      }
    }
  }
  return 0;
}

/*------------------------------------------------------------------*/

/* ccp4_banner

   Write CCP4 banner to standard output.
*/
void ccp4_banner(void) {

  int diag=0,i,npad;
  char date[11],time[9],prog_vers_str[19],infoline[100];
  char prog_vers_full[MAXLEN_PROGVERSION];

  if (diag) printf("Entering ccp4_banner \n");

  /* Program version number */
  strcpy(prog_vers_full,ccp4_prog_vers(NULL));
  if (strlen(prog_vers_full)) {
    strcpy(prog_vers_str,"version ");
    strncpy(prog_vers_str+8,prog_vers_full,10);
    prog_vers_str[18] = '\0';
  } else {
    /* If no program version available then use the major library
       version number */
    sprintf(prog_vers_str,"version %-10s",ccp4_vers_no());
  }
  /* Trim back the spaces in this string */
  i = strlen(prog_vers_str);
  while (prog_vers_str[--i] == ' ') {
    prog_vers_str[i] = '\0';
  }
  
  printf(" \n");
  printf(" ###############################################################\n");
  printf(" ###############################################################\n");
  printf(" ###############################################################\n");
  /* Information line
     This has information on the version numbers, names and RCS date
     Originally it was printed using the line:

     printf(" ### CCP4 %3s: %-17s  %-18s: %-8s##\n",
     ccp4_vers_no(),ccp4ProgramName(NULL),prog_vers_str,ccp4RCSDate(NULL));

     If the CCP4 version number exceeded three characters this would lead
     to the tail of the line being misaligned.

     This version tries to account for components of the line being
     longer than expected (nb it is still possible to run out of space).
  */
  sprintf(infoline," ### CCP4 %3s: %-17s",ccp4_vers_no(),ccp4ProgramName(NULL));
  /* Trim back spaces in this string */
  i = strlen(infoline);
  while ( i != 0 && infoline[--i] == ' ') {
    infoline[i] = '\0';
  }
  /* Determine how much padding we need based on length of the
     program version number plus what's already been printed*/
  npad = 51 - strlen(infoline) - strlen(prog_vers_str);
  i = strlen(infoline);
  while (npad > 0) {
    infoline[i++] = ' ';
    infoline[i] = '\0';
    --npad;
  }
  sprintf(infoline+i,"%s : %-8s##",prog_vers_str,ccp4RCSDate(NULL));
  printf("%s\n",infoline);
  /* Rest of the banner */
  printf(" ###############################################################\n");
  printf(" User: %s  Run date: %s Run time: %s \n\n\n",
	 ccp4_utils_username(),ccp4_utils_date(date),ccp4_utils_time(time)); 
  printf(" Please reference: Collaborative Computational Project, Number 4. 2011.\n");
  printf(" \"Overview of the CCP4 suite and current developments\". Acta Cryst. D67, 235-242.\n");
  printf(" as well as any specific reference in the program write-up.\n\n");

  if (diag) printf("Leaving ccp4_banner \n");

}
