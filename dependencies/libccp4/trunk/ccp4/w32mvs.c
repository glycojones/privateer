/*
     w32mvs.c: functions required by MVS
     Copyright (C) 2003  Alun Ashton

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

#ifdef _MSC_VER

#include "w32mvs.h"


/* For communication from `getopt' to the caller.
   When `getopt' finds an option that takes an argument,
   the argument value is returned here.
   Also, when `ordering' is RETURN_IN_ORDER,
   each non-option ARGV-element is returned here.  */

char *optarg = 0;

/* Index in ARGV of the next element to be scanned.
   This is used for communication to and from the caller
   and for communication between successive calls to `getopt'.

   On entry to `getopt', zero means this is the first call; initialize.

   When `getopt' returns EOF, this is the index of the first of the
   non-option elements that the caller should itself scan.

   Otherwise, `optind' communicates from one call to the next
   how much of ARGV has been scanned so far.  */

int optind = 0;

/* The next char to be scanned in the option-element
   in which the last option character we returned was found.
   This allows us to pick up the scan where we left off.

   If this is zero, or a null string, it means resume the scan
   by advancing to the next ARGV-element.  */

static char *nextchar;

/* Callers store zero here to inhibit the error message
   for unrecognized options.  */

int opterr = 1;

/* Describe how to deal with options that follow non-option ARGV-elements.

   UNSPECIFIED means the caller did not specify anything;
   the default is then REQUIRE_ORDER if the environment variable
   _OPTIONS_FIRST is defined, PERMUTE otherwise.

   REQUIRE_ORDER means don't recognize them as options.
   Stop option processing when the first non-option is seen.
   This is what Unix does.

   PERMUTE is the default.  We permute the contents of `argv' as we scan,
   so that eventually all the options are at the end.  This allows options
   to be given in any order, even with programs that were not written to
   expect this.

   RETURN_IN_ORDER is an option available to programs that were written
   to expect options and other ARGV-elements in any order and that care about
   the ordering of the two.  We describe each non-option ARGV-element
   as if it were the argument of an option with character code zero.
   Using `-' as the first character of the list of option characters
   requests this mode of operation.

   The special argument `--' forces an end of option-scanning regardless
   of the value of `ordering'.  In the case of RETURN_IN_ORDER, only
   `--' can cause `getopt' to return EOF with `optind' != ARGC.  */

static enum { REQUIRE_ORDER, PERMUTE, RETURN_IN_ORDER } ordering;

static void exchange (char **argv)
{
  int nonopts_size
    = (last_nonopt - first_nonopt) * sizeof (char *);
  char **temp = (char **) alloca (nonopts_size);

  /* Interchange the two blocks of data in argv.  */

  bcopy (&argv[first_nonopt], temp, nonopts_size);
  bcopy (&argv[last_nonopt], &argv[first_nonopt],
	 (optind - last_nonopt) * sizeof (char *));
  bcopy (temp, &argv[first_nonopt + optind - last_nonopt],
	 nonopts_size);

  /* Update records for the slots the non-options now occupy.  */

  first_nonopt += (optind - last_nonopt);
  last_nonopt = optind;
}

int getopt (int argc,char **argv,char *optstring)
{
  /* Initialize the internal data when the first call is made.
     Start processing options with ARGV-element 1 (since ARGV-element 0
     is the program name); the sequence of previously skipped
     non-option ARGV-elements is empty.  */

  if (optind == 0)
    {
      first_nonopt = last_nonopt = optind = 1;

      nextchar = 0;

      /* Determine how to handle the ordering of options and nonoptions.  */

      if (optstring[0] == '-')
	ordering = RETURN_IN_ORDER;
      else if (getenv ("_POSIX_OPTION_ORDER") != 0)
	ordering = REQUIRE_ORDER;
      else
	ordering = PERMUTE;
    }

  if (nextchar == 0 || *nextchar == 0)
    {
      if (ordering == PERMUTE)
	{
	  /* If we have just processed some options following some non-options,
	     exchange them so that the options come first.  */

	  if (first_nonopt != last_nonopt && last_nonopt != optind)
	    exchange (argv);
	  else if (last_nonopt != optind)
	    first_nonopt = optind;

	  /* Now skip any additional non-options
	     and extend the range of non-options previously skipped.  */

	  while (optind < argc
		 && (argv[optind][0] != '-'
		     || argv[optind][1] == 0))
	    optind++;
	  last_nonopt = optind;
	}

      /* Special ARGV-element `--' means premature end of options.
	 Skip it like a null option,
	 then exchange with previous non-options as if it were an option,
	 then skip everything else like a non-option.  */

      if (optind != argc && !strcmp (argv[optind], "--"))
	{
	  optind++;

	  if (first_nonopt != last_nonopt && last_nonopt != optind)
	    exchange (argv);
	  else if (first_nonopt == last_nonopt)
	    first_nonopt = optind;
	  last_nonopt = argc;

	  optind = argc;
	}

      /* If we have done all the ARGV-elements, stop the scan
	 and back over any non-options that we skipped and permuted.  */

      if (optind == argc)
	{
	  /* Set the next-arg-index to point at the non-options
	     that we previously skipped, so the caller will digest them.  */
	  if (first_nonopt != last_nonopt)
	    optind = first_nonopt;
	  return EOF;
	}
	 
      /* If we have come to a non-option and did not permute it,
	 either stop the scan or describe it to the caller and pass it by.  */

      if (argv[optind][0] != '-' || argv[optind][1] == 0)
	{
	  if (ordering == REQUIRE_ORDER)
	    return EOF;
	  optarg = argv[optind++];
	  return 0;
	}

      /* We have found another option-ARGV-element.
	 Start decoding its characters.  */

      nextchar = argv[optind] + 1;
    }

  /* Look at and handle the next option-character.  */

  {
    char c = *nextchar++;
    char *temp = (char *) strchr (optstring, c);

    /* Increment `optind' when we start to process its last character.  */
    if (*nextchar == 0)
      optind++;

    if (temp == 0 || c == ':')
      {
	if (opterr != 0)
	  {
	    if (c < 040 || c >= 0177)
	      fprintf (stderr, "%s: unrecognized option, character code 0%o\n",
		       argv[0], c);
	    else
	      fprintf (stderr, "%s: unrecognized option `-%c'\n",
		       argv[0], c);
	  }
	return '?';
      }
    if (temp[1] == ':')
      {
	if (temp[2] == ':')
	  {
	    /* This is an option that accepts an argument optionally.  */
	    if (*nextchar != 0)
	      {
	        optarg = nextchar;
		optind++;
	      }
	    else
	      optarg = 0;
	    nextchar = 0;
	  }
	else
	  {
	    /* This is an option that requires an argument.  */
	    if (*nextchar != 0)
	      {
		optarg = nextchar;
		/* If we end this ARGV-element by taking the rest as an arg,
		   we must advance to the next element now.  */
		optind++;
	      }
	    else if (optind == argc)
	      {
		if (opterr != 0)
		  fprintf (stderr, "%s: no argument for `-%c' option\n",
			   argv[0], c);
		optarg = 0;
	      }
	    else
	      /* We already incremented `optind' once;
		 increment it again when taking next ARGV-elt as argument.  */
	      optarg = argv[optind++];
	    nextchar = 0;
	  }
      }
    return c;
  }
}
#endif
