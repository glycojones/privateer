/*
     ccp4_parser.c: Functions to read in and "parse" CCP4 keyworded input. 
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


/** @file ccp4_parser.c
 *
 *  @brief Functions to read in and "parse" CCP4-style keyworded input.
 *
 *  @author Peter Briggs
 *  @date April 2001
 */

/*   ccp4_parser.c
     Peter Briggs CCP4 April 2001

     Functions to read in and "parse" (scan, in reality) CCP4-style
     "keyworded" input, plus utility routines.

     ccp4_parse_start(maxtokens)
     initialise a CCP4PARSERARRAY to be used in subsequent calls to
     ccp4_parser routines.

     ccp4_parse_end(parser)
     clean up CCP4PARSERARRAY parser after use

     ccp4_parse_init_token(parser,itok)
     initialise a single token in CCP4PARSERARRAY before use

     ccp4_parse_reset(parser)
     initialise CCP4PARSERARRAY before use (includes calls to
     ccp4_parse_init_token to initialise all tokens)

     ccp4_parse_delimiters(parser,delimiters,nulldelimiters)
     set up or restore non-default delimiters

     int ccp4_parse_comments(parser,comment_chars)
     set up or restore non-default comment characters

     ccp4_parse_maxmin(parser,max_exp,min_exp)
     set non-default maximum and minimum values for numerical tokens

     ccp4_parse(line,parser)
     given a string "line", break up into tokens and store in
     CCP4PARSERARRAY "parser"

     ccp4_parser(line,parser,print)
     read input from stdin or external file, break up into tokens
     and store in CCP4PARSERARRAY "parser"

     ccp4_keymatch(keyin1,keyin2)
     compare input strings to see if they match as CCP4-style keywords

     strtoupper(str1,str2)
     convert string to uppercase

     strmatch(str1,str2)
     check if strings are identical

     charmatch(character,charlist)
     check if character appears in a list of possibilities

     doublefromstr(str,....)
     convert a string representation of a number into the number, also
     checks for exponent over/underflow, returns numbers of digits and
     values of "components" (integer and fractional parts, and base-10
     exponent)

*/

/* Header files */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef __FLOAT_H__
#  include <float.h>
#endif
#include <math.h>
#include "ccp4_parser.h"
#include "ccp4_errno.h"
#include "ccp4_sysdep.h"
/* rcsid[] = "$Id$" */

/* stuff for error reporting */
#define CPARSER_ERRNO(n) (CCP4_ERR_PARS | (n))

/* error defs */
#define  CPARSERR_Ok                  0
#define  CPARSERR_MaxTokExceeded      1
#define  CPARSERR_AllocFail           2
#define  CPARSERR_NullPointer         3
#define  CPARSERR_LongLine            4
#define  CPARSERR_CantOpenFile        5
#define  CPARSERR_NoName              6
#define  CPARSERR_ExpOverflow         7
#define  CPARSERR_ExpUnderflow        8
#define  CPARSERR_MatToSymop          9
#define  CPARSERR_SymopToMat         10

/*------------------------------------------------------------------*/

/* Parser Functions */

/*------------------------------------------------------------------*/

/* ccp4_parse_start

   This initialises a CCP4PARSERARRAY to be used with the ccp4_parse/
   ccp4_parser functions.
   The calling function must supply the maximum number of tokens.
*/
CCP4PARSERARRAY* ccp4_parse_start(const int maxtokens)
{
  int itok,diag=0;
  CCP4PARSERARRAY *parsePtr;

  if (diag) printf("CCP4_PARSE_START: ccp4_parse_start starting\n");

  /* Initial check for sensible values */
  if (maxtokens < 1) return NULL;

  /* Return a pointer to a CCP4PARSERARRAY */
  parsePtr = (CCP4PARSERARRAY *) ccp4_utils_malloc(sizeof(CCP4PARSERARRAY));
  if (parsePtr) {
    if (diag) printf("CCP4_PARSE_START: allocated parsePtr\n");
    parsePtr->token = (CCP4PARSERTOKEN *) ccp4_utils_malloc(sizeof(CCP4PARSERTOKEN)*maxtokens);
    if (!parsePtr->token) {
      free(parsePtr);
      parsePtr = NULL;
    } else {
      if (diag) printf("CCP4_PARSE_START: allocated parsePtr->token\n");
      parsePtr->maxtokens = maxtokens;
      parsePtr->fp = NULL;

      /* Explicitly ensure that each token's fullstring member is
	 set to NULL before calling ccp4_parse_reset (which tries to
	 free memory associated with non-NULL fullstrings) since
	 we can't rely on them being created with a NULL value */
      for (itok = 0; itok < maxtokens; itok++)
	parsePtr->token[itok].fullstring = NULL;
      ccp4_parse_reset(parsePtr);
      if (diag) printf("CCP4_PARSE_START: fullstring set to NULL\n");

      /* Initialise the default maximum and minimum allowed
	 exponents for numerical tokens */
      ccp4_parse_maxmin(parsePtr,DBL_MAX_10_EXP,DBL_MIN_10_EXP);
      if (diag) printf("CCP4_PARSE_START: max and min set\n");

      /* Initialise the default delimiter and null delimiter characters
	 with a call to ccp4_parse_delimiters */
      parsePtr->delim=NULL;
      parsePtr->nulldelim=NULL;
      if (!ccp4_parse_delimiters(parsePtr,NULL,NULL)) {
	ccp4_parse_end(parsePtr);
	parsePtr = NULL;
      }
      if (diag) printf("CCP4_PARSE_START: delimiters set\n");

      /* Initialise the default comment characters with a call
	 to ccp4_parse_comments */
      parsePtr->comment=NULL;
      if (!ccp4_parse_comments(parsePtr,NULL)) {
	ccp4_parse_end(parsePtr);
	parsePtr = NULL;
      }
      if (diag) printf("CCP4_PARSE_START: comments set\n");
    }
  }
  
  if (diag) printf("CCP4_PARSE_START: Returning from ccp4_parse_start\n");
  return parsePtr;
}

/*------------------------------------------------------------------*/

/* ccp4_parse_end

   This cleans up a CCP4PARSEARRAY after being used by ccp4_parse/
   ccp4_parser functions.
*/
int ccp4_parse_end(CCP4PARSERARRAY *parsePtr)
{
  int i,maxtokens;

  /* Anything to do? */
  if (parsePtr) {
    /* Free memory for each token */
    maxtokens = parsePtr->maxtokens;
    if (parsePtr->token && parsePtr->maxtokens > 0) {
      for (i=0; i<maxtokens; i++)
	if(parsePtr->token[i].fullstring) free(parsePtr->token[i].fullstring); 
      free(parsePtr->token);
    }
    /* Free memory for lists of comments and delimiters */
    if (parsePtr->comment) free(parsePtr->comment);
    if (parsePtr->delim) free(parsePtr->delim);
    if (parsePtr->nulldelim) free(parsePtr->nulldelim);
    /* Free memory for rest of parserarray structure */
    free(parsePtr);
  }
  /* All done */
  return 0;
}

/*------------------------------------------------------------------*/

/* ccp4_parse_init_token

   Initialise a token in a parser array
   This sets all string members of the specified token to NULL and
   all numerical values (including flags) to zero
*/

int ccp4_parse_init_token(const CCP4PARSERARRAY *parsePtr, const int itok)
{
  if (parsePtr) {
    if (itok < parsePtr->maxtokens) {
      /* Full string is dynamically allocated - free the
	 associated memory, if assigned */
      if (parsePtr->token[itok].fullstring) {
	  free(parsePtr->token[itok].fullstring);
	  parsePtr->token[itok].fullstring = NULL;
      }
      /* Set fixed string tokens to empty string */
      strcpy(parsePtr->token[itok].word,"");
      /* Set numerical value to zero */
      parsePtr->token[itok].value = 0.0;
      /* Set flags to zero */
      parsePtr->token[itok].isstring = 0;
      parsePtr->token[itok].isnumber = 0;
      parsePtr->token[itok].isquoted = 0;
      parsePtr->token[itok].isnull   = 0;
      /* Set start and end positions to zero */
      parsePtr->token[itok].ibeg = 0;
      parsePtr->token[itok].iend = 0;
    }
  }
  return 0;
}

/*------------------------------------------------------------------*/

/* ccp4_parse_reset

   Reinitialise a parser array before calling ccp4_parse

   An application using ccp4_parse (rather than ccp4_parser, which
   also calls this function) can call this function to reset the
   parser array, rather than reinitialising the structure members
   explicitly
*/
int ccp4_parse_reset(CCP4PARSERARRAY *parsePtr)
{
  int itok;

  if (parsePtr) {
    /* Initialise the tokens to have null values */
    for (itok=0; itok<parsePtr->maxtokens; itok++)
      ccp4_parse_init_token(parsePtr,itok);
    /* Initialise number of tokens to zero */
    parsePtr->ntokens = 0;
  }
  return 0;
}

/*------------------------------------------------------------------*/

/* ccp4_parse_delimiters

   This allows the application to set its own delimiter characters
   to be used in the ccp4_parser routines.

   If a NULL pointer is supplied for either of the two lists then
   then the default delimiters are (re)set.

   Returns 1 on success, 0 if there was an error. In the event of
   an error the delimiter lists will be unchanged.
*/

int ccp4_parse_delimiters(CCP4PARSERARRAY *parsePtr,
			  const char *delim, const char *nulldelim)
{
  const char defdelim[]=" \t,=\r",defnulldelim[]=",=";
  char *delimPtr=NULL,*nulldelimPtr=NULL;
  int  ldelim,lnulldelim,istatus=1;

  if (parsePtr) {

    /* If supplied delim is NULL then set to the default */
    if (!delim) {
      ldelim = strlen(defdelim) + 1;
    } else {
      ldelim = strlen(delim) + 1;
    }
    delimPtr = (char *) ccp4_utils_malloc(sizeof(char)*ldelim);
    if (delimPtr) {
      ldelim--;
      if (!delim) {
	strncpy(delimPtr,defdelim,ldelim+1);
      } else {
	strncpy(delimPtr,delim,ldelim+1);
      }
      delimPtr[ldelim] = '\0';
    }

    /* If supplied nulldelim is NULL then set to the default */
    if (!nulldelim) {
      lnulldelim = strlen(defnulldelim) + 1;
    } else {
      lnulldelim = strlen(nulldelim) + 1;
    }
    nulldelimPtr = (char *) ccp4_utils_malloc(sizeof(char)*lnulldelim);
    if (nulldelimPtr) {
      lnulldelim--;
      if (!nulldelim) {
	strncpy(nulldelimPtr,defnulldelim,lnulldelim+1);
      } else {
	strncpy(nulldelimPtr,nulldelim,lnulldelim+1);
      }
      nulldelimPtr[lnulldelim] = '\0';
    }

    /* Assign new delimiters in parser array */
    if (delimPtr && nulldelimPtr) {
      if (parsePtr->delim) free(parsePtr->delim);
      parsePtr->delim = delimPtr;
      if (parsePtr->nulldelim) free(parsePtr->nulldelim);
      parsePtr->nulldelim = nulldelimPtr;
    } else {
      /* There is an error - don't reset the parser array */
      if (delimPtr) free(delimPtr);
      if (nulldelimPtr) free(nulldelimPtr);
      istatus = 0;
    }
  } else {
    istatus = 0;
  }
  return istatus;
}

/*------------------------------------------------------------------*/

/* ccp4_parse_comments

   This allows the application to set its own comment characters
   to be used in the ccp4_parser routines.

   If a NULL pointer is supplied for the list of comment characters
   then the default comment characters are (re)set.

   Returns 1 on success, 0 if there was an error. In the event of
   an error the comment lists will be unchanged.
*/

int ccp4_parse_comments(CCP4PARSERARRAY *parsePtr, const char *comment_chars)
{
  const char def_comment_chars[]="#!";
  char  *commentPtr=NULL;
  int   lcomment,istatus=1;

  if (parsePtr) {

    /* If the supplied comment list is NULL then set to the default */
    if (!comment_chars) {
      lcomment = strlen(def_comment_chars) + 1;
    } else {
      lcomment = strlen(comment_chars) + 1;
    }
    commentPtr = (char *) ccp4_utils_malloc(sizeof(char)*lcomment);
    if (commentPtr) {
      if (!comment_chars) {
	strncpy(commentPtr,def_comment_chars,lcomment);
      } else {
	strncpy(commentPtr,comment_chars,lcomment);
      }
      lcomment--;
      commentPtr[lcomment] = '\0';
    }

    /* Assign the new comments in the parser array */
    if (commentPtr) {
      if (parsePtr->comment) free(parsePtr->comment);
      parsePtr->comment = commentPtr;
    } else {
      /* There was an error - don't reset the parser array */
      istatus = 0;
    }
  } else {
    /* The parser was unset on entry - also an error */
    istatus = 0;
  }
  return istatus;
}

/*------------------------------------------------------------------*/

/* ccp4_parse_maxmin

   This allows the application to set its own maximum and minimum
   exponent values, which are used as limits when evaluating the values
   of numerical tokens in order to avoid over/underflow.
*/

int ccp4_parse_maxmin(CCP4PARSERARRAY *parsePtr, const double max_exponent,
		      const double min_exponent)
{
  if (parsePtr) {
    parsePtr->max_exponent = (double) max_exponent;
    parsePtr->min_exponent = (double) min_exponent;
  }
  return 1;
}

/*------------------------------------------------------------------*/

/* ccp4_parse

   This is a scanner based on the old CCP4 Fortranic PARSE routine.
 
   It takes an input string ("line") and returns the number of
   tokens ("ntokens") which are delimited by certain characters
   (defaulted to space, tab, comma, equals - these can be changed by
   the application using a call to ccp4_parse_delimiters).
   Information about the tokens themselves is returned as members of
   elements in an array ("tokens") of type CCP4PARSERTOKEN (see header
   file for definition and members).

   Substrings can be delimited by single- or double-quotes but must 
   be surrounded by delimiters to be recognised.

   An unquoted comment character (defaulted to ! or #) in the input line
   introduces a trailing comment which is ignored. The comment
   characters can be changed using a call to ccp4_parse_comments.

   Null fields are denoted by two adjacent null delimiters (defaulted
   to comma and equals - these can be changed by the application
   using a call to ccp4_parse_delimiters).

   ccp4_parse returns the number of tokens found in the line. The
   tokens are returned via the CCP4PARSERARRAY parser.

   Arguments:

   line   = pointer to a null-terminated string of characters,
            forming the input to be processed. Unaltered on
	    output.
   parser = pointer to a CCP4PARSERARRAY structure which will
            be used to hold the results of processing the input
	    line.
*/

int ccp4_parse(const char *line, CCP4PARSERARRAY *parser)
{
  int quotedstring,starttoken,endtoken;
  char this_char,next_char,matchquote;

  int llen,ich,lword,diag=0;
  int token,nulltoken,isquote,iscommt=0,isdelim;
  double value;
  char *delim,*nulldelim,*comm;
  char quot[]="\"\'";
  int  ibeg,iend,start;

  double intvalue,frcvalue,expvalue;
  int    intdigits,frcdigits,expdigits;

  /* Local parser variables */
  int ntok,maxtok;
  CCP4PARSERTOKEN *tokenarray;

  /* Begin */

  /* Set diag = 1 and recompile to switch on diagnostic output */
  if (diag) printf("CCP4_PARSE: ccp4_parse starting\n");

  maxtok = parser->maxtokens;
  ntok = parser->ntokens;
  if (ntok < 0) ntok = 0;

  /* Initialise pointer for local version of the token array */
  tokenarray = parser->token;

  /* Initialise pointers for lists of comments and delimiters */
  comm = parser->comment;
  delim = parser->delim;
  nulldelim = parser->nulldelim;

  /* Don't process any tokens already dealt with */
  if (ntok > 0) {
    start = tokenarray[ntok-1].iend + 1;
    /* Special case: if the last token was quoted
       then in fact the character at this position will be
       a quote - if so then step on another character */
    if (charmatch(line[start],quot)) {
      if (diag) printf("CCP4_PARSE: start character is a quote, stepping on\n");
      start++;
    }
  } else {
    start = 0;
  }

  /* Don't process empty lines */
  llen = strlen(line);
  if (diag) printf("CCP4_PARSE: Line is: \"%s\"\nLength of line is %d\n",line,llen);
  if (llen > 0) {

    /* Initialise flags and counters */
    quotedstring = 0;
    token        = 0;
    nulltoken    = 0;

    /* Process the line two characters at a time */

    if (diag) printf("CCP4_PARSE: Start position for parsing line is %d\n",start);
    for (ich=start-1; ich<llen; ich++) {

      /* Examine characters in pairs */

      /* Deal with special case, start and end of line
	 This is to avoid accessing non-existant parts of the string */
      if (ich < start) {
	this_char = delim[0];
      } else {
	this_char = line[ich];
      }
      if (ich == llen-1) {
	next_char = delim[0];
      } else {
	next_char = line[ich+1];
      }
      if (diag) printf("CCP4_PARSE: %d: Current char = \"%c\" : \n",ich,this_char);

      /* Set flags for this round
	 The pairs of characters are analysed and flags set
	 accordingly to signal actions after the analysis */
      starttoken = 0;
      endtoken   = 0;

      if (!quotedstring) {
	/* Not in a quoted string so look for opening quotes
	   This is delimiter followed by quote character */

	/* Is current character a delimiter character? */
	isdelim = charmatch(this_char,delim);
	
	/* Is next charcter a quote character? */
	isquote = charmatch(next_char,quot);

	/* Start of quoted string? */
	if (isdelim && isquote) {

	  if (diag) printf("CCP4_PARSE: start of a quoted string\n");
	  quotedstring = 1;
	  starttoken = 1;
	  matchquote = next_char;
	  /* Advance to the next character immediately */
	  ich++;

	} else {
	  /* Not the start of a quoted string */

	  /* Is the current character the start of a comment? */
	  iscommt = charmatch(this_char,comm);
	  if (diag)
	    printf("CCP4_PARSE: character = %c comments = %s iscommt = %d\n",
		   this_char,comm,iscommt);
	  if (iscommt) {
	    if (diag) printf("CCP4_PARSE: start of a comment\n");
	  }

	  /* Is the current character the start or end of a token? */
	  if (!isdelim) {

	    /* End of the current token?
	       Only if we are in a token now, and the next
	       character is a delimiter or a comment */
	    if (token && (charmatch(next_char,delim) || charmatch(next_char,comm))) {
	      if (diag) printf("CCP4_PARSE: end of a token\n");
	      endtoken = 1;
	    }

	  } else {

	    /* Start of a token?
	       Only if we are not in a token now, and the next character
	       is not a delimiter (or a comment) too */
	    if (!token && !(charmatch(next_char,delim) || charmatch(next_char,comm))) {
		if (diag) printf("CCP4_PARSE: start of a token\n");
		starttoken = 1;
	    }
	    
	    /* Or - could be a null token?
	       This is a pair of null delimiters together */
	    if (!token && charmatch(this_char,nulldelim)
		&& charmatch(next_char,nulldelim)) {
	      if (diag) printf("CCP4_PARSE: null token\n");
	      starttoken = 1;
	      nulltoken  = 1;
	    }

	  }
	  /* End of token identification */
	}
	
      } else {
	/* Inside a quoted string so look for closing quotes
	   This is a matching quote followed by a delimiter */
	
	/* Is current character a matching quote? */
	isquote = (this_char == matchquote);

	/* Is next charcter a delimiter character? */
	isdelim = charmatch(next_char,delim);

	/* End of quoted string */
	if (isdelim && isquote) {
	  if (diag) printf("CCP4_PARSE: end of a quoted string\n");
	  quotedstring = 0;
	  endtoken = 1;
	  matchquote = ' ';
	}
	/* End of character analyses */
      }

      /* Start of new token */
      if (starttoken) {
	/* Check we don't have maximum number of tokens already */
	if (ntok < maxtok) {
	  /* Set flags */
	  token = 1;
	  if (quotedstring || !nulltoken) {
	    ibeg = ich + 1;
	  } else {
	    ibeg = ich;
	  }
	  if (diag) printf("CCP4_PARSE: Start of a new token... ibeg = %d\n",ibeg);
	  /* Initialise values */
	  tokenarray[ntok].fullstring = NULL;
	  tokenarray[ntok].value    = 0.0;
	  tokenarray[ntok].isstring = 0;
	  tokenarray[ntok].strlength = 0;
	  tokenarray[ntok].isnumber = 0;
	  tokenarray[ntok].intdigits= 0;
	  tokenarray[ntok].frcdigits= 0;
	  tokenarray[ntok].isquoted = 0;
	  tokenarray[ntok].isnull   = 0;
	  /* Flag start of quoted string */
	  if (quotedstring) { tokenarray[ntok].isquoted = 1; }
	  /* Flag null token */
	  if (nulltoken) {
	    if (diag) printf("CCP4_PARSE: Null token\n");
	    tokenarray[ntok].isnull = 1;
	    tokenarray[ntok].iend   = ich;
	    token = 0;
	    nulltoken = 0;
	  }
	} else {
	  /* Maximum number of tokens already found */
          ccp4_signal(CPARSER_ERRNO(CPARSERR_MaxTokExceeded),"ccp4_parse",NULL);
	  parser->ntokens = ntok;
	  return ntok;
	}
	if (diag) printf("CCP4_PARSE: This is the start of token %d\n",ntok);
      }
      /* End of new token */

      /* End of current token */
      if (endtoken) {
	token = 0;
	/* Exclude trailing quote from the token? */
	if (tokenarray[ntok].isquoted) {
	  iend = ich - 1;
	} else {
	  iend = ich;
	}
	if (diag) printf("CCP4_PARSE: End of a token... iend = %d\n",iend);
	/* Store the full token in the array */
	lword = iend - ibeg + 1;
	if (diag) printf("CCP4_PARSE: lword = %d - start char = %c, end char = %c\n",
			 lword,line[ibeg],line[iend]);
	tokenarray[ntok].fullstring = (char *) ccp4_utils_malloc(sizeof(char)*(lword+1));
	if (tokenarray[ntok].fullstring) {
	  strncpy(tokenarray[ntok].fullstring,&line[ibeg],lword);
	  tokenarray[ntok].fullstring[lword] = '\0';
	  if (diag) printf("CCP4_PARSE: Token is \"%s\"\n",tokenarray[ntok].fullstring);
	} else {
          ccp4_signal(CPARSER_ERRNO(CPARSERR_AllocFail),"ccp4_parse",NULL);
	}
	tokenarray[ntok].ibeg = ibeg;
	tokenarray[ntok].iend = iend;
	/* Store the 4 character token in the array */
	if (lword > 4) lword = 4; 
	strncpy(tokenarray[ntok].word,&line[ibeg],lword);
	tokenarray[ntok].word[lword] = '\0';
	/* Determine numerical value (if any) */
	if (doublefromstr(tokenarray[ntok].fullstring,parser->max_exponent,
			  parser->min_exponent,&value,&intvalue,&intdigits,
			  &frcvalue,&frcdigits,&expvalue,&expdigits)) {
	  if (diag) printf("CCP4_PARSE: This has a numerical value of %lf\n",value);
	  tokenarray[ntok].value     = value;
	  tokenarray[ntok].isnumber  = 1;
	  tokenarray[ntok].intdigits = intdigits;
	  tokenarray[ntok].frcdigits = frcdigits;
	} else {
	  if (diag) printf("CCP4_PARSE: There is no numerical value for this token\n");
	  tokenarray[ntok].isstring = 1;
	  tokenarray[ntok].strlength   = strlen(tokenarray[ntok].fullstring);
	}
	/* Reset flags etc ready for next token*/
	token  = 0;
	value  = 0.0;
	/* Increment number of tokens */
	ntok++;
	if (diag) printf("CCP4_PARSE: This is the end of a token\n");
      }

      /* Don't do any more processing after a comment */
      
      if (iscommt) {
	parser->ntokens = ntok;
	if (diag) printf("CCP4_PARSE: returning after a comment\n");
	return ntok;
      }
      /* Check the next pair of characters */
    }
    /* Reset the number of tokens in the parser array */
    parser->ntokens = ntok;
    if (diag) printf("CCP4_PARSE: ntokens = %d, and ntok = %d\n",parser->ntokens,ntok);
  }
  if (diag) printf("CCP4_PARSE: returning at function end\n");
  return ntok;
}

/*------------------------------------------------------------------*/

/* ccp4_parser

   This is based on the old CCP4 Fortranic PARSER routine.

   The normal behaviour is to read "keyworded" data from the input
   stream, and interpret it. Stdin is the default, but a line
   starting with @<name> starts reading from file <name> until eof.

   Each input line may be continued on the next line by the continuation
   characters `&', `-' or `\' at the end of the input line. This
   character is dropped from the list returned to the calling application.

   Pass in a zero length line to force reading from the command line.
   nchars is the maximum number of characters which will be read into the line.
   (If line is not blank then it will be processed and more input
   read in if it ends in a continuation character, or forces reading from
   an external file.)

   The "print" argument should be supplied as 0 to suppress echoing of the
   input lines to standard output.

   ccp4_parser returns the number of tokens parsed in the input line. The
   results of the parsing are stored as members of the CCP4PARSEARRAY
   structure "parser" and can be accessed by the application program.

   The function returns the number of tokens, or 0 on reaching end of file.
   On encountering an unrecoverable error ccp4_parser returns -1. 

   Arguments:

   line   = pointer to a null-terminated string of characters,
            forming the input to be processed.
	    On input can either be an empty string ("") or
	    contain characters to be processed (see above for
	    description).
	    On output "line" will be overwritten with the actual
	    input line, up to nchar characters.
   nchars = maximum number of characters that can be read into
            "line" i.e. the size of "line" in memory.
   parser = pointer to a CCP4PARSERARRAY structure which will
            be used to hold the results of processing the input
	    line.
   print  = flag controlling echoing of input lines to stdout.
            print=0: suppress echoing of lines to stdout
	    Otherwise echoing is turned on.
*/

int ccp4_parser(char *line, const int nchars, CCP4PARSERARRAY *parser,
		const int print)
{
  int fromstdin=0,fromfile=0,fromapp=0,diag=0;
  int nch,nold,continuation,first,trunc,llen,buflen;
  char *linein=NULL,filename[200];

  /* Local parser variables */
  int  ntok;
  FILE *filein;
  CCP4PARSERTOKEN *tokenarray;

  /* Undocumented feature - if print < 0 then also print out
     diagnostic info */
  if (print < 0) {
    /*print = 1;*/
    diag  = 1;
  }

  /* Begin */
  if (diag) printf("CCP4_PARSER: ccp4_parser starting\n");

  /* Abort if parser is a NULL pointer */
  if (!parser) {
    ccp4_signal(CPARSER_ERRNO(CPARSERR_NullPointer),"ccp4_parser",NULL);
    return -1;
  }

  /* Abort if line is NULL pointer */
  if (!line) {
    ccp4_signal(CPARSER_ERRNO(CPARSERR_NullPointer),"ccp4_parser",NULL);
    return -1;
  }

  /* Reset the parser array for this sweep
     This will prevent phantom values from an earlier call
     to ccp4_parser persisting in the parser array */
  ccp4_parse_reset(parser);

  /* Blank the keyword */
  strcpy(parser->keyword,"");

  /* Initialise local variables and pointers */
  tokenarray = parser->token;
  ntok       = parser->ntokens;
  filein     = parser->fp;
  if (diag) printf("CCP4_PARSER: parser->ntokens = %d, ntok = %d\n",
		   parser->ntokens,ntok);

  /* Set up an internal buffer for input
     The buffer is over-allocated (twice as long as the max string
     length allocated for line by the calling application)
  */
  buflen = (nchars*2)+1;
  linein = (char *) ccp4_utils_malloc(buflen*sizeof(char));

  if (!linein) {
    ccp4_signal(CPARSER_ERRNO(CPARSERR_AllocFail),"ccp4_parser",NULL);
    return 0;
  }

  /* Use nch as a count of the number of remaining characters
     in line */
  nch = nchars;

  /* If line is empty then read from the standard input
     Otherwise process the line from the application first */
  if (strlen(line)==0) {
    if (!filein) {
      if (diag) printf("CCP4_PARSER: Reading from stdin\n");
      fromstdin = 1;
    } else {
      if (diag) printf("CCP4_PARSER: Reading from file\n");
      fromfile = 1;
    }
  } else {
    if (diag) printf("CCP4_PARSER: Reading line supplied by the application program\n");
    fromapp = 1;
  }

  /* Set flag for first line of input */
  first = 1;
  
  /* Set flag for line continuation */
  continuation = 1;
  
  /* Start the input loop */
  while (continuation) {

    if (diag) printf("CCP4_PARSER: starting loop\n");
    /* Read input from stdin a line at a time */
    if (fromstdin) {
      if (diag) printf("CCP4_PARSER: reading from stdin...\n");
      if (!fgets(linein,buflen,stdin)) {
	/* Jump out at this point if eof is reached from
	   stdin */
	return 0;
      }
    } else if (fromfile) {
      if (diag) printf("CCP4_PARSER: reading from external file...\n");
      if (!fgets(linein,buflen,filein)) {
	/* Return to input from stdin if eof is read from
	   the external file */
	if (diag) printf("CCP4_PARSER: End of external file reached\n");
	fclose(filein);
	filein = NULL;
	fromfile  = 0;
	fromstdin = 1;
	/* Blank the line and reset the first flag to
	   force reading from standard input immediately */
	linein[0] = '\0';
	ntok      = 0;
	parser->ntokens = ntok;
	first     = 1;
      }
    } else if (fromapp) {
      if (diag) printf("CCP4_PARSER: reading from application...\n");
      /* If this line contains a continuation mark then
	 read from stdin next time around */
      strncpy(linein,line,nchars);
      linein[nchars]='\0';
    }

    /* Strip any trailing newline e.g. from fgets */
    llen = strlen(linein);
    if (llen > 0)
      if (linein[llen-1] == '\n') {
	linein[llen-1] = '\0';
	llen--;
      }
    
    /* If previous line ended with a continuation character
       then append this one to it
       Check that we don't overflow the number of characters
       specified by the application */

    if (llen > nch) {
      ccp4_signal(CPARSER_ERRNO(CPARSERR_LongLine),"ccp4_parser",NULL);
    }
    if (first) {
      strncpy(line,linein,nch);
      first = 0;
    } else {
      strncat(line,linein,nch);
    }
    nch = nchars - llen;
    if (diag) {
      printf("CCP4_PARSER: line = \"%s\"\n",line);
      printf("CCP4_PARSER: remaining available characters = %d\n",nch);
    }
      
    /* Use ccp4_parse to break the input line up into tokens
       Only parse the latest chunk - ccp4_parse will append
       new tokens onto the tokenarray */
    nold = ntok;
    ntok = ccp4_parse(line,parser);

    if (diag) printf("CCP4_PARSER: ntok = %d, nold = %d\n",ntok,nold);

    /* Have we found more tokens since last time? */
    if (ntok != nold) {
      /* Check first token to see if it is an instruction
	 to read from an external file */
      if (!fromfile && tokenarray[0].word[0] == '@') {
	if (diag) printf("CCP4_PARSER: Instruction to read from external file\n");
	/* Get filename and attempt to open file */
	if (tokenarray[0].fullstring) {
	  llen = strlen(tokenarray[0].fullstring);
	  strncpy(filename,&tokenarray[0].fullstring[1],llen);
	  if (diag) printf("CCP4_PARSER: External file name is \"%s\"\n",filename);
	  /* Open the external file as read-only */
	  filein = fopen(filename,"r");
	  if (!filein) {
            ccp4_signal(CPARSER_ERRNO(CPARSERR_CantOpenFile),"ccp4_parser",NULL);
	  } else {
	    fromstdin = 0;
	    fromfile  = 1;
	  }
	} else {
	  /* Token with file name is null */
          ccp4_signal(CPARSER_ERRNO(CPARSERR_NoName),"ccp4_parser",NULL);
	}
	/* Blank the line and reset the number of tokens
	   to force reading from the external file immediately */
	line[0] = '\0';
	ntok    = 0;
	parser->ntokens = ntok;

      /* Check last token to see if it is continuation
	 character */
      } else if (ntok > 0 &&
	  (strmatch("&",tokenarray[ntok-1].word) ||
	   strmatch("\\",tokenarray[ntok-1].word) ||
	   strmatch("-",tokenarray[ntok-1].word))) {
	if (diag) printf("CCP4_PARSER: Detected continuation character\n");
	/* It's a continuation mark 
	   Set flag to indicate this fact in later rounds */
	continuation = 1;
	/* Truncate the line to remove the continuation
	   character */
	if (ntok > 1)
	  trunc = tokenarray[ntok-1].ibeg;
	else
	  trunc = 0;
	if (diag) printf("CCP4_PARSER: Continuation character should be at position %d\n\"%c\" is the character at this position\n",trunc,line[trunc]);
	line[trunc] = '\0';
	/* Lose the last token */ 
	ntok--;
	parser->ntokens = ntok;
      } else {
	/* Not a continuation character */
	continuation = 0;
      }
      
    } else {
      /* Didn't get any more tokens from the last pass
	 Check if it is a blank line or comment line */
      if (ntok == 0) {
	/* Echo comment line to stdout and blank
	   the line */
	if (strlen(line) > 0) {
	  if (print) printf(" Comment line--- %s\n",line);
	  line[0] = '\0';
          nch = nchars;
	}
        if (fromapp) continuation = 0;
      }
    }
    if (diag) printf("CCP4_PARSER: Continuation = %d\n",continuation);

    /* If the line was supplied by the application but is now being continued
       then make sure we read from stdin next time */
    if (continuation && fromapp) {
      if (diag) printf("CCP4_PARSER: switching to stdin\n");
      fromapp   = 0;
      fromstdin = 1;
    }
  }

  /* Fetch and uppercase keyword */
    if (ntok > 0) {
      strtoupper(parser->keyword,tokenarray[0].word);
      parser->keyword[strlen(tokenarray[0].word)] = '\0';
      if (diag) printf("CCP4_PARSER: Keyword is %s\n",parser->keyword);
      /*Echo the line to standard output */ 
      if (print) printf(" Data line--- %s\n",line); 
    } else {
      parser->keyword[0] = '\0';
    }

  free(linein);
  /* Update the returned variables */
  parser->fp = filein;
  
  if (diag) printf("CCP4_PARSER: Returning from ccp4_parser\n");
  return ntok;
}

/*------------------------------------------------------------------*/

/* ccp4_keymatch

   Returns 1 if keywords keyin1 and keyin2 are "identical", 0 otherwise.

   Keywords are identical if they are the same up to the first four
   characters, independent of case.
*/
int ccp4_keymatch(const char *keyin1, const char *keyin2)
{
  int  len1,len2;
  char key1[5],key2[5],keyup1[5],keyup2[5];

  /* Initial check */
  if (!keyin1 || !keyin2) return 0;

  /* Compare truncated lengths */
  len1 = strlen(keyin1);
  if (len1 > 4) len1 = 4;
 
  len2 = strlen(keyin2);
  if (len2 > 4) len2 = 4;

  /* If lengths don't match then keywords can't be identical */
  if (len1 != len2) return 0;

  /* If supplied words are longer than four characters then
     truncate them after the fourth character */
  strncpy(key1,keyin1,len1);
  key1[len1] = '\0';

  strncpy(key2,keyin2,4);
  key2[len2] = '\0';

  /* Convert strings to uppercase */
  strtoupper(keyup1,key1);
  keyup1[len1] = '\0';
  strtoupper(keyup2,key2);
  keyup2[len2] = '\0';

  /* Compare using strmatch */
  return strmatch(keyup1,keyup2);
}

char *strtoupper (char *str1, const char *str2)
{
  int len2,i;

  if (!str2) return NULL;
  
  len2 = strlen(str2);
  if (len2 > 0) for (i=0; i<len2 ; i++) str1[i] = toupper(str2[i]);
  str1[len2] = '\0';
  return str1;
}

char *strtolower (char *str1, const char *str2)
{
  int len2,i;

  if (!str2) return NULL;
  
  len2 = strlen(str2);
  if (len2 > 0) for (i=0; i<len2 ; i++) str1[i] = tolower(str2[i]);
  str1[len2] = '\0';
  return str1;
}

/*------------------------------------------------------------------*/

/* strmatch

   Compare two strings.

   Returns 1 if strings are identical, 0 if they differ.
*/
int strmatch (const char *str1, const char *str2)
{
  int len1,len2,i;

  /* Don't process null strings */
  if (!str1 || !str2) return 0;

  len1 = strlen(str1);
  len2 = strlen(str2);

  /* Cannot be identical if lengths differ */
  if (len1 != len2) return 0;

  /* Check character by character
     If any character differs then strings cannot be identical */
  for (i=0; i<len1; i++)
    if (str1[i] != str2[i]) return 0;

  /* Must be a match */
  return 1;
}

/*------------------------------------------------------------------*/

/* charmatch

   Returns 1 if character matches one of the characters in the string,
   0 otherwise.
*/
int charmatch(const char character, const char *charlist)
{
  int jdo,ismatch=0;

  if (charlist) {
    jdo = 0;
    while (jdo<strlen(charlist) && !ismatch) {
      if (charlist[jdo] == character) ismatch = 1;
      jdo++;
    }
  }

  return ismatch;
}

/*------------------------------------------------------------------*/

/* doublefromstr

   Determine whether string represents a valid number, and return the
   numerical value if it is.
   The function returns 1 for a valid number, 0 otherwise.

   Valid numbers are represented by:
   ^[+-]?[0-9]*\.?[0-9]*[Ee]?[+-]?[0-9]*$
   (I think this is the correct regular expression...? -pjx)

   The component parts are: an integer part representing the value
   of the number upto the decimal point; a decimal fraction
   representing the value of the number upto the base-10 exponent;
   a base-10 exponent. The full value of the string is then
   int_part + frac_part * 10^(exp_part)

   The value of each component is returned separately via the argument
   list. The number of "digits" (= number of characters in the string)
   of each component are also returned - note that these include non-
   significant digits such as leading zeroes in the integer part,
   but non-numeric characters (+-eE) are not counted.

   The function can also trap for large or small exponents, to
   avoid overflow or underflow.
   The maximum and minimum allowed exponents are supplied by the calling
   application as the arguments max_exp and min_exp. (Values can be
   be taken from float.h).
   In the event of the limits being exceeded, only the values of the
   components are evaluated; the total value of the number is returned as
   zero.
*/
int doublefromstr(const char *str, const double max_exp, const double min_exp,
		  double *valuePtr, double *intvaluePtr, int *intdigitsPtr,
		  double *frcvaluePtr, int *frcdigitsPtr,
		  double *expvaluePtr, int *expdigitsPtr)
{
  int    lstr,ichar=0,char_value,diag=0;
  int    sign,expsign,point,exponent,is_int,is_frc,is_exp;
  int    n_int_digits=0,n_frc_digits=0,n_exp_digits=0;
  char   this_char, this_str[2];
  double int_value=0,frc_value=0,exp_value=0;

  if (diag) printf("DOUBLEFROMSTR: starting\n");
  if (diag) printf("DOUBLEFROMSTR: Biggest exponent is %lf, smallest is %lf\n",max_exp,min_exp);

  /* Initialise exported variables */
  *valuePtr = 0.0;
  *intvaluePtr = 0.0;
  *frcvaluePtr = 0.0;
  *expvaluePtr = 0.0;
  *intdigitsPtr = 0;
  *frcdigitsPtr = 0;
  *expdigitsPtr = 0;

  /* Nothing to do for an empty string */
  lstr = strlen(str);
  if (!lstr) return 0;

  /* Initialise internal flags and variables */
  sign     = 1;
  expsign  = 1;
  point    = -1;
  exponent = -1;

  is_int   = 1;
  is_frc   = 0;
  is_exp   = 0;

  /* Process the string character by character */
  while (ichar < lstr) {
    /* Get current character */
    this_char = str[ichar];
    if (diag) printf("DOUBLEFROMSTR: This char = %c",this_char);

    /* Check: is this character a digit? */
    if (!isdigit(this_char)) {

      /* Not a digit */
      if (diag) printf(" is not a digit...\n");

      if (this_char == '+' || this_char == '-') {
	/* Sign? i.e. + or -
	   This can only occur in two places: immediately at the
	   start of the putative number, or immediately at the start
	   of the putative exponent */
	if (ichar == 0 && this_char == '-') {
	  sign = -1;
	} else if (ichar == exponent + 1 && this_char == '-') {
	  expsign = -1;
	} else if (this_char != '+') {
	  return 0;
	}

      } else if (this_char == '.') {
	/* Decimal point? i.e. .
	   There can only be one decimal point */
	if (point > -1) return 0;
	point = ichar;
	is_int = 0;
	is_frc = 1;

      } else if (toupper(this_char) == 'E') {
        char next_char = (ichar+1 < lstr ) ? str[ichar+1] : '\0';
        if ( next_char == '+' || next_char == '-')
           next_char = (ichar+2 < lstr ) ? str[ichar+2] : '\0';
        /* require the next active character after E to be a digit */
        if ( !isdigit(next_char) )  return 0;
	/* Exponent? i.e. e or E
	   There can only be one exponent */
        if (exponent > -1) return 0;
        exponent = ichar;
        is_int = 0;
        is_frc = 0;
        is_exp = 1;
      } else {
	/* Not a permissible character
	   This is not a number so get out now */
	if (diag) printf("DOUBLEFROMSTR: Not a permitted character - exiting\n");
	return 0;
      }
    } else {
      /* It is a digit
	 Build up the value of each component */

      if (diag) printf(" is a digit ...\n");

      this_str[0] = this_char;
      this_str[1] = '\0';
      char_value = atoi(this_str);

      if (is_int) {
	/* Integer part of the number */
	n_int_digits++;
	int_value = int_value * 10.0 + (double) char_value;
	if (diag) printf("DOUBLEFROMSTR: Processing integer component: value = %lf, #digits = %d\n",int_value,n_int_digits);
      } else if (is_frc) {
	/* Decimal part of the number */
	n_frc_digits++;
	frc_value = frc_value + ((double) char_value)/pow(10.0,(double) n_frc_digits);
	if (diag) printf("DOUBLEFROMSTR: Processing decimal component: value = %lf, #digits = %d\n",frc_value,n_frc_digits);
      } else if (is_exp) {
	/* Exponent */
	n_exp_digits++;
	exp_value = exp_value * 10.0 + (double) char_value;
	if (diag) printf("DOUBLEFROMSTR: Processing exponential component: value = %lf, #digits = %d\n",exp_value,n_exp_digits);
      }

    }
    /* Next character */
    ichar++;
  }
  /*
    Done loop over characters - if we have got this far then it
    must be a number
  */

  /* Set component values */
  int_value = int_value * (double) sign;
  frc_value = frc_value * (double) sign;
  exp_value = exp_value * (double) expsign;

  /* Export component values */
  *intvaluePtr = int_value;
  *frcvaluePtr = frc_value;
  *expvaluePtr = exp_value;

  /* Export numbers of 'digits' */
  *intdigitsPtr = n_int_digits;
  *frcdigitsPtr = n_frc_digits;
  *expdigitsPtr = n_exp_digits;

  /* Is the exponent out-of-range? */

  /* There are two considerations:
     (i) can pow(10.0,exp_value) actually be evaluated?
     (ii) can int_part * pow(10.0,exp_value) be evaluated?
     This second is an issue for numbers with int_part > 0.
  */
  if ( (exp_value + (double) (n_int_digits - 1) > max_exp) && 
        (n_int_digits || n_frc_digits) ) {
    ccp4_signal(CPARSER_ERRNO(CPARSERR_ExpOverflow),"doublefromstr",NULL);
    printf("DOUBLEFROMSTR: Token is \"%s\"\n",str);
    *valuePtr = 0.0;
  } else if ( (exp_value < min_exp) && 
        (n_int_digits || n_frc_digits) ) {
    ccp4_signal(CPARSER_ERRNO(CPARSERR_ExpUnderflow),"doublefromstr",NULL);
    printf("DOUBLEFROMSTR: Token is \"%s\"\n",str);
    *valuePtr = 0.0;
  } else {
    /* Evaluate the number to get a value */
    *valuePtr = int_value + frc_value;
    if (is_exp) *valuePtr = (*valuePtr)*pow(10.0,exp_value);
  }

  if (diag) printf("DOUBLEFROMSTR: Integer component = %lf, (%d digits)\n",
		   *intvaluePtr,*intdigitsPtr);
  if (diag) printf("DOUBLEFROMSTR: Decimal component = %lf, (%d digits)\n",
		   *frcvaluePtr,*frcdigitsPtr);
  if (diag) printf("DOUBLEFROMSTR: Exponent component = %lf, (%d digits)\n",
		   *expvaluePtr,*expdigitsPtr);
  if (diag) printf("DOUBLEFROMSTR: Finished - value is determined to be %lf\n",*valuePtr);

  return 1;
}

ccp4_symop symop_to_rotandtrn(const char *symchs_begin, const char *symchs_end) {

  float rsm[4][4];

  symop_to_mat4(symchs_begin, symchs_end, rsm[0]);
  return (mat4_to_rotandtrn((const float (*)[4])rsm));

}

/*------------------------------------------------------------------*/

/* symop_to_mat4

   Translates a single symmetry operator string into a 4x4 quine
   matrix representation
   NB Uses a utility function (symop_to_mat4_err) when reporting
   failures.

   Syntax of possible symop strings:

   real space symmetry operations, e.g. X+1/2,Y-X,Z
   reciprocal space operations,    e.g. h,l-h,-k
   reciprocal axis vectors,        e.g. a*+c*,c*,-b*
   real space axis vectors,        e.g. a,c-a,-b

   The strings can contain spaces, and the coordinate and translation
   parts may be in either order.

   The function returns 1 on success, 0 if there was a failure to
   generate a matrix representation.
*/
const char *symop_to_mat4(const char *symchs_begin, const char *symchs_end, float *rot)
{
  int no_real =0, no_recip = 0, no_axis = 0;          /* counters */
  int col = 3, nops = 0;
  int nsym = 0, init_array = 1;
  float sign = 1.0f, value = 0.0f, value2;
  char *cp, ch;
  const char *ptr_symchs = symchs_begin;
  int j,k;                                 /* loop variables */
  int Isep = 0;                             /* parsed seperator? */

  while (ptr_symchs < symchs_end) {
    ch = *ptr_symchs;

    /* Parse symop */
    if (isspace(ch)) {
      /* Have to allow symop strings to contain spaces for
	 compatibility with older MTZ files
	 Ignore and step on to next character */
      ++ptr_symchs;
      continue;
    } else if (ch == ',' || ch == '*') {           
      ++ptr_symchs;
      if (value == 0.0f && col == 3) {
        /* nothing set, this is a problem */
        ccp4_signal(CPARSER_ERRNO(CPARSERR_SymopToMat),"symop_to_mat4",NULL);
        return NULL ;
      } else {
        Isep = 1;     /* drop through to evaluation*/
      }
    } else if (ch == 'X' || ch == 'x') {
      no_real++, col = 0;
      if (value == 0.0f) value = sign * 1.0f;
      ++ptr_symchs;
      continue;
    } else if (ch == 'Y' || ch == 'y') {
      no_real++, col = 1;
      if (value == 0.0f) value = sign * 1.0f;
      ++ptr_symchs;
      continue;
    } else if (ch == 'Z' || ch == 'z') {
      no_real++, col = 2;
      if (value == 0.0f) value = sign * 1.0f;
      ++ptr_symchs;
      continue;
    } else if (ch == 'H' || ch == 'h') {
      no_recip++, col = 0;
      if (value == 0.0f) value = sign * 1.0f;
      ++ptr_symchs;
      continue;
    } else if (ch == 'K' || ch == 'k') {
      no_recip++, col = 1;
      if (value == 0.0f) value = sign * 1.0f;
      ++ptr_symchs;
      continue;
    } else if (ch == 'L' || ch == 'l') {
      no_recip++, col = 2;
      if (value == 0.0f) value = sign * 1.0f;
      ++ptr_symchs;
      continue;
    } else if (ch == 'A' || ch == 'a') {
      no_axis++, col = 0;
      if (value == 0.0f) value = sign * 1.0f;
      if (*++ptr_symchs == '*' && ( no_axis != 3 || no_recip )) ++ptr_symchs;
      continue;
    } else if (ch == 'B' || ch == 'b') {
      no_axis++, col = 1;
      if (value == 0.0f) value = sign * 1.0f;
      if (*++ptr_symchs == '*' && ( no_axis != 3 || no_recip )) ++ptr_symchs;
      continue;
    } else if (ch == 'C' || ch == 'c') {
      no_axis++, col = 2;
      if (value == 0.0f) value = sign * 1.0f;
      if (*++ptr_symchs == '*' && ( no_axis != 3 || no_recip )) ++ptr_symchs;
      continue;
    } else if (ch == '+' || ch == '-') {
      sign = ((ch == '+')? 1.0f : -1.0f) ;
      ++ptr_symchs;
      if ( value == 0.0f && col == 3) 
	continue;
      /* drop through to evaluation */
    } else if ( ch == '/') {
      ++ptr_symchs;
      if (value == 0.0f) {
	/* error */
	symop_to_mat4_err(symchs_begin);
	return NULL;
      }
      value2 = strtod(ptr_symchs, &cp);
      if (!value2) {
	/* error */
	symop_to_mat4_err(symchs_begin);
	return NULL;
      }
      /* Nb don't apply the sign to value here
	 It will already have been applied in the previous round */
      value = (float) value/value2;
      ptr_symchs = cp;
      continue;
    } else if ( isdigit(ch) || ch == '.') {
      value = sign*strtod(ptr_symchs, &cp);
      ptr_symchs = cp;
      continue;
    } else {
      ++ptr_symchs;
      continue;
    }
   
  /* initialise and clear the relevant array  (init_array == 1)*/
  /* use knowledge that we are using a [4][4] for rot */
  if (init_array) {
    init_array = 0;       
    for (j = 0 ; j !=4 ; ++j)
      for (k = 0; k !=4 ; ++k)
        rot[(((nsym << 2) + k ) << 2) +j] = 0.0f;
     rot[(((nsym << 2 ) + 3 )  << 2) +3] = 1.0f;
  }
 
    /* value to be entered in rot */
    rot[(((nsym << 2) + nops) << 2) + col] = value;

    /* have we passed a operator seperator */
    if (Isep) {
      Isep = 0;
      ++nops;
      sign = 1.0f;
      if (nops == 3 ) { ++nsym; nops=0 ; init_array = 1; }
    }

    /* reset for next cycle */
    col = 3;
    value = 0.0f;
    no_recip = 0, no_axis = 0, no_real = 0;
  }

  /* Tidy up last value */
  if (value) rot[(((nsym << 2) + nops) << 2) + col] = value;

  if (nops<2) {
    /* Processed fewer than 3 operators, raise an error */
    symop_to_mat4_err(symchs_begin);
    return NULL;
  }

  /* Return with success */
  return ptr_symchs;
}

/* Internal function: report error from symop_to_mat4_err */
int symop_to_mat4_err(const char *symop)
{
  printf("\n **SYMMETRY OPERATOR ERROR**\n\n Error in interpreting symop \"%s\"\n\n",
	 symop);
  ccp4_signal(CPARSER_ERRNO(CPARSERR_SymopToMat),"symop_to_mat4",NULL);
  return 1;
}

ccp4_symop mat4_to_rotandtrn(const float rsm[4][4]) {

  int i,j;
  ccp4_symop symop;

  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) 
      symop.rot[i][j]=rsm[i][j];
    symop.trn[i]=rsm[i][3];
  }

  return (symop);
}

char *rotandtrn_to_symop(char *symchs_begin, char *symchs_end, const ccp4_symop symop)
{
  float rsm[4][4];

  rotandtrn_to_mat4(rsm,symop);
  return(mat4_to_symop(symchs_begin,symchs_end,(const float (*)[4])rsm));
}

void rotandtrn_to_mat4(float rsm[4][4], const ccp4_symop symop) {

  int i,j;

  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) 
      rsm[i][j]=symop.rot[i][j];
    rsm[i][3]=symop.trn[i];
    rsm[3][i]=0.0;
  }
  rsm[3][3]=1.0;
}

char *mat4_to_symop(char *symchs_begin, char *symchs_end, const float rsm[4][4])
{
  static char axiscr[] = {'X','Y','Z'};
  static char numb[] = {'0','1','2','3','4','5','6','7','8','9'};
  static int npntr1[12] = { 0,1,1,1,0,1,0,2,3,5,0,0 };
  static int npntr2[12] = { 0,6,4,3,0,2,0,3,4,6,0,0 };

  int jdo10, jdo20, irsm, itr, ist;
  register char *ich;
  int debug=0;

  if (debug)
    for (jdo20 = 0; jdo20 != 4; ++jdo20) 
      printf("Input matrix: %f %f %f %f \n",rsm[jdo20][0],rsm[jdo20][1],
          rsm[jdo20][2],rsm[jdo20][3]);

  /* blank output string */
  for (ich = symchs_begin; ich < symchs_end; ++ich) 
    *ich = ' ';
  ich = symchs_begin;

  for (jdo20 = 0; jdo20 != 3; ++jdo20) {
    *ich = '0';
    ist = 0;    /* ---- Ist is flag for first character of operator */
    for (jdo10 = 0; jdo10 != 4; ++jdo10) {
      
      if (rsm[jdo20][jdo10] != 0.f) {
	irsm = (int) rint(fabs(rsm[jdo20][jdo10]));
	
	if ( rsm[jdo20][jdo10] > 0. && ist) {
	  if (ich >= symchs_end) {
	    ccp4_signal(CCP4_ERRLEVEL(3) | CPARSER_ERRNO(CPARSERR_MatToSymop), 
			"mat4_to_symop 1", NULL);
	    return NULL; }
	  *ich++ = '+';
	} else if ( rsm[jdo20][jdo10] < 0.f ) {
	  if (ich >= symchs_end) {
	    ccp4_signal(CCP4_ERRLEVEL(3) | CPARSER_ERRNO(CPARSERR_MatToSymop), 
			"mat4_to_symop 2", NULL);
	    return NULL; }
	  if (jdo10 != 3) {
	    *ich++ = '-';
	  } else {
	    /* translation part is forced to be positive, see below */
	    *ich++ = '+';
	  }
          ist = 1;
	}
      
	if (jdo10 != 3) {
	  /* rotation part */
	  if (ich+1 >= symchs_end) {
	    ccp4_signal(CCP4_ERRLEVEL(3) | CPARSER_ERRNO(CPARSERR_MatToSymop), 
		      "mat4_to_symop 3", NULL);
	    return NULL; }
	  if (irsm != 1) {
	    *ich++ = numb[irsm];
	    *ich++ = axiscr[jdo10];
	  } else {
	    *ich++ = axiscr[jdo10];
	  }
	  ist = 1;
        } else {
	  /* translation part */
	  itr = (int) rint(rsm[jdo20][3]*12.0);
          while (itr < 0) itr += 12;
          itr = (itr - 1) % 12;
          if (npntr1[itr] > 0) {
	   if (ich+2 >= symchs_end) {
	    ccp4_signal(CCP4_ERRLEVEL(3) | CPARSER_ERRNO(CPARSERR_MatToSymop), 
		      "mat4_to_symop 4", NULL);
	    return NULL; }
	   *ich++ = numb[npntr1[itr]];
	   *ich++ = '/';
           *ich++ = numb[npntr2[itr]];
	  } else {
           *--ich = ' ';
	  }
	}
      }
    }
    if (jdo20 != 2) {
      if (*ich == '0')
        ++ich;
      if (ich+2 >= symchs_end) {
        ccp4_signal(CCP4_ERRLEVEL(3) | CPARSER_ERRNO(CPARSERR_MatToSymop), 
                     "mat4_to_symop 5", NULL);
        return NULL; }
      *ich++ = ',';
      *ich++ = ' ';
      *ich++ = ' ';
    }
  }
  return symchs_begin;
}

char *mat4_to_recip_symop(char *symchs_begin, char *symchs_end, const float rsm[4][4])
{
  char *symop;
  size_t lsymop;
  register char *ich, *ich_out;

  lsymop = symchs_end-symchs_begin;
  symop = (char *) ccp4_utils_malloc(lsymop*sizeof(char));

  mat4_to_symop(symop, symop+lsymop, rsm);
  ich_out = symchs_begin;
  for (ich = symop; ich < symop+lsymop; ++ich) {
    if (*ich == 'X') {
      if (ich_out == symchs_begin || (ich_out > symchs_begin && 
          *(ich_out-1) != '-' && *(ich_out-1) != '+')) *ich_out++ = '+';
      *ich_out++ = 'h';
    } else if (*ich == 'Y') {
      if (ich_out == symchs_begin || (ich_out > symchs_begin && 
          *(ich_out-1) != '-' && *(ich_out-1) != '+')) *ich_out++ = '+';
      *ich_out++ = 'k';
    } else if (*ich == 'Z') {
      if (ich_out == symchs_begin || (ich_out > symchs_begin && 
          *(ich_out-1) != '-' && *(ich_out-1) != '+')) *ich_out++ = '+';
      *ich_out++ = 'l';
    } else if (*ich == ' ') {
      /* skip */
    } else {
      *ich_out++ = *ich;
    }
  }
  while (ich_out < symchs_end) *ich_out++ = ' ';

  free (symop);
  return symchs_begin;
}
