/*
     ccp4_parser.h: Headers for functions to read in and "parse" CCP4 keyworded input.
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

/** @page cparser_page CParser library
 *
 * @verbatim

<!-- ::INDEX_INFO::CParser library::Library::::C Software Library for CCP4-style parsing:::::::: -->

   @endverbatim
 *
 *  @section cparser_file_list File list

<ul>
<li>ccp4_parser.h - contains details of the C/C++ API
</ul>

 *  @section cparser_overview Overview
 
These functions do CCP4-style parsing, as used for processing keywords
of CCP4 programs, MTZ header records, etc.

 *  @section cparser_usage Usage

The following code snippets illustrate how the functions might be used
to read from stdin:
<pre>

int           ntok=0;
char          line[201],*key;
CCP4PARSERTOKEN * token=NULL;
CCP4PARSERARRAY * parser;

  parser = (CCP4PARSERARRAY *) ccp4_parse_start(20);
  key   = parser->keyword;
  token = parser->token;

  RC   = 0;
  while (!RC) {

    line[0] = '\0';
    ntok = ccp4_parser(line,200,parser,1);

    if (ntok < 1) {

      RC = 111;

    } else {      

      if (ccp4_keymatch("MINDIST",key))  {
	if (ntok != 2) {
	  ccperror ( 1,"MINDIST requires a single numerical argument" );
	  RC = -100;
	} else {
	  minDist = token[1].value;
        }
      }	else  {
	printf ( "Unrecognised keyword \"%s\"\n",token[0].fullstring );
	RC = -118;
      }
    }
  }

  ccp4_parse_end ( parser );

</pre>

 *  @section cparser_examples Examples

See the distributed programs <a href="../ncont.html">NCONT</a> and
<a href="../pdbcur.html">PDBCUR</a>.

 */

/** @file ccp4_parser.h
 *
 *  @brief Functions to read in and "parse" CCP4-style keyworded input.
 *
 *  @author Peter Briggs
 *  @date April 2001
 */

/*------------------------------------------------------------------*/

/* Macro definitions */

/*------------------------------------------------------------------*/

#ifndef __CCP4_Parser__
#define __CCP4_Parser__

/* rcsidhhh[] = "$Id$" */

/* note order: these must be outside CCP4 namespace */
#include <stdio.h>
#include"ccp4_utils.h"
#include"ccp4_spg.h"

/* Macro to make C functions callable from C++ */
#ifdef  __cplusplus
namespace CCP4 {
extern "C" {
typedef CSym::ccp4_symop ccp4_symop;
#endif

/*------------------------------------------------------------------*/

/* Structures and typedefs */

/*------------------------------------------------------------------*/

/* CCP4 Parser token
   Construct to hold the information about a single token */

typedef struct {
  char   *fullstring;   /* Full string containing all of token */
  char   word[5];       /* First four characters of token */
  double value;         /* Equivalent numerical value */
  int    isstring;      /* Flag: true if token is character string */
  int    strlength;     /* Number of characters in whole token (strings only) */
  int    isnumber;      /* Flag: true if token is number */
  int    intdigits;     /* Number of 'digits' preceeding the decimal point
			   (numbers only) */
  int    frcdigits;     /* Number of 'digits' after the decimal point (numbers
			   only) */
  int    isquoted;      /* Flag: true if token is contained in quotes */
  int    isnull;        /* Flag: true if token is null field */
  int    ibeg,iend;     /* Begin and end character positions of token
			   in input line */
} CCP4PARSERTOKEN;

/* CCP4 Parser array
   Construct to hold the information about a parsed line */

typedef struct {
  /* "Public" members */
  char   keyword[5];      /* Keyword (=token[1].token, uppercased) */
  int    ntokens;         /* Number of tokens */
  CCP4PARSERTOKEN *token; /* Array of tokens */
  /* "Private" members */
  FILE   *fp;             /* Pointer to an external command file */
  int    maxtokens;       /* Maximum number of tokens allowed */
  char   *delim;          /* List of delimiter characters */
  char   *nulldelim;      /* List of null delimiter characters */
  char   *comment;        /* List of comment characters */
  double max_exponent;    /* Largest allowed exponent for numerical tokens */
  double min_exponent;    /* Smallest allowed exponent for numerical tokens */
} CCP4PARSERARRAY;     

/*------------------------------------------------------------------*/

/* Function Prototypes */

/*------------------------------------------------------------------*/

/* Core cparser functions */

/** Initialise a CCP4PARSERARRAY to be used in subsequent calls to
 *  ccp4_parser routines. The calling function must supply the maximum 
 *  number of tokens on a line (including continuation lines).
 * @param maxtokens maximum number of tokens on a line
 * @return pointer to a new CCP4PARSERARRAY structure
 */
CCP4PARSERARRAY* ccp4_parse_start(const int maxtokens);

/** Cleans up a CCP4PARSEARRAY after being used by ccp4_parse/
   ccp4_parser functions.
 * @param parsePtr pointer to a CCP4PARSERARRAY structure
 * @return 0 on completion
 */
int ccp4_parse_end(CCP4PARSERARRAY *parsePtr);

int ccp4_parse_init_token(const CCP4PARSERARRAY *parsePtr, const int itok);

int ccp4_parse_delimiters(CCP4PARSERARRAY *parsePtr, const char *delim,
				  const char *nulldelim);

int ccp4_parse_comments(CCP4PARSERARRAY *parsePtr, const char *comment_chars);

int ccp4_parse_maxmin(CCP4PARSERARRAY *parsePtr, const double max_exponent,
			      const double min_exponent);

int ccp4_parse_reset(CCP4PARSERARRAY *parsePtr);

int ccp4_parse(const char *line, CCP4PARSERARRAY *parser);

/** The main function for parsing lines, either supplied or read
 * from stdin.
 * @param line pointer to a null-terminated string of characters,
 * forming the input to be processed. On input can either be an empty 
 * string ("") which forces reading from stdin, or contain characters 
 * to be processed. On output "line" will be overwritten with the actual
 * input line.
 * @param n maximum number of characters that can be read into
 * "line" i.e. the size of "line" in memory.
 * @param parser pointer to a CCP4PARSERARRAY structure which will
 * be used to hold the results of processing the input line.
 * @param print flag controlling echoing of input lines to stdout.
 * print=0: suppress echoing of lines to stdout. Otherwise echoing is 
 * turned on.
 * @return Number of tokens found.
 */
int ccp4_parser(char *line, const int n, CCP4PARSERARRAY *parser,
			const int print);

/* External utility functions */

/** Test whether two keywords are identical. Keywords are identical if 
 * they are the same up to the first four characters, independent of case.
 * @param keyin1 keyword 1.
 * @param keyin2 keyword 2.
 * @return 1 if keywords keyin1 and keyin2 are "identical", 0 otherwise.
 */
int ccp4_keymatch(const char *keyin1, const char *keyin2);

/* Internal utility functions */

/** Convert string to uppercase.
 * @param str1 On exit str1 will contain uppercased copy of str2
 * @param str2 Input string
 * @return str1
 */
char *strtoupper (char *str1, const char *str2);

/** Convert string to lowercase.
 * @param str1 On exit str1 will contain lowercased copy of str2
 * @param str2 Input string
 * @return str1
 */
char *strtolower (char *str1, const char *str2);

int strmatch (const char *str1, const char *str2);

int charmatch(const char character, const char *charlist);

int doublefromstr(const char *str, const double max_exp, const double min_exp,
			  double *valuePtr, double *intvaluePtr, int *intdigitsPtr,
			  double *frcvaluePtr, int *frcdigitsPtr,
			  double *expvaluePtr, int *expdigitsPtr);

/** Convert symmetry operator as string to ccp4_symop struct.
 * @param symchs_begin pointer to beginning of string
 * @param symchs_end pointer to end of string (i.e. last character
 *   is *(symchs_end-1) )
 * @return pointer to ccp4_symop struct
 */
ccp4_symop symop_to_rotandtrn(const char *symchs_begin, const char *symchs_end);

/** Convert symmetry operator as string to matrix.
 * This is Charles' version of symfr. Note that translations
 * are held in elements [*][3] and [3][3] is set to 1.0
 * @param symchs_begin pointer to beginning of string
 * @param symchs_end pointer to end of string (i.e. last character
 *   is *(symchs_end-1) )
 * @param rot 4 x 4 matrix operator
 * @return  NULL on error, final position pointer on success
 */
const char * symop_to_mat4(const char *symchs_begin, const char *symchs_end, float *rot);
int symop_to_mat4_err(const char *symop);
ccp4_symop mat4_to_rotandtrn(const float rsm[4][4]);
/* This is Charles' version of symtr */
char *rotandtrn_to_symop(char *symchs_begin, char *symchs_end, const ccp4_symop symop);
void rotandtrn_to_mat4(float rsm[4][4], const ccp4_symop symop);

/** Convert symmetry operator as matrix to string.
 * This is Charles' version of symtr. Note that translations
 * are held in elements [*][3] and [3][3] is set to 1.0
 * @param symchs_begin pointer to beginning of string
 * @param symchs_end pointer to end of string (i.e. last character
 *   is *(symchs_end-1) )
 * @param rsm 4 x 4 matrix operator
 * @return pointer to beginning of string
 */
char *mat4_to_symop(char *symchs_begin, char *symchs_end, const float rsm[4][4]);

/** Convert symmetry operator as matrix to string in reciprocal space notation.
 * This is Charles' version of symtr. Note that translations
 * are held in elements [*][3] and [3][3] is set to 1.0
 * @param symchs_begin pointer to beginning of string
 * @param symchs_end pointer to end of string (i.e. last character
 *   is *(symchs_end-1) )
 * @param rsm 4 x 4 matrix operator
 * @return pointer to beginning of string
 */
char *mat4_to_recip_symop(char *symchs_begin, char *symchs_end, const float rsm[4][4]);

#ifdef __cplusplus
}
}
#endif

#endif  /* __CCP4_Parser__ */
