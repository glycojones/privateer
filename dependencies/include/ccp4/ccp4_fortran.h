/*
     ccp4_fortran.h: header file for Fortran APIs
     Copyright (C) 2001   Eugene Krissinel

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

/** @file ccp4_fortran.h
 *  header file for Fortran APIs
 *  Eugene Krissinel
 */

#ifndef __CCP4_FORTRAN
#define __CCP4_FORTRAN

#include "ccp4_types.h"
/* rcsidhh[] = "$Id$" */

/* stardent is now obsolete, but we retain this category in case it is useful later */
#ifdef CALL_LIKE_STARDENT
  /* SStrParam is used in Ardent-like machines' fortran calls */
  /* for passing a string parameter */
  DefineStructure(SStrPar)
  struct SStrPar  {
    pstr S;
    int  len;
    int  id;
  };
#endif

#define _LVTOB(l) ((long) ((l) == 0 ? 0 : 1))  
#define _BTOLV(l) ((int) ((l) == 0 ? 0 : 1))
#if defined (__OSF1__) || defined (__osf__)
#undef _BTOLV
#define _BTOLV(l) ((int) ((l) == 0 ? 0 : -1))  
#endif    

/*
     Macro  FORTRAN_SUBR(NAME,name,p_send,p_sstruct,p_sflw)
   makes function header statements that allow for linking with
   programs written in FORTRAN.
  
     Parameters:
  
     NAME      name of the FORTRAN subroutine in capital letters
     name      name of the FORTRAN subroutine in small letters
     p_send    parameter list (in brackets) with string lengths
               attached to the end of it (see below)
     p_sstruct parameter list (in brackets) with strings passed
               as complex parameters, or structures
     p_sflw    parameter list (in brackets) with string lengths
               following immediately the string parameters
               (see below)
  
     All non-string parameters must be passed as pointers, in
   the same order as they enter the FORTRAN call. Rules for
   the string parameters are as follows.
  
     1. All strings should be specified as of 'fpstr' type.
        The 'fpstr' type is defined below and depends on the
        platform:
  
          a) whenever length of string is passed as a separate
             parameter ( CALL_LIKE_SUN, CALL_LIKE_HPUX,
             CALL_LIKE_MVS )  'fpstr' is identical to 'pstr'.
             You may choose arbitrary name for the string,
             but you MUST use the same name, appended with
             suffix '_len', for its length (see example below).
  
          b) whenever string and its length are passed as
             complex parameter, 'fpstr' is identical to the
             pointer on the corresponding structure:
               CALL_LIKE_STARDENT :
                   'fpstr' is identical to 'PSStrPar'
               CALL_LIKE_VMS      :
                   'fpstr' is identical to 'dsc$descriptor_s *'
  
        With 'fpstr' type, two important macro definition come:
  
          i)  FTN_STR(s)  - returns pointer to fortran-passed
                            string s. This pointer is always
                            of 'pstr' type
          ii) FTN_LEN(s)  - returns integer length of fortran-
                            passed string s. For this macro to
                            work properly with SUN- and MVS-like
                            machines, always use suffix '_len' 
                            for the string length parameters as
                            described in a) above.
  
     2. Three parameter lists, each enclosed in brackets, should
        be given. These lists retain the general order of
        parameters in the corresponding fortran call. Non-string
        parameters are passed as pointers. String parameters
        and their lengths are passed differently in different
        lists:
  
         p_send    strings enter their place in the list as in
                   the corresponding FORTRAN call, having 'fpstr'
                   parameter type. Their lengths are appended as
                   'int' to the end of the list. They should
                   retain the order in which the strings appear
                   in the list.
  
         p_sstruct strings enter their place in the list as in
                   the corresponding FORTRAN call, having 'fpstr'
                   parameter type.
  
         p_sflw    strings enter their place in the list as in
                   the corresponding FORTRAN call, having 'fpstr'
                   type and being immediately followed by their
                   lengths as 'int' parameters.
  
  
  
   Example:
  
     FORTRAN statement
  
       subroutine  SomeSub ( k,s1,a,s2,m )
       integer       k,m
       real          a
       character*(*) s1,s2
  
     is translated to
  
       FORTRAN_SUBR ( SOMESUB, somesub,
         ( int * k, fpstr s1, float * a, fpstr s2, int * m,
           int s1_len, int s2_len ),
         ( int * k, fpstr s1, float * a, fpstr s2, int * m ),
         ( int * k, fpstr s1, int s1_len, float * a,
           fpstr s2, int s2_len, int * m ) )
  
  
     The macro should replace ordinary function header
   statements to assure compatibility with FORTRAN links.
   In header files, do not forget to add semicolumn:
  
     FORTRAN_SUBR ( .... );
  
   while in source files use simply
  
     FORTRAN_SUBR ( .... )  {
      <source body, operators>
     }
  
  
  
     Macro  FORTRAN_CALL(NAME,name,p_send,p_sstruct,p_sflw)
   calls function defined with macro FORTRAN_SUBR(...), from
   a C/C++ application. Its parameters and their meaning are
   exactly identical to those of FORTRAN_SUBR(...).
   FORTRAN_CALL(...) should be followed by semicolon.                    */


#if  defined(CALL_LIKE_SUN)

  typedef pstr fpstr;

#define FTN_STR(s)  s
#define FTN_LEN(s)  s##_len

#define char_struct(s)           \
    pstr  s;                       \
    int   s##_len;
#define fill_char_struct(s,str)  \
    s  = str;                      \
    s##_len = strlen(str);
#define init_char_struct(s,str,size)  \
    s  = str;                      \
    s##_len = size;

#define FORTRAN_SUBR(NAME,name,p_sun,p_stardent,p_mvs) \
    void name##_ p_sun
#define FORTRAN_CALL(NAME,name,p_sun,p_stardent,p_mvs) \
    name##_ p_sun
#define FORTRAN_FUN(val,NAME,name,p_sun,p_stardent,p_mvs) \
    val name##_ p_sun
#elif defined(CALL_LIKE_HPUX)

  typedef pstr fpstr;

#  define FTN_STR(s)  s
#  define FTN_LEN(s)  s##_len

#  define char_struct(s)  \
    pstr  s;              \
    int   s##_len;
#  define fill_char_struct(s,str)  \
    s  = str;                      \
    s##_len = strlen(str);
#  define init_char_struct(s,str,size)  \
    s  = str;                      \
    s##_len = size;

#  define FORTRAN_SUBR(NAME,name,p_sun,p_stardent,p_mvs) \
    void name p_sun
#  define FORTRAN_CALL(NAME,name,p_sun,p_stardent,p_mvs) \
    name p_sun
#  define FORTRAN_FUN(val,NAME,name,p_sun,p_stardent,p_mvs) \
    val name p_sun
#elif defined(CALL_LIKE_STARDENT)

  typedef PStrPar fpstr;

#  define FTN_STR(s)  s->Str_pointer
#  define FTN_LEN(s)  s->Str_length

#  define char_struct(s)           \
    SStrPar s;
#  define fill_char_struct(s,str)  \
    s.S   = str;                   \
    s.len = strlen(FName);         \
    s.id  = 0;
#  define init_char_struct(s,str,size)  \
    s.S   = str;                   \
    s.len = size;         \
    s.id  = 0;

#  define FORTRAN_SUBR(NAME,name,p_send,p_sstruct,p_sflw) \
    void NAME p_stardent
#  define FORTRAN_CALL(NAME,name,p_send,p_sstruct,p_sflw) \
    NAME p_stardent
#  define FORTRAN_FUN(val,NAME,name,p_send,p_sstruct,p_sflw) \
    val NAME p_stardent

#elif defined(CALL_LIKE_VMS)

  typedef dsc$descriptor_s * fpstr;

#  define FTN_STR(s)  s->dsc$a_pointer;
#  define FTN_LEN(s)  s->dsc$w_length;

#  define char_struct(s)                \
    dsc$descriptor_s s;
#  define fill_char_struct(s,str)     \
    s.dsc$a_pointer = str;            \
    s.dsc$w_length  = strlen(str);    \
    s.dsc$b_dtype   = DSC$K_DTYPE_T;  \
    s.dsc$b_class   = DSC$K_CLASS_S;
#  define init_char_struct(s,str,size)     \
    s.dsc$a_pointer = str;            \
    s.dsc$w_length  = size;    \
    s.dsc$b_dtype   = DSC$K_DTYPE_T;  \
    s.dsc$b_class   = DSC$K_CLASS_S;

#  define FORTRAN_SUBR(NAME,name,p_sun,p_stardent,p_mvs) \
    void NAME p_stardent
#  define FORTRAN_CALL(NAME,name,p_sun,p_stardent,p_mvs) \
    NAME p_stardent
#  define FORTRAN_FUN(val,NAME,name,p_sun,p_stardent,p_mvs) \
    val NAME p_stardent

#elif defined(CALL_LIKE_MVS)

#if (CALL_LIKE_MVS == 2)

  typedef pstr fpstr;

#define FTN_STR(s)  s
#define FTN_LEN(s)  s##_len

#define char_struct(s)           \
    pstr  s;                       \
    int   s##_len;
#define fill_char_struct(s,str)  \
    s  = str;          \
	s##_len = strlen(str);
#define init_char_struct(s,str,size)  \
    s  = str;			\
	s##_len = size;

#define FORTRAN_SUBR(NAME,name,p_sun,p_stardent,p_mvs) \
    void NAME p_sun
#define FORTRAN_CALL(NAME,name,p_sun,p_stardent,p_mvs) \
    NAME p_sun
#define FORTRAN_FUN(val,NAME,name,p_sun,p_stardent,p_mvs) \
    val NAME p_sun

#else

  typedef pstr fpstr;

#  define FTN_STR(s)  s
#  define FTN_LEN(s)  s##_len

#  define char_struct(s)  \
    pstr  s;              \
    int   s##_len;
#  define fill_char_struct(s,str)  \
    s  = str;                      \
    s##_len = strlen(str);
#  define init_char_struct(s,str,size)  \
    s  = str;                      \
    s##_len = size;

#  define FORTRAN_SUBR(NAME,name,p_sun,p_stardent,p_mvs) \
    void __stdcall NAME p_mvs
#  define FORTRAN_CALL(NAME,name,p_sun,p_stardent,p_mvs) \
    NAME p_mvs
#  define FORTRAN_FUN(val,NAME,name,p_sun,p_stardent,p_mvs) \
    val __stdcall NAME p_mvs

# endif

#else

#  error  Unknown machine!!!

  typedef pstr fpstr;

#  define FTN_STR(s)  s
#  define FTN_LEN(s)  s##_len

#  define char_struct(s)  \
    pstr  s;              \
    int   s##_len;
#  define fill_char_struct(s,str)  \
    s  = str;                      \
    s##_len = strlen(str);
#  define init_char_struct(s,str,size)  \
    s  = str;                      \
    s##_len = size;

/** Macro to define a function such that it is callable as
 * a Fortran subroutine.
 * @param NAME Subroutine name in upper case
 * @param name Subroutine name in lower case
 * @param p_sun Argument list in Sun style
 * @param p_stardent Argument list in Stardent style
 * @param p_mvs Argument list in MVS style
*/
#  define FORTRAN_SUBR(NAME,name,p_sun,p_stardent,p_mvs) \
    void name##_ p_sun

/** Macro to call a Fortran subroutine from a C function.
 * @param NAME Subroutine name in upper case
 * @param name Subroutine name in lower case
 * @param p_sun Argument list in Sun style
 * @param p_stardent Argument list in Stardent style
 * @param p_mvs Argument list in MVS style
*/
#  define FORTRAN_CALL(NAME,name,p_sun,p_stardent,p_mvs) \
    name##_ p_sun

/** Macro to define a function such that it is callable as
 * a Fortran function.
 * @param val Data type of return value.
 * @param NAME Function name in upper case
 * @param name Function name in lower case
 * @param p_sun Argument list in Sun style
 * @param p_stardent Argument list in Stardent style
 * @param p_mvs Argument list in MVS style
*/
#  define FORTRAN_FUN(val,NAME,name,p_sun,p_stardent,p_mvs) \
    val name##_ p_sun

#endif

/* Define Fortran logical */
typedef unsigned int ftn_logical;
#define FORTRAN_LOGICAL_TRUE  1
#define FORTRAN_LOGICAL_FALSE 0
#if defined (__OSF1__) || defined (__osf__)
#  undef FORTRAN_LOGICAL_TRUE
#  define FORTRAN_LOGICAL_TRUE  -1
#endif     

char *ccp4_FtoCString(fpstr str1, int str1_len);
void ccp4_CtoFString(fpstr str1, int str1_len, const char *cstring);

#endif    /* __CCP4_FORTRAN */
