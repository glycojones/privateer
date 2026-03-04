/*
     ccp4_types.h: CCP4 library.c macro definitions etc
     Copyright (C) 2001  CCLRC

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
#ifndef __CCP4_TYPES
#define __CCP4_TYPES

#include "ccp4_sysdep.h"

typedef unsigned short uint16;
#ifdef SIXTEENBIT
typedef unsigned long uint32;
#else
typedef unsigned int uint32;
#endif
typedef float float32;
typedef unsigned char uint8;
union float_uint_uchar {
    float32 f;
    uint32 i;
    uint8 c[4];
  };

typedef   char     *        pstr;

/* CCP4 library.c macro definitions */

#ifndef FALSE
#define FALSE 0
#define TRUE 1
#endif

typedef struct { double r;             /* radial and */
                 double phi;           /* angular component of */
               } POLAR;                /* a complex number */

/* some simple macros, which may exist anyway */
#ifndef SQR
#define SQR(x) ((x)*(x))
#endif
#ifndef DEGREE
#define DEGREE(x) ((((x < 0)?(x)+2*M_PI:(x))*360)/(2*M_PI))
#endif
#ifndef RADIAN
#define RADIAN(x) ((((x<0)?(x)+360:(x))*2*M_PI)/360)
#endif
#ifndef MAX
#define MAX(x, y) (((x)>(y))?(x):(y))
#endif
#ifndef MIN
#define MIN(x, y) (((x)<(y))?(x):(y))
#endif
#ifndef ABS
#define ABS(x) (((x)<0)?-(x):(x))
#endif
#ifndef SIGN
#define SIGN(x) (((x)<0)?-1:1)
#endif

#endif   /* __CCP4_TYPES */
