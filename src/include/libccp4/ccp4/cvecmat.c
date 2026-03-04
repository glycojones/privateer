/*
     cvecmat.c: C library for vector and matrix manipulations
     Copyright (C) 2001  CCLRC, Martyn Winn

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

/** @file cvecmat.c
 *  C library for vector and matrix manipulations.
 *  Martyn Winn 
 */

#include <math.h>
#include "cvecmat.h"
/* rcsid[] = "$Id$" */

/*  c = a X b  */

void ccp4_dcross(const double a[3], const double b[3], double c[3])
{
  c[0] = a[1]*b[2] - b[1]*a[2];
  c[1] = a[2]*b[0] - b[2]*a[0];
  c[2] = a[0]*b[1] - b[0]*a[1];
}

void ccp4_3matmul(double c[3][3], const double a[3][3], const double b[3][3])
{
  int i,j,k;

  for ( i = 0; i < 3; i++ ) 
    for ( j = 0; j < 3; j++ ) {
      c[i][j] = 0.0;
      for ( k = 0; k < 3; k++ ) 
        c[i][j] += a[i][k]*b[k][j];
    }
}

void ccp4_4matmul( float c[4][4], const float  a[4][4], const float b[4][4])
{
  int i,j,k;

  for ( i = 0; i < 4; i++ ) 
    for ( j = 0; j < 4; j++ ) {
      c[i][j] = 0.0;
      for ( k = 0; k < 4; k++ ) 
        c[i][j] += a[i][k]*b[k][j];
    }
}

/*   A (I)   3*3 matrix to be inverted */
/*   AI (O)   inverse matrix */
/*   returns determinant */

double invert3matrix(const double a[3][3], double ai[3][3])

{ int i,j;
  double c[3][3],d;

  ccp4_dcross(a[1],a[2],c[0]);
  ccp4_dcross(a[2],a[0],c[1]);
  ccp4_dcross(a[0],a[1],c[2]);

  d = a[0][0]*c[0][0] + a[0][1]*c[0][1] + a[0][2]*c[0][2];

  if (fabs(d) > 1.0e-30) {
    for ( i = 0; i < 3; i++ ) 
      for ( j = 0; j < 3; j++ ) 
        ai[i][j] = c[j][i] / d;
  } else {
    return 0.0;
  }
  return d;
}

/*   A (I)   4*4 matrix to be inverted */
/*   AI (O)   inverse matrix */
/*   returns determinant */

float invert4matrix(const float a[4][4], float ai[4][4])

{
    double c[4][4], d;
    int i, j;
    double x[3][3];
    int i1, j1, i2 ;
    double am, q;
    int ii, jj;

    /* Function Body */
    for (ii = 0; ii < 4; ++ii) {
	for (jj = 0; jj < 4; ++jj) {
            ai[ii][jj] = 0.0;
	    i = -1;
	    for (i1 = 0; i1 < 4; ++i1) {
		if (i1 != ii) {
		    ++i;
		    j = -1;
		    for (j1 = 0; j1 < 4; ++j1) {
			if (j1 != jj) {
			    ++j;
			    x[i][j] = a[i1][j1];
			}
		    }
		}
	    }

	    am = x[0][0]*x[1][1]*x[2][2] - x[0][0]*x[1][2]*x[2][1] +
                 x[0][1]*x[1][2]*x[2][0] - x[0][1]*x[1][0]*x[2][2] +
                 x[0][2]*x[1][0]*x[2][1] - x[0][2]*x[1][1]*x[2][0];
	    i2 = ii + jj;
	    c[ii][jj] = ccp4_pow_ii(-1.0, i2) * am;
	}
    }

/* ---- Calculate determinant */

    d = 0.0;

    for (i = 0; i < 4; ++i) {
	d = a[i][0] * c[i][0] + d;
    }

/* ---- Get inverse matrix */


  if (fabs(d) > 1.0e-30) {
    q = 1.0/d;
    for (i = 0; i < 4; ++i) {
	for (j = 0; j < 4; ++j) {
	  ai[i][j] = (float) (c[j][i] * q);
	}
    }
  } else {
    return 0.0;
  }

  return ((float) d);
} 

float ccp4_pow_ii(const float base, const int power) {

  int i = 0;
  float pow = 1;

  while (++i <= power)
    pow *= base;

  return pow;
}

