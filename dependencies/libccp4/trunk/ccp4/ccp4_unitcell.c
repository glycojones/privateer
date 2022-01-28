/*
     ccp4_unitcell.c: C library for manipulations based on cell parameters.
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

/** @file ccp4_unitcell.c
 *  C library for manipulations based on cell parameters.
 *  Martyn Winn 
 */

#include <stdio.h>

#include "ccp4_unitcell.h"
#include "cvecmat.h"
#include "ccp4_errno.h"
/* rcsid[] = "$Id$" */

/* from input cell and orthogonalisation code, find orthogonalisation
   and fractionalisation matrices. Returns cell volume. */

double ccp4uc_frac_orth_mat(const double cell[6], const int ncode, 
			 double ro[3][3], double rf[3][3])
{
  int i,j;
  double conv,alph,bet,gamm,sina,cosa,sinb,cosb,sing,cosg,
    sinas,cosas,sinbs,cosbs,sings,cosgs,a,b,c;

  conv = atan(1.0)*4.0/180.0;
  alph = cell[3]*conv;
  bet = cell[4]*conv;
  gamm = cell[5]*conv;
  sina = sin(alph);
  cosa = cos(alph);
  sinb = sin(bet);
  cosb = cos(bet);
  sing = sin(gamm);
  cosg = cos(gamm);
  cosas = (cosg*cosb-cosa)/ (sinb*sing);
  sinas = sqrt(1.0-cosas*cosas);
  cosbs = (cosa*cosg-cosb)/ (sina*sing);
  sinbs = sqrt(1.0-cosbs*cosbs);
  cosgs = (cosa*cosb-cosg)/ (sina*sinb);
  sings = sqrt(1.0-cosgs*cosgs);
  a = cell[0];
  b = cell[1];
  c = cell[2];
	
  /* calculate ro */
  for ( i = 0; i < 3; i++ ) 
    for ( j = 0; j < 3; j++ ) 
      ro[i][j] = 0.0;

  /* ncode 1 -  xo along a  zo along c* */

  switch (ncode) {
  case 1:
    ro[0][0] = a;
    ro[0][1] = b*cosg;
    ro[0][2] = c*cosb;
    ro[1][1] = b*sing;
    ro[1][2] = -c*sinb*cosas;
    ro[2][2] = c*sinb*sinas;
    break;

 /* ncode 2 -  xo along b  zo along a* */

  case 2:
    ro[0][0] = a*cosg;
    ro[0][1] = b;
    ro[0][2] = c*cosa;
    ro[1][0] = -a*sing*cosbs;
    ro[1][2] = c*sina;
    ro[2][0] = a*sing*sinbs;
    break;

 /* ncode 3 -  xo along c  zo along b* */

  case 3:
    ro[0][0] = a*cosb;
    ro[0][1] = b*cosa;
    ro[0][2] = c;
    ro[1][0] = a*sinb;
    ro[1][1] = -b*sina*cosgs;
    ro[2][1] = b*sina*sings;
    break;

 /* ncode 4 -  trigonal only - xo along a+b  yo alon a-b  zo along c* */

  case 4:
    ro[0][0] = a/2.0;
    ro[0][1] = a/2.0;
    ro[1][0] = -a*sing;
    ro[1][1] = a*sing;
    ro[2][2] = c;
    break;

 /* ncode 5 -  xo along a*   zo along c */

  case 5:
    ro[0][0] = a*sinb*sings;
    ro[1][0] = -a*sinb*cosgs;
    ro[1][1] = b*sina;
    ro[2][0] = a*cosb;
    ro[2][1] = b*cosa;
    ro[2][2] = c;
    break;

 /* ncode 6 -  grr*! to  gerard bricogne - his setting for p1 in skew.
     xo along a  yo along b* */

  case 6:
    ro[0][0] = a;
    ro[0][1] = b*cosg;
    ro[0][2] = c*cosb;
    ro[1][1] = b*sing*sinas;
    ro[2][1] = -b*sing*cosas;
    ro[2][2] = c*sinb;
    break;
  }

   /* now calculate rf from ro, determinant gives cell volume */

  return invert3matrix((const double (*)[3]) ro, rf);

}

/* from input cell, find dimensions of reciprocal cell. 
   Returns reciprocal cell volume. */

double ccp4uc_calc_rcell(const double cell[6], double rcell[6])
{
  double conv,alph,bet,gamm,vol,sina,cosa,sinb,cosb,sing,cosg,
    sinas,cosas,sinbs,cosbs,sings,cosgs,a,b,c;

  conv = 3.14159/180.0;
  alph = cell[3]*conv;
  bet = cell[4]*conv;
  gamm = cell[5]*conv;
  vol = ccp4uc_calc_cell_volume(cell);
  sina = sin(alph);
  cosa = cos(alph);
  sinb = sin(bet);
  cosb = cos(bet);
  sing = sin(gamm);
  cosg = cos(gamm);
  cosas = (cosg*cosb-cosa)/ (sinb*sing);
  sinas = sqrt(1.0-cosas*cosas);
  cosbs = (cosa*cosg-cosb)/ (sina*sing);
  sinbs = sqrt(1.0-cosbs*cosbs);
  cosgs = (cosa*cosb-cosg)/ (sina*sinb);
  sings = sqrt(1.0-cosgs*cosgs);
  a = cell[0];
  b = cell[1];
  c = cell[2];
  rcell[0] = b*c*sina/vol;
  rcell[1] = c*a*sinb/vol;
  rcell[2] = a*b*sing/vol;
  rcell[3] = atan2(sinas,cosas)/conv;
  rcell[4] = atan2(sinbs,cosbs)/conv;
  rcell[5] = atan2(sings,cosgs)/conv;

  return (1.0 / vol);
}

/* Convert orthogonal to fractional coordinates. Translation only if
   deliberate origin shift - does this ever happen? Leave it to the
   application. */

void ccp4uc_orth_to_frac(const double rf[3][3], const double xo[3], double xf[3])
{
  xf[0] = rf[0][0]*xo[0] + rf[0][1]*xo[1] + rf[0][2]*xo[2];
  xf[1] = rf[1][0]*xo[0] + rf[1][1]*xo[1] + rf[1][2]*xo[2];
  xf[2] = rf[2][0]*xo[0] + rf[2][1]*xo[1] + rf[2][2]*xo[2];
}

/* Convert fractional to orthogonal coordinates. */

void ccp4uc_frac_to_orth(const double ro[3][3], const double xf[3], double xo[3])
{
  xo[0] = ro[0][0]*xf[0] + ro[0][1]*xf[1] + ro[0][2]*xf[2];
  xo[1] = ro[1][0]*xf[0] + ro[1][1]*xf[1] + ro[1][2]*xf[2];
  xo[2] = ro[2][0]*xf[0] + ro[2][1]*xf[1] + ro[2][2]*xf[2];
}

/* Convert orthogonal to fractional u matrix. */

void ccp4uc_orthu_to_fracu(const double rf[3][3], const double uo[6], double uf[6])
{
  int i,j;
  double uomat[3][3], ufmat[3][3], rft[3][3], temp[3][3];

  uomat[0][0] = uo[0]; uomat[0][1] = uo[3]; uomat[0][2] = uo[4]; 
  uomat[1][0] = uo[3]; uomat[1][1] = uo[1]; uomat[1][2] = uo[5]; 
  uomat[2][0] = uo[4]; uomat[2][1] = uo[5]; uomat[2][2] = uo[2]; 
  for ( i = 0; i < 3; i++ ) 
    for ( j = 0; j < 3; j++ ) 
      rft[i][j] = rf[j][i];

  ccp4_3matmul(temp,(const double (*)[3]) uomat,(const double (*)[3]) rft);
  ccp4_3matmul(ufmat,rf,(const double (*)[3]) temp);

  uf[0] = ufmat[0][0]; uf[1] = ufmat[1][1]; uf[2] = ufmat[2][2]; 
  uf[3] = ufmat[0][1]; uf[4] = ufmat[0][2]; uf[5] = ufmat[1][2]; 
}

/* Convert fractional to orthogonal u matrix. */

void ccp4uc_fracu_to_orthu(const double ro[3][3], const double uf[6], double uo[6])
{
  int i,j;
  double uomat[3][3], ufmat[3][3], rot[3][3], temp[3][3];

  ufmat[0][0] = uf[0]; ufmat[0][1] = uf[3]; ufmat[0][2] = uf[4]; 
  ufmat[1][0] = uf[3]; ufmat[1][1] = uf[1]; ufmat[1][2] = uf[5]; 
  ufmat[2][0] = uf[4]; ufmat[2][1] = uf[5]; ufmat[2][2] = uf[2]; 
  for ( i = 0; i < 3; i++ ) 
    for ( j = 0; j < 3; j++ ) 
      rot[i][j] = ro[j][i];

  ccp4_3matmul(temp,(const double (*)[3]) ufmat,(const double (*)[3]) rot);
  ccp4_3matmul(uomat,ro,(const double (*)[3]) temp);

  uo[0] = uomat[0][0]; uo[1] = uomat[1][1]; uo[2] = uomat[2][2]; 
  uo[3] = uomat[0][1]; uo[4] = uomat[0][2]; uo[5] = uomat[1][2]; 
}

/* Calculate cell volume from cell parameters */

double ccp4uc_calc_cell_volume(const double cell[6])
{
  double conv,alph,bet,gamm,sum,v;

  conv = 3.14159/180.0;
  alph = cell[3]*conv;
  bet = cell[4]*conv;
  gamm = cell[5]*conv;
  sum = (alph+bet+gamm)*0.5;
  v = sqrt(sin(sum-alph)*sin(sum-bet)*sin(sum-gamm)*sin(sum));
  return (2.0*cell[0]*cell[1]*cell[2]*v);
}

/* Check cells agree within tolerance */

int ccp4uc_cells_differ(const double cell1[6], const double cell2[6], const double tolerance)
{
  int i;
  double vol1, vol2, acheck;

  vol1 = ccp4uc_calc_cell_volume(cell1);  
  vol2 = ccp4uc_calc_cell_volume(cell2);

  /* check cell volumes */
  acheck = fabs(0.5*(vol1 - vol2))/(vol1 + vol2);
  if (acheck > tolerance) {
    if (ccp4_liberr_verbosity(-1)) {
      printf("Difference in cell volumes detected.\n");
      printf("  vol1 = %lf  vol2 = %lf \n",vol1,vol2);
    }
    return 1;
  }

  /* check cell parameters */
  acheck = 0.0;
  for ( i = 0; i < 6; i++ ) 
    acheck += fabs(0.5*(cell2[i]-cell1[i]))/(cell2[i]+cell1[i]);
  if (acheck > 3.0*tolerance) {
    if (ccp4_liberr_verbosity(-1)) {
      printf("Large difference in cell parameters detected.\n");
      printf("  cell1 = %lf %lf %lf %lf %lf %lf \n",
	   cell1[0],cell1[1],cell1[2],cell1[3],cell1[4],cell1[5]);
      printf("  cell2 = %lf %lf %lf %lf %lf %lf \n",
	   cell2[0],cell2[1],cell2[2],cell2[3],cell2[4],cell2[5]);
    }
    return 1;
  } else if (acheck > tolerance) {
    if (ccp4_liberr_verbosity(-1)) {
      printf("Small difference in cell parameters detected.\n");
      printf("  cell1 = %lf %lf %lf %lf %lf %lf \n",
	   cell1[0],cell1[1],cell1[2],cell1[3],cell1[4],cell1[5]);
      printf("  cell2 = %lf %lf %lf %lf %lf %lf \n",
	   cell2[0],cell2[1],cell2[2],cell2[3],cell2[4],cell2[5]);
    }
    return 1;
  }
  return 0;
}

int ccp4uc_is_rhombohedral(const float cell[6], const float tolerance) {

  double acheck;

  acheck = fabs(cell[0]-cell[1]);
  acheck += fabs(cell[1]-cell[2]);
  acheck += fabs(cell[0]-cell[2]);
  acheck += fabs(cell[3]-cell[4]);
  acheck += fabs(cell[3]-cell[5]);
  acheck += fabs(cell[4]-cell[5]);
  if (acheck > (double) tolerance) return 0;
  return 1;

}

int ccp4uc_is_hexagonal(const float cell[6], const float tolerance) {

  double acheck;

  acheck = fabs(cell[0]-cell[1]);
  acheck += fabs(cell[3]-90.0);
  acheck += fabs(cell[4]-90.0);
  acheck += fabs(cell[5]-120.0);
  if (acheck > (double) tolerance) return 0;
  return 1;

}
