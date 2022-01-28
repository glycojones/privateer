/*
     ccp4_spg.h: Data structure for symmetry information
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

/** @file ccp4_spg.h
 *
 *  @brief Data structure for symmetry information.
 *
 *  A data structure for spacegroup information and related
 *  quantities. Some items are loaded from syminfo.lib while
 *  others are calculated on the fly.
 *
 *  There is also a struct for symop rotation matrices and translation
 *  vectors. This can be converted to other representations using
 *  functions in ccp4_parser. 
 *
 *  @author Martyn Winn 
 */

#ifndef __CCP4_SPG__
#define __CCP4_SPG__
/* rcsidhsp[] = "$Id$" */

#ifdef  __cplusplus
namespace CSym {
extern "C" {
#endif

/* Kevin's symmetry operator */

typedef struct ccp4_symop_
{
  float rot[3][3];
  float trn[3];
} ccp4_symop;

typedef struct ccp4_spacegroup_
{
  int spg_num;            /* true spacegroup number */
  int spg_ccp4_num;       /* CCP4 spacegroup number */
  char symbol_Hall[40];   /* Hall symbol */
  char symbol_xHM[20];    /* Extended Hermann Mauguin symbol  */
  char symbol_old[20];    /* old spacegroup name */

  char point_group[20];   /* point group name */
  char crystal[20];       /* crystal system e.g. MONOCLINIC */

  int nlaue;              /* CCP4 Laue class number, inferred from asu_descr */
  char laue_name[20];     /* Laue class name */
  int laue_sampling[3];   /* sampling factors for FFT */

  int npatt;              /* Patterson spacegroup number, inferred from asu_descr */
  char patt_name[40];     /* Patterson spacegroup name */

  int nsymop;             /* total number of symmetry operations */
  int nsymop_prim;        /* number of primitive symmetry operations */
  ccp4_symop *symop;      /* symmetry matrices */
  ccp4_symop *invsymop;   /* inverse symmetry matrices */

  float chb[3][3];        /* change of basis matrix from file */

  char asu_descr[80];     /* asu description from file */
  int (*asufn)(const int, const int, const int); /* pointer to ASU function */

  int centrics[12];       /* symop which generates centric zone, 0 if none */
  int epsilon[13];        /* flag which epsilon zones are applicable */

  char mapasu_zero_descr[80];  /* origin-based map asu: description from file */
  float mapasu_zero[3];   /* origin-based map asu: upper limits */

  char mapasu_ccp4_descr[80];  /* CCP4 map asu: defaults to mapasu_zero */
  float mapasu_ccp4[3];   /* CCP4 map asu: upper limits */

} CCP4SPG;

#ifdef __cplusplus
} }
#endif

#endif  /*!__CCP4_SPG__ */

