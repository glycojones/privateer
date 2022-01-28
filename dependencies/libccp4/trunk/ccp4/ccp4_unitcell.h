/*
     ccp4_unitcell.h: headers for C library for ccp4_unitcell.c
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

/** @file ccp4_unitcell.h
 *  C library for manipulations based on cell parameters.
 *  Martyn Winn 
 */

#ifndef __CCP4_UNITCELL
#define __CCP4_UNITCELL

#ifdef  __cplusplus
namespace CCP4uc {
extern "C" {
#endif

#include <math.h>

/** From input cell and orthogonalisation code, find orthogonalisation
   and fractionalisation matrices.
 * @param cell
 * @param ncode
 * @param ro
 * @param rf
 * @return Cell volume
 */
double ccp4uc_frac_orth_mat(const double cell[6], const int ncode, 
			   double ro[3][3], double rf[3][3]);

/** From input cell, find dimensions of reciprocal cell.
 * @param cell
 * @param rcell
 * @return Reciprocal cell volume
 */
double ccp4uc_calc_rcell(const double cell[6], double rcell[6]);

/** Convert orthogonal to fractional coordinates. Translation only if
   deliberate origin shift - does this ever happen? Leave it to the
   application.
 * @param rf
 * @param xo
 * @param xf
 * @return void
 */
void ccp4uc_orth_to_frac(const double rf[3][3], const double xo[3], double xf[3]);

/** Convert fractional to orthogonal coordinates.
 * @param ro
 * @param xf
 * @param xo
 * @return void
 */
void ccp4uc_frac_to_orth(const double ro[3][3], const double xf[3], double xo[3]);

/** Convert orthogonal to fractional u matrix.
 * @param rf
 * @param uo
 * @param uf
 * @return void
 */
void ccp4uc_orthu_to_fracu(const double rf[3][3], const double uo[6], double uf[6]);

/** Convert fractional to orthogonal u matrix.
 * @param ro
 * @param uf
 * @param uo
 * @return void
 */
void ccp4uc_fracu_to_orthu(const double ro[3][3], const double uf[6], double uo[6]);

/** Calculate cell volume from cell parameters.
 * @param cell
 * @return Cell volume.
 */
double ccp4uc_calc_cell_volume(const double cell[6]);

/** Check cells agree within tolerance.
 * @param cell1 First cell.
 * @param cell2 Second cell.
 * @param tolerance A tolerance for agreement.
 * @return 1 if cells differ by more than tolerance, 0 otherwise.
 */
int ccp4uc_cells_differ(const double cell1[6], const double cell2[6], const double tolerance);

/** Check if cell parameters conform to a rhombohedral setting.
 * @param cell Cell parameters. Angles are assumed to be in degrees.
 * @param tolerance A tolerance for agreement.
 * @return 1 if cell parameters conform, 0 otherwise.
 */
int ccp4uc_is_rhombohedral(const float cell[6], const float tolerance);

/** Check if cell parameters conform to a hexagonal setting.
 * @param cell Cell parameters. Angles are assumed to be in degrees.
 * @param tolerance A tolerance for agreement.
 * @return 1 if cell parameters conform, 0 otherwise. 
 */
int ccp4uc_is_hexagonal(const float cell[6], const float tolerance);

#ifdef __cplusplus
} }
#endif

#endif  /*!CCP4_UNITCELL */
