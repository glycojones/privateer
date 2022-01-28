/*
     csymlib.h: header file for csymlib.c
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

/** @page csym_page CSYM library
 *
 * @verbatim

<!-- ::INDEX_INFO::CSYM library::Library::::C/C++ Software Library for symmetry information:::::::: -->

   @endverbatim
 *
 *  @section csym_file_list File list

<ul>
<li>csymlib.h - contains details of the C/C++ API
<li>ccp4_spg.h - contains details of the spacegroup data structure
</ul>
 *
 *  @section csym_overview Overview
 
The CSYM library is centred around a data file <tt>syminfo.lib</tt> which is
auto-generated from sgtbx (the Space Group Toolbox of 
<a href="http://cctbx.sourceforge.net/">cctbx</a>). A description of
the contents of this file is given in the <a href="../symlib.html">
documentation</a> of the Fortran API.

<p>A particular spacegroup in a particular setting
is loaded into an in-memory data structure by requesting a particular
spacegroup name, number, or set of operators. See the functions
<tt>ccp4spg_load_by_standard_num</tt>, <tt>ccp4spg_load_by_ccp4_num</tt>, 
<tt>ccp4spg_load_by_spgname</tt>,  <tt>ccp4spg_load_by_ccp4_spgname</tt>
and <tt>ccp4_spgrp_reverse_lookup</tt>. Information on the in-memory 
data structure is given in ccp4_spg.h The memory can be freed by the
function <tt>ccp4spg_free</tt>.

<p>Functions are provided to:
<ul>
<li>Query the data structure, e.g. <tt>ccp4spg_symbol_Hall</tt>, etc. (members
of the structure can of course be obtained directly)
<li>Check reciprocal space indices for a particular spacegroup, 
e.g. <tt>ccp4spg_is_in_asu</tt>, <tt>ccp4spg_is_centric</tt>, 
<tt>ccp4spg_get_multiplicity</tt>, <tt>ccp4spg_is_sysabs</tt>, etc.
<li>Set appropriate grids for FFT, e.g. <tt>set_fft_grid</tt>
</ul>
 *
 *  @section csym_operators Symmetry operators

Symmetry operators are expressed in a variety of ways:
<ul>
<li>Using the struct <tt>ccp4_symop</tt>, which consists of a 3 x 3 rotation
matrix and a translation vector.
<li>As a 4 x 4 matrix, in which the rotation matrix is in the top-left-hand
corner and the translation vector is in elements [*][3]. Element [3][3] is
set to 1.0
<li>As a string, such as "-x+1/2,-y,z+1/2"
</ul>
Check the function description for which form is expected. Often, there 
are alternative functions if you wish to supply the operators in a
different form. There are also the following conversion functions:
<ul>
<li>rotandtrn_to_mat4
<li>rotandtrn_to_symop
<li>mat4_to_rotandtrn
<li>mat4_to_symop
<li>mat4_to_recip_symop
<li>symop_to_rotandtrn
<li>symop_to_mat4
</ul>
Note that the order of symmetry operators may be important in some cases, for
example in MTZ files with a M/ISYM column where ISYM encodes the symmetry operation
used.

 *  @section csym_examples Examples

See examples on <a href="ftp://ftp.ccp4.ac.uk/mdw/csym">ftp area</a>

*/

/** @file csymlib.h
 *
 *  @brief C-level library for symmetry information.
 *
 *  Functions defining the C-level API for accessing spacegroup properties.
 *  The primary spacegroup information comes from the data file syminfo.lib
 *
 *  @author Martyn Winn 
 */

#ifndef __CSymLib__
#define __CSymLib__

/* rcsidhs[] = "$Id$" */

/* note that definitions in ccp4_spg.h are within the CSym namespace */
#include "ccp4_spg.h"

#ifdef  __cplusplus
namespace CSym {
extern "C" {
#endif

/** Look up spacegroup in standard setting by number, and load properties.
 * Allocates memory for the spacegroup structure. This can be freed
 * later by ccp4spg_free().
 * @param numspg spacegroup number
 * @return pointer to spacegroup
 */
CCP4SPG *ccp4spg_load_by_standard_num(const int numspg); 

/** Look up spacegroup by CCP4 number, and load properties.
 * Allocates memory for the spacegroup structure. This can be freed
 * later by ccp4spg_free().
 * @param ccp4numspg CCP4 spacegroup number
 * @return pointer to spacegroup
 */
CCP4SPG *ccp4spg_load_by_ccp4_num(const int ccp4numspg); 

/** Look up spacegroup by the extended Hermann Mauguin symbol.
 * Allocates memory for the spacegroup structure. This can be freed
 * later by ccp4spg_free().
 * @param spgname Spacegroup name in form of extended Hermann Mauguin symbol.
 * @return pointer to spacegroup
 */
CCP4SPG *ccp4spg_load_by_spgname(const char *spgname);

/** Look up spacegroup by name. This is for use by CCP4 programs
 * and is more complicated than ccp4spg_load_by_spgname. For each
 * spacegroup in syminfo.lib it checks the CCP4 spacegroup name
 * first, and then the extended Hermann Mauguin symbol.
 * Allocates memory for the spacegroup structure. This can be freed
 * later by ccp4spg_free().
 * @param ccp4spgname Spacegroup name.
 * @return pointer to spacegroup
 */
CCP4SPG *ccp4spg_load_by_ccp4_spgname(const char *ccp4spgname);

/** Look up spacegroup by symmetry operators and load properties.
 * Allocates memory for the spacegroup structure. This can be freed
 * later by ccp4spg_free().
 * @param nsym1 number of operators (including non-primitive)
 * @param op1 pointer to array of operators
 * @return pointer to spacegroup
 */
CCP4SPG * ccp4_spgrp_reverse_lookup(const int nsym1, const ccp4_symop *op1);

/** Look up spacegroup from SYMOP.
 *  This would not normally be called directly, but via one of
 *  the wrapping functions. 
 *  Allocates memory for the spacegroup structure. This can be freed
 *  later by ccp4spg_free().
 * @param numspg spacegroup number
 * @param ccp4numspg CCP4 spacegroup number
 * @param spgname Spacegroup name.
 * @param ccp4spgname Spacegroup name.
 * @param nsym1 number of operators (including non-primitive)
 * @param op1 pointer to array of operators
 * @return pointer to spacegroup
 */
CCP4SPG *ccp4spg_load_spacegroup(const int numspg, const int ccp4numspg,
        const char *spgname, const char *ccp4spgname, 
        const int nsym1, const ccp4_symop *op1); 

/** Free all memory malloc'd from static pointers.
 * To be called before program exit. The function can be
 * registered with atexit.
 */
void ccp4spg_mem_tidy(void);

/** Generate symop matrices from description strings
 *  This would not normally be called directly, but via one of
 *  the wrapping functions SYMFR2 and SYMFR3 in the Fortran API.
 * @param line null-terminated string containing symop descriptions
 * @param rot array of 4x4 matrices
 * @return number of symops read, or -1 on failure
 */
int symfr_driver (const char *line, float rot[][4][4]);

/** Free memory associated with spacegroup.
 * @param sp pointer to spacegroup
 */
void ccp4spg_free(CCP4SPG **sp);

/** Look up spacegroup in standard setting by number and load into
 * static storage of csymlib_f.
 * @param numspg spacegroup number
 * @return void
 */
void ccp4spg_register_by_ccp4_num(int numspg);

/** Look up spacegroup by set of symmetry operators and load into
 * static storage of csymlib_f.
 * @param nops number of symops
 * @param rsm symmetry operators
 * @return void
 */
void ccp4spg_register_by_symops(int nops, float rsm[][4][4]);

/** Derive centering operators from Hall symbol (deprecated).
 * Centering operators are now read from syminfo.lib
 * @param symbol_Hall Hall symbol for spacegroup
 * @param cent_ops centering operators
 * @return number of centering operators (0 if none found)
 */
int ccp4_spg_get_centering(const char *symbol_Hall, float cent_ops[4][3]);

/** Load Laue data into spacegroup structure.
 * @param nlaue CCP4 code for Laue group
 * @param spacegroup Pointer to CCP4 spacegroup structure
 * @return 0 on success, 1 on failure to load Laue data
 */
int ccp4spg_load_laue(CCP4SPG* spacegroup, const int nlaue);

/** Test if reflection is in asu of Laue group 1bar.
 * @return 1 if in asu else 0
 */
int ASU_1b   (const int h, const int k, const int l);

/** Test if reflection is in asu of Laue group 2/m.
 * @return 1 if in asu else 0
 */
int ASU_2_m  (const int h, const int k, const int l);

/** Test if reflection is in asu of Laue group mmm.
 * @return 1 if in asu else 0
 */
int ASU_mmm  (const int h, const int k, const int l);

/** Test if reflection is in asu of Laue group 4/m.
 * @return 1 if in asu else 0
 */
int ASU_4_m  (const int h, const int k, const int l);

/** Test if reflection is in asu of Laue group 4/mmm.
 * @return 1 if in asu else 0
 */
int ASU_4_mmm(const int h, const int k, const int l);

/** Test if reflection is in asu of Laue group 3bar.
 * @return 1 if in asu else 0
 */
int ASU_3b   (const int h, const int k, const int l);

/** Test if reflection is in asu of Laue group 3bar1m.
 * @return 1 if in asu else 0
 */
int ASU_3bm  (const int h, const int k, const int l);

/** Test if reflection is in asu of Laue group 3barm.
 * @return 1 if in asu else 0
 */
int ASU_3bmx (const int h, const int k, const int l);

/** Test if reflection is in asu of Laue group 6/m.
 * @return 1 if in asu else 0
 */
int ASU_6_m  (const int h, const int k, const int l);

/** Test if reflection is in asu of Laue group 6/mmm.
 * @return 1 if in asu else 0
 */
int ASU_6_mmm(const int h, const int k, const int l);

/** Test if reflection is in asu of Laue group m3bar.
 * @return 1 if in asu else 0
 */
int ASU_m3b  (const int h, const int k, const int l);

/** Test if reflection is in asu of Laue group m3barm.
 * @return 1 if in asu else 0
 */
int ASU_m3bm  (const int h, const int k, const int l);

/** Function to return Hall symbol for spacegroup.
 * @param sp pointer to spacegroup
 * @return pointer to Hall symbol for spacegroup
 */
char *ccp4spg_symbol_Hall(CCP4SPG* sp);

/** inverts a symmetry operator. The input operator is
 * converted to a 4 x 4 matrix, inverted, and converted back.
 * @param ccp4_symop input symmetry operator
 * @return inverted symmetry operator
 */
ccp4_symop ccp4_symop_invert( const ccp4_symop op1 );

/** Compare two spacegroup names. Strings are converted to upper
 * case before making the comparison, but otherwise match must be
 * exact.
 * @param spgname1 First spacegroup name.
 * @param spgname2 Second spacegroup name.
 * @return 1 if they are equal else 0.
*/
int ccp4spg_name_equal(const char *spgname1, const char *spgname2);

/** Try to match a spacegroup name to one from SYMINFO. Blanks are 
 * removed when making the comparison. Strings are converted to upper
 * case before making the comparison. If spgname_lib has " 1 " and 
 * spgname_match doesn't, then strip out " 1" to do "short" comparison.
 * @param spgname1 First spacegroup name, assumed to be a standard one
 *  obtained at some point from SYMINFO
 * @param spgname2 Second spacegroup name that you are trying to match
 *  to a standard SYMINFO one. E.g. it might have been provided by the
 *  user.
 * @return 1 if they are equal else 0.
*/
int ccp4spg_name_equal_to_lib(const char *spgname_lib, const char *spgname_match);

/** Function to create "short" name of spacegroup. Blanks
 * are removed, as are " 1" elements (except for the special case
 * of "P 1").
 * @param shortname String long enough to hold short name.
 * @param longname Long version of spacegroup name.
 * @return Pointer to shortname.
*/
char *ccp4spg_to_shortname(char *shortname, const char *longname);

/** Function to deal with colon-specified spacegroup settings.
 * E.g. 'R 3 :H' is converted to 'H 3   '. Note that spaces are
 * returned and should be dealt with by the calling function.
 * @param name Spacegroup name.
 * @return void
*/
void ccp4spg_name_de_colon(char *name);

/** Compare two point group names. Blanks are removed when
 * making the comparison. Strings are converted to upper
 * case before making the comparison. Any initial "PG" is ignored.
 * @param pgname1 First point group name.
 * @param pgname2 Second point group name.
 * @return 1 if they are equal else 0.
*/
int ccp4spg_pgname_equal(const char *pgname1, const char *pgname2);

/** Function to normalise translations of a symmetry operator,
 * i.e. to ensure 0.0 <= op.trn[i] < 1.0.
 * @param op pointer to symmetry operator.
 * @return Pointer to normalised symmetry operator.
*/
ccp4_symop *ccp4spg_norm_trans(ccp4_symop *op);

/** Sort and compare two symmetry operator lists.
 * Kevin's code. The lists are coded as ints, which are then sorted and compared.
 * Note that no changes are made to the input operators, so that operators
 * differing by an integral number of unit cell translations are considered
 * unequal. If this is not what you want, normalise the operators with 
 * ccp4spg_norm_trans first.
 * @param nsym1 number of symmetry operators in first list
 * @param op1 first list of symmetry operators
 * @param nsym2 number of symmetry operators in second list
 * @param op2 second list of symmetry operators
 * @return 1 if they are equal else 0.
*/
int ccp4_spgrp_equal( int nsym1, const ccp4_symop *op1, int nsym2, const ccp4_symop *op2);

/** Compare two symmetry operator lists.
 * Kevin's code. The lists are coded as ints, which are compared.
 * Unlike ccp4_spgrp_equal, the lists are not sorted, so the same operators
 * in a different order will be considered unequal.
 * @param nsym1 number of symmetry operators in first list
 * @param op1 first list of symmetry operators
 * @param nsym2 number of symmetry operators in second list
 * @param op2 second list of symmetry operators
 * @return 1 if they are equal else 0.
*/
int ccp4_spgrp_equal_order( int nsym1, const ccp4_symop *op1, int nsym2, const ccp4_symop *op2);

/** Make an integer coding of a symmetry operator.
 * The coding takes 30 bits: 18 for the rotation and 12 for the translation.
 * @param op symmetry operator
 * @return int code.
 */
int ccp4_symop_code(ccp4_symop op);

/** Comparison of symmetry operators encoded as integers.
 * In ccp4_spgrp_equal, this is passed to the stdlib qsort.
 * @param p1 pointer to first integer
 * @param p1 pointer to second integer
 * @return difference between integers
*/
int ccp4_int_compare( const void *p1, const void *p2 );

/** Test whether reflection or it's Friedel mate is in asu.
 * @param sp pointer to spacegroup
 * @param h reflection index
 * @param k reflection index
 * @param l reflection index
 * @return 1 if in asu, -1 if -h -k -l is in asu, 0 otherwise
 */
int ccp4spg_is_in_pm_asu(const CCP4SPG* sp, const int h, const int k, const int l);

/** Test whether reflection is in asu.
 * @param sp pointer to spacegroup
 * @param h reflection index
 * @param k reflection index
 * @param l reflection index
 * @return 1 if in asu, 0 otherwise
 */
int ccp4spg_is_in_asu(const CCP4SPG* sp, const int h, const int k, const int l);

/** Place reflection (hin,kin,lin) in the asymmetric unit of spacegroup "sp".
 * Resultant indices are placed in (hout,kout,lout).
 * @param sp pointer to spacegroup
 * @param hin input reflection index
 * @param kin input reflection index
 * @param lin input reflection index
 * @param hout output reflection index
 * @param kout output reflection index
 * @param lout output reflection index
 * @return "isym" if successful, 0 otherwise. "isym" = 2*isymop - 1 for 
 * reflections placed in the positive asu, i.e. I+ of a Friedel pair, and
 * "isym" = 2*isymop for reflections placed in the negative asu, i.e. I- of 
 * a Friedel pair. Here "isymop" is the number of the symmetry operator used.
 */
int ccp4spg_put_in_asu(const CCP4SPG* sp, const int hin, const int kin, const int lin,
		       int *hout, int *kout, int *lout );

/** Transform reflection (hin,kin,lin) according to spacegroup "sp" and
 * operation "isym". Resultant indices are placed in (hout,kout,lout).
 * @param sp pointer to spacegroup
 * @param isym required operation, see ccp4spg_put_in_asu
 * @param hin input reflection index
 * @param kin input reflection index
 * @param lin input reflection index
 * @param hout output reflection index
 * @param kout output reflection index
 * @param lout output reflection index
 * @return void
 */
void ccp4spg_generate_indices(const CCP4SPG* sp, const int isym,
                  const int hin, const int kin, const int lin,
			      int *hout, int *kout, int *lout );

/** Shift phase value associated with hin,kin,lin according to translation 
and optional sign change. Return in range 0,360.
 * @param hin reflection index
 * @param kin reflection index
 * @param lin reflection index
 * @param phasin Input phase.
 * @param trans Requested translation
 * @param isign If -1, change sign of phase
 * @return shifted phase
 */
float ccp4spg_phase_shift(const int hin, const int kin, const int lin,
			  const float phasin, const float trans[3], const int isign);

/** Check whether change of basis is necessary, i.e. whether the
 * change of basis matrix is not the identity.
 * @param chb change of basis matrix
 * @return 1 if change of basis is necessary, 0 otherwise
 */
int ccp4spg_do_chb(const float chb[3][3]);

/** Set up centric zones for a given spacegroup. This is called
 * upon loading a spacegroup.
 * @param sp pointer to spacegroup
 * @return void
 */
void ccp4spg_set_centric_zones(CCP4SPG* sp);

/** Function to determine whether or not h,k,l is a centric reflection
 * in spacegroup "sp".
 * @param sp pointer to spacegroup
 * @param h input reflection index
 * @param k input reflection index
 * @param l input reflection index
 * @return 1 if h,k,l is centric, 0 if not centric, and -1 if there is
 *  an error.
 */
int ccp4spg_is_centric(const CCP4SPG* sp, const int h, const int k, const int l);

/** Check indices against a centric zone for a given spacegroup.
 * @param  nzone index of centric zone
 * @param h reflection index
 * @param k reflection index
 * @param l reflection index
 * @return 0 if in zone "nzone", non-zero otherwise
 */
int ccp4spg_check_centric_zone(const int nzone, const int h, const int k, const int l);

/** Return phase of a centric reflection in the range 0.0 <= phase < 180.0.
 * You should first check that reflection really is centric.
 * @param sp pointer to spacegroup
 * @param h reflection index
 * @param k reflection index
 * @param l reflection index
 * @return phase of a centric reflection
 */
float ccp4spg_centric_phase(const CCP4SPG* sp, const int h, const int k, const int l);

/** Print a summary of the centric zones of a spacegroup.
 * @param sp pointer to spacegroup
 * @return void
 */
void ccp4spg_print_centric_zones(const CCP4SPG* sp);

/** Obtain string description of centric zone.
 * @param nzone index of centric zone
 * @param centric_zone string description of centric zone
 * @return string description of centric zone
 */
char *ccp4spg_describe_centric_zone(const int nzone, char *centric_zone);

/** Set up epsilon zones for a given spacegroup. This is called
 * upon loading a spacegroup.
 * @param sp pointer to spacegroup
 * @return void
 */
void ccp4spg_set_epsilon_zones(CCP4SPG* sp);

/** Return reflection multiplicity factor for a given hkl in a given
 * spacegroup.
 * @param sp pointer to spacegroup
 * @param h reflection index
 * @param k reflection index
 * @param l reflection index
 * @return reflection multiplicity factor
 */
int ccp4spg_get_multiplicity(const CCP4SPG* sp, const int h, const int k, const int l);

/** Check indices against an epsilon zone for a given spacegroup.
 * @param nzone index of epsilon zone (runs from 1 to 13)
 * @param h reflection index
 * @param k reflection index
 * @param l reflection index
 * @return 0 if in zone "nzone", non-zero otherwise
 */
int ccp4spg_check_epsilon_zone(const int nzone, const int h, const int k, const int l);

/** Print a summary of the epsilon zones of a spacegroup.
 * @param sp pointer to spacegroup
 * @return void
 */
void ccp4spg_print_epsilon_zones(const CCP4SPG* sp);

/** Obtain string description of epsilon zone.
 * @param nzone index of epsilon zone
 * @param epsilon_zone string description of epsilon zone
 * @return string description of epsilon zone
 */
char *ccp4spg_describe_epsilon_zone(const int nzone, char *epsilon_zone);


/** Check if reflection is a systematic absence.
 * @param sp pointer to spacegroup
 * @param h reflection index
 * @param k reflection index
 * @param l reflection index
 * @return 1 if reflection is a systematic absence, 0 otherwise.
 */
int ccp4spg_is_sysabs(const CCP4SPG* sp, const int h, const int k, const int l);

/** Translated from Alexei Vagin's CALC_ORIG_PS.
 * @param namspg Spacegroup name for printing only.
 * @param nsym Input number of symmetry operators.
 * @param rsym Input symmetry operators.
 * @param origins Array containing alternative origins on output.
 * @param polarx Return whether polar along x axis.
 * @param polary Return whether polar along y axis.
 * @param polarz Return whether polar along z axis.
 * @param iprint If true, print out list of alternative origins.
 * @return Number of alternate origins for spacegroup.
 */
int ccp4spg_generate_origins(const char *namspg, const int nsym, const float rsym[][4][4],
			     float origins[][3], int *polarx, int *polary, int *polarz,
			     const int iprint);

/** Print details on reciprocal spacegroup.
 * @param sp pointer to spacegroup
 * @return void
 */
void ccp4spg_print_recip_spgrp(const CCP4SPG* sp);

/** Print reciprocal symops.
 * @param sp pointer to spacegroup
 * @return void
 */
void ccp4spg_print_recip_ops(const CCP4SPG* sp);

/** Convert string of type 0<=y<=1/4 to 0.0-delta, 0.25+delta, where
 * delta is set to 0.00001 Makes many assumptions about string.
 * @param range input string.
 * @param limits output range limits.
 * @return 0 on success
 */
int range_to_limits(const char *range, float limits[2]);

/** Sets an FFT grid for a spacegroup.
 * @param sp pointer to spacegroup
 * @param nxmin minimum sampling on x
 * @param nymin minimum sampling on y
 * @param nzmin minimum sampling on z
 * @param sample default fineness of sample
 * @param nx returns sampling intervals along x
 * @param ny returns sampling intervals along y
 * @param nz returns sampling intervals along z
 * @return void
 */
void set_fft_grid(CCP4SPG* sp, const int nxmin, const int nymin, const int nzmin, 
		  const float sample, int *nx, int *ny, int *nz);

/** Checks whether all factors of a number n are less than or
 * equal to 19.
 * @param n Number to be tested.
 * @return 1 on success, O on failure.
 */
int all_factors_le_19(const int n);

/** Sets a grid sample greater than minsmp, which has no prime
 * factors greater than 19, and contains the factor nmul.
 * @param minsmp
 * @param nmul
 * @param sample
 * @return Grid sample or -1 on failure.
 */
int get_grid_sample(const int minsmp, const int nmul, const float sample);

/** Check for consistency between cell dimensions and spacegroup. Latter
 * is identified from symmetry operators.
 * @param nsym No. of symmetry operators.
 * @param rsym Symmetry operators.
 * @param cell Cell dimensions.
 * @return 1 if they are consistent, 0 if there is a problem. 
 */
int ccp4spg_check_symm_cell(int nsym, float rsym[][4][4], float cell[6]);

#ifdef __cplusplus
} }
#endif
#endif
