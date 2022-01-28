/*
     cmtzlib.h: header file for cmtzlib.c
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

/** @page cmtz_page CMTZ library
 *
 * @verbatim

<!-- ::INDEX_INFO::CMTZ library::Library::::C/C++ Software Library for MTZ files:::::::: -->

   @endverbatim
 *
 *  @section cmtz_file_list File list

<ul>
<li>cmtzlib.h - contains details of the C/C++ API
<li>mtzdata.h - contains details of the MTZ data structure
</ul>

 *  @section cmtz_overview Overview
 
The CMTZ library is centred around a data structure, defined in
mtzdata.h Reading an MTZ file causes a data structure to be
created and populated from the MTZ header. There are a variety of functions 
for manipulating the data structure, for example adding crystals, datasets
or columns. The data structure can be dumped to an output MTZ data file.
<p>
The library operates in two modes with respect to the reflection data.
If <tt>mtz->refs_in_memory = 1</tt> (as set e.g. by the argument to 
<tt>MtzGet</tt>), then all reflection data is read into memory from
the file, and can be manipulated in memory. Else if 
<tt>mtz->refs_in_memory = 0</tt> then explicit calls to <tt>ccp4_lrrefl</tt>,
<tt>ccp4_lrreff</tt> and <tt>ccp4_lwrefl</tt> are required to read
and write reflections from/to disk.
<p>
Information on the data structure is given in mtzdata.h

 *  @section cmtz_princ_func Principal Functions
 *  @subsection cmtz_reading Reading MTZ files

Start by looking at <tt>MtzGet</tt> and <tt>ccp4_lrrefl</tt> / <tt>ccp4_lrreff</tt>.

 *  @subsection cmtz_writing Writing MTZ files

If you have a structure in memory already, use <tt>MtzPut</tt> followed
by <tt>MtzFree</tt> to release the memory. 
<p>
If you need to create a structure from scratch (i.e. without reading
from an input file) then use <tt>MtzMalloc</tt>, <tt>MtzOpenForWrite</tt>, 
<tt>ccp4_lwsymm</tt>, <tt>MtzAddXtal</tt>, <tt>MtzAddDataset</tt>, 
<tt>MtzAddColumn</tt> and <tt>ccp4_lwrefl</tt>.

 * @section cmtz_symmetry Symmetry Information

All reflection data in an MTZ file is assumed to belong to the same
spacegroup. The spacegroup is identified in the MTZ file by SYMINF
and SYMM records in the file header. This information is copied 
into the in-memory data structure. The list of symmetry operators
(copied from the SYMM header records) is taken to be the definitive
indicator of the spacegroup.
<p>
The functions <tt>ccp4_lrsymi</tt>, <tt>ccp4_lrsymm</tt> and 
<tt>ccp4_lwsymm</tt> read from and write to the symmetry sections
of the data structure. No symmetry manipulations are done within
the CMTZ library itself. Within CCP4, the CSYM library provides 
appropriate functions, but other symmetry libraries could be used.

 *  @section cmtz_examples Examples

See examples on <a href="ftp://ftp.ccp4.ac.uk/mdw/cmtz">ftp area</a>

 */

/** @file cmtzlib.h
 *
 *  @brief C-level library for input, output and manipulation of MTZ files.
 *
 *  Functions defining the C-level API for accessing MTZ files.
 *  MtzGet and MtzPut read and write MTZ files to/from a data
 *  structure defined in mtzdata.h  Other functions allow one
 *  to access data structure members, and to manipulate the structure
 *  and the values of structure members. Functions with names
 *  beginning <tt>ccp4_lr</tt> or <tt>ccp4_lw</tt> are primarily 
 *  to support the Fortran API.
 *
 *  @author Martyn Winn 
 */

#ifndef __CMTZLib__
#define __CMTZLib__

/* rcsidhm[100] = "$Id$" */

/* defines CCP4::CCP4File */
#include "ccp4_utils.h"

#ifdef  __cplusplus
namespace CMtz {
extern "C" {
typedef CCP4::CCP4File CCP4File;
#endif

/* needs to be here for C++ to use CCP4File typedef */
#include "mtzdata.h"

/**** MTZ i/o ****/

/** Reads the contents of the MTZ file into an MTZ structure.
 * @param logname (I) Logical name of MTZ file
 * @param read_refs (I) Whether to read reflections into memory (non-zero) or
 *        to read later from file (zero)
 * @return Pointer to MTZ struct
 */
MTZ *MtzGet(const char *logname, int read_refs);

/** Reads the contents of the MTZ file into an MTZ structure. As for function
 * MtzGet except for extra argument cell_tolerance.
 * @param logname (I) Logical name of MTZ file
 * @param read_refs (I) Whether to read reflections into memory (non-zero) or
 *        to read later from file (zero)
 * @param cell_tolerance (I) User-defined tolerance for ccp4uc_cells_differ.
 *        Setting this to zero means that a new crystal is generated whenever
 *        dataset cell dimensions are different. MtzGet allows for a certain
 *        variation within a single crystal.
 * @return Pointer to MTZ struct
 */
MTZ *MtzGetUserCellTolerance(const char *logname, int read_refs, const double cell_tolerance);

/** Reads reflection data from MTZ file.
 * @param filein pointer to input file
 * @param ncol number of columns to read
 * @param refldata array of reflection data
 * @return istat from ccp4_file_read
 */
int MtzRrefl(CCP4File *filein, int ncol, float *refldata);

/** Writes an MTZ data structure to disk. If file is already open, MtzPut
 * uses file pointer in mtz struct, else uses logical name of file.
 * @param mtz pointer to MTZ struct.
 * @param logname logical name for output file or blank.
 * @return 1 on success, 0 on failure
 */
int MtzPut(MTZ *mtz, const char *logname);

/** Opens a new MTZ file for writing. The output file can be specified
 * either with a true filename, or more likely as a logical name
 * corresponding to an environment variable or a CCP4 command line
 * argument such as HKLOUT.
 * @param logname logical name or filename for output file.
 * @return pointer to file or NULL on failure
 */
CCP4File *MtzOpenForWrite(const char *logname);

/** Write header record to fileout. Record is filled from
 * buffer and padded by blanks to a total length of MTZRECORDLENGTH.
 * @param fileout Pointer to output file.
 * @param nitems Number of characters in buffer.
 * @param buffer Character buffer containing MTZ header line.
 * @return Number of bytes written on success, EOF on failure.
 */
int MtzWhdrLine(CCP4File *fileout, int nitems, char buffer[]);

/** Write ncol column entries to fileout from refldata.
 * @param fileout pointer to output MTZ file.
 * @param ncol number of reflection data items to write.
 * @param refldata array of reflection data items.
 * @return Number of items written. If this is less than ncol, then
 * that indicates a write error.
 */
int MtzWrefl(CCP4File *fileout, int ncol, float *refldata);

/** Delete a reflection from the data structure. Beware, there
 * is no going back!
 * @param mtz pointer to MTZ struct.
 * @param iref index of reflection
 * @return 0 if successful, 1 otherwise
 */
int MtzDeleteRefl(MTZ *mtz, int iref);

/** Position input file at start of reflections. Useful if you need
 * to read the reflections a second time.
 * @param mtz pointer to MTZ struct.
 */
void MtzRewdInput(MTZ *mtz);

/**** Memory allocation ****/

/** Allocates memory for an MTZ header structure. The structure can
 * contain 0, 1 or more crystals, and for each crystal 0, 1 or more 
 * datasets. Crystals have a name based on the time "NULL_xnameHHMMSS"
 * to ensure uniqueness (compared to crystals defined elsewhere - all
 * new crystals created here will (probably) have the same name).
 * Crystals have the project name "NULL_pname", and datasets have the
 * name "NULL_dname".
 * @param nxtal Number of crystals to allocate.
 * @param nset Number of datasets for each crystal to allocate.
 * @return pointer to MTZ header struct
 */
MTZ *MtzMalloc(int nxtal, int nset[]);

/** Frees the memory reserved for the MTZ header struct.
 * @param mtz pointer to MTZ header struct. 
 * @return 1 on success, 0 on failure
 */
int MtzFree(MTZ *mtz);

/** Allocates memory for an MTZ column. Space is allocated for
 * the reflection data if and only if mtz->refs_in_memory is set.
 * @param mtz pointer to MTZ header struct. 
 * @param nref number of reflections in column.
 * @return pointer to MTZ column.
 */
MTZCOL *MtzMallocCol(MTZ *mtz, int nref);

/** Frees the memory reserved for 'col'
 * @param col pointer to MTZ column.
 * @return 1 on success, 0 on failure
 */
int MtzFreeCol(MTZCOL *col);

/** Allocates memory for a single batch header.
 * @return pointer to batch
 */
MTZBAT *MtzMallocBatch(void);

/** Frees the memory reserved for 'batch'.
 * @param batch
 * @return 1 on success, 0 on failure
 */
int MtzFreeBatch(MTZBAT *batch);
  
/** Allocates memory for the mtz history with 'nhist' lines.
 * @param nhist
 * @return pointer to history
 */
char *MtzCallocHist(int nhist);

/** Frees the memory reserved for 'hist'.
 * @param hist
 * @return 1 on success, 0 on failure
 */
int MtzFreeHist(char *hist);

/** Free all memory malloc'd from static pointers.
 * To be called before program exit. The function can be
 * registered with atexit.
 * @return void
 */
void MtzMemTidy(void);

/**** Header operations ****/
  
/** Get the number of batches in the mtz.
 * @param mtz pointer to MTZ struct
 * @return Number of batches.
 */
int MtzNbat(const MTZ *mtz);

/** Get the number of reflections in the mtz.
 * @param mtz pointer to MTZ struct
 * @return Number of reflections.
 */
int MtzNref(const MTZ *mtz);

/** Get the spacegroup number (likely CCP4 convention).
 * @param mtz pointer to MTZ struct
 * @return Spacegroup number.
 */
int MtzSpacegroupNumber(const MTZ *mtz);

/** Return the overall resolution limits of the MTZ structure.
 * These are the widest limits over all crystals present.
 * @param mtz pointer to MTZ struct
 * @param minres minimum resolution
 * @param maxres maximum resolution
 * @return 1 on success, 0 on failure 
 */
int MtzResLimits(const MTZ *mtz, float *minres, float *maxres);

/**** Crystal operations ****/

/** Get the total number of crystals in the MTZ structure
 * @param mtz pointer to MTZ struct
 * @return number of active crystals
 */
int MtzNxtal(const MTZ *mtz);

/** Get the number of active crystals in the MTZ structure
 * @param mtz pointer to MTZ struct
 * @return number of active crystals
 */
int MtzNumActiveXtal(const MTZ *mtz);

/** Return array of pointers to crystals.
 * @param mtz pointer to MTZ struct
 * @return array of pointers to crystals
 */
MTZXTAL **MtzXtals(MTZ *mtz);

/** Return pointer to the ixtal'th crystal. 
 * @param mtz pointer to MTZ struct
 * @param ixtal number of the particular crystal (ixtal = 0 ... MtzNxtal(xtal) -1 
 * @return pointer to the specified crystal
 */
MTZXTAL *MtzIxtal(const MTZ *mtz, const int ixtal);

/** Return the full path name of a crystal as "/xname" 
 * Memory for the path name is assigned with malloc, and can
 * be free'd by the calling function.
 * @param xtal pointer to the crystal struct
 * @return pointer to string containing path name
 */
char *MtzXtalPath(const MTZXTAL *xtal);

/** Returns a pointer to the crystal of mtz with the given `label`, or NULL.
 * @param mtz pointer to MTZ struct
 * @param label
 * @return pointer to crystal
 */
MTZXTAL *MtzXtalLookup(const MTZ *mtz, const char *label);

/** Add a crystal to header mtz.
 * @param mtz pointer to MTZ struct
 * @param xname Crystal name.
 * @param pname Name of associated project.
 * @param cell Cell dimensions of crystal.
 * @return Pointer to crystal.
 */
MTZXTAL *MtzAddXtal(MTZ *mtz, const char *xname, const char *pname,
                  const float cell[6]);

/** For a given crystal, return number of datasets in that crystal.
 * @param xtal pointer to the crystal struct
 * @return number of datasets
 */
int MtzNsetsInXtal(const MTZXTAL *xtal);

/** For a given crystal, return number of active datasets in that crystal.
 * @param xtal pointer to the crystal struct
 * @return number of datasets
 */
int MtzNumActiveSetsInXtal(const MTZ *mtz, const MTZXTAL *xtal);

/** For a given crystal, return array of pointers to datasets 
 * in that crystal.
 * @param xtal pointer to the crystal struct
 * @return array of pointers to datasets 
 */
MTZSET **MtzSetsInXtal(MTZXTAL *xtal);

/** For a given crystal, return pointer to the iset'th dataset
 * in that crystal.
 * @param xtal pointer to the crystal struct
 * @param iset number of the particular dataset (iset = 0 ... MtzNsetsInXtal(xtal) -1
 * @return pointer to specified dataset
 */
MTZSET *MtzIsetInXtal(const MTZXTAL *xtal, const int iset);

/**** Dataset operations ****/

/** Get the number of datasets in the MTZ structure
 * @param mtz pointer to MTZ struct
 * @return total number of datasets
 */
int MtzNset(const MTZ *mtz);

/** Get the number of active datasets in the MTZ structure
 * @param mtz pointer to MTZ struct
 * @return total number of datasets
 */
int MtzNumActiveSet(const MTZ *mtz);

/** Get the crystal associated with a dataset 
 * The pointer to MTZ is required to do reverse lookup of xname. 
 * @param mtz pointer to MTZ struct
 * @param set pointer to dataset
 * @return pointer to parent crystal, or NULL if "set" is
 *    not present in "mtz".
 */
MTZXTAL *MtzSetXtal(const MTZ *mtz, const MTZSET *set);

/** Return the full path name of a dataset as "/xname/dname" 
 * The pointer to MTZ is required to do reverse lookup of xname.
 * Memory for the path name is assigned with malloc, and can
 * be free'd by the calling function.
 * @param mtz pointer to MTZ struct
 * @param set pointer to dataset
 * @return pointer to string containing path name
 */
char *MtzSetPath(const MTZ *mtz, const MTZSET *set);

/** Returns a pointer to the dataset of MTZ with the given label.
 * @param mtz pointer to MTZ struct.
 * @param label Label of desired set. This could be <dname> or
 *   <xname>/<dname>.
 * @return pointer to set or NULL if not found
 */
MTZSET *MtzSetLookup(const MTZ *mtz, const char *label);

/** Add a dataset to crystal xtl
 * @param mtz pointer to MTZ struct.
 * @param xtl pointer to crystal struct.
 * @param dname Dataset name
 * @param wavelength X-ray wavelength of dataset
 * @return pointer to set
 */
MTZSET *MtzAddDataset(MTZ *mtz, MTZXTAL *xtl, const char *dname,
                    const float wavelength);

/** For a given dataset, return number of columns in that dataset.
 * This is simply set->ncol and so includes all columns irrespective
 * of col->active
 * @param set pointer to dataset
 * @return number of columns
 */
int MtzNcolsInSet(const MTZSET *set);

/** For a given dataset, return number of active columns in that dataset.
 * @param set pointer to dataset
 * @return number of active columns
 */
int MtzNumActiveColsInSet(const MTZSET *set);

/** For a given dataset, return number of columns in that dataset
 * which have a source in an input file (i.e. non-zero source attribute).
 * @param set pointer to dataset
 * @return number of source columns
 */
int MtzNumSourceColsInSet(const MTZSET *set);

/** For a given dataset, return number of batches in that dataset.
 * @param mtz pointer to MTZ struct
 * @param set pointer to dataset
 * @return number of batches
 */
int MtzNbatchesInSet(const MTZ *mtz, const MTZSET *set);

/** For a given dataset, return array of pointers to columns 
 * in that dataset.
 * @param set pointer to dataset
 * @return array of pointers to columns 
 */
MTZCOL **MtzColsInSet(MTZSET *set);

/** For a given dataset, return pointer to the icol'th column
 * in that dataset.
 * @param set pointer to dataset
 * @param icol number of 
 * the particular column (icol = 0 ... MtzNcolsInSet(set) -1
 * @return pointer to specified column
 */
MTZCOL *MtzIcolInSet(const MTZSET *set, const int icol);

/**** Column operations ****/

/** Add a column to dataset set and create + fill with NAN
 * @param mtz pointer to MTZ struct
 * @param set pointer to dataset
 * @param label Column label
 * @param type Column type
 * @return pointer to column
 */
MTZCOL *MtzAddColumn(MTZ *mtz, MTZSET *set, const char *label,
                   const char *type);

/** Assigns HKL columns to the base dataset.
 * @param mtz pointer to MTZ struct
 * @return 1 on success, 0 on failure
 */
int MtzAssignHKLtoBase(MTZ *mtz);

/** Assigns a column to a dataset identified by crystal_name and
 * dataset_name. First, the function checks whether the
 * column already belongs to this dataset, in which case it does nothing.
 * Then it checks if the requested dataset exists. If not, it is created,
 * though it is better to explicitly create it beforehand. Finally, the column
 * is assigned to the dataset.
 * @param mtz pointer to MTZ struct
 * @param col pointer to column
 * @param crystal_name name of crystal containing dataset
 * @param dataset_name name of dataset
 * @return 1 on success, 0 on failure
 */
int MtzAssignColumn(MTZ *mtz, MTZCOL *col, const char crystal_name[],  
	     const char dataset_name[]);

/** Toggle active flag of column. A value of 0 means inactive and will
 * not be written out, whereas a value of 1 means active and will
 * be written out.
 * @param col pointer to column
 * @return New value of col->active.
 */
int MtzToggleColumn(MTZCOL *col);

/** Get the dataset associated with a column.
 * @param mtz pointer to MTZ struct
 * @param col pointer to column of interest
 * @return pointer to set containing column of interest, or NULL
 *   if "col" is not contained in "mtz".
 */
MTZSET  *MtzColSet(const MTZ *mtz, const MTZCOL *col);

/** Get the number of columns in the MTZ data structure.
 * @param mtz pointer to MTZ struct
 * @return number of columns
 */
int MtzNcol(const MTZ *mtz);
  
/** Get the number of active columns in the MTZ data structure.
 * @param mtz pointer to MTZ struct
 * @return number of columns
 */
int MtzNumActiveCol(const MTZ *mtz);
  
/** Get the number of columns in the MTZ data structure which have 
 * a source in an input file (i.e. non-zero source attribute).
 * @param mtz pointer to MTZ struct
 * @return number of columns
 */
int MtzNumSourceCol(const MTZ *mtz);

/** Return the full path name of a column as "/xname/dname/label" 
 * Memory for the path name is assigned with malloc, and can
 * be free'd by the calling function.
 * @param mtz pointer to MTZ struct
 * @param col pointer to MTZ column.
 * @return full path name of column
 */
char *MtzColPath(const MTZ *mtz, const MTZCOL *col);

/** Complete a right-justified path by prefixing with wildcards
 * @param path Completed path.
 * @param partial Partial right-justified path
 * @param njust
 * @return 1 on success, 0 on failure.
 */
int MtzRJustPath(char *path, const char *partial, const int njust);

/** Test for match between two paths, including wildcards
 * @param path1 First path
 * @param path2 Second path
 * @return 1 if paths match, else 0.
 */
int MtzPathMatch(const char *path1, const char *path2);

/** Returns a pointer to the column of mtz with the given `label`, or NULL
 * @param mtz pointer to MTZ struct
 * @param label Column label.
 * @return pointer to column
 */
MTZCOL *MtzColLookup(const MTZ *mtz, const char *label);

/** Get the MTZ column type of a particular column.
 * @param col pointer to MTZ column.
 * @return column type
 */
char *MtzColType(MTZCOL *col);

/** Print summary of current crystal/dataset/column hierarchy. This
 * is designed for debugging purposes rather than for the user.
 * @param mtz pointer to MTZ struct
 * @return void
 */
void MtzDebugHierarchy(const MTZ *mtz);

/** List of column information: label, type, dataset.
 * @param mtz pointer to MTZ struct
 * @param clabs List of labels (output).
 * @param ctyps List of column types (output).
 * @param csetid List of dataset IDs (output).
 * @return number of columns in current structure.
 */
int MtzListColumn(const MTZ *mtz, char clabs[][31], char ctyps[][3], int csetid[]);

/** List of column information from input file: label, type, dataset.
 * @param mtz pointer to MTZ struct
 * @param clabs List of labels (output).
 * @param ctyps List of column types (output).
 * @param csetid List of dataset IDs (output).
 * @return number of columns in input file.
 */
int MtzListInputColumn(const MTZ *mtz, char clabs[][31], char ctyps[][3], int csetid[]);

/**** helper functions ****/

/** Find where indices h, k, l are in MTZ structure. Usually, they
 * will be first 3 columns of 1st dataset, but safest not to assume this.
 * @param mtz pointer to MTZ struct
 * @param ind_xtal crystal containing indices
 * @param ind_set dataset containing indices
 * @param ind_col 3 columns containing indices
 * @return 1 on success, 0 on failure
 */
int MtzFindInd(const MTZ *mtz, int *ind_xtal, int *ind_set, int ind_col[3]);

/** Calculate resolution from indices and coefhkl.
 * coefhkl is obtained from MtzHklcoeffs.
 * @param in integer array of 3 indices
 * @param coefhkl double array of 6 coefficients
 * @return resolution
 */
float MtzInd2reso(const int in[3], const double coefhkl[6]);

/** Generate coefhkl coefficients from given cell parameters.
 * @param cell cell dimensions to be used for resolution calculation.
 * @param coefhkl double array of 6 coefficients
 * @return 1 on success, 0 on failure
 */
int MtzHklcoeffs(const float cell[6], double coefhkl[6]);

/** Reads batch arrays into data structure.
 * @param intbuf pointer to integer batch array
 * @param fltbuf pointer to float batch array
 * @param batch pointer to batch structure
 * @return 1 on success, 0 on failure
 */
int MtzArrayToBatch(const int *intbuf, const float *fltbuf, MTZBAT *batch);

/** Writes data structure to batch arrays.
 * @param batch pointer to batch structure
 * @param intbuf pointer to integer batch array
 * @param fltbuf pointer to float batch array
 * @return 1 on success, 0 on failure
 */
int MtzBatchToArray(MTZBAT *batch, int *intbuf, float *fltbuf);

/** Sort a linked list of batches on batch number. The function
 * first checks whether batches are already in order, as they
 * will usually be.
 * @param batch pointer to beginning of batch list
 * @param numbat number of batches to be sorted
 * @return Sorted list if batches sorted, original list if batches 
 *  already in order
 */
MTZBAT *sort_batches(MTZBAT *batch, int numbat);

/**** pseudo-mtzlib routines ****/

/** Returns title of MTZ structure 'mtz'
 * @param mtz pointer to MTZ struct
 * @param title MTZ title as character string
 * @return length of title excluding trailing blanks and terminating null
 *  character.
 */
int ccp4_lrtitl(const MTZ *mtz, char *title);

/** Get history lines from MTZ structure.
 * @param mtz Pointer to MTZ struct.
 * @param history Returned history lines.
 * @param nlines Requested number of history lines.
 * @return Actual number of history lines returned.
 */
int ccp4_lrhist(const MTZ *mtz, char history[][MTZRECORDLENGTH], int nlines);

/** Get sort order from MTZ structure.
 * @param mtz Pointer to MTZ struct.
 * @param isort Returned sort order.
 * @return 1 on success, 0 on failure
 */
int ccp4_lrsort(const MTZ *mtz, int isort[5]);

/** Get batch numbers from MTZ structure.
 * @param mtz Pointer to MTZ struct.
 * @param nbatx Number of batches in input file.
 * @param batchx Returned array of batch numbers.
 * @return Number of batches.
 */
int ccp4_lrbats(const MTZ *mtz, int *nbatx, int batchx[]);

/** Get cell dimensions for a particular crystal.
 * @param xtl Pointer to crystal.
 * @param cell Output cell dimensions.
 * @return 1 on success, 0 on failure.
 */
int ccp4_lrcell(const MTZXTAL *xtl, float cell[]);

/** Get spacegroup information as held in MTZ header.
 * @param mtz Pointer to MTZ struct.
 * @param nsympx Number of primitive symmetry operators.
 * @param ltypex Lattice type (P,A,B,C,I,F,R).
 * @param nspgrx Spacegroup number.
 * @param spgrnx Spacegroup name.
 * @param pgnamx Pointgroup name.
 * @return Spacegroup number.
 */
int ccp4_lrsymi(const MTZ *mtz, int *nsympx, char *ltypex, int *nspgrx, 
       char *spgrnx, char *pgnamx);

/** Get spacegroup information as held in MTZ header. Extended version of
 * ccp4_lrsymi to return confidence flag as well.
 * @param mtz Pointer to MTZ struct.
 * @param nsympx Number of primitive symmetry operators.
 * @param ltypex Lattice type (P,A,B,C,I,F,R).
 * @param nspgrx Spacegroup number.
 * @param spgrnx Spacegroup name.
 * @param pgnamx Pointgroup name.
 * @param spgconf One-character flag indicating confidence in nominal spacegroup.
 * @return Spacegroup number.
 */
int ccp4_lrsymi_c(const MTZ *mtz, int *nsympx, char *ltypex, int *nspgrx, 
       char *spgrnx, char *pgnamx, char *spgconf);

/** Get symmetry matrices from MTZ structure. Note: ordering of matrices
 * in rsymx was changed in April 2004.
 * @param mtz Pointer to MTZ struct.
 * @param nsymx Number of symmetry operators held in MTZ header.
 * @param rsymx Symmetry operators as 4 x 4 matrices, in the order they
 *   are held in the MTZ header. Each matrix has translations in
 *   elements [*][3].
 * @return Number of symmetry operators.
 */
int ccp4_lrsymm(const MTZ *mtz, int *nsymx, float rsymx[192][4][4]);

/** Uses LABIN or LABOUT line to convert program labels to user labels.
 * This is a helper function, but does not access reflection structure at all.
 * @param labin_line (I) LABIN/LABOUT line from Parser.
 * @param prog_labels (I) Progam labels.
 * @param nlprgi (I) Number of program labels.
 * @param user_labels (O) On output, user-supplied file labels in corresponding
 *  positions. For unassigned program labels, user_labels is empty string.
 * @return Number of program labels matched, or -1 if there was an error.
 */
int MtzParseLabin(char *labin_line, const char prog_labels[][31], 
	   const int nlprgi, char user_labels[][2][31]);

/** Finds columns in an MTZ struct according to column labels. Column types
 *  are checked for agreement between requested type (in argument 'types')
 *  and file type. If requested type is blank, file type is returned in 
 *  argument 'types'. Note, this function is different from Fortranic LRASSN, 
 *  in that any conversion from program labels to user labels should have been done
 *  previously.
 * @param mtz Pointer to MTZ struct.
 * @param labels Input array of column labels to be found in MTZ struct.
 * @param nlabels Number of columns to be found.
 * @param types Input array of column types of columns to be found.
 * @return Array of pointers to columns in MTZ struct. Array element is
 *  NULL if column not found.
 */
MTZCOL **ccp4_lrassn(const MTZ *mtz, const char labels[][31], const int nlabels, 
		char types[][3]);

/** Report information on a particular dataset. This represents the
 * collection of data held in one series of dataset records in the MTZ header.
 * It is mainly useful for supporting old Fortran calls.
 * @param mtz pointer to MTZ struct
 * @param set pointer to dataset
 * @param crystal_name Crystal name
 * @param dataset_name Dataset name
 * @param project_name Project name
 * @param isets Dataset ID.
 * @param datcell Cell dimensions of crystal that dataset belongs to.
 * @param datwave X-ray wavelength associated with dataset.
 * @return 1 on success, 0 on failure
 */
int ccp4_lridx(const MTZ *mtz, const MTZSET *set, char crystal_name[64], 
            char dataset_name[64], char project_name[64], int *isets, 
            float datcell[6], float *datwave);

/** Returns iref'th reflection from file held in MTZ struct mtz. Returns
 * data for all columns held in input file, in the order that they are 
 * held in the source file. The value of col->source can be used to find
 * particular columns.
 * In "in-memory" mode, reflection data is taken from arrays held in
 * memory. In the traditional file-based mode, a reflection record is
 * read from the input file.
 * @param mtz pointer to MTZ struct
 * @param resol resolution of reflection (output).
 * @param adata array of requested values (output).
 * @param logmss array of flags for missing data (output). A value of 1 indicates
 * a Missing Number Flag, and a value of 0 indicates usable data.
 * @param iref index of requested reflection (starting at 1).
 * @return 1 if past last reflection, else 0
 */
int ccp4_lrrefl(const MTZ *mtz, float *resol, float adata[], int logmss[], int iref);

/** Returns iref'th reflection from file held in MTZ struct mtz. Returns
 * data for certain columns held in input file, as specified by the
 * column pointers held in lookup.
 * In "in-memory" mode, reflection data is taken from arrays held in
 * memory. In the traditional file-based mode, a reflection record is
 * read from the input file.
 * @param mtz pointer to MTZ struct
 * @param resol resolution of reflection (output).
 * @param adata array of requested values (output).
 * @param logmss array of flags for missing data (output). A value of 1 indicates
 * a Missing Number Flag, and a value of 0 indicates usable data.
 * @param lookup array of pointers to requested columns
 * @param ncols number of requested columns
 * @param iref index of requested reflection (starting at 1).
 * @return 1 if past last reflection, else 0
 */
int ccp4_lrreff(const MTZ *mtz, float *resol, float adata[], int logmss[],
		const MTZCOL *lookup[], const int ncols, const int iref);

/** Checks whether a particular reflection value represents missing data.
 * @param mtz Pointer to the MTZ struct, which holds the value representing
 * missing data (the Missing Number Flag) against which the input datum
 * is compared.
 * @param datum Reflection value to be checked.
 * @return Returns 1 is datum is the MNF, and 0 otherwise.
 */
int ccp4_ismnf(const MTZ *mtz, const float datum);

/** Function to print header information in traditional format.
 * @param mtz Pointer to MTZ struct
 * @param iprint Print level
 * @return 1 on success, 0 on failure 
 */
int ccp4_lhprt(const MTZ *mtz, int iprint);

/** Function to print header information in format appropriate
 * to data structure hierarchy.
 * @param mtz Pointer to MTZ struct
 * @param iprint Print level
 * @return 1 on success, 0 on failure 
 */
int ccp4_lhprt_adv(const MTZ *mtz, int iprint);

/** Function to return batch header data for a specified batch.
 * @param batch Pointer to requested batch.
 * @param buf On return, real and integer batch data.
 * @param charbuf On return, character batch data (title and axes names).
 * @param iprint =0 no printing, =1 print title only, >1 print full header.
 * @return 1 on success, 0 on failure
 */
int ccp4_lrbat(MTZBAT *batch, float *buf, char *charbuf, int iprint);

/** Function to print batch header data for a specified batch to stdout.
 * @param batch Pointer to requested batch.
 * @return 1 on success, 0 on failure 
 */
int MtzPrintBatchHeader(const MTZBAT *batch);

/** Write header title for later output to file.
 * @param mtz Pointer to MTZ struct.
 * @param ftitle Title string.
 * @param flag If 0 overwrite existing title, else append to
 * existing title.
 * @return 1 on success, 0 on failure 
 */
int ccp4_lwtitl(MTZ *mtz, const char *ftitle, int flag);

/** Sets the sort order in the MTZ header. The sort order is
 * a list of columns to be used for sorting, to be applied in
 * the order they appear in the list, i.e. sort first on
 * colsort[0], then on colsort[1], etc. If there are less than
 * 5 columns to be used for sorting, some of colsort[] may be NULL.
 * @param mtz Pointer to MTZ struct
 * @param colsort Array of pointers to columns.
 * @return 1 on success, 0 on failure
 */
int MtzSetSortOrder(MTZ *mtz, MTZCOL *colsort[5]);

/** Adds history lines to the MTZ header in front of existing history lines.
 * @param mtz pointer to MTZ struct
 * @param history lines to be added
 * @param nlines number of lines to be added
 * @return total number of history lines
 */
int MtzAddHistory(MTZ *mtz, const char history[][MTZRECORDLENGTH], const int nlines);

/** Write or update symmetry information for MTZ header. This provides support
 * for the Fortran API, and is not particularly convenient for C programs. 
 * Note: ordering of matrices in rsymx was changed in November 2003.
 * @param mtz Pointer to MTZ struct
 * @param nsymx Number of symmetry operators. If zero, symmetry operations
 *   are not updated.
 * @param nsympx Number of primitive symmetry operators.
 * @param rsymx Array of symmetry operators (dimensions ordered in C convention,
 *   with translations in elements [*][3])
 * @param ltypex Lattice type. If blank, not updated.
 * @param nspgrx Spacegroup number. If zero, not updated.
 * @param spgrnx Spacegroup name. If blank, not updated.
 * @param pgnamx Point group name. If blank, not updated.
 * @return 1 on success, 0 on failure 
 */
int ccp4_lwsymm(MTZ *mtz, int nsymx, int nsympx, float rsymx[192][4][4], 
		char ltypex[], int nspgrx, char spgrnx[], char pgnamx[]);

/** Write or update symmetry information for MTZ header. This provides support
 * for the Fortran API, and is not particularly convenient for C programs. Extended
 * version of ccp4_lwsymm to set confidence flag as well.
 * Note: ordering of matrices in rsymx was changed in November 2003.
 * @param mtz Pointer to MTZ struct
 * @param nsymx Number of symmetry operators. If zero, symmetry operations
 *   are not updated.
 * @param nsympx Number of primitive symmetry operators.
 * @param rsymx Array of symmetry operators (dimensions ordered in C convention,
 *   with translations in elements [*][3])
 * @param ltypex Lattice type. If blank, not updated.
 * @param nspgrx Spacegroup number. If zero, not updated.
 * @param spgrnx Spacegroup name. If blank, not updated.
 * @param pgnamx Point group name. If blank, not updated.
 * @param spgconf One-character flag indicating confidence in nominal spacegroup.
 * @return 1 on success, 0 on failure 
 */
int ccp4_lwsymm_c(MTZ *mtz, int nsymx, int nsympx, float rsymx[192][4][4], 
		  char ltypex[], int nspgrx, char spgrnx[], char pgnamx[], 
                  char spgconf[]);

/** Write or update symmetry confidence information for MTZ header.
 * @param mtz Pointer to MTZ struct
 * @param spgconf One-character flag indicating confidence in nominal spacegroup.
 * @return 1 on success, 0 on failure 
 */
int ccp4_lwsymconf(MTZ *mtz, char spgconf[]);

/* Assign columns for writing. Check to see if columns already exist,
 * else create them. New columns are assigned to the base dataset if it 
 * exists, else the first dataset.
 * @param mtz pointer to MTZ struct
 * @param labels Input array of column labels to be assigned.
 * @param nlabels Number of columns.
 * @param types Input array of column types of columns.
 * @param iappnd If iappnd = 0, then deactivate columns which are
 *   not selected (allows one to write out a subset of columns). Else
 *   if iappnd = 1, append these columns to existing ones.
 * @return Array of pointers to columns in MTZ struct.
 */
MTZCOL **ccp4_lwassn(MTZ *mtz, const char labels[][31], const int nlabels, 
             const char types[][3], const int iappnd);

/* Add or update a dataset in the MTZ structure. If the crystal name is
 * not recognised, then a new crystal is created containing a single
 * dataset. If the crystal name is recognised, then the parameters of
 * the crystal are updated. The child dataset is then created if necessary
 * or an existing dataset updated.
 * Note that this function is used to update crystal parameters, as well
 * as add crystals. If a duplicate crystal name is used by mistake, then 
 * the old crystal parameters are lost.
 * @param mtz pointer to MTZ struct
 * @param crystal_name Crystal name
 * @param dataset_name Dataset name
 * @param project_name Project name
 * @param datcell Cell dimensions of crystal that dataset belongs to.
 * @param datwave X-ray wavelength associated with dataset.
 * @return 1 on success, 0 on failure
 */
int ccp4_lwidx(MTZ *mtz, const char crystal_name[],  const char dataset_name[],
	const char project_name[], const float datcell[6], const float *datwave);


/** Function to output reflection values for iref'th reflection.
 * The reflection values are provided in an array "adata". The value
 * in adata[i] will be assigned to the MTZ column pointed to by lookup[i].
 * Typically, lookup[i] is a pointer returned from an earlier call to
 * MtzAddColumn().
 * In "in-memory" mode, values are added/updated to column data held
 * in memory. In the traditional file-based mode, a reflection record is
 * written to file.
 * This function will also update the column ranges and the crystal/file
 * resolution limits.
 * @param mtz pointer to MTZ struct
 * @param adata array of reflection column values.
 * @param lookup array of pointers to columns.
 * @param ncol number of columns.
 * @param iref Reflection number such that 1st reflection is iref=1.
 * @return 1 on success, 0 on failure
 */
int ccp4_lwrefl(MTZ *mtz, const float adata[], MTZCOL *lookup[], 
		 const int ncol, const int iref);

/** Write new batch information to 'batch' or if 'batch' is NULL create 
 * new batch header with batch number 'batno'. If you try to create more 
 * than one new batch header with the same batch number, the function will
 * complain and return. It is OK to create a new batch with the same
 * number as one read from file - this is the mechanism for changing
 * batch headers.
 * @param mtz pointer to MTZ struct
 * @param batch pointer to batch
 * @param batno batch number
 * @param buf pointer to batch array
 * @param charbuf pointer to character batch array
 * @return 1 on success, 0 on failure
 */
int ccp4_lwbat(MTZ *mtz, MTZBAT *batch, const int batno, const float *buf, const char *charbuf);

int ccp4_lwbsetid(MTZ *mtz, MTZBAT *batch, const char xname[], const char dname[]);

/* -- Below here there are no implementations -- */

/* COMPLEX HLToSF(float hla, float hlb, float hlc, float hld, BOOLEAN centric); */
/* Returns the mean structure factor as a complex number from a structure
   factor probability distribution described by Hendrickson/Lattman
   coefficients. If `centric == TRUE`, the coefficients describe a centric
   distribution. */

/* MTZ *MtzSort(MTZ *mtz, char *ident); */
/* Sorts `mtz` using the identifiers (separated by spaces) in `ident` as
   keys. Sorting can be done on up to 200 columns. A pointer to `*mtz` is
   returned. */

/* MTZ *HLCombine (MTZ *to, float toscale, MTZ *frm, float frmscale); */
/* Combines extra phase information for common reflections between 'frm'
   and 'to' into the phase information of 'to'. The phase probabilities
   are described by Hendrickson Lattman coefficients, with the identifiers
   "HLA", "HLB", HLC", and "HLD", the indices are identified by "H", "K" 
   and "L". HL-coeffs from 'to' are scaled by 'toscale', the coeffs from
   'frm' are scaled by 'frmscale'. A pointer to `to` is returned. */

/* void MtzPhiFom(MTZ *mtz); */
/* Calculates the best phase and the figure of from Hendrickson Lattman
   coefficients. The following columns should be present in `mtz`:
   "H", "K" & "L" (indices); "PHIB" & "FOM" (best phase (degrees) and figure of
   merit); "HLA", "HLB", "HLC" & "HLD" (Hendrickson Lattman coefficients). */

/* MTZ *MtzCopy(MTZ *frm); */
/* Produces an exact copy of `frm` and returns a pointer to it. */

/* MTZ *MtzColAppend(MTZ *mtz, char *ident, char type); */
/* Appends a column to `*mtz` with identity `ident` and type `type`, provided
   no column with identity `ident` exists. */

/* MTZ *MtzColRemove(MTZ *mtz, char *ident); */
/* Removes a column from `*mtz` with identity `ident`. */

/* MTZ *MtzUpdateRanges(MTZ *mtz); */
/* Updates ranges of all columns in `mtz` and returns `mtz`. */

/* MTZCOL *MtzColNewRange(MTZCOL *col, int nref); */
/* Updates the minimum & maximum values in `col` and returns `col`. */

/* int *MtzUnique(MTZ *mtz, char *ident); */
/* Returns an array (k) of indices: k[j] returns the first occurrence of a set
   of idents, eg. MtzUnique(mtz, "H K L") returns an array from which all the 
   unique reflections can be determined by indexing with k: k[i] is the index
   of the last relection of a set which all have the same hkl indices. It is
   assumed that `mtz` is sorted, using `ident` as keys. */

/* float PhaseProb(float phase, float hla, float hlb, float hlc, float hld,
		BOOLEAN centric); */
/* Returns the probability of `phase` (expressed in radians) as determined by
   the Hendrickson-Lattman coefficients `hla`, `hlb`, `hlc` and `hld`. If
   `centric == TRUE`, the coefficients describe a centric distribution. */

#ifdef __cplusplus
} }
#endif
#endif
