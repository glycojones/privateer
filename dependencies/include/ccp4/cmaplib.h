/*
     cmaplib.h: C/C++ level API for accessing CCP4 map files
     Copyright (C) 2001  CCLRC, Charles Ballard

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

/** @page cmap_page CMAP library
 *
 * @verbatim

<!-- ::INDEX_INFO::CMAP library::Library::::C/C++ Software Library for CCP4 map files:::::::: -->

   @endverbatim
 *
 *  @section cmap_file_list File list

<ul>
<li>cmaplib.h - contains details of the C/C++ API
<li>cmap_data.h     
<li>cmap_header.h   
<li>cmap_skew.h     
<li>cmap_errno.h    
<li>cmap_labels.h   
<li>cmap_stats.h    
</ul>

 *  @section cmap_overview Overview
 
  Functions defining the C-level API for accessing CCP4 map files.

 */

/** @file cmaplib.h
 *
 *  @brief ccp4 map i/o user-level library header file
 *
 *  Functions defining the C-level API for accessing CCP4 map files.
 *
 *  @author Charles Ballard
 */

#ifndef __GUARD_MAPLIB
#define __GUARD_MAPLIB

#include "ccp4_utils.h"

#ifdef __cplusplus
namespace CMap_io {
typedef CCP4::CCP4File CCP4File;
extern "C" {
#endif

typedef struct _CMMFile_Skew CMMFile_Skew;
typedef struct _CMMFile_Labels CMMFile_Labels;
typedef struct _CMMFile_Symop CMMFile_Symop;
typedef struct _CMMFile_Data CMMFile_Data;
typedef struct _CMMFile_Stats CMMFile_Stats;
typedef struct _CMMFile CMMFile;

struct _CMMFile_Labels {
  unsigned int number;
  char *labels[10];
};

struct _CMMFile_Skew {
  float rotation[3][3];
  float translation[3];
};

struct _CMMFile_Symop {
unsigned int offset;
unsigned int size;
unsigned int number;
};

struct _CMMFile_Data {
  size_t offset;
  size_t section_size;
  size_t header_size;
  size_t block_size;
  unsigned int number;
};

struct _CMMFile_Stats {
  float offset;                /* pseudo zero value */
  float min;                   /* minimum density value */
  float max;                   /* maximum density value */
  double mean;               /* sum of densities (less offset) */
  double rms;              /* sum of square of densities (less offset) */
  int total;                    /* number of summed densities */
};

struct _CMMFile {
CCP4File *stream;
char *file_name;
unsigned int data_mode;
unsigned int close_mode;
float cell[6];
int spacegroup;
int EM_spacegroup;
char EM_exthead_type[5];
char EM_contents[5];
int map_dim[3];
int origin[3];
int cell_grid[3];
int axes_order[3];
CMMFile_Symop symop;
CMMFile_Data data;
CMMFile_Stats stats;
CMMFile_Labels labels;
CMMFile_Skew skew;
int reserved[8];
char user_access[28];
};

/* open a file for read/write */
void *ccp4_cmap_open(const char *filename, int mode);

/* close a file for read/write (dumping the header if write) */
void ccp4_cmap_close(CMMFile *mfile);

/* set the close mode (calculation of map statistics) */
void ccp4_cmap_closemode(CMMFile *mfile, unsigned int closemode);

/* seek to a section in the map (read mode only)*/
int ccp4_cmap_seek_section(CMMFile *mfile, int offset, unsigned int seek_mode);

/* seek to a row in a section (read mode only)*/
int ccp4_cmap_seek_row(CMMFile *, int offset, unsigned int seek_mode);

/* raw seek (read mode only)*/
int ccp4_cmap_seek_data(CMMFile *, int offset, unsigned int seek_mode);

/* read a map section from file to memory */
int ccp4_cmap_read_section(CMMFile *mfile, void *section);

/* read a row from file to memory */
int ccp4_cmap_read_row(CMMFile *mfile, void *row);

/* read n_items from file to memory (item determined by data mode) */
int ccp4_cmap_read_data(const CMMFile *mfile, void *items, int n_items);

/* write a map section from memory to file */
int ccp4_cmap_write_section(CMMFile *mfile, const void *section);

/* write a map row from memory to file */
int ccp4_cmap_write_row(CMMFile *mfile, const void *row);

/* write n_items from memory to file (item determined by data mode) */
int ccp4_cmap_write_data(CMMFile *mfile, const void *items, int n_items);

/* read the section header corresponding to the current section */
int ccp4_cmap_read_section_header(const CMMFile *mfile, char *header);

/* write the section header corresponding to the current section */
int ccp4_cmap_write_section_header(CMMFile *mfile, const char *header);

/* get the header parameters */
void ccp4_cmap_get_cell(const CMMFile *mfile, float *cell);
void ccp4_cmap_get_grid(const CMMFile *mfile, int *grid);
void ccp4_cmap_get_origin(const CMMFile *mfile, int *origin);
void ccp4_cmap_get_order(const CMMFile *mfile, int *axes_order);
void ccp4_cmap_get_dim(const CMMFile *mfile, int *map_dim);
int ccp4_cmap_get_spacegroup(const CMMFile *mfile);
void ccp4_cmap_get_mapstats(const CMMFile *mfile, float *min, float* max, 
                           double *mean, double *rms);

/* set the header parameters */
void ccp4_cmap_set_cell(CMMFile *mfile, const float *cell);
void ccp4_cmap_set_grid(CMMFile *mfile, const int *grid);
void ccp4_cmap_set_origin(CMMFile *mfile, const int *origin);
void ccp4_cmap_set_order(CMMFile *mfile, const int *axes_order);
void ccp4_cmap_set_dim(CMMFile *mfile, const int *map_dim);
void ccp4_cmap_set_spacegroup(CMMFile *mfile, int spacegroup);
void ccp4_cmap_set_mapstats(CMMFile *mfile, const float min, const float max,
                           const double mean, const double rms);

/* get map file datamode */
unsigned int ccp4_cmap_get_datamode(const CMMFile *mfile);

/* set map file datamode */
void ccp4_cmap_set_datamode(CMMFile *mfile, unsigned int datamode);

/* get the local header size */
size_t ccp4_cmap_get_local_header(CMMFile *mfile);

/* set the local header size (before data writing begins) */
void ccp4_cmap_set_local_header(CMMFile *mfile, size_t size);

/* get the number of symops in the file */
int ccp4_cmap_num_symop(const CMMFile *mfile);

/* seek among the symops strings */
int ccp4_cmap_seek_symop(CMMFile *mfile, int isymop, unsigned int whence);

/* read a symop string of 80 characters */
int ccp4_cmap_get_symop(CMMFile *mfile, char *buffer);

/* write a symop string of 80 characters */
int ccp4_cmap_set_symop(CMMFile *mfile, const char *buffer);

/* get the mask */
int ccp4_cmap_get_mask(const CMMFile *mfile, float *skew_mat, float *skew_trans);

/* set the mask */
int ccp4_cmap_set_mask(CMMFile *mfile, const float *skew_mat, const float *skew_trans);

/* the number of labels used */
int ccp4_cmap_number_label(const CMMFile *mfile);

/* set label at posn from C-string */
int ccp4_cmap_set_label(CMMFile *mfile, const char *label, int posn);

/* return label at posn as C-string */
char *ccp4_cmap_get_label(const CMMFile *mfile, int posn);

/* set title (label=0) */
int ccp4_cmap_set_title(CMMFile *mfile, const char *label);

/* get title (label=0) */
char *ccp4_cmap_get_title(const CMMFile *mfile);

#ifdef __cplusplus
}
}
#endif

#endif  /* __GUARD_MAPLIB */
