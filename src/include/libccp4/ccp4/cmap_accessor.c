/*
     cmap_accessor.c: get and set map header information
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
#include <math.h>
#include "cmaplib.h"
#include "cmap_errno.h"

/* accessors */

/*! Get the cell parameters
 \param mfile (const CMMFile *)
 \param cell (float *) contains the cell parameter on exit (dim 6) */
void ccp4_cmap_get_cell(const CMMFile *mfile, float *cell)
{
  cell[0] = mfile->cell[0];
  cell[1] = mfile->cell[1];
  cell[2] = mfile->cell[2];
  cell[3] = mfile->cell[3];
  cell[4] = mfile->cell[4];
  cell[5] = mfile->cell[5];
}

/*! Set the cell parameters.
  Only allowed when file is opened in write mode.
 \param mfile (CMMFile *)
 \param cell (const float *) the cell parameters */
void ccp4_cmap_set_cell(CMMFile *mfile, const float *cell)
{
  if (ccp4_file_is_write(mfile->stream)) {
    mfile->cell[0] = cell[0];
    mfile->cell[1] = cell[1];
    mfile->cell[2] = cell[2];
    mfile->cell[3] = cell[3];
    mfile->cell[4] = cell[4];
    mfile->cell[5] = cell[5];
  }
}

/*! Get the grid for the complete cell (X,Y,Z) ordering
 \param mfile (const CMMFile *)
 \param grid (int *) contains the grid dimension on exit (dim 3) */
void ccp4_cmap_get_grid(const CMMFile *mfile, int *grid)
{
  grid[0] = mfile->cell_grid[0];
  grid[1] = mfile->cell_grid[1];
  grid[2] = mfile->cell_grid[2];
}

/*! Set the cell grid dimension.
  Only allowed when file is opened in write mode.
 \param mfile (CMMFile *)
 \param grid (const int *) the cell grid dimension (X,Y,Z) */
void ccp4_cmap_set_grid(CMMFile *mfile, const int *grid)
{
  if (ccp4_file_is_write(mfile->stream)) {
    mfile->cell_grid[0] = grid[0];
    mfile->cell_grid[1] = grid[1];
    mfile->cell_grid[2] = grid[2];
  }
}

/*! Get the stored map origin (rows,sections,columns)
 \param mfile (const CMMFile *)
 \param origin (int *) contains the origin on exit (dim 3) */
void ccp4_cmap_get_origin(const CMMFile *mfile, int *origin)
{
  origin[0] = mfile->origin[0];
  origin[1] = mfile->origin[1];
  origin[2] = mfile->origin[2];
}

/*! Set the stored map origin (rows,sections,columns)
  Only allowed when file is opened in write mode.
 \param mfile (CMMFile *)
 \param origin (const int *) the origin */
void ccp4_cmap_set_origin(CMMFile *mfile, const int *origin)
{
  if (ccp4_file_is_write(mfile->stream)) {
    mfile->origin[0] = origin[0];
    mfile->origin[1] = origin[1];
    mfile->origin[2] = origin[2];
  }
}

/*! Get the stored map axes order (rows,sections,columns)
 where 1=X, 2=Y, 3=Z
 \param mfile (const CMMFile *)
 \param axes_order (float *) contains the ordering on exit (dim 3) */
void ccp4_cmap_get_order(const CMMFile *mfile, int *axes_order)
{
  axes_order[0] = mfile->axes_order[0];
  axes_order[1] = mfile->axes_order[1];
  axes_order[2] = mfile->axes_order[2];
}

/*! Set the stored map axes order (rows,sections,columns)
 where 1=X, 2=Y, 3=Z.
 Only allowed when file is opened in write mode.
 \param mfile (CMMFile *)
 \param axes_order (const float *) the axes ordering */
void ccp4_cmap_set_order(CMMFile *mfile, const int *axes_order)
{
  if (ccp4_file_is_write(mfile->stream)) {
    mfile->axes_order[0] = axes_order[0];
    mfile->axes_order[1] = axes_order[1];
    mfile->axes_order[2] = axes_order[2];
  }
}

/*! Get the stored map dimension (rows,sections,columns)
 \param mfile (const CMMFile *)
 \param map_dim (int *) contains the map dimension on exit (dim 3) */
void ccp4_cmap_get_dim(const CMMFile *mfile, int *map_dim)
{
  map_dim[0] = mfile->map_dim[0];
  map_dim[1] = mfile->map_dim[1];
  map_dim[2] = mfile->map_dim[2];
}

/*! Set the stored map dimension (rows,sections,columns)
 Only allowed when file is opened in write mode before any data
 is written.  
 Note: the row dimension will be overridden during writing
 \param mfile (CMMFile *)
 \param map_dim (const int *) the map dimension */
void ccp4_cmap_set_dim(CMMFile *mfile, const int *map_dim)
{
  if (ccp4_file_is_write(mfile->stream) && !mfile->data.number) {
    mfile->map_dim[0] = map_dim[0];
    mfile->map_dim[1] = map_dim[1];
    mfile->map_dim[2] = map_dim[2]; 
    mfile->data.section_size = map_dim[0]*map_dim[1]*
      ccp4_file_itemsize(mfile->stream);
    mfile->data.block_size = mfile->data.section_size +
                             mfile->data.header_size;
  }
}
/*! Return the spacegroup listed in the map header.
 This is overriden by the symops.
 \param mfile (CMMFile *)
 \return spacegroup number */
int ccp4_cmap_get_spacegroup(const CMMFile *mfile)
{
  return mfile->spacegroup;
}

/*! Set the spacegroup listed in the map header.
 Only allowed when file is opened in write mode.
 \param mfile (CMMFile *) 
 \param spacegroup (int) spacegroup number */
void ccp4_cmap_set_spacegroup(CMMFile *mfile, int spacegroup)
{
  if (ccp4_file_is_write(mfile->stream))
    mfile->spacegroup = spacegroup;
}

/*! Return the datamode
  \param mfile (const CMMFile *)
  \return datamode */
unsigned int ccp4_cmap_get_datamode(const CMMFile *mfile)
{
  return mfile->data_mode;
}

/*! Set the datamode.
  This is only allowed if the file is opened in write mode, and
  no data has been written.
  \param mfile (CMMFile *)
  \param datamode (unsigned int) major mode of map */
void ccp4_cmap_set_datamode(CMMFile *mfile, unsigned int datamode)
{
  if (ccp4_file_is_write(mfile->stream) && !mfile->data.number && 
      datamode <= 6 && datamode != 5) {
    mfile->data_mode = datamode;
    ccp4_file_setmode(mfile->stream, datamode);
    mfile->data.section_size = mfile->map_dim[0]*mfile->map_dim[1]*
      ccp4_file_itemsize(mfile->stream);
    mfile->data.block_size = mfile->data.section_size +
                             mfile->data.header_size;
  }
}

/*! Get the map statistics, including maximum, minimum, mean and standard 
  deviation.  This is only meaningful for datamode FLOAT32.
  \param mfile (const CMMFile *)
  \param min (float *)
  \param max (float *)
  \param mean (double *)
  \param rms (double *) */
void ccp4_cmap_get_mapstats(const CMMFile *mfile, float *min, float* max, 
                           double *mean, double *rms)
{
  double f1,f2,f3;
  *min = mfile->stats.min;
  *max = mfile->stats.max;
  if (ccp4_file_is_write(mfile->stream)  && mfile->close_mode == 0) {
    f1 = (mfile->stats.total != 0) ? mfile->stats.mean / mfile->stats.total : 0;
    f2 = (mfile->stats.total != 0) ? mfile->stats.rms / mfile->stats.total : 0;
    f3 = f2 - f1*f1; 
    *rms = (f3 > 0) ? sqrt(f3) : 0;
    *mean = f1 - (double) mfile->stats.offset;
  } else {
    *mean = mfile->stats.mean;
    *rms = mfile->stats.rms;
  }
}

/*! Set the map statistics, including maximum, minimum, mean and standard 
  deviation.  This is only meaningful for datamode FLOAT32 and the file
  open in write mode.
  \param mfile (CMMFile *)
  \param min (float)
  \param max (float)
  \param mean (double)
  \param rms (double) */
void ccp4_cmap_set_mapstats(CMMFile *mfile, const float min, const float max,
                           const double mean, const double rms)
{
  if (ccp4_file_is_write(mfile->stream)) {
    mfile->stats.min = min;
    mfile->stats.max = max;
    mfile->stats.mean = mean;
    mfile->stats.rms = rms;
  }
}

/*! Set the local header size (in bytes)
  \param mfile (CMMFile *)
  \param size (size_t) header size associated with each section (in bytes) */
void ccp4_cmap_set_local_header(CMMFile *mfile, size_t size)
{
  if (ccp4_file_is_write(mfile->stream) && mfile->data.number == 0) {
    mfile->data.header_size = size;
    mfile->data.block_size = mfile->data.section_size + mfile->data.header_size;
  }
  return;
}

/*! Return the local header size
 \param mfile (CMMFile *)
 \return header size associated with each section (in bytes) */
size_t ccp4_cmap_get_local_header(CMMFile *mfile)
{
  return mfile->data.header_size;
}
    
