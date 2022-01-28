/*
     cmap_close.c: close map file
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
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include "cmaplib.h"
#include "cmap_header.h"
#include "cmap_labels.h"
#include "cmap_errno.h"

/*! Close the file.
 In write mode the header is output, along with the machine
 stamp.  In read mode the file is just closed. 
 Write mode supports ways of updating the map statistics ( 
 only active for FLOAT32).
 /param mfile (CMMFile *)
 /return void */
void ccp4_cmap_close(CMMFile *mfile)
{
  int i;
  
  if ( mfile == NULL) 
    return;

  if (ccp4_file_is_write(mfile->stream) ) {
    if ( mfile->data_mode == FLOAT32) {
      switch (mfile->close_mode) {
      case 1:
        break;
      case 2:
	mfile->stats.offset = 0.0f;
      case 0:
      default:
        if (mfile->stats.total != 0) {
          mfile->stats.mean /= mfile->stats.total;
          mfile->stats.rms /= mfile->stats.total;
          mfile->stats.rms -= mfile->stats.mean*mfile->stats.mean;
          mfile->stats.rms = (mfile->stats.rms > 0) ? sqrt(mfile->stats.rms) : 0;
          mfile->stats.mean += (double) mfile->stats.offset;
        }
        break;
      }
    }
    write_mapheader(mfile);
    write_maplabels(mfile);
    ccp4_file_warch(mfile->stream);
  }
  ccp4_file_close(mfile->stream);
  for (i=0 ; i != mfile->labels.number ; i++)
    if (mfile->labels.labels[i] != NULL)
      free(mfile->labels.labels[i]);
  free(mfile);
}

/*! Set the close mode:
    0: calculate based on stored values (default)
    1: just dump the current values
 /param mfile (CMMFile *)
 /param mask (unsigned int) close mode
 /return void */
void ccp4_cmap_closemode(CMMFile *mfile, unsigned int mask)
{
  mfile->close_mode = mask;
}
