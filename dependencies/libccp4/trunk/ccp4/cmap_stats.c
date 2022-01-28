/*
     cmap_stats.c: deal with map statistics.
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
#include "cmap_stats.h"
#include "cmap_errno.h"

/*! Internal: use floats in range section_begin to section_end
  to update the map statistics.
  \param stats (CMMFile_Stats *)
  \param section_begin (void *) start of section
  \param section_end (void *) one past end-of-section
  \return total of map elements so far */
int stats_update(CMMFile_Stats *stats, void *section_begin,
                         void *section_end)
{        
  float *ufp = (float *) section_begin;
  double val;
      
  if (stats->total == 0 && *ufp < -1.0e10 ) {         
    stats->offset = *ufp;
  } 
  while (ufp < (float *) section_end) {
    val = (double) (*ufp - stats->offset);
    stats->mean += val;
    stats->rms += val * val;
    stats->min = MIN( stats->min, *ufp);
    stats->max = MAX( stats->max, *ufp);
    
    ufp++;
  }
  
  stats->total += (float *)section_end - (float *)section_begin;
                
  return (stats->total);
}

