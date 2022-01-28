/*
     vmslibrary.c
     Copyright (C) 1999  Martyn Winn

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
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <types.h>
#include <descrip.h>
#include <stat.h>

#define MAXFLEN       500       /* the maximum length of a filename in CCP4 */

/* prototype definitions */
void HGETLIMITS (int *IValueNotDet, float *ValueNotDet);
static size_t flength (char *s, int len);
void CMKDIR (struct dsc$descriptor_s *path, struct dsc$descriptor_s *cmode, 
          int *result);
void CCHMOD (struct dsc$descriptor_s *path, struct dsc$descriptor_s *cmode, 
          int *result);

void HGETLIMITS (int *IValueNotDet, float *ValueNotDet)
{
  *IValueNotDet = INT_MAX;
  *ValueNotDet  = FLT_MAX;
}

static size_t flength (char *s, int len)
{
  while (s[--len] == ' ');
  return (++len);
}

/* Wrap-around for mkdir function. Returns 0 if successful, 1 if directory already exists,  */
/* and -1 if other error.                                                                   */
void CMKDIR (struct dsc$descriptor_s *path, struct dsc$descriptor_s *cmode, 
      int *result)
{ size_t Length;
  char name[MAXFLEN];
  mode_t mode;

  /* truncate path to MAXFLEN - 1 characters, MAXFLEN defined in library.h */
  Length = flength (path->dsc$a_pointer, path->dsc$w_length);
  if (Length > (size_t) MAXFLEN) Length = (size_t) MAXFLEN - 1;
  (void) strncpy (name, path->dsc$a_pointer, Length);
  name[Length] = '\0'; 

/* Possible modes (see stat.h)
  Currently pass 3-character string and interpret as octal.
  Try also S_IRWXU, S_IRWXG, etc. */
  sscanf(cmode->dsc$a_pointer,"%o",&mode);
   
  *result = mkdir(name,mode); 

  if (*result == -1) {
/* Distinguish directory-exists error from others, since usually not a problem. */
    if (errno == EEXIST) {
      *result = 1;
    }
  }
    
}

void CCHMOD (struct dsc$descriptor_s *path, struct dsc$descriptor_s *cmode, 
     int *result)
{ size_t Length;
  char name[MAXFLEN];
  mode_t mode;

  /* truncate path to MAXFLEN - 1 characters, MAXFLEN defined in library.h */
  Length = flength (path->dsc$a_pointer, path->dsc$w_length);
  if (Length > (size_t) MAXFLEN) Length = (size_t) MAXFLEN - 1;
  (void) strncpy (name, path->dsc$a_pointer, Length);
  name[Length] = '\0'; 

/* Possible modes (see stat.h)
  Currently pass 3-character string and interpret as octal.
  Try also S_IRWXU, S_IRWXG, etc. */
  sscanf(cmode->dsc$a_pointer,"%o",&mode);

  *result = chmod(name,mode); 
}
