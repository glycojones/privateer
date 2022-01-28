/*
     cmap_header.h: header file for cmap_header.c
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
#ifndef __GUARD_MAPLIB_HEADER
#define __GUARD_MAPLIB_HEADER

#ifdef __cplusplus
extern "C" {
#endif

int parse_mapheader(CMMFile *mfile);

int write_mapheader(CMMFile *mfile);

#ifdef __cplusplus
}
#endif

#endif  /* __GUARD_MAPLIB_HEADER */
