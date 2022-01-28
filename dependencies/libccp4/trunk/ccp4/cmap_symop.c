/*
     cmap_symop.c: set and fetch symmetry operations in map header
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
#include <stdlib.h>
#include "cmaplib.h"
#include "cmap_errno.h"

/*! Return the number of symops (estimated as the size/80)
 \param mfile (const CMMFile *)
 \return number of symops */
int ccp4_cmap_num_symop(const CMMFile *mfile)
{
  if (mfile == NULL)
    return 0;
  return (mfile->symop.number);
}

/*! navigate around the symops, seeking in 80 byte units
 The result must lie within the symop strings in the file.
 \param mfile (CMMFile *)
 \param isymop (int) the number of the symop "string" of interest
 \param whence (unsigned int) mode of seek
 \return symop string number or EOF */
int ccp4_cmap_seek_symop(CMMFile *mfile, int isymop, unsigned int whence)
{
  const int n_byt_symop = 80;
  div_t symops;
  int result = EOF;
  
  if (ccp4_file_is_read(mfile->stream) == 0)
    return EOF;
    
  switch (whence) {
  case SEEK_SET:
    if (isymop < 0 || isymop > mfile->symop.number) 
      ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_ParamError), 
		   "ccp4_cmap_seek_symop",
		   NULL);
    else 
      result = ccp4_file_raw_seek(mfile->stream, mfile->symop.offset +
                              isymop*n_byt_symop, whence);
    break;
  case SEEK_END:
    if (isymop > 0 || abs(isymop) > mfile->symop.number )
      ccp4_signal( CCP4_ERRLEVEL(2) |  CMAP_ERRNO(CMERR_ParamError), 
		   "ccp4_cmap_seek_symop",
		   NULL);
    else
      result = ccp4_file_raw_seek(mfile->stream, mfile->symop.offset + 
                              mfile->symop.size + isymop*n_byt_symop, 
                              SEEK_SET);
    break;
  case SEEK_CUR:
    symops = div(ccp4_file_tell(mfile->stream) - mfile->symop.offset,n_byt_symop);
    if (symops.quot < 0 || symops.quot >= mfile->symop.number ||
        symops.quot + isymop < 0 || symops.quot + isymop >= mfile->symop.number)
      ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_ParamError), 
		   "ccp4_cmap_seek_symop",
		   NULL);
    else
      result = ccp4_file_raw_seek(mfile->stream,(isymop > 0) ?
                              (n_byt_symop - symops.rem + n_byt_symop * (isymop-1)) :
                              (n_byt_symop * isymop -symops.rem), SEEK_CUR);
  }
  return (result == EOF) ? EOF : (result - mfile->symop.offset)/n_byt_symop;
}

/*! get a symop string of 80 characters 
 \param mfile (CMMFile *)
 \param buffer (char *) array of bytes which will contain the symop string.
 This must be at least 81 characters long (including space for null terminator).
 \return 1 on success, 0 if no symops, EOF on failure */
int ccp4_cmap_get_symop(CMMFile *mfile, char *buffer)
{
  const int n_byt_symop = 80;
  off_t file_posn;

  if ( mfile->symop.size == 0)  {
    ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_SymErr),
		 "cmap_get_symop",NULL);
    return (0);}

  file_posn = ccp4_file_tell(mfile->stream);
  if (file_posn < mfile->symop.offset ||
      file_posn > mfile->symop.offset + mfile->symop.size) {
    ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_SymErr),
		 "cmap_get_symop",NULL);
    return (EOF);}

  if (ccp4_file_readchar(mfile->stream, (uint8 *) buffer, n_byt_symop) != n_byt_symop) {
    ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_ReadFail),
		 "cmap_get_symop",NULL);
    return (EOF);
  }
  buffer[n_byt_symop] = '\0';
  return (1);
}

/*! write symops to file.
 This wraps a raw write.  It is up to the calling program to
 ensure the positioning (effectively assume appends).  Writing
 is blocked if data has alread been written to the file.  80
 bytes of continuous memory is written to the file.
 \param mfile (CMMFile *)
 \param symop (const char *) character array containing the 
 symop string (at least 80 characters in length
 \return 1 on success, EOF on failure */
int ccp4_cmap_set_symop(CMMFile *mfile, const char *symop)
{
  const int n_byt_symop = 80;
  char buffer[80];
  memset(buffer,' ',80U);
  memcpy(buffer, symop, (strlen(symop) > n_byt_symop) ?
         n_byt_symop : strlen(symop) );   
  if (ccp4_file_is_write(mfile->stream) && mfile->data.number == 0) {
    if (ccp4_file_writechar(mfile->stream, (uint8 *) buffer, n_byt_symop) != n_byt_symop) {
      ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_WriteFail),
		   "cmap_set_symop",NULL);
      return (EOF);
    }
  mfile->symop.number++;
  mfile->symop.size += n_byt_symop;
  mfile->data.offset = mfile->symop.offset + mfile->symop.size;
}
  return (1);
}
