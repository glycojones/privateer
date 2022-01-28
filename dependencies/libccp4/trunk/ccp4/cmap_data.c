/*
     cmap_data.c: read and write map sections.
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
#include <string.h>
#include <stdlib.h>
#include "cmaplib.h"
#include "cmap_data.h"
#include "cmap_stats.h"
#include "cmap_errno.h"

/*! Internal: return the an estimate of the number of sections in the map
  file based upon the length.
  Update mfile->data.number as a side effect.
  \param mfile (CMMFile *)
  \return number of sections according to length-data/section_size */
int number_sections(CMMFile *mfile)
{
  div_t sections;
  
  sections = div(ccp4_file_length(mfile->stream)-mfile->data.offset, 
                 mfile->data.block_size);
  
  return mfile->data.number = sections.quot;
}

/*! seek among the map sections.  The units are of size block_size.
 \param mfile (CMMFile *)
 \param sec (int) section number
 \param whence (unsigned int) SEEK_SET, SEEK_CUR or SEEK_END
 \return offset in file, or EOF */
int ccp4_cmap_seek_section(CMMFile *mfile, int sec, unsigned int whence)
{
  size_t curr_posn;
  div_t secs;
  int result = EOF;
  
  if ( mfile == NULL ) {
    ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_NoChannel),
		 "ccp4_cmap_seekdata",NULL);
    return EOF; }

  switch (whence) {
  case SEEK_SET:
    if ( ccp4_file_is_read(mfile->stream) && 
         ( sec < 0 || sec > mfile->data.number) )
      ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_ParamError),
		 "ccp4_cmap_seek_section",NULL);
    else 
      result = ccp4_file_raw_seek(mfile->stream, mfile->data.offset +
                              sec * mfile->data.block_size, SEEK_SET);
    break;
  case SEEK_END:
    if ( ccp4_file_is_read(mfile->stream) &&
         ( sec > 0 || abs(sec) > mfile->data.number) )
      ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_ParamError),
		 "ccp4_cmap_seek_section",NULL);
    else 
      result = ccp4_file_raw_seek(mfile->stream, sec * mfile->data.block_size,
                              SEEK_END);
    break;
  case SEEK_CUR:
    curr_posn = ccp4_file_tell(mfile->stream);
    secs = div(curr_posn - mfile->data.offset,mfile->data.block_size);
    if ( ccp4_file_is_read(mfile->stream) &&
       ( (secs.quot + sec) < 0 || (secs.quot + sec) >= mfile->data.number) )
      ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_ParamError),
		 "ccp4_cmap_seek_section",NULL);
    else 
      result = ccp4_file_raw_seek(mfile->stream,
                              (sec > 0) ? (mfile->data.block_size - secs.rem +
                                         (sec - 1)*mfile->data.block_size) :
                                        (sec*mfile->data.block_size - secs.rem),
                              SEEK_CUR);
  }
  return (result == EOF) ? EOF :
    ((result - mfile->data.offset)/mfile->data.block_size);
}

/*! write map section to file.  
  Note:  this wraps a raw write, with no location checking.  It is
  therefore the responsibility of the calling program to ensure that
  everything is correct.  Effectively assume appending to file.
  \param mfile (CMMFile *)
  \param section (const void *)
  \return 1 on success, 0 on failure */
int ccp4_cmap_write_section(CMMFile *mfile, const void *section)
{
  int result=0;
  size_t write_dim;
  
  if (mfile == NULL || section == NULL) {
    ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_NoChannel),
		 "ccp4_cmap_write_section",NULL);
    return 0; }

  if (!ccp4_file_is_write(mfile->stream)) {
    ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_WriteFail),
		 "ccp4_cmap_write_section",NULL);
    return 0; }
    
  write_dim = mfile->map_dim[0] * mfile->map_dim[1];
  result = ccp4_file_write(mfile->stream, section, write_dim);

/* note that we have started writing */
  mfile->data.number++;
  
  if (result != write_dim)
        ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_WriteFail),
		 "ccp4_cmap_write_section",NULL);
  else
    if (mfile->data_mode == FLOAT32)
      stats_update(&mfile->stats, (float *)section, 
                   (float *)section+write_dim);
  
  return (result == write_dim) ? 1 : 0;
}

/*! read current map section from file to section.
  Some checking is performed to ensure we are at the start of a
  legitimate map section.
  \param mfile (CMMFile *)
  \param section (void *) array large enough to hold the map section
  \return 1 on success, 0 on failure */
int ccp4_cmap_read_section(CMMFile *mfile, void *section)
{
  int result = 0;
  div_t secs;
  off_t curr_posn;
  const off_t data_offset = 0;
  size_t read_dim;
  
  if (mfile == NULL ) {
    ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_NoChannel),
		 "ccp4_cmap_read_section",NULL);
    return 0; }
    
  if (!ccp4_file_is_read(mfile->stream)) {
    ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_ReadFail),
		 "ccp4_cmap_read_section",NULL);
    return 0; }

  curr_posn = ccp4_file_tell(mfile->stream);
              
  secs = div(curr_posn - mfile->data.offset, 
             mfile->data.block_size);

/* ensure legit section (although rely upon EOF ) */  
  if (secs.quot < 0 || secs.rem < 0) {
    ccp4_file_raw_seek(mfile->stream, mfile->data.offset, SEEK_SET);
    secs.quot = 0;
  } else if( secs.rem > data_offset && secs.rem < mfile->data.section_size ) 
    ccp4_file_raw_seek(mfile->stream, - secs.rem, SEEK_CUR);
  else if ( secs.rem >= mfile->data.section_size ) {
    ccp4_file_raw_seek(mfile->stream, (mfile->data.block_size-secs.rem), SEEK_CUR);
    secs.quot++; }

  read_dim = mfile->map_dim[0] * mfile->map_dim[1];
/* do not read if at end */
  if (secs.quot < 0 || secs.quot < mfile->data.number)
    result = ccp4_file_read(mfile->stream, section, read_dim);

  if (result != read_dim)
    ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_ReadFail),
		 "ccp4_cmap_read_section",NULL);

  return (result == read_dim) ? 1 : 0;
}

/*! read current section header (character array)
  After reading we are at the end of the local header
  \param mfile (CMMFile *)
  \param header (char *) character array large enough to hold
  the local header (raw read so not string)
  \return 1 on success, 0 on failure */
int ccp4_cmap_read_section_header(const CMMFile *mfile, char *header)
{
  int result;
  div_t secs;
  
  if (mfile == NULL || header == NULL) {
    ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_NoChannel),
		 "ccp4_cmap_read_section_header",NULL);
    return EOF; }

  if (!ccp4_file_is_read(mfile->stream)) {
    ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_ReadFail),
		 "ccp4_cmap_read_section header",NULL);
    return EOF; }
    
  if ( mfile->data.header_size == 0) return (0);

  result = ccp4_file_tell(mfile->stream);
  secs = div(result - mfile->data.offset, mfile->data.block_size);
  
  if ( secs.quot < 0 || secs.quot >= mfile->data.number ) return (0);
 
/* navigate to nearest header */ 
  if ( secs.rem != mfile->data.section_size)
    ccp4_file_raw_seek(mfile->stream,(mfile->data.section_size 
				      - secs.rem), SEEK_CUR);     
  
  if ( (result = ccp4_file_readchar( mfile->stream, (uint8 *) header, 
        mfile->data.header_size)) != mfile->data.header_size) 
    ccp4_signal(ccp4_errno,
		"ccp4_cmap_read_section_header",
		NULL);
            
  return (result == mfile->data.header_size) ? 1 : 0;
}

/*! write the local section header to the file.  This must be of
  size mfile->data.header.size.  
  Note: no checking is done so it is up to the calling program
  to ensure that the file is in the correct location.  As seeking
  is turned off, this assumes we are appending to the file.
 \param mfile (CMMFile *)
 \param header (const char *) the local header character array
 (not necessarily a string)
 \return number of bytes written or EOF */
int ccp4_cmap_write_section_header(CMMFile *mfile, const char *header)
{
  char *output;
  int result;
  
  if (mfile == NULL ) {
    ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_NoChannel),
		 "ccp4_cmap_write_section_header",NULL);
    return EOF; }
    
  if (!ccp4_file_is_write(mfile->stream)) {
    ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_NoChannel),
		 "ccp4_cmap_write_section_header",NULL);
    return EOF; }
   
  if ( mfile->data.header_size == 0)  return (0);

  output = (char *) malloc(mfile->data.header_size);
  memset(output,' ', mfile->data.header_size);
  if (header) memcpy(output, header,mfile->data.header_size);

  if ( (result = ccp4_file_writechar( mfile->stream, (uint8 *) output, 
                                     mfile->data.header_size)) 
	 != mfile->data.header_size) 
    ccp4_signal(ccp4_errno,
		"ccp4_cmap_write_section_header",
		NULL);
            
  return (result == mfile->data.header_size) ? 1 : 0;
}

/*! seek a row within a map section
 \param mfile (CMMFile *)
 \param row (int) 
 \param whence (unsigned int) SEEK_SET, SEEK_END, SEEK_CUR
 \return offset in file or EOF */
int ccp4_cmap_seek_row(CMMFile *mfile, int row, unsigned int whence)
{
  size_t curr_posn;
  div_t secs, rows;
  int result = EOF;
  size_t item_size;
  
  if ( mfile == NULL ) {
    ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_NoChannel),
		 "ccp4_cmap_seek_row",NULL);
    return EOF; }

  item_size = ccp4_file_itemsize(mfile->stream);
  curr_posn = ccp4_file_tell(mfile->stream);
  secs = div(curr_posn - mfile->data.offset,mfile->data.block_size);
  
  switch (whence) {
  case SEEK_SET:
    if ( row < 0 || row >= mfile->map_dim[1]) 
      ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_ParamError),
		 "ccp4_cmap_seek_row",NULL);
    else 
      result = ccp4_file_raw_seek(mfile->stream, mfile->data.offset +
                              (secs.quot * mfile->data.block_size +
                              row * mfile->map_dim[0]*item_size),
                              SEEK_SET);
    break;
  case SEEK_END:
    if ( row >= 0 || abs(row) > mfile->map_dim[1]) 
      ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_ParamError),
		 "ccp4_cmap_seek_row",NULL);
    else 
      result = ccp4_file_raw_seek(mfile->stream, mfile->data.offset +
                              (secs.quot * mfile->data.block_size +
                              mfile->data.section_size +
                              row * mfile->map_dim[0]*item_size),
                              SEEK_SET);
     break;
  case SEEK_CUR:
    rows = div(secs.rem,mfile->map_dim[0]*item_size);
    if ( (rows.quot + row) < 0 || (rows.quot + row) >= mfile->data.number)
       ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_ParamError),
		   "ccp4_cmap_seek_row",NULL);
    else 
      result = ccp4_file_raw_seek(mfile->stream, 
                              ( row > 0) ? (mfile->map_dim[0]*item_size - rows.rem
                               + (row-1)*mfile->map_dim[0]*item_size) :
                               ( row*mfile->map_dim[0]*item_size - rows.rem),
                              SEEK_CUR);
  }
  return (result);
}

/*! write map row to file.  
  Note:  this wraps a raw write, with no location checking.  It is
  therefore the responsibility of the calling program to ensure that
  everything is correct.  Effectively assume appending to file.
  \param mfile (CMMFile *)
  \param row (const void *) data to be written
  \return 1 on success, 0 on failure */
int ccp4_cmap_write_row(CMMFile *mfile, const void *row)
{
  int result=0;
  
  if (mfile == NULL || row == NULL) {
    ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_NoChannel),
		 "ccp4_cmap_write_row",NULL);
    return EOF; }

  if (!ccp4_file_is_write(mfile->stream)) {
    ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_WriteFail),
		 "ccp4_cmap_write_row",NULL);
    return EOF; }
   
  result = ccp4_file_write(mfile->stream, row, mfile->map_dim[0]);

/* note that we have started writing */
  mfile->data.number++;

  if (result != mfile->map_dim[0])
        ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_WriteFail),
		 "ccp4_cmap_write_row",NULL);
  else
    if (mfile->data_mode == FLOAT32)
      stats_update(&mfile->stats, (float *)row, (float *)row+mfile->map_dim[0]);
  
  return (result == mfile->map_dim[0]) ? 1 : 0;
}

/*! read current map section from file to section.
  Some checking is performed to ensure we are at the start of a
  legitimate map row.
  \param mfile (CMMFile *)
  \param row (void *) array large enough to hold the map row
  \return 1 on success, 0 on failure */
int ccp4_cmap_read_row(CMMFile *mfile, void *row)
{
  int result = 0, item_size;
  div_t secs, rows;
  off_t curr_posn;
  
  if (mfile == NULL) {
    ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_NoChannel),
		 "ccp4_cmap_read_row",NULL);
    return EOF; }
    
  if (!ccp4_file_is_read(mfile->stream) || row == NULL) {
    ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_ReadFail),
		 "ccp4_cmap_read_row",NULL);
    return EOF; }

  item_size = ccp4_file_itemsize(mfile->stream); 
  curr_posn = ccp4_file_tell(mfile->stream);
              
  secs = div(curr_posn - mfile->data.offset, 
             mfile->data.block_size);
  rows = div(secs.rem, mfile->map_dim[0]*item_size);
    
  if (secs.quot < 0 || secs.rem < 0)
    ccp4_file_raw_seek(mfile->stream, mfile->data.offset, SEEK_SET);
  else if(rows.quot >= mfile->map_dim[1] )
    ccp4_file_raw_seek(mfile->stream, (mfile->data.block_size
                       - secs.rem), SEEK_CUR);
  else if( rows.rem != 0) 
    ccp4_file_raw_seek(mfile->stream, ( - secs.rem), SEEK_CUR);

  result = ccp4_file_read(mfile->stream, row, mfile->map_dim[0]);

  if (result != mfile->map_dim[0])
    ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_ReadFail),
		 "ccp4_cmap_read_row",NULL);

  return (result == mfile->map_dim[0]) ? 1 : 0;
}

/*! raw seek in items
  \param mfile (CMMFile *)
  \param offset (int) number of items
  \param whence (unsigned int) SEEK_SET, SEEK_CUR, SEEK_END;
  \return 0 on success, EOF on failure */
int ccp4_cmap_seek_data(CMMFile *mfile, int offset, unsigned int whence)
{
  int result = EOF;
  
  if ( mfile == NULL ) {
    ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_NoChannel),
		 "ccp4_cmap_seekdata",NULL);
    return (result); }

  if ((result = ccp4_file_seek( mfile->stream, offset, whence)) == -1) 
    ccp4_signal(ccp4_errno, "ccp4_cmap_seek_data",NULL);

  return (result);
}

/*! raw write of nelements items to file, according to the datamode,
   at current location
 \param mfile (const CMMFile *)
 \param section (void *) values written, should contain at least 
 nelements items
 \param n_items (int) number of items to be written
 \return number of items written or EOF */
int ccp4_cmap_write_data(CMMFile *mfile, const void *items, int n_items)
{
  int result=0;

  if (mfile == NULL || items == NULL) {
    ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_NoChannel),
		 "ccp4_cmap_write_data",NULL);
    return EOF; }
  
  if (ccp4_file_is_write(mfile->stream)) {
    result = ccp4_file_write(mfile->stream, (uint8 *) items, n_items);
    if (result != n_items)
      ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_WriteFail),
		 "ccp4_cmap_write_data",NULL);
     else if (mfile->data_mode == FLOAT32)
      stats_update(&mfile->stats, (float *)items, (float *)items+result);
  }
  
  return (result);
}

/*! raw read of nelements items from file according to the datamode
   at current location
 \param mfile (const CMMFile *)
 \param items (void *) values read to here, so should have enough space
 for nelements items
 \param n_items (int) number of items to be read
 \return number of items read or EOF */
int ccp4_cmap_read_data(const CMMFile *mfile, void *items, int n_items)
{
  int result=0;

  if (mfile == NULL || items == NULL) {
    ccp4_signal( CCP4_ERRLEVEL(2) | CMAP_ERRNO(CMERR_NoChannel),
		 "ccp4_cmap_read_data",NULL);
    return EOF; }
      
  if (ccp4_file_is_read(mfile->stream)) 
    result = ccp4_file_read(mfile->stream, (uint8 *) items, n_items);

  return (result);
}

