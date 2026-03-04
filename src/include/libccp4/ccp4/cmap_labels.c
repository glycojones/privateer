/*
     cmap_labels.c: read and write map header labels
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
#include "cmaplib.h"
#include "cmap_labels.h"
#include "cmap_errno.h"

/*! Internal: read the labels from file header and copy into char * array
  Called when the file is opened in read mode.
  \param mfile (CMMFile *)
  \return 1 on succes */
int parse_maplabels(CMMFile *mfile)
{
  char buffer[81], *cptr;
  const unsigned int n_byt_label = 80U, max_label = 10U;
/*  const unsigned int labels_offset = 224U; */
  int i;
  
/*  ccp4_file_seek(mfile->stream,labels_offset,SEEK_SET); */
  for (i=0 ; i!=mfile->labels.number ; i++) {
    ccp4_file_readchar(mfile->stream,(uint8 *) buffer,n_byt_label);
    cptr = buffer+n_byt_label;
    while (cptr> buffer && *--cptr == ' ');
    *(++cptr) = '\0';
    mfile->labels.labels[i] = strdup(buffer);
  }
  ccp4_file_raw_seek(mfile->stream,(max_label-mfile->labels.number)
		     *n_byt_label,
		     SEEK_CUR);
  return 1;
}

/*! Internal: dump the labels char * array to file, offset at 224 bytes.
  Called when the file is opened or closed in write mode, immediately after the
  header is written.
 \param mfile (const CMMFile *)
 \return 1 on success, 0 on failure */
int write_maplabels(const CMMFile *mfile)
{
  char buffer[80];
/*  const unsigned int labels_offset = 224U; */
  int i, result = 0;
  size_t slen;

/*  ccp4_file_seek(mfile->stream,labels_offset,SEEK_SET);  */
  for (i=0 ; i != mfile->labels.number ; i++) {
    memset(buffer,' ',80U);
    slen = strlen(mfile->labels.labels[i]);
    if (slen > 80U) slen = 80U;
    strncpy(buffer,mfile->labels.labels[i],slen);
    result += ccp4_file_writechar(mfile->stream,(uint8 *) buffer,80U);
  }
  memset(buffer,' ',80U);
  while(i != 10) {
    result += ccp4_file_writechar(mfile->stream,(uint8 *) buffer,80U);
    i++;
  }
  return (result == 800) ? 1 : 0 ;
}

/*! Set the label in the map header.  Headers are 80 characters long.
  The labels are written to the file when it is closed. Therefore,
  the file must be in write mode.
  If label == NULL the element corresponding to posn is removed.
  The number of labels is recalculated on each call.
 \param mfile (CMMFile *)
 \param label (const char *) the C-style character array
 \param posn (int) the label number (C-style, 0 -> 9) 
 \return number of label effected, or EOF */
int ccp4_cmap_set_label(CMMFile *mfile, const char *label, int posn)
{  
  int i,j;
  
  if (mfile == NULL) {
    ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_NoChannel),
                "ccp4_cmap_set_label",NULL);
    return (EOF);}

  if (ccp4_file_is_write(mfile->stream) == 0) {
    ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_WriteFail),
                "ccp4_cmap_label_set",NULL);
    return (EOF);}    
  
/*posn must be between 0 and 9 */
  if (posn < 0) {
    posn = 0;
  } else if (posn > mfile->labels.number) {
    posn = mfile->labels.number;
  }

  if (mfile->labels.labels[posn] != NULL) 
    free(mfile->labels.labels[posn]);

/* if label == NULL reset the value and compress set */    
  if (label == NULL) {
    mfile->labels.labels[posn] = NULL;
    for ( i=posn ; i!=10 ; i++)
      if (mfile->labels.labels[i] == NULL)
        for ( j=i+1 ; j!=10; j++)
          if (mfile->labels.labels[j] != NULL) {
            mfile->labels.labels[i] = mfile->labels.labels[j];
            mfile->labels.labels[j] = NULL;
            break;
          }
    }
  else
    mfile->labels.labels[posn] = strdup(label);

/* recalculate number */
  for ( i=0 ; i!=10 ; i++)
    if (mfile->labels.labels[i] == NULL)
      break;
  mfile->labels.number = i;
  
  return posn;
}

/*! Get the label corresponding to position posn
  \param mfile (const CMMFile *)
  \param posn (int) desired label number
  \return pointer to label posn */
char *ccp4_cmap_get_label(const CMMFile *mfile, int posn)
{
  char *label;
  if (mfile == NULL) {
    ccp4_signal( CCP4_ERRLEVEL(3) | CMAP_ERRNO(CMERR_NoChannel),
                "ccp4_cmap_get_label",NULL);
    return (NULL);}

  if (posn < 0 || posn >= mfile->labels.number) 
    label = NULL;
  else 
    label = mfile->labels.labels[posn];

  return label;
}

/*! Return the number of labels.
 \param mfile (CMMFile *)
 \return the number of labels */
int ccp4_cmap_number_label(const CMMFile *mfile)
{
  return mfile->labels.number;
}

/*! Get the label corresponding to the title
    wrapping ccp4_cmap_get_label.
 \param mfile (const CMMFile *)
 \return pointer to label 0, or NULL */
char *ccp4_cmap_get_title(const CMMFile *mfile)
{
  return ccp4_cmap_get_label(mfile, 0);
}

/*! Set the label corresponding to the title,
    wrapping ccp4_cmap_set_label
 \param mfile (CMMFile *)
 \param label
 \return 0 or EOF on failure */
int ccp4_cmap_set_title(CMMFile *mfile, const char *title)
{
  return ccp4_cmap_set_label(mfile, title, 0);
}
