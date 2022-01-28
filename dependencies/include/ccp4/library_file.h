/*
     library_file.h: header file for library_file.c
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

/** @file library_file.h
 *  Functions for file i/o.
 *  Charles Ballard
 */

#ifndef __CCP4_LIB_FILE
#define __CCP4_LIB_FILE

#include "ccp4_sysdep.h"
#include "ccp4_types.h"

#ifdef __cplusplus
namespace CCP4 {
extern "C" {
#endif

/** Generic CCP4 file. */
typedef struct _CFileStruct CCP4File;

struct _CFileStruct {
  char *name;
  FILE *stream;
  int fd;
  unsigned int read : 1;
  unsigned int write : 1;
  unsigned int append : 1;
  unsigned int binary : 1;
  unsigned int scratch : 1 , : 3;
  unsigned int buffered : 1;
  unsigned int sync : 1, : 6;
  unsigned int direct : 1, : 7;
  unsigned int open : 1;
  unsigned int own : 1;
  unsigned int last_op : 2;
  unsigned int getbuff : 1, : 4;
  int iostat;
  unsigned int mode : 8;
  unsigned int itemsize : 8;
  unsigned int iconvert : 8;
  unsigned int fconvert: 8;
  off_t length;
  off_t loc;
  size_t stamp_loc;
  int (*_read) (CCP4File *, uint8 *, size_t);
  int (*_write) (CCP4File *, const uint8 *, size_t);
  char buff[8];
  void *priv;
};


CCP4File *ccp4_file_open (const char *, const int);

CCP4File *ccp4_file_open_file (const FILE *, const int);

CCP4File *ccp4_file_open_fd (const int, const int);

int ccp4_file_rarch ( CCP4File*);

int ccp4_file_warch ( CCP4File*);

int ccp4_file_close ( CCP4File*);

int ccp4_file_mode ( const CCP4File*);

int ccp4_file_setmode ( CCP4File*, const int);

int ccp4_file_setstamp( CCP4File *, const size_t);

int ccp4_file_itemsize( const CCP4File*);

int ccp4_file_setbyte( CCP4File *, const int);

int ccp4_file_byteorder( CCP4File *);

int ccp4_file_is_write(const CCP4File *);

int ccp4_file_is_read(const CCP4File *);

int ccp4_file_is_append(const CCP4File *);

int ccp4_file_is_scratch(const CCP4File *);

int ccp4_file_is_buffered(const CCP4File *);

int ccp4_file_status(const CCP4File *);

char *ccp4_file_name( CCP4File *);

int ccp4_file_read ( CCP4File*, uint8 *, size_t);

int ccp4_file_readcomp ( CCP4File*, uint8 *, size_t);

int ccp4_file_readshortcomp ( CCP4File*, uint8 *, size_t);

int ccp4_file_readfloat ( CCP4File*, uint8 *, size_t);

int ccp4_file_readint ( CCP4File*, uint8 *, size_t);

int ccp4_file_readshort ( CCP4File*, uint8 *, size_t);

int ccp4_file_readchar ( CCP4File*, uint8 *, size_t);

int ccp4_file_write ( CCP4File*, const uint8 *, size_t);

int ccp4_file_writecomp ( CCP4File*, const uint8 *, size_t);

int ccp4_file_writeshortcomp ( CCP4File*, const uint8 *, size_t);

int ccp4_file_writefloat ( CCP4File*, const uint8 *, size_t);

int ccp4_file_writeint ( CCP4File*, const uint8 *, size_t);

int ccp4_file_writeshort ( CCP4File*, const uint8 *, size_t);

int ccp4_file_writechar ( CCP4File*, const uint8 *, size_t);

int ccp4_file_seek ( CCP4File*, long, int);

void ccp4_file_rewind ( CCP4File*);

void ccp4_file_flush (CCP4File *);

long ccp4_file_length ( CCP4File*);

long ccp4_file_tell ( CCP4File*);

int ccp4_file_feof(CCP4File *);

void ccp4_file_clearerr(CCP4File *);

void ccp4_file_fatal (CCP4File *, char *);

char *ccp4_file_print(CCP4File *, char *, char *);

int ccp4_file_raw_seek( CCP4File *, long, int);
int ccp4_file_raw_read ( CCP4File*, char *, size_t);
int ccp4_file_raw_write ( CCP4File*, const char *, size_t);
int ccp4_file_raw_setstamp( CCP4File *, const size_t);
#ifdef __cplusplus
}
}
#endif

#endif  /* __CCP4_LIB_FILE */
