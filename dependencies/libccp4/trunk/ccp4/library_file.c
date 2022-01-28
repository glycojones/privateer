/*
     library_file.c: functions for file i/o
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

/** @file library_file.c
 *  Functions for file i/o.
 *  Charles Ballard
 */

#include<string.h>
#include <limits.h>
#include <fcntl.h>
#if defined _MSC_VER
#include <io.h>
#endif
#include "library_file.h"
#include "ccp4_errno.h"
#include "ccp4_file_err.h"
/* rcsid[] = "$Id: library_file.c,v 1.26 2012/08/20 12:21:16 gxg60988 Exp $" */
                                                        
static uint16 nativeIT = NATIVEIT; /* machine integer type */ 
static uint16 nativeFT = NATIVEFT; /* machine float type */

static int _item_sizes[] = {
  (int) sizeof (char),          /* 0: bytes */
  (int) sizeof (short int),     /* 1: (integer) half words */
  (int) sizeof (float),         /* 2: reals/words */
  (int) sizeof (int),           /* 3: `short complex' (pairs of half words).
                                   NB int rather than 2*short since must fit
                                   into fortran integer */
  (int) 2*sizeof (float),        /* 4: complex (pairs of words) */
  (int) sizeof (int),           /* 5: not used */
  (int) sizeof (int)            /* 6: integers */
};

static int (*_read_mode[])(CCP4File *, uint8 *, size_t) = {
  ccp4_file_readchar,
  ccp4_file_readshort,
  ccp4_file_readfloat,
  ccp4_file_readshortcomp,
  ccp4_file_readcomp,
  NULL,
  ccp4_file_readint
};

static int (*_write_mode[])(CCP4File *, const uint8 *, size_t) = {
  ccp4_file_writechar,
  ccp4_file_writeshort,
  ccp4_file_writefloat,
  ccp4_file_writeshortcomp,
  ccp4_file_writecomp,
  NULL,
  ccp4_file_writeint    
};

/**
 * vaxF2ieeeF:
 * @param buffer (float_uint_uchar *) vax order float array
 * @param size (unsigned int) number of items
 *
 * Translation from Vax floating point format to ieee.
 *
 */
static void vaxF2ieeeF(union float_uint_uchar *buffer, const unsigned int size)
{
  union float_uint_uchar out;
  unsigned char exp;
  unsigned int i;

  if ( buffer == NULL || size == 0) return;
  
  for (i = 0; i < size; i++) {
    exp = (buffer[i].c[1] << 1) | (buffer[i].c[0] >> 7); /* extract exponent */
    if (!exp && !buffer[i].c[1])        /* zero value */
      out.c[0] = out.c[1] = out.c[2] = out.c[3] = 0;
    else if (exp > 2) {         /* normal value */
      out.c[0] = buffer[i].c[1] - (uint8)1; /* subtracts 2 from exponent */
      /* copy mantissa, LSB of exponent */
      out.c[1] = buffer[i].c[0];
      out.c[2] = buffer[i].c[3];
      out.c[3] = buffer[i].c[2];
    } else if (exp) {           /* denormalized number */
      int shft;

      out.c[0] = buffer[i].c[1] & 0x80; /* keep sign, zero exponent */
      shft = 3 - exp;
      /* shift original mant by 1 or 2 to get denormalized mant */
      /* prefix mantissa with '1'b or '01'b as appropriate */
      out.c[1] = (uint8)((buffer[i].c[0] & 0x7f) >> shft) |
        (uint8)(0x10 << exp);
      out.c[2] = (uint8)(buffer[i].c[0] << (8-shft)) |
        (uint8)(buffer[i].c[3] >> shft);
      out.c[3] = (uint8)(buffer[i].c[3] << (8-shft)) |
        (uint8)(buffer[i].c[2] >> shft);
    } else {                    /* sign=1 -> infinity or NaN */
      out.c[0] = 0xff;          /* set exp to 255 */
      /* copy mantissa */
      out.c[1] = buffer[i].c[0] | (uint8)0x80; /* LSB of exp = 1 */
      out.c[2] = buffer[i].c[3];
      out.c[3] = buffer[i].c[2];
    }
    buffer[i] = out;            /* copy back result */
  }
}

/**
 * ieeeF2vaxF:
 * @param buffer (float_uint_uchar *) big endian order float array
 * @param size (unsigned int) number of items
 *
 * Translation from ieee floating point format to vax.
 *
 */
static void ieeeF2vaxF(union float_uint_uchar *buffer, const unsigned int size)
{
  union float_uint_uchar out;
  unsigned char exp;
  unsigned int i;

  if ( buffer == NULL || size == 0) return;

  for (i=0; i<size; i++) {
    exp = (buffer[i].c[0]<<1) | (buffer[i].c[1]>>7); /* extract exponent */
    if (exp) {                  /* non-zero exponent */
      /* copy mantissa, last bit of exponent */
      out.c[0] = buffer[i].c[1];
      out.c[2] = buffer[i].c[3];
      out.c[3] = buffer[i].c[2];
      if (exp < 254)            /* normal value */
        out.c[1] = buffer[i].c[0] + (uint8)1; /* actually adds two to exp */
      else {                    /* infinity or NaN */
        if (exp == 254)         /* unrepresentable - OFL */
          /* set mant=0 for overflow */
          out.c[0] = out.c[1] = out.c[2] = out.c[3] = 0; 
        out.c[0] &= 0x7f;       /* set last bit of exp to 0 */
        out.c[1] = 0x80;        /* sign=1 exp=0 -> OFL or NaN.  this will raise
                                   a reserved operand exception if used. */
      }
    } else if (buffer[i].c[1] & 0x60) { /* denormalized value */
      int shft;
      
      shft = (buffer[i].c[1] & 0x40) ? 1 : 2; /* shift needed to normalize */
      /* shift mantissa */
      /* note last bit of exp set to 1 implicitly */
      out.c[0] = (uint8)(buffer[i].c[1] << shft) |
        (uint8)(buffer[i].c[2] >> (8-shft));
      out.c[3] = (uint8)(buffer[i].c[2] << shft) |
        (uint8)(buffer[i].c[3] >> (8-shft));
      out.c[2] = (uint8)(buffer[i].c[3] << shft);
      out.c[1] = (uint8)(buffer[i].c[0] & 0x80); /* sign */
      if (shft==1) {            /* set exp to 2 */
        out.c[1] |= 0x01;
        out.c[0] &= 0x7f;       /* set LSB of exp to 0 */
      }
    } else                      /* zero */
      out.c[0] = out.c[1] = out.c[2] = out.c[3] = 0;
    buffer[i] = out;            /* copy back the result */
  }
}

/**
 * convexF2ieeeF:
 * @param buffer (float_uint_uchar *) float array with convex byte ordering
 * @param size (unsigned int) number of items
 *
 * Translation from convex floating point format to ieee.
 *
 */
static void convexF2ieeeF(union float_uint_uchar *buffer, const unsigned int size)
{
  union float_uint_uchar out;
  unsigned char exp;
  unsigned int i;

  if ( buffer == NULL || size == 0) return;
  
  for (i = 0; i < size; i++) {
    exp = (buffer[i].c[0]<<1) | (buffer[i].c[1]>>7); /* extract exponent */
    if (!exp && !buffer[i].c[0])        /* zero value */
      out.c[0] = out.c[1] = out.c[2] = out.c[3] = 0;
    else if (exp > 2) {         /* normal value */
      out.c[0] = buffer[i].c[0] - (uint8)1; /* subtracts 2 from exponent */
      /* copy mantissa, LSB of exponent */
      out.c[1] = buffer[i].c[1];
      out.c[2] = buffer[i].c[2];
      out.c[3] = buffer[i].c[3];
    } else if (exp) {           /* denormalized number */
      int shft;
      
      out.c[0] = buffer[i].c[0] & 0x80; /* keep sign, zero exponent */
      shft = 3 - exp;
      /* shift original mant by 1 or 2 to get denormalized mant */
      /* prefix mantissa with '1'b or '01'b as appropriate */
      out.c[1] = (uint8)((buffer[i].c[1] & 0x7f) >> shft) |
        (uint8)(0x10 << exp);
      out.c[2] = (uint8)(buffer[i].c[1] << (8-shft)) |
        (uint8)(buffer[i].c[2] >> shft);
      out.c[3] = (uint8)(buffer[i].c[2] << (8-shft)) |
        (uint8)(buffer[i].c[3] >> shft);
    } else {                    /* sign=1 -> infinity or NaN */
      out.c[0] = 0xff;          /* set exp to 255 */
      /* copy mantissa */
      out.c[1] = buffer[i].c[1] | (uint8)0x80; /* LSB of exp = 1 */
      out.c[2] = buffer[i].c[2];
      out.c[3] = buffer[i].c[3];
    }
    buffer[i] = out;            /* copy back result */
  }
}

/**
 * ieeeF2convexF:
 * @param buffer (float_uint_uchar *) float array with big endian byte ordering
 * @param size (const unsigned int) numnber of items
 *
 * Translation from ieee floating point format to convex.
 *
 */
static void ieeeF2convexF(union float_uint_uchar *buffer, const unsigned int size)
{
  union float_uint_uchar out;
  unsigned char exp;
  unsigned int i;

  if ( buffer == NULL || size == 0) return;

  for (i=0; i < size; i++) {
    exp = (uint8)(buffer[i].c[0] << 1) |
      (uint8)(buffer[i].c[1] >> 7); /* extract exponent */
    if (exp) {                  /* non-zero exponent */
      /* copy mantissa, last bit of exponent */
      out.c[1] = buffer[i].c[1];
      out.c[3] = buffer[i].c[3];
      out.c[2] = buffer[i].c[2];
      if (exp < 254)            /* normal value */
        out.c[0] = buffer[i].c[0] + (uint8)1; /* actually adds two to exp */
      else {                    /* infinity or NaN */
        if (exp == 254)         /* unrepresentable - OFL */
          /* set mant=0 for overflow */
          out.c[0] = out.c[1] = out.c[2] = out.c[3] = 0; 
        out.c[1] &= 0x7f;       /* set last bit of exp to 0 */
        out.c[0] = 0x80;        /* sign=1 exp=0 -> OFL or NaN.  this will raise
                                   a reserved operand exception if used. */
      }
    } else if (buffer[i].c[1] & 0x60) { /* denormalized value */
      int shft;
      
      shft = (buffer[i].c[1] & 0x40) ? 1 : 2; /* shift needed to normalize */
      /* shift mantissa */
      /* note last bit of exp set to 1 implicitly */
      out.c[1] = (uint8)(buffer[i].c[1] << shft) |
        (uint8)(buffer[i].c[2] >> (8-shft));
      out.c[2] = (uint8)(buffer[i].c[2] << shft) |
        (uint8)(buffer[i].c[3] >> (8-shft));
      out.c[3] = (uint8)(buffer[i].c[3] << shft);
      out.c[0] = (uint8)(buffer[i].c[0] & 0x80); /* sign */
      if (shft==1) {            /* set exp to 2 */
        out.c[0] |= 0x01;
        out.c[1] &= 0x7f;       /* set LSB of exp to 0 */
      }
    } else                      /* zero */
      out.c[0] = out.c[1] = out.c[2] = out.c[3] = 0;
    buffer[i] = out;            /* copy back the result */
  }
}

/**
 * ccp4_file_raw_read:
 * @param cfile *  (CCP4File *)
 * @param buffer * (char *) input array
 * @param n_items (size_t) number of items
 *
 * reads block of n_items bytes from cfile to buffer via
 * FILE struct cfile->stream(fread) or file desc cfile->fd
 * read/_read).  Increments location value cfile->loc. The
 * cfile->iostat flag is set on failure.
 * @return number of bytes read.
 */
int ccp4_file_raw_read(CCP4File *cfile, char *buffer, size_t n_items)
{
  int result;
  
  if (cfile->buffered && cfile->stream) {
    result = fread (buffer, (size_t) sizeof(char), n_items,
                    cfile->stream);
    if (result != n_items && feof(cfile->stream)) {
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_EOF), 
		  "ccp4_file_raw_read", NULL); 
      cfile->iostat = CIO_EOF;
      result = EOF;
    } else if (result != n_items && ferror(cfile->stream)) {
      ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_raw_read", NULL);
      cfile->iostat = CIO_ReadFail; 
      result = 0; }
  } else {
#if defined _MSC_VER
    result = _read (cfile->fd, buffer, n_items);
#else
    result = read (cfile->fd, buffer, n_items);
#endif
    if (n_items && result <= 0) {
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_ReadFail),
		  "ccp4_file_raw_read", NULL);
      cfile->iostat = CIO_ReadFail;
      result = 0; } 
  }
  cfile->last_op = READ_OP;

  cfile->loc += result; 
  cfile->getbuff = 0;

  return result;
}

/**
 * ccp4_file_raw_write:
 * @param cfile  (CCP4File *)
 * @param buffer (char *) output array
 * @param n_items (size_t) number of items
 *
 * writes block of @n_items bytes from @buffer to @cfile via FILE struct
 * @cfile->stream(fwrite) or file desc @cfile->fd(write/_write).  Increments
 * @cfile->loc on success, or resets on failure, which is then used to
 * determine the file length.  On failure @cfile->iostat is set.
 * @return number of bytes written.
 */
int ccp4_file_raw_write(CCP4File *cfile, const char *buffer, size_t n_items)
{
  int result;
  
  if (cfile->buffered && cfile->stream)
    result = fwrite (buffer, (size_t) sizeof(char), n_items,
                    cfile->stream);
  else 
#if defined _MSC_VER
    result = _write (cfile->fd, buffer, n_items);
#else
    result = write (cfile->fd, buffer, n_items);
#endif
        
  cfile->last_op = WRITE_OP;
    
  if (result == n_items) 
    cfile->loc += n_items; 
  else {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_WriteFail),
		"ccp4_file_raw_write", NULL);
    cfile->iostat = CIO_WriteFail; 
    ccp4_file_tell(cfile); }
  cfile->length = MAX(cfile->loc,cfile->length);
  cfile->getbuff = 0;

  return result;
}

/**
 * ccp4_file_raw_seek:
 * @param cfile  (CCP4File *)
 * @param offset (long) offset in bytes
 * @param whence (int) SEEK_SET, SEEK_CUR, or SEEK_END
 *
 * if the file is "seekable" (not stdin) the function
 * seeks on @cfile by offset bytes using fseek/ftell (@cfile->stream)
 * or lseek (@cfile->fd).  %SEEK_SET is relative
 * to start of file, %SEEK_CUR to current, %SEEK_END to
 * end.
 * @return offset in bytes on success, -1 on failure.
 */
int ccp4_file_raw_seek(CCP4File *cfile, long offset, int whence)
{
  int result = -1;
      
  if (!cfile->direct)  {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode),
		"ccp4_file_raw_seek", NULL);
    return result; }
  
  if (cfile->buffered) {
#if defined (__alpha) && defined (vms)
    (void) fflush (cfile->stream);
#endif
    if (!(result = (fseek (cfile->stream,offset,whence))))
      result = ftell(cfile->stream);   
  } else {
#if defined _MSC_VER
     result = _lseek(cfile->fd,offset,whence);
#else
     result = lseek(cfile->fd, offset, whence);
#endif
  } 
  
  cfile->last_op = IRRELEVANT_OP;
  
  if (result ==  -1) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_SeekFail),
		"ccp4_file_raw_seek", NULL);
    cfile->iostat = CIO_SeekFail;
  } else 
    cfile->loc = result;
  cfile->getbuff = 0;
  
  return (result);
}

/**
 * _file_free:
 * @param cfile  (CCP4File *)
 *
 * free up @cfile
 */
static void _file_free(CCP4File *cfile)
{
  if(cfile->name) free(cfile->name);
  cfile->name = NULL;
  free(cfile);
}

/**
 * _file_init:
 * @param cfile ()
 * 
 * initialise cfile struct
 * @return cfile struct
 */
static CCP4File *_file_init()
{
  CCP4File *cfile = (CCP4File *) malloc(sizeof(CCP4File));
  char *foreign = getenv ("CONVERT_FROM");
  char *native = getenv ("NATIVEMTZ");
  
  memset(cfile,'\0',sizeof(CCP4File));
  cfile->fd = -1;
  cfile->buffered = 1;
  cfile->binary = 1;
  cfile->last_op = IRRELEVANT_OP;
  cfile->mode = DEFMODE;
  cfile->itemsize = _item_sizes[DEFMODE]; 
  if (native == NULL && foreign != NULL) {
    if (strcmp (foreign, "BEIEEE") == 0) {
      cfile->fconvert = DFNTF_BEIEEE ;
      cfile->iconvert = DFNTI_MBO ; }
    else if (strcmp (foreign, "LEIEEE") == 0) {
      cfile->fconvert = DFNTF_LEIEEE;
      cfile->iconvert = DFNTI_IBO ; }
    else if (strcmp (foreign, "VAX") == 0) {
      cfile->fconvert = DFNTF_VAX ;
      cfile->iconvert = DFNTI_IBO ; }
    else if (strcmp (foreign, "CONVEXNATIVE") == 0) {
      cfile->fconvert = DFNTF_CONVEXNATIVE ;
      cfile->iconvert = DFNTI_MBO ; }  
  } else {
    cfile->fconvert = nativeFT;
    cfile->iconvert = nativeIT;
  }
  cfile->_read=_read_mode[DEFMODE];
  cfile->_write=_write_mode[DEFMODE]; 
  return (cfile);
}

/** 
 * _file_open_mode:
 * @param cfile (CCP4File *)
 * @param flag (const int) mode flag
 *
 * set file open mode elements of @cfile (see ccp4_sysdep.h)
 *  O_TMP    = 0x0010
 *  O_RDONLY = 0x0000
 *  O_WRONLY = 0x0001
 *  O_RDWR   = 0x0002
 *  O_APPEND = 0x0008
 */
static void _file_open_mode(CCP4File * cfile, const int flag)
{
  if (flag & O_TMP)
    cfile->scratch = 1;
  if (flag & (O_WRONLY | O_RDWR | O_APPEND) ) {
    cfile->write = 1; 
    if (flag & O_RDWR)
      cfile->read = 1;
    if (flag & O_APPEND)
      cfile->append = 1; 
  } else 
    cfile->read = 1; 
}

/**
 * _file_close:
 * @param cfile  (CCP4File *)
 *
 * close @cfile if it is "owned" (@cfile->own) using fclose or close.
 * Reset @cfile to some safe value.  Note: flush anyway.
 * @return 0 on success, -1 on failure.
 */
static int _file_close (CCP4File *cfile)
{
  int result = 0;
  
  if(cfile->buffered && cfile->stream) {
    if (cfile->own)
      result = fclose (cfile->stream); 
    else
      result = fflush(cfile->stream);
  } else {
    if (cfile->own)
#if defined _MSC_VER
      result = _close(cfile->fd);
#else
      result = close (cfile->fd);
#endif
  }
  
  if (result == EOF) 
    cfile->iostat = CIO_CloseFail;
  else 
    cfile->stream = NULL; 

  return (result);
}

/**
 * ccp4_file_is_write:
 * @param cfile  (CCP4File *)
 *
 * is the @cfile writeable
 * @return 1 if true 
 */
int ccp4_file_is_write(const CCP4File *cfile)
{
  return (cfile->write);
}

/**
 * ccp4_file_is_read:
 * @param cfile  (CCP4File *)
 *
 * is the @cfile readable
 * @return 1 if true.
 */
int ccp4_file_is_read(const CCP4File *cfile)
{
  return (cfile->read);
}

/**
 * ccp4_file_is_append:
 * @param cfile  (CCP4File *)
 *
 * is the @cfile in append mode
 * @return 1 if true.
 */
int ccp4_file_is_append(const CCP4File *cfile)
{
  return (cfile->append);
}

/**
 * ccp4_file_is_scratch:
 * @param cfile  (CCP4File *)
 *
 * is scratch file
 * @return 1 if true.
 */
int ccp4_file_is_scratch(const CCP4File *cfile)
{
  return (cfile->scratch);
}

/**
 * ccp4_file_is_buffered:
 * @param cfile  (CCP4File *)
 *
 * is the file buffered
 * @return 1 if true
 */
int ccp4_file_is_buffered(const CCP4File *cfile)
{
  return (cfile->buffered);
}

/**
 * ccp4_file_status:
 * @param cfile  (CCP4File *)
 *
 * @return @cfile error status
 */
int ccp4_file_status(const CCP4File *cfile)
{
  return (cfile->iostat);
}

int ccp4_file_raw_setstamp(CCP4File *cfile, const size_t offset)
{
  cfile->stamp_loc = offset;
  return 0;
}

/**
 * ccp4_file_setstamp:
 * @param cfile  (CCP4File *)
 * @param stamp_loc (size_t) offset in items
 *
 * set the machine stamp offset in CCP4 items determined
 * by the mode of @cfile.  See ccp4_file_setmode().
 * @return 0 on success, %EOF on failure
 */
int ccp4_file_setstamp(CCP4File *cfile, const size_t offset)
{
  if (!cfile) {
    ccp4_signal(CCP4_ERRLEVEL(3)| CCP4_ERRNO(CIO_NullPtr),
		"ccp4_file_setstamp", NULL);
    return EOF; }

  return ccp4_file_raw_setstamp(cfile, offset*cfile->itemsize);
}

/**
 * ccp4_file_setmode:
 * @param cfile  (CCP4File *)
 * @param mode (int) io_mode
 *
 * set the data mode of cfile to mode
 * (CCP4_BYTE (8 bit) = 0, 
 *  CCP4_INT16 (16 bit) = 1, 
 *  CCP4_INT32 (32 bit) = 6, 
 *  FLOAT32 (32 bit) = 2, 
 *  COMP32 (2*16 bit) = 3, 
 *  COMP64 (2*32 bit) = 4).  
 * @return 0 on success, EOF on failure.
 */
int ccp4_file_setmode (CCP4File *cfile, const int mode)
{
  if (!cfile) {
    ccp4_signal(CCP4_ERRLEVEL(3)| CCP4_ERRNO(CIO_NullPtr),
		"ccp4_file_mode", NULL);
    return EOF; }

  if (mode >= 0 && mode <= 6 && mode != 5) {
    cfile->mode = mode;
    cfile->itemsize = _item_sizes[mode];
    cfile->_read=_read_mode[mode];
    cfile->_write=_write_mode[mode];
  } else {
    ccp4_signal(CCP4_ERRLEVEL(3)| CCP4_ERRNO(CIO_BadMode),
		"ccp4_file_mode", NULL);
    return EOF; }

  return  0;
}

/**
 * ccp4_file_mode:
 * @param cfile  (CCP4File *)
 *
 * get data mode of @cfile (CCP4_BYTE =0, CCP4_INT16 =1, CCP4_INT32=6,
 * FLOAT32 =2, COMP32 =3, COMP64 =4)
 * @return %mode
 */
int ccp4_file_mode (const CCP4File *cfile)
{
  if (!cfile) {
    ccp4_signal(CCP4_ERRLEVEL(3)| CCP4_ERRNO(CIO_NullPtr),
		"ccp4_file_mode", NULL);
    return EOF; }
    
  return (cfile->mode);
}

/**
 * ccp4_file_itemsize:
 * @param cfile  (CCP4File *)
 *
 * @return %itemsize of @cfile.
 */
int ccp4_file_itemsize(const CCP4File *cfile)
{
  if (!cfile) {
    ccp4_signal(CCP4_ERRLEVEL(3)| CCP4_ERRNO(CIO_BadMode),
		"ccp4_file_itemsize", NULL);
    return EOF; }

  return (cfile->itemsize);
}

/**
 * ccp4_file_name:
 * @param cfile (CCP4File *)
 *
 * strdup @cfile->name
 * @return name of file as char *
 */
char *ccp4_file_name( CCP4File *cfile)
{
#if defined _MSC_VER
  return ( cfile == NULL ? NULL : _strdup(cfile->name));
#else
  return ( cfile == NULL ? NULL : strdup(cfile->name));
#endif
}

/**
 * ccp4_file_setbyte:
 * @param cfile (CCP4File *)
 * @param byte_order (int)
 *
 * set byte ordering for file
 * Return:
  */
int ccp4_file_setbyte(CCP4File *cfile, const int byte_order)
{
  int result = (cfile->fconvert | (cfile->iconvert<<8));
  
  switch (byte_order) {
  case DFNTF_BEIEEE:
    cfile->fconvert = DFNTF_BEIEEE;
    cfile->iconvert = DFNTI_MBO;
    break;
  case DFNTF_LEIEEE:
    cfile->fconvert = DFNTF_LEIEEE;
    cfile->iconvert = DFNTI_IBO ; 
    break;
  case DFNTF_VAX:
    cfile->fconvert = DFNTF_VAX ;
    cfile->iconvert = DFNTI_IBO ; 
    break;
  case DFNTF_CONVEXNATIVE:
    cfile->fconvert = DFNTF_CONVEXNATIVE ;
    cfile->iconvert = DFNTI_MBO ;
    break;
  default:
    ccp4_signal(CCP4_ERRLEVEL(3)| CCP4_ERRNO(CIO_BadMode),
		"ccp4_file_setbyte", NULL);
    result = 0;
  }
  return result;
}

/**
 * ccp4_file_byte:
 * @param cfile (CCP4File *)
 *
 * get byte ordering for file
 * @return byte ordering information
  */
int ccp4_file_byte(CCP4File *cfile)
{
  return (cfile->fconvert | (cfile->iconvert<<8));
}

/**
 * ccp4_file_open_file:
 * @param file (const FILE *) FILE struct
 * @param flag (const int) io mode (O_RDONLY =0, O_WRONLY =1, O_RDWR =2,
 *        O_TMP =, O_APPEND =)
 *
 * open @cfile with existing handle FILE struct file and mode @flag.
 * The struct stat is check to determine if file is a regular file,
 * if it is, and is not stdin, it is assumed to be direct access. 
 * @return (CCP4File *) on success, NULL on failure
 */
CCP4File *ccp4_file_open_file (const FILE *file, const int flag)
{
  CCP4File *cfile;
#if defined _MSC_VER
  struct _stat st;
#else
  struct stat st;
#endif

  if (!file) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr),
                "ccp4_file_open_file", NULL);
    return NULL; }

  if (!(cfile = _file_init() ) ) {
      ccp4_signal(CCP4_ERRLEVEL(3),
                  "ccp4_file_open_file", NULL);    
      return NULL; }
          
  /* set values in structure */
  _file_open_mode(cfile,flag);
  cfile->stream = (FILE *) file;
  cfile->buffered = 1;
  cfile->open = 1;   

#if defined _MSC_VER
  _fstat(_fileno(cfile->stream), &st);
  if ( !(st.st_mode & S_IFREG) || file == stdin) {
#else
  fstat(fileno(cfile->stream), &st);
  if ( !S_ISREG(st.st_mode) || file == stdin ) {
#endif
    cfile->length = INT_MAX;
    cfile->direct = 0;
  } else {
    cfile->length = st.st_size;
    cfile->direct = 1;
  }
  cfile->loc = ftell( (FILE *)file);
  
  return cfile;
}

/**
 * ccp4_file_open_fd:
 * @param fd (const int) file descriptor
 * @param flag (const int) io mode (O_RDONLY =0, O_WRONLY =1, O_RDWR =2,
 *        O_TMP =, O_APPEND =)
 *
 * initialise CCP4File struct with file descriptor @fd and mode @flag 
 * The struct stat is check to determine if file is a regular file,
 * if it is, and is not stdin, it is assumed to be direct access. 
 * @return (CCP4File *) on success, NULL on failure
 */
CCP4File *ccp4_file_open_fd (const int fd, const int flag)
{
  CCP4File * cfile;
#if defined _MSC_VER
  struct _stat st;
#else
  struct stat st;
#endif

  if (!(cfile = _file_init() ) ) {
      ccp4_signal(CCP4_ERRLEVEL(3),
                  "ccp4_file_open_fd", NULL);    
      return NULL; }
      
  _file_open_mode(cfile, flag);      
  cfile->fd = fd;
  cfile->open = 1;    
  cfile->buffered = 0;

#if defined _MSC_VER
  _fstat(fd, &st);
  if ( !(st.st_mode & S_IFREG) || fd == 0) {
#else
  fstat(fd, &st);
  if ( !S_ISREG(st.st_mode) || fd == 0 ) {
#endif
    cfile->length = INT_MAX;
    cfile->direct = 0;
    cfile->loc = 0;
  } else {
    cfile->length = st.st_size;
    cfile->direct = 1;
#if defined _MSC_VER
    cfile->loc = _lseek(fd, 0L, SEEK_CUR);
#else
    cfile->loc = lseek(fd, 0L, SEEK_CUR);
#endif
  }

  return cfile;
}

/**
 * ccp4_file_open:
 * @param filename (const char *) filename
 * @param flag (const int) i/o mode, possible values are O_RDONLY, O_WRONLY, 
 *      O_RDWR, O_APPEND, O_TMP, O_CREAT, O_TRUNC - see ccp4_sysdep.h
 *
 * initialise CCP4File struct for file filename with mode @flag.
 * If !buffered use open(), otherwise fopen()
  * The struct stat is check to determine if file is a regular file,
 * if it is, and is not stdin, it is assumed to be direct access. 
 *
 * @return (CCP4File *) on success, NULL on failure
 */
CCP4File *ccp4_file_open (const char *filename, const int flag)
{
  CCP4File *cfile; 
  int openflags = O_RDONLY;
  char fmode[5];
#if defined _MSC_VER
  struct _stat st;
#else
  struct stat st;
#endif

  if (!(cfile = _file_init())) {
    ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_open", NULL);
    return NULL; }
  _file_open_mode(cfile, flag);

  if (!cfile->buffered) {
    if (cfile->read && cfile->write) openflags = (O_RDWR | O_CREAT);
    else if (cfile->write)  openflags = (O_WRONLY | O_CREAT); 
    if (cfile->append) openflags |= O_APPEND;
    if (flag & O_TRUNC) openflags |= O_TRUNC;
#if defined _MSC_VER
    if (cfile->scratch) openflags |= O_TEMPORARY;
#endif
#if defined(__DECC) && defined(VMS) || defined (_MSC_VER)
    openflags |= O_BINARY;
#endif
#if defined _MSC_VER
    cfile->fd = _open(filename, openflags);
#else
    cfile->fd = open(filename, openflags);
#endif
    if (cfile->fd == -1) {
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_CantOpenFile),
                  "ccp4_file_open1", NULL);
      return NULL; 
    } else {
#if defined _MSC_VER
      _fstat(cfile->fd, &st); }
#else
      fstat(cfile->fd, &st); }
#endif
  } else {
    char *mptr = fmode;
    if (cfile->append) {
      *mptr++ = 'a';
      if (cfile->read) *mptr++ = '+';
    } else {
      if (cfile->read && cfile->write) {
        if (flag & O_TRUNC) {*mptr++ = 'w'; }
        else *mptr++ = 'r';
        *mptr++ = '+';
      } else if (cfile->write) 
        *mptr++ = 'w';
      else 
        *mptr++ = 'r';
    }
#if defined(__DECC) && defined(VMS) || defined (_WIN32)
    *mptr++ = 'b';
#endif
    *mptr++ = '\0';

#ifdef VMS
    if (cfile->scratch)
      cfile->stream = fopen (filename, fmode,
                             "mbc=16",        /* bigger blocksize */
                             "fop=tmd");      /* temporary, delete on close */
    else
      cfile->stream = fopen (filename, fmode,
                             "mbc=16",        /* bigger blocksize */
                             "ctx=stm", "mrs=0", "rat=cr", "rfm=stmlf");
#elif defined(_WIN32)
    if (cfile->scratch) {
      cfile->stream = tmpfile();
      if (!cfile->stream) {
        ccp4_signal(CCP4_ERRLEVEL(2) | CCP4_ERRNO(CIO_CantOpenFile),
                    "tmpfile() failed, opening normal file instead.", NULL);
        cfile->stream = fopen (filename, fmode);
      }
    }
    else
      cfile->stream = fopen (filename, fmode);
#else
    cfile->stream = fopen (filename, fmode);
    if (cfile->stream) 
      if (cfile->scratch && unlink (filename)!=0) {
        ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_UnlinkFail),
                    "ccp4_file_open(unlink)", NULL);
        cfile->iostat = CIO_UnlinkFail; return NULL; }
#endif
    if (!cfile->stream) {
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_CantOpenFile),
                  "ccp4_file_open2", NULL);
      cfile->iostat = CIO_CantOpenFile;
      free(cfile);
      return NULL; } 
#if defined (__alpha) && defined (vms)
(void) fflush (cfile->stream);
#endif
#if defined _MSC_VER
  _fstat(_fileno(cfile->stream), &st);
#else
  fstat(fileno(cfile->stream), &st);
#endif
  }
#if defined _MSC_VER
  cfile->name = _strdup(filename);
#else
  cfile->name = strdup(filename);
#endif
  cfile->open = 1;  
  cfile->own = 1;  
#if defined _MSC_VER
  if ( !(st.st_mode & S_IFREG) ) {
#else
  if ( !cfile->scratch && !S_ISREG(st.st_mode) ) {
#endif
    cfile->length = INT_MAX;
    cfile->direct = 0;
  } else {
    cfile->length = st.st_size;
    cfile->direct = 1;
  }
  cfile->loc = cfile->append ? cfile->length : 0;
  
  return cfile;
}

/**
 * ccp4_file_close:
 * @param cfile (CCP4File *) 
 *
 * close @cfile if owned, close (non-buffered) or 
 * fclose (buffered),  or fflush if stream not owned.
 * Free resources.
 * @return 0 on success, EOF on failure
 */
int ccp4_file_close (CCP4File *cfile)
{

  if (!cfile) {
    ccp4_signal(CCP4_ERRLEVEL(3)| CCP4_ERRNO(CIO_NullPtr),
		"ccp4_file_close", NULL);
    return EOF; }

  if (_file_close (cfile) == EOF) {
    ccp4_signal(CCP4_ERRLEVEL(3),"ccp4_file_close", NULL);
    return EOF; }
 
  _file_free(cfile);

  return (0);
}

/**
 * ccp4_file_rarch:
 * @param cfile (CCP4File *) 
 *
 * read machine stamp from file @cfile->stream.
 * The machine stamp is at @cfile->stamp_loc items, set
 * by ccp4_file_setstamp() (default 0).
 * NB. these values may be overrriden with the environmental
 * variable CONVERT_FROM.
 * @return fileFT | (fileIT<<8)
 */
int ccp4_file_rarch (CCP4File *cfile)
{
  unsigned char mtstring[4];    /* machine stamp */
  char *native = getenv ("NATIVEMTZ");
  char *foreign = getenv ("CONVERT_FROM");

  if (!cfile) {
    ccp4_signal(CCP4_ERRLEVEL(3)| CCP4_ERRNO(CIO_NullPtr),
		"ccp4_file_rarch", NULL);
    return EOF; }

  if (native != NULL) return (nativeFT | (nativeIT<<8)); 
  if (foreign == NULL) {
    if (ccp4_file_raw_seek(cfile, cfile->stamp_loc, SEEK_SET) == -1) { 
      ccp4_signal(CCP4_ERRLEVEL(3),"ccp4_file_rarch", NULL);
      return EOF; }
    
    if (ccp4_file_raw_read(cfile, (char *) mtstring, 4UL) != 4) { 
        ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_rarch", NULL);
        return EOF; }
        
    cfile->iconvert = (mtstring[1]>>4) & 0x0f;
    cfile->fconvert = (mtstring[0]>>4) & 0x0f;

    /* iconvert and fconvert should be one of the DFNTI/DFNTF values listed
       in ccp4_sysdep.h and hence non-zero. Some machine stamps can be corrupted
       (e.g. mrc files from chimera). We try to trap for this, and revert to
       native. */
    if (cfile->iconvert == 0 || cfile->fconvert == 0) {
       if (ccp4_liberr_verbosity(-1))
          printf("Warning: Machine stamp corrupted? Assuming native format. \n");
       cfile->iconvert = nativeIT;
       cfile->fconvert = nativeFT;       
    } 
  }
  
  return (cfile->fconvert | (cfile->iconvert<<8));
}

/**
 * ccp4_file_warch:
 * @param cfile (CCP4File *) 
 *
 * write machine stamp to file @cfile->stream.
 * The machine stamp is placed at @cfile->stamp_loc items,
 * set by ccp4_file_setstamp() (defaults to 0).
 *
 * @return 0 on success, EOF on failure
 */
int ccp4_file_warch (CCP4File *cfile)
{
  unsigned char mtstring[4];    /* machine stamp */

  if (!cfile) {
    ccp4_signal(CCP4_ERRLEVEL(3)| CCP4_ERRNO(CIO_NullPtr),
		"ccp4_file_warch", NULL);
    return EOF; }

  if (ccp4_file_raw_seek(cfile, cfile->stamp_loc, SEEK_SET) == -1) { 
      ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_warch", NULL);
      return EOF; }

  mtstring[0] = cfile->fconvert | (cfile->fconvert << 4);
  mtstring[1] = 1 | (cfile->iconvert << 4);
  mtstring[2] = mtstring[3] = 0;

  if (ccp4_file_raw_write(cfile, (const char *) mtstring, 4) != 4) { 
      ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_warch", NULL);
    return EOF; }

  return 0;
}

/**
 * ccp4_file_read:
 * @param cfile (CCP4File *)
 * @param buffer (uint8 *) buffer
 * @param nitems (size_t) number of items 
 *
 * mode dependent read function.  Reads @nitems items from stream
 * @cfile->stream to @buffer as determined by cfile->mode.
 *
 * @return number of items read on success, EOF on failure
 */
int ccp4_file_read (CCP4File *cfile, uint8 *buffer, size_t nitems)
{
  int result;
 
  result = cfile->_read(cfile,(uint8 *) buffer,nitems);
  
  if (result != nitems) 
    ccp4_signal(CCP4_ERRLEVEL(3), 
	       "ccp4_file_read", NULL);
  return (result);
}

/**
 * ccp4_file_readcomp:
 * @param cfile (CCP4File *)
 * @param buffer (uint8 *) buffer
 * @param nitems (size_t) number of items 
 *
 * float complex {float,float} read function.  Reads @nitems complex from stream
 * @cfile->stream to @buffer.  Allows short count when eof is detected (
 * buffered input only).
 *
 * @return number of complex read on success, EOF on failure
 */
int ccp4_file_readcomp (CCP4File *cfile, uint8 *buffer, size_t nitems)
{
  int i, n, result;

  if (!cfile)  {
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr), 
		  "ccp4_file_readcomp", NULL);
    return EOF; }

  if ( !cfile->read || cfile->iostat) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		"ccp4_file_readcomp", NULL);
    return EOF; }

  if (cfile->last_op == WRITE_OP)
    if (ccp4_file_raw_seek(cfile,0L,SEEK_CUR) == -1) {
      ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_readcomp", NULL);
      return EOF; }   

  n = _item_sizes[COMP64] * nitems;
  if ( (result = ccp4_file_raw_read (cfile, (char *) buffer, n)) != n) {
    ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_readcomp", NULL);
    if (cfile->stream && !feof(cfile->stream)) 
      return EOF;  }  /* short count on stream is OK if EOF */

  result /= _item_sizes[COMP64];
  n = result;
  if (cfile->fconvert != nativeFT) {
    n *= 2;                        /* pairs of reals */
    switch (cfile->fconvert) {     /* get to BE IEEE */
    case DFNTF_VAX :
        vaxF2ieeeF((union float_uint_uchar *) buffer, n);
        break;   
    case DFNTF_CONVEXNATIVE :
        convexF2ieeeF((union float_uint_uchar *) buffer, n);
        break;
    case DFNTF_BEIEEE :
        break;
    case DFNTF_LEIEEE :
        {
            char j;
            for (i=0; i < n*4; i+=4) {
                j = buffer[i];
                buffer[i] = buffer[i+3];
                buffer[i+3] = j;
                j = buffer[i+1];
                buffer[i+1] = buffer[i+2];
                buffer[i+2] =j; }
        }
        break;
    default :
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		  "ccp4_file_readcomp", NULL);
      return EOF; 
    }
    switch (nativeFT) {              /* get to Native if not BE IEEE */
    case DFNTF_BEIEEE :
        break;                      
    case DFNTF_LEIEEE :
        {
            char j;
            for (i=0; i < n*4; i+=4) {
                j = buffer[i];
                buffer[i] = buffer[i+3];
                buffer[i+3] = j;
                j = buffer[i+1];
                buffer[i+1] = buffer[i+2];
                buffer[i+2] =j; }
        }
        break;
    case DFNTF_CONVEXNATIVE :
        ieeeF2convexF((union float_uint_uchar *) buffer, n);
        break;
    case DFNTF_VAX :
        ieeeF2vaxF((union float_uint_uchar *) buffer, n);
        break;
    default :
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		  "ccp4_file_readcomp", NULL);
      return EOF; 
    }
  }
  return (result);
}

/**
 * ccp4_file_readshortcomp:
 * @param cfile (CCP4File *)
 * @param buffer (uint8 *) buffer
 * @param nitems (size_t) number of items 
 *
 * short complex {short,short} read function.  Reads @nitems complex from stream
 * @cfile->stream to @buffer. Allows short count when eof is detected (
 * buffered input only).
 *
 * @return number of complex read on success, EOF on failure
 */
int ccp4_file_readshortcomp (CCP4File *cfile, uint8 *buffer, size_t nitems)
{
  int i, n, result;

  if (!cfile)  {
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr), 
		  "ccp4_file_readshortcomp", NULL);
    return EOF; }

  if ( !cfile->read || cfile->iostat) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		"ccp4_file_readshortcomp", NULL);
    return EOF; }

  if (cfile->last_op == WRITE_OP)
    if (ccp4_file_raw_seek(cfile,0L,SEEK_CUR) == -1) {
      ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_readshortcomp", NULL);
      return EOF; }   

  n = _item_sizes[COMP32] * nitems;
  if ( (result = ccp4_file_raw_read (cfile, (char *) buffer, n)) != n) {
      ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_readshortcomp", NULL);
    if (cfile->stream && !feof(cfile->stream))
      return EOF; } 

  result /= _item_sizes[COMP32];

  n = result;
  if (cfile->iconvert != nativeIT) {
    n *= 2;                  /* pairs of ints */
    {
      if ((cfile->iconvert==DFNTI_MBO && nativeIT==DFNTI_IBO) ||
          (cfile->iconvert==DFNTI_IBO && nativeIT==DFNTI_MBO)) {
        char j;
        for (i=0; i < n*2; i+=2) {
          j = buffer[i];
          buffer[i] = buffer[i+1];
          buffer[i+1] = j; } }
      else {
        ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
                    "ccp4_file_readshortcomp", NULL);
        return EOF; }
    }
  }
  return (result);
}

/**
 * ccp4_file_readfloat:
 * @param cfile (CCP4File *)
 * @param buffer (uint8 *) buffer
 * @param nitems (size_t) number of items 
 *
 * float read function.  Reads @nitems floats from stream
 * @cfile->stream to @buffer.
 *
 * @return number of floats read on success, EOF on failure
 */
int ccp4_file_readfloat (CCP4File *cfile, uint8 *buffer, size_t nitems)
{
  int i, n, result;

  if (!cfile)  {
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr), 
		  "ccp4_file_readfloat", NULL);
    return EOF; }

  if (!cfile->read || cfile->iostat) {
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		  "ccp4_file_readfloat", NULL);
      return EOF; }

  if (cfile->last_op == WRITE_OP)
    if (ccp4_file_raw_seek(cfile,0L,SEEK_CUR) == -1) {
      ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_readfloat", NULL);
      return EOF; }   

  n = _item_sizes[FLOAT32] * nitems;
  if ( (result = ccp4_file_raw_read (cfile, (char *) buffer, n)) != n) { 
    ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_readfloat", NULL);
    if (cfile->stream && !feof(cfile->stream)) 
      return EOF; } 

  result /= _item_sizes[FLOAT32];

  n = result;
  if (cfile->fconvert != nativeFT) {
    switch (cfile->fconvert) {     /* get to BE IEEE */
    case DFNTF_VAX :
      vaxF2ieeeF((union float_uint_uchar *) buffer, n);
      break;   
    case DFNTF_CONVEXNATIVE :
      convexF2ieeeF((union float_uint_uchar *) buffer, n);
      break;
    case DFNTF_BEIEEE :
      break;
    case DFNTF_LEIEEE :
      {
	char j;
	for (i=0; i < n*4; i+=4) {
	  j = buffer[i];
	  buffer[i] = buffer[i+3];
	  buffer[i+3] = j;
	  j = buffer[i+1];
	  buffer[i+1] = buffer[i+2];
	  buffer[i+2] =j; }
      }
      break;
    default :
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		  "ccp4_file_readfloat", NULL);
      return EOF; 
    }
    switch (nativeFT) {
    case DFNTF_BEIEEE :
      break;                      /* done enough */
    case DFNTF_LEIEEE :
      {
	char j;
	for (i=0; i < n*4; i+=4) {
	  j = buffer[i];
	  buffer[i] = buffer[i+3];
	  buffer[i+3] = j;
	  j = buffer[i+1];
	  buffer[i+1] = buffer[i+2];
	  buffer[i+2] =j; }
      }
      break;
    case DFNTF_CONVEXNATIVE :
      ieeeF2convexF((union float_uint_uchar *) buffer, n);
      break;
    case DFNTF_VAX :
      ieeeF2vaxF((union float_uint_uchar *) buffer, n);
      break;
    default :
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		  "ccp4_file_readfloat", NULL);
      return EOF; 
    }
  }
  return (result);
}

/**
 * ccp4_file_readint:
 * @param cfile (CCP4File *)
 * @param buffer (uint8 *) buffer
 * @param nitems (size_t) number of items 
 *
 * integer read function.  Reads @nitems int from stream
 * @cfile->stream to @buffer.
 *
 * @return number of int read on success, EOF on failure
 */
int ccp4_file_readint (CCP4File *cfile, uint8 *buffer, size_t nitems)
{
  int n, result;

  if (!cfile)  {
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr), 
		  "ccp4_file_readint", NULL);
    return EOF; }

  if ( !cfile->read || cfile->iostat) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		"ccp4_file_readint", NULL);
    return EOF; }

  if (cfile->last_op == WRITE_OP)
    if (ccp4_file_raw_seek(cfile,0L,SEEK_CUR) == -1) {
      ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_readint", NULL);
      return EOF; }   

  n = _item_sizes[CCP4_INT32] * nitems;
  if ( (result = ccp4_file_raw_read (cfile, (char *) buffer, n)) != n) {
    ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_readint", NULL);
    if (cfile->stream && !feof(cfile->stream)) 
      return EOF; } 

  result /= _item_sizes[CCP4_INT32];

  n = result;

  if (cfile->iconvert != nativeIT) {
    if ((cfile->iconvert==DFNTI_MBO && nativeIT==DFNTI_IBO) ||
	(cfile->iconvert==DFNTI_IBO && nativeIT==DFNTI_MBO)) {
      char j;
      int i;
      for (i=0; i < n*4; i+=4) {
	j = buffer[i];
	buffer[i] = buffer[i+3];
	buffer[i+3] = j;
	j = buffer[i+1];
	buffer[i+1] = buffer[i+2];
	buffer[i+2] =j; }
    } else { 
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		"ccp4_file_readint", NULL);
      return EOF; }
  }
  return (result);
}

/**
 * ccp4_file_readshort:
 * @param cfile (CCP4File *)
 * @param buffer (uint8 *) buffer
 * @param nitems (size_t) number of items 
 *
 * short read function.  Reads @nitems shorts from stream
 * @cfile->stream to @buffer.
 *
 * @return number of shorts read on success, EOF on failure
 */
int ccp4_file_readshort (CCP4File *cfile, uint8 *buffer, size_t nitems)
{
  int n, result;

  if (!cfile)  {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr), 
		"ccp4_file_readshort", NULL);
    return EOF; }

  if (!cfile->read || cfile->iostat) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		"ccp4_file_readshort", NULL);
    return EOF; }

  if (cfile->last_op == WRITE_OP)
    if (ccp4_file_raw_seek(cfile,0L,SEEK_CUR) == -1) {
      ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_readshort", NULL);
      return EOF; }   

  n = _item_sizes[CCP4_INT16] * nitems;
  if ( (result = ccp4_file_raw_read (cfile, (char *) buffer, n)) != n) {
    ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_readshort", NULL);
    if (cfile->stream && !feof(cfile->stream)) 
      return EOF; } 

  result /= _item_sizes[CCP4_INT16];

  n = result;
  if (cfile->iconvert != nativeIT) {
    if ((cfile->iconvert==DFNTI_MBO && nativeIT==DFNTI_IBO) ||
	(cfile->iconvert==DFNTI_IBO && nativeIT==DFNTI_MBO)) {
      char j;
      int i;
      for (i=0; i < n*2; i+=2) {
	j = buffer[i];
	buffer[i] = buffer[i+1];
	buffer[i+1] = j; } 
    } else {
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		  "ccp4_file_readshort", NULL);
      return EOF; }
  }
  return (result);
}

/**
 * ccp4_file_readchar:
 * @param cfile (CCP4File *)
 * @param buffer (uint8 *) buffer
 * @param nitems (size_t) number of items 
 *
 * character read function.  Reads @nitems characters from stream
 * @cfile->stream to @buffer.
 *
 * @return number of characters read on success, EOF on failure
 */
int ccp4_file_readchar (CCP4File *cfile, uint8 *buffer, size_t nitems)
{
  size_t result;

  if (! cfile) {
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr), 
		  "ccp4_file_readchar", NULL);
      return EOF; }

  if (!cfile->read || cfile->iostat) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		"ccp4_file_readchar", NULL);
    return EOF; }

  if (cfile->last_op == WRITE_OP)
    if (ccp4_file_raw_seek(cfile,0L,SEEK_CUR) == -1) {
      ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_readchar", NULL);
      return EOF; }   

  if ( (result = ccp4_file_raw_read (cfile, (char *) buffer, nitems)) != nitems) {
    ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_readchar", NULL);
    if (cfile->stream && !feof(cfile->stream))
      return EOF; } 

  return (result);
}

/**
 * ccp4_file_write:
 * @param cfile (CCP4File *)
 * @param buffer (uint8 *) buffer
 * @param nitems (size_t) number of items 
 *
 * mode dependent write function.  Write @nitems items from @buffer
 * to @cfile->stream as determined by cfile->mode.
 *
 * @return number of items written on success, EOF on failure
 */
int ccp4_file_write (CCP4File *cfile, const uint8 *buffer, size_t nitems)
{
  size_t result;

  result = cfile->_write(cfile, buffer, nitems);

  if ( result != nitems) 
    ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_write", NULL);
  
  return (result);
}

/**
 * ccp4_file_writecomp:
 * @param cfile (CCP4File *)
 * @param buffer (uint8 *) buffer
 * @param nitems (size_t) number of items 
 *
 * complex {float,float} write function.  Write @nitems items from @buffer
 * to @cfile->stream.
 *
 * @return number of complex items written on success, EOF on failure
 */
int ccp4_file_writecomp (CCP4File *cfile, const uint8 *buffer, size_t nitems)
{
  size_t result = 0, n;

  if (!cfile) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr), 
		"ccp4_file_writecomp", NULL);
    return EOF; }

  if (!cfile->write ||cfile->iostat) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		"ccp4_file_writecomp", NULL);
    return EOF;}

  if (cfile->last_op == READ_OP) 
    if (ccp4_file_raw_seek(cfile,0L,SEEK_CUR) == -1) {
      ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_writecomp", NULL);
      return EOF; }   
      
  n = nitems * _item_sizes[COMP64];
  
  if (cfile->fconvert != nativeFT) {
    char out_buffer[8];
    const char *out_ptr = (char *) buffer;
    size_t i;
    for (i = 0; i != nitems; i++) {
      switch (nativeFT) {
      case DFNTF_BEIEEE :
        memcpy(out_buffer, out_ptr, _item_sizes[COMP64]);
        out_ptr += _item_sizes[COMP64];
        break;                      
      case DFNTF_LEIEEE : 
        out_buffer[3] = *out_ptr++;
        out_buffer[2] = *out_ptr++;
        out_buffer[1] = *out_ptr++;
        out_buffer[0] = *out_ptr++;
        out_buffer[7] = *out_ptr++;
        out_buffer[6] = *out_ptr++;
        out_buffer[5] = *out_ptr++;
        out_buffer[4] = *out_ptr++;
        break;
    case DFNTF_CONVEXNATIVE :
        memcpy(out_buffer, out_ptr, _item_sizes[COMP64]);
        out_ptr += _item_sizes[COMP64];
        ieeeF2convexF((union float_uint_uchar *) out_buffer, 2);
        break;
    case DFNTF_VAX :
        memcpy(out_buffer, out_ptr, _item_sizes[COMP64]);
        out_ptr += _item_sizes[COMP64];
        ieeeF2vaxF((union float_uint_uchar *) out_buffer, 2);
        break;
    default :
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		  "CCP4_File::writecomp", NULL);
    }
    switch (cfile->fconvert) {    
    case DFNTF_VAX :
        vaxF2ieeeF((union float_uint_uchar *) out_buffer, 2);
        break;   
    case DFNTF_CONVEXNATIVE :
        convexF2ieeeF((union float_uint_uchar *) out_buffer, 2);
        break;
    case DFNTF_BEIEEE :
        break;
    case DFNTF_LEIEEE :
        {
          char j;
          j = out_buffer[0];
          out_buffer[0] = out_buffer[3];
          out_buffer[3] = j;
          j = out_buffer[1];
          out_buffer[1] = out_buffer[2];
          out_buffer[2] =j; 
          j = out_buffer[4];
          out_buffer[4] = out_buffer[7];
          out_buffer[7] = j;
          j = out_buffer[5];
          out_buffer[5] = out_buffer[6];
          out_buffer[6] =j; 
        }
        break;
    default :
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		  "ccp4_file_writecomp", NULL);
      return EOF; 
    }
    result += ccp4_file_raw_write (cfile, out_buffer, _item_sizes[COMP64]);
    }
  } else {
    result = ccp4_file_raw_write (cfile, (char *) buffer, n);
  }

  if (result != n)
    ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_writecomp", NULL);
    
  return (result / _item_sizes[COMP64]);
}

/**
 * ccp4_file_writeshortcomp:
 * @param cfile (CCP4File *)
 * @param buffer (uint8 *) buffer
 * @param nitems (size_t) number of items 
 *
 * short complex {short,short} write function.  Write @nitems items from @buffer
 * to @cfile->stream.
 *
 * @return number of complex items written on success, EOF on failure
 */
int ccp4_file_writeshortcomp (CCP4File *cfile, const uint8 *buffer, size_t nitems)
{
  size_t result = 0, n;

  if (!cfile) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr), 
		"ccp4_file_writeshortcomp", NULL);
    return EOF; }

  if (!cfile->write ||cfile->iostat) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		"ccp4_file_writeshortcomp", NULL);
    return EOF;}

  if (cfile->last_op == READ_OP) 
    if (ccp4_file_raw_seek(cfile,0L,SEEK_CUR) == -1) {
      ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_writeshortcomp", NULL);
      return EOF; }   
      
  n = nitems * _item_sizes[COMP32];

  if (cfile->iconvert != nativeIT) {
    char out_buffer[4];
    const char *out_ptr = (char *) buffer;
    size_t i;
    for (i = 0; i != nitems; i++) {
      if ((cfile->iconvert==DFNTI_MBO && nativeIT==DFNTI_IBO) ||
	  (cfile->iconvert==DFNTI_IBO && nativeIT==DFNTI_MBO)) {
        out_buffer[1] = *out_ptr++;
        out_buffer[0] = *out_ptr++;
        out_buffer[3] = *out_ptr++;
        out_buffer[2] = *out_ptr++; 
      } else { 
        ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		"ccp4_file_writeshortcomp", NULL);
        return EOF; }
    result += ccp4_file_raw_write (cfile, out_buffer, _item_sizes[COMP32]);
    }
   } else {
    result = ccp4_file_raw_write (cfile, (char *) buffer, n);
  }
  
  if ( result != n) 
    ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_writeshortcomp", NULL);

  return (result / _item_sizes[COMP32]);
}

/**
 * ccp4_file_writefloat:
 * @param cfile (CCP4File *)
 * @param buffer (uint8 *) buffer
 * @param nitems (size_t) number of items 
 *
 * float write function.  Write @nitems items from @buffer
 * to @cfile->stream.
 *
 * Returns number of floats written on success, EOF on failure
 */
int ccp4_file_writefloat (CCP4File *cfile, const uint8 *buffer, size_t nitems)
{
  size_t result = 0, n;

  if (!cfile) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr), 
		"ccp4_file_writefloat", NULL);
    return EOF; }

  if (!cfile->write ||cfile->iostat) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		"ccp4_file_writefloat", NULL);
    return EOF;}

  if (cfile->last_op == READ_OP) 
    if (ccp4_file_raw_seek(cfile,0L,SEEK_CUR) == -1) {
      ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_writefloat", NULL);
      return EOF; }   
      
  n = nitems * _item_sizes[FLOAT32];

  if (cfile->fconvert != nativeFT) {
    char out_buffer[4];
    const char *out_ptr = (char *) buffer;
    size_t i;
    for (i = 0; i != nitems; i++) {
      switch (nativeFT) {
      case DFNTF_BEIEEE :
        memcpy(out_buffer, out_ptr, _item_sizes[FLOAT32]);
        out_ptr += _item_sizes[FLOAT32];
        break;                      
      case DFNTF_LEIEEE : 
        out_buffer[3] = *out_ptr++;
        out_buffer[2] = *out_ptr++;
        out_buffer[1] = *out_ptr++;
        out_buffer[0] = *out_ptr++;
        break;
    case DFNTF_CONVEXNATIVE :
        memcpy(out_buffer, out_ptr, _item_sizes[FLOAT32]);
        out_ptr += _item_sizes[FLOAT32];
        ieeeF2convexF((union float_uint_uchar *) out_buffer, 1);
        break;
    case DFNTF_VAX :
        memcpy(out_buffer, out_ptr, _item_sizes[FLOAT32]);
        out_ptr += _item_sizes[FLOAT32];
        ieeeF2vaxF((union float_uint_uchar *) out_buffer, 1);
        break;
    default :
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		  "CCP4_File::writefloat", NULL);
    }
    switch (cfile->fconvert) {    
    case DFNTF_VAX :
        vaxF2ieeeF((union float_uint_uchar *) out_buffer, 1);
        break;   
    case DFNTF_CONVEXNATIVE :
        convexF2ieeeF((union float_uint_uchar *) out_buffer, 1);
        break;
    case DFNTF_BEIEEE :
        break;
    case DFNTF_LEIEEE :
        {
          char j;
          j = out_buffer[0];
          out_buffer[0] = out_buffer[3];
          out_buffer[3] = j;
          j = out_buffer[1];
          out_buffer[1] = out_buffer[2];
          out_buffer[2] =j; 
        }
        break;
    default :
      ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		  "ccp4_file_writefloat", NULL);
      return EOF; 
    }
    result += ccp4_file_raw_write (cfile, out_buffer, _item_sizes[FLOAT32]);
    }
  } else {
    result = ccp4_file_raw_write (cfile, (char *) buffer, n);
  }

  if (result != n)
    ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_writefloat", NULL);
    
  return (result / _item_sizes[FLOAT32]);
}
  
/**
 * ccp4_file_writeint:
 * @param cfile (CCP4File *)
 * @param buffer (uint8 *) buffer
 * @param nitems (size_t) number of items 
 *
 * int write function.  Write @nitems items from @buffer
 * to @cfile->stream.
 *
 * @return number of int written on success, EOF on failure
 */
int ccp4_file_writeint (CCP4File *cfile, const uint8 *buffer, size_t nitems)
{
  size_t result = 0, n;

  if (!cfile) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr), 
		"ccp4_file_writeint", NULL);
    return EOF; }

  if (!cfile->write ||cfile->iostat) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		"ccp4_file_writeint", NULL);
    return EOF;}

  if (cfile->last_op == READ_OP) 
    if (ccp4_file_raw_seek(cfile,0L,SEEK_CUR) == -1) {
      ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_writeint", NULL);
      return EOF; }   
      
  n = nitems * _item_sizes[CCP4_INT32];

  if (cfile->iconvert != nativeIT) {
    char out_buffer[4];
    const char *out_ptr = (char *) buffer;
    size_t i;
    for (i = 0; i != nitems; i++) {
      if ((cfile->iconvert==DFNTI_MBO && nativeIT==DFNTI_IBO) ||
	  (cfile->iconvert==DFNTI_IBO && nativeIT==DFNTI_MBO)) {
        out_buffer[3] = *out_ptr++;
        out_buffer[2] = *out_ptr++;
        out_buffer[1] = *out_ptr++;
        out_buffer[0] = *out_ptr++; 
      } else { 
        ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		"ccp4_file_writeint", NULL);
        return EOF; }
    result += ccp4_file_raw_write (cfile, out_buffer, _item_sizes[CCP4_INT32]);
    }
   } else {
    result = ccp4_file_raw_write (cfile, (char *) buffer, n);
  }

  if ( result != n) 
    ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_writeint", NULL);

  return (result / _item_sizes[CCP4_INT32]);
}

/**
 * ccp4_file_writeshort:
 * @param cfile (CCP4File *)
 * @param buffer (uint8 *) buffer
 * @param nitems (size_t) number of items 
 *
 * short write function.  Write @nitems items from @buffer
 * to @cfile->stream.
 *
 * @return number of short written on success, EOF on failure
 */
int ccp4_file_writeshort (CCP4File *cfile, const uint8 *buffer, size_t nitems)
{
  size_t result = 0, n;

  if (!cfile) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr), 
		"ccp4_file_writeshort", NULL);
    return EOF; }

  if (!cfile->write ||cfile->iostat) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		"ccp4_file_writeshort", NULL);
    return EOF;}

  if (cfile->last_op == READ_OP) 
    if (ccp4_file_raw_seek(cfile,0L,SEEK_CUR) == -1) {
      ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_writeshort", NULL);
      return EOF; }   
      
  n = nitems * _item_sizes[CCP4_INT16];
  
  if (cfile->iconvert != nativeIT) {
    char out_buffer[2];
    const char *out_ptr = (char *) buffer;
    size_t i;
    for (i = 0; i != nitems; i++) {
      if ((cfile->iconvert==DFNTI_MBO && nativeIT==DFNTI_IBO) ||
	  (cfile->iconvert==DFNTI_IBO && nativeIT==DFNTI_MBO)) {
        out_buffer[1] = *out_ptr++;
        out_buffer[0] = *out_ptr++;
      } else { 
        ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		"ccp4_file_readint", NULL);
        return EOF; }
    result += ccp4_file_raw_write (cfile, out_buffer, _item_sizes[CCP4_INT16]);
    }
   } else {
    result = ccp4_file_raw_write (cfile, (char *) buffer, n);
  }
  
  if ( result != n) 
    ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_writeshort", NULL);

  return (result / _item_sizes[CCP4_INT16]);
}

/**
 * ccp4_file_writechar:
 * @param cfile (CCP4File *)
 * @param buffer (uint8 *) buffer
 * @param nitems (size_t) number of items 
 *
 * char write function.  Write @nitems items from @buffer
 * to @cfile->stream.
 *
 * @return number of bytes written on success, EOF on failure
 */
int ccp4_file_writechar (CCP4File *cfile, const uint8 *buffer, size_t nitems)
{
  size_t result;

  if (!cfile) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr), 
		"ccp4_file_writechar", NULL);
    return EOF; }

  if (!cfile->write ||cfile->iostat) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_BadMode), 
		"ccp4_file_writechar", NULL);
    return EOF;}

  if (cfile->last_op == READ_OP) 
    if (ccp4_file_raw_seek(cfile,0L,SEEK_CUR) == -1) {
      ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_writechar", NULL);
      return EOF; }   
  
  if ( (result = ccp4_file_raw_write (cfile, (char *) buffer, nitems)) != nitems) 
    ccp4_signal(CCP4_ERRLEVEL(3), "ccp4_file_writechar", NULL);

  return (result);
}

/**
 * ccp4_file_seek:
 * @param cfile (CCP4File *)
 * @param offset (long) offset in items
 * @param whence (int) SEEK_SET, SEEK_CUR, or SEEK_END
 *
 * seeks on file by offset items.  SEEK_SET is relative
 * to start of file, SEEK_CUR to current, SEEK_END to
 * end.
 *
 * @return 0 on success, -1 on failure
 */
int ccp4_file_seek (CCP4File *cfile, long offset, int whence)
{
  int result;

  if (!cfile) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr), 
		"ccp4_file_seek", NULL);
    return -1; }
        
  result = ccp4_file_raw_seek(cfile, offset*cfile->itemsize, whence);

  if (result != -1) ccp4_file_clearerr(cfile);

  return ((result == -1) ? -1 : 0);
}

/**
 * ccp4_file_rewind:
 * @param cfile (CCP4File *)
 *
 * Seek to start of file.  Clear error status.
 *
 * @return 0 on success, EOF on failure
 */
void ccp4_file_rewind (CCP4File *cfile)
{
  if (!cfile) {   
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr), 
		"ccp4_file_rewind", NULL);
    return ; }
   
  if (ccp4_file_raw_seek(cfile, 0, SEEK_SET)) {
      ccp4_signal(CCP4_ERRLEVEL(3), 
	  	  "ccp4_file_rewind", NULL);
  } else {
    ccp4_file_clearerr(cfile);
  }
}

/**
 * ccp4_file_length:
 * @param cfile (CCP4File *)
 * 
 * Length of file on disk.
 * @return length of @cfile on success, EOF on failure
 */
long ccp4_file_length (CCP4File *cfile)
{
#if defined _MSC_VER
  struct _stat st;
#else
  struct stat st;
#endif

  if (!cfile) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr), 
		"ccp4_file_length", NULL);
    return EOF; } 
        
  cfile->last_op = IRRELEVANT_OP;
  
  if (cfile->buffered && cfile->stream)
      fflush (cfile->stream);
#if defined _MSC_VER
    _fstat(cfile->stream ? _fileno(cfile->stream) : cfile->fd, &st);
#else
    fstat(cfile->stream ? fileno(cfile->stream) : cfile->fd, &st);
#endif
    cfile->length = st.st_size;

  return (st.st_size);
}

/**
 * ccp4_file_tell:
 * @param cfile (CCP4File *)
 *
 * Current location in file, uses either ftell or lseek.
 * @return current offset of @cfile in bytes.
 */
long ccp4_file_tell (CCP4File *cfile)
{
  long result;
  
  if (! cfile) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr), 
		"ccp4_file_tell", NULL);
    return EOF; } 

  cfile->last_op = IRRELEVANT_OP;

  if (cfile->buffered && cfile->stream) {
#if !defined (_MSC_VER)
    if ( cfile->last_op == WRITE_OP ) fflush (cfile->stream);
#endif	
    result = (long) ftell(cfile->stream);
  } else
#if defined _MSC_VER
    result = _lseek(cfile->fd, 0L, SEEK_CUR);
#else
    result = lseek(cfile->fd, 0L, SEEK_CUR);
#endif
  
  cfile->loc = result;
    
  return (result);
}

/**
 * ccp4_file_feof:
 * @param cfile (CCP4File *)
 *
 * @return true if @cfile is at EoF.
 *
 */
int ccp4_file_feof(CCP4File *cfile)
{
  if (!cfile) {
    ccp4_signal(CCP4_ERRLEVEL(3) | CCP4_ERRNO(CIO_NullPtr), 
		"ccp4_file_feof", NULL);
    return EOF; } 

  return (cfile->stream) ? feof(cfile->stream) : cfile->loc >= cfile->length;
}

/**
 * ccp4_file_clearerr:
 * @param cfile (CCP4File *)
 *
 * Clears error status of @cfile.
 *
 */
void ccp4_file_clearerr(CCP4File *cfile)
{
  if (!cfile || !cfile->stream)
    return;
  cfile->iostat = CIO_Ok;
  clearerr(cfile->stream);
}

/**
 * ccp4_file_fatal:
 * @param cfile (CCP4File *)
 *
 * Die with error message based on @cfile error status.
 */
void ccp4_file_fatal (CCP4File *cfile, char *message)
{
  char *buff;
  size_t l;

  if (!cfile)
    ccp4_signal(CCP4_ERRLEVEL(4) | CCP4_ERRNO(CIO_NullPtr), "ccp4_file_fatal",
                NULL);
    
  l = strlen (message) + strlen (cfile->name) + 1;
  if ( !(buff = malloc(l)))
      ccp4_signal(CCP4_ERRLEVEL(4), "ccp4_file_fatal", NULL);
  buff[0] = '\0';
  strcat (buff, message);
  strcat (buff, cfile->name);
  ccp4_fatal(buff);
}

/**
 * ccp4_file_error:
 * @param cfile (CCP4File *)
 *
 * print error mesage.
 * @return associated error code
 */
int ccp4_file_error(CCP4File *cfile)
{  
  if (!cfile->iostat)
    return 0;
  fprintf(stderr,"%s %s \n",
            cfile->name,ccp4_strerror(cfile->iostat));
  return CCP4_ERRGETCODE(cfile->iostat); 
}

/**
 * ccp4_file_flush:
 * @param cfile (CCP4File *)
 *
 * flush buffer contents of @cfile
 */
void ccp4_file_flush(CCP4File *cfile)
{
   if (cfile && cfile->stream && cfile->buffered)
     fflush(cfile->stream);
}

/**
 * ccp4_file_print:
 * @param cfile (CCP4File *)
 *
 * @return @cfile information in char array for printing.
 */
char *ccp4_file_print(CCP4File *cfile, char *msg_start, char *msg_end)
{
  char *msg_curr = msg_start;
  
  if (!cfile)
    return msg_start;
    
  if (cfile->name)
    if ((msg_end - msg_curr) > strlen(cfile->name)) {
      strcpy(msg_curr,cfile->name);
      msg_curr = strrchr(msg_curr,'\0'); }

  if (cfile->open) {
    if ((msg_end - msg_curr) > 6 ) {
      strcat(msg_start, " opened");
      msg_curr = strrchr(msg_curr,'\0'); }
  } else {
    if ((msg_end - msg_curr) > 7 ) {
      strcat(msg_start, " closed");
      msg_curr = strrchr(msg_curr,'\0'); }
  }

  if (cfile->append) {
    if ((msg_end - msg_curr) > 13 ) {
      strcat(msg_start, ", append mode");
      msg_curr = strrchr(msg_curr,'\0'); }
  } else if (cfile->read && cfile->write) {
    if ((msg_end - msg_curr) > 17 ) {
      strcat(msg_start, ", read-write mode");
      msg_curr = strrchr(msg_curr,'\0'); }
  } else if (cfile->write) {
    if ((msg_end - msg_curr) > 12 ) {
      strcat(msg_start, ", write mode");
      msg_curr = strrchr(msg_curr,'\0'); }
  } else {
    if ((msg_end - msg_curr) > 11 ) {
      strcat(msg_start, ", read mode");
      msg_curr = strrchr(msg_curr,'\0'); }
  }

  return msg_curr;
}

