//  $Id: memio_.h $
//  =================================================================
//
//   CCP4 SRS Library: Storage, Retrievak and Search support for
//   CCP4 ligand data.
//
//   Copyright (C) Eugene Krissinel 2010-2013.
//
//   This library is free software: you can redistribute it and/or
//   modify it under the terms of the GNU Lesser General Public
//   License version 3, modified in accordance with the provisions
//   of the license to address the requirements of UK law.
//
//   You should have received a copy of the modified GNU Lesser
//   General Public License along with this library. If not, copies
//   may be downloaded from http://www.ccp4.ac.uk/ccp4license.php
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU Lesser General Public License for more details.
//
//  =================================================================
//
//    03.02.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  memio_  <interface>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::MemIO  - buffered I/O
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2014
//
//  =================================================================
//

#ifndef MEMIO_H
#define MEMIO_H

#include "mmdb2/mmdb_io_file.h"

namespace ccp4srs  {

  DefineClass(MemIO);

  class MemIO  {

    public:

      enum { Ok,CantOpenFile,EmptyFile,ReadError,WriteError,MemReadError,
             CompressionError,UncompressionError,SizeError };

      MemIO();
      virtual ~MemIO();

      void  setCompressionLevel ( int compressionLevel ); // 0,1,2,3

      int   read  ( mmdb::cpstr  fileName,
                    mmdb::io::GZ_MODE gzMode=mmdb::io::GZM_NONE );
      int   read  ( mmdb::io::RFile f );

      int   write ( mmdb::cpstr  fileName,
                    mmdb::io::GZ_MODE gzMode=mmdb::io::GZM_NONE );
      int   write ( mmdb::io::RFile f );

      int   length() { return buffer_length; }

      void  reset();
      void  free ();

      void * get_buffer        ( const int length );
      void   set_buffer_length ( const int len )  { buffer_length = len; }
      void   get_buffer        ( mmdb::pstr * buf, int * length );

      void  write_buffer ( const void * src,  const int length  );
      bool  read_buffer  ( void       * dest, const int length,
                           bool * Ok = NULL);

      void  put_integer   ( int        I );
      void  put_ishort    ( int        I );  // saves as short
      void  put_ibyte     ( int        I );  // saves as byte
      void  put_byte      ( mmdb::byte B );
      void  put_short     ( short      I );
      void  put_long      ( long       I );
      void  put_word      ( mmdb::word W );
      void  put_real      ( mmdb::realtype  R );
      void  put_float     ( mmdb::realtype  R );  // stores single precision
      void  put_shortreal ( mmdb::shortreal R );
      void  put_line      ( mmdb::cpstr     L );  // for static strings
      void  put_string    ( mmdb::cpstr     L );  // for dynamic strings
      void  put_bool      ( bool   B );

      bool  get_integer   ( int        & I, bool * Ok = NULL );
      bool  get_ishort    ( int        & I, bool * Ok = NULL );
      bool  get_ibyte     ( int        & I, bool * Ok = NULL );
      bool  get_byte      ( mmdb::byte & B, bool * Ok = NULL );
      bool  get_short     ( short      & I, bool * Ok = NULL );
      bool  get_long      ( long       & I, bool * Ok = NULL );
      bool  get_word      ( mmdb::word      & W, bool * Ok = NULL );
      bool  get_real      ( mmdb::realtype  & R, bool * Ok = NULL );
      bool  get_float     ( mmdb::realtype  & R, bool * Ok = NULL );
      bool  get_shortreal ( mmdb::shortreal & R, bool * Ok = NULL );
      bool  get_line      ( mmdb::pstr        L, bool * Ok = NULL );
      bool  get_string    ( mmdb::pstr      & L, bool * Ok = NULL );
      bool  get_bool      ( bool   & B, bool * Ok = NULL );

    private:
      int           compression_level;
      mmdb::bvector buffer,compression_buffer;
      int           buffer_length,buffer_pos,alloc_length,chunk_size;
      int           compression_alloc;

  };

}  // namespace ccp4srs

#endif // MEMIO_H
