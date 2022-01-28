//  $Id: mpfile_.h $
//  =================================================================
//
//   CCP4 SRS Library: Storage, Retrieval and Search support for
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
//  **** Module  :  mpfile_  <interface>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::MPFile  - multi-part file support
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2014
//
//  =================================================================
//

#ifndef  MPFILE_H
#define  MPFILE_H

#include "mmdb2/mmdb_io_file.h"

namespace ccp4srs  {

  //  ===========================  MPFile  ===========================

  DefineStructure(MPFilePos);

  struct MPFilePos  {
    mmdb::byte chunkNo;   // file part number (1-254)
    long       offset;    // offset in chunk
    void Copy  ( RMPFilePos dbf_pos );
    void write ( mmdb::io::RFile f );
    void read  ( mmdb::io::RFile f );
  };

  DefineClass(MPFile)

  class MPFile  {

    public :
      MPFile ();
      ~MPFile();

      void  assign  ( mmdb::cpstr FName, bool UniB=true );
      bool  reset   ( bool rdOnly=false );
      bool  rewrite ( long chunk_size=1000000000  );
      bool  append  ( long chunk_size=1000000000  );

      mmdb::pstr FileName    () { return fpath; }
      void  getPosition ( RMPFilePos MPFilePos );
      bool  seek        ( RMPFilePos MPFilePos );
      bool  MPFileEnd   ();
      bool  Success     ();
      void  shut        ();

      mmdb::io::PFile   getChunkRead ();
      mmdb::io::PFile   getChunkWrite();

    protected :
      mmdb::io::PPFile f;
      mmdb::pstr       fpath;
      long             chunkSize;
      mmdb::byte       chunkNo;
      bool             uniBinary,readOnly;

      mmdb::pstr getChunkPath  ();
      bool getChunk      ( int toChunkNo, bool onRead );
      bool checkChunkEnd ( bool onRead );

    private :
      mmdb::pstr chunkPath;
      int        nchunks;

  };

}  // namespace ccp4srs

#endif  //  MPFILE_H
