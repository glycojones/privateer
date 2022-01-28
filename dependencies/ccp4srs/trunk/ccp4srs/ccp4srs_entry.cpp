//  $Id: ccp4srs_entry.cpp $
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
//    18.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  ccp4srs_entry  <implementation>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Entry  - base class for all SRS entries
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#include "ccp4srs_entry.h"
#include "ccp4srs_defs.h"

namespace ccp4srs  {

  Entry::Entry()  {}

  Entry::~Entry() {}

  int Entry::write ( mmdb::io::RFile f, int version, PMemIO memIO ) {
  PMemIO mem_io;
  int    rc;

    if (memIO)  mem_io = memIO;
          else  mem_io = new MemIO();

    write_mem ( mem_io,version );
    rc = mem_io->write ( f );

    if (!memIO)
      delete mem_io;

    return rc;

  }

  int Entry::read ( mmdb::io::RFile f, int version, PMemIO memIO )  {
  PMemIO mem_io;
  int    rc;

    if (memIO)  mem_io = memIO;
          else  mem_io = new MemIO();

    rc = mem_io->read ( f );
    if (rc==MemIO::Ok)  {
      if (!read_mem(mem_io,version))
        rc = MemIO::MemReadError;
    }

    if (!memIO)
      delete mem_io;

    return rc;

  }

  void Entry::write_mem ( PMemIO memIO, int version )  {
  UNUSED_ARGUMENT(memIO);
  UNUSED_ARGUMENT(version);
  }

  bool Entry::read_mem ( PMemIO memIO, int version, bool * Ok )  {
  UNUSED_ARGUMENT(memIO);
  UNUSED_ARGUMENT(version);
  UNUSED_ARGUMENT(Ok);
    return true;
  }

}  // namespace ccp4srs

