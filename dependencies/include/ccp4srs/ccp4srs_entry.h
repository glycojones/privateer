//  $Id: ccp4srs_entry.h $
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
//  **** Module  :  ccp4srs_entry  <interface>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Entry  - base class for all SRS entries
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#ifndef CCP4SRS_ENTRY_H
#define CCP4SRS_ENTRY_H

#include "memio_.h"

namespace ccp4srs  {

  DefineClass(Entry);

  class Entry  {

    public:

      // entry types
      enum {None,Monomer,Link};

      Entry();
      virtual ~Entry();

      virtual int type()  { return Entry::None; }

      int write ( mmdb::io::RFile f, int version, PMemIO memIO=NULL );
      int read  ( mmdb::io::RFile f, int version, PMemIO memIO=NULL );

    protected:

      virtual void write_mem ( PMemIO memIO, int version );
      virtual bool read_mem  ( PMemIO memIO, int version,
                               bool * Ok = NULL );

  };

}  // namespace ccp4srs


#endif // CCP4SRS_ENTRY_H
