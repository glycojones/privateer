//  $Id: ccp4srs_index.h $
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
//  **** Module  :  ccp4srs_index  <interface>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Index - SRS data index class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#ifndef CCP4SRS_INDEX_H
#define CCP4SRS_INDEX_H

#include "ccp4srs_monomer.h"

namespace ccp4srs  {

  DefineClass(Index);

  class Index  {

    public :
      mmdb::ResName entryID;    //!< monomer or link ID
      int           nAtoms;     //!< number of atoms in the entry
      int           nNonHAtoms; //!< number of non-hydrogen atoms in the entry
      int           nBonds;     //!< number of bonds in the entry
      long          fGraphPos;  //!< file offset for CGraph
      long          fStructPos; //!< file offset for CCP4SRSMonomer
      mmdb::pstr    Comp1;      //!< composition string
      mmdb::pstr    Comp2;      //!< composition string with leaving atom

      Index ();
      virtual ~Index();

      int   MakeCompositions ( PMonomer monomer );

      int write ( mmdb::io::RFile f, int version, PMemIO memIO=NULL );
      int read  ( mmdb::io::RFile f, int version, PMemIO memIO=NULL );

    protected:

      void  IndexInit();

      virtual void write_mem ( PMemIO memIO, int version );
      virtual bool read_mem  ( PMemIO memIO, int version,
                               bool * Ok = NULL );

  };

}  // namespace ccp4srs


#endif // CCP4SRS_INDEX_H
