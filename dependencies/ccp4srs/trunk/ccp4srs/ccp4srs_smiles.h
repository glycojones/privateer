//  $Id: ccp4srs_smiles.h $
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
//    12.04.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  ccp4srs_smiles  <interface>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Smiles  - atom description class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2014
//
//  =================================================================
//

#ifndef CCP4SRS_SMILES_H
#define CCP4SRS_SMILES_H

#include "memio_.h"

#include "mmdb2/mmdb_mattype.h"
#include "mmdb2/mmdb_mmcif_.h"

namespace ccp4srs  {

  DefineClass(Smiles);

  class Smiles  {

    public:

      Smiles();
      virtual ~Smiles();

      void freeSmiles();

      int  readFromCIF_rcsb ( mmdb::mmcif::PLoop mmCIFLoop,
                              mmdb::cpstr program );

      bool   isEmpty();

      inline mmdb::cpstr getVersion  ()  { return ver;       }
      inline mmdb::cpstr getSmiles   ()  { return plain;     }
      inline mmdb::cpstr getCanonical()  { return canonical; }
      inline mmdb::cpstr getInChI    ()  { return plain;     }
      inline mmdb::cpstr getInChIKey ()  { return canonical; }

      void copy ( PSmiles smiles );

/*
      static void makeCIFTags ( mmdb::mmcif::PLoop mmCIFLoop );
      int  readFromCIF_ccp4 ( mmdb::mmcif::PLoop mmCIFLoop, int bondNo,
                              Container<Atom> & atoms );
      int  readFromCIF_rcsb ( mmdb::mmcif::PLoop mmCIFLoop, int bondNo,
                              Container<Atom> & atoms );
      int  readFromCIF ( mmdb::mmcif::PStruct mmCIFStruct,
                         Container<Atom> & atoms );
      void writeToCIF  ( mmdb::mmcif::PLoop mmCIFLoop, mmdb::cpstr monID,
                         Container<Atom> & atoms );

      bool compare ( Container<Atom> & atoms, PBond bond,
                     Container<Atom> & bond_atoms );

      bool copy    ( Container<Atom> & atoms, PBond bond,
                     Container<Atom> & bond_atoms );
*/

      void write_mem ( PMemIO memIO, int version );
      bool read_mem  ( PMemIO memIO, int version,
                          bool * Ok = NULL );

    protected:
      mmdb::pstr  ver;       //!< software version
      mmdb::pstr  plain;     //!< plain smiles
      mmdb::pstr  canonical; //!< canonical smiles

  };

}  // namespace ccp4srs

#endif // CCP4SRS_SMILES_H
