//  $Id: ccp4srs_bond.h $
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
//  **** Module  :  ccp4srs_bond  <interface>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Bond  - atom description class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#ifndef CCP4SRS_BOND_H
#define CCP4SRS_BOND_H

#include "ccp4srs_atom.h"

namespace ccp4srs  {

  DefineClass(Bond);

  class Bond  {

    public:

      //  bond orders
      enum {noOrder,Single,Aromatic,Double,Triple,Deloc,Covalent,Metal};

      Bond();
      virtual ~Bond();

      inline int      atom1     () { return atom_1;          }
      inline int      atom2     () { return atom_2;          }
      inline int      order     () { return bond_order;      }
      inline mmdb::realtype length    () { return bond_length;     }
      inline mmdb::realtype length_esd() { return bond_length_esd; }

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

      void write_mem ( PMemIO memIO, int version );
      bool read_mem  ( PMemIO memIO, int version,
                          bool * Ok = NULL );

      void copy  ( PBond bond, mmdb::ivector anmatch=NULL );

    protected:
      int            atom_1,atom_2;   //!< ordinal numbers of bonded atoms [0,1,.]
      int            bond_order;      //!< bond order
      mmdb::realtype bond_length;     //!< bond length in A; set to 0.0 if
                                      ///< not provided
      mmdb::realtype bond_length_esd; //!< bond length esd in A; set to 0.0 if
                                      /// not provided

  };

}  // namespace ccp4srs

#endif // CCP4SRS_BOND_H
