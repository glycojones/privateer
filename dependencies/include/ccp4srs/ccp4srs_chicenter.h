//  $Id: ccp4srs_chicenter.h $
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
//  **** Module  :  ccp4srs_chicenter  <interface>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::ChiCenter - chiral center class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#ifndef CCP4SRS_CHICENTER_H
#define CCP4SRS_CHICENTER_H

#include "ccp4srs_atom.h"

namespace ccp4srs  {

  DefineClass(ChiCenter);

  typedef char ChiCenterID[20];

  class ChiCenter  {

    // Chiral center signs
    enum  {noSign,Positive,Negative,Both};

    public:
      ChiCenter();
      virtual ~ChiCenter();

      inline mmdb::cpstr id    () { return chicenter_id; }
      inline int         center() { return atom_center;  }
      inline int         atom1 () { return atom_1;       }
      inline int         atom2 () { return atom_2;       }
      inline int         atom3 () { return atom_3;       }
      inline int         sign  () { return volume_sign;  }

      char get_chirality();

      static void makeCIFTags ( mmdb::mmcif::PLoop mmCIFLoop );
      int  readFromCIF ( mmdb::mmcif::PLoop mmCIFLoop, int chiNo,
                         Container<Atom> & atoms );
      void writeToCIF  ( mmdb::mmcif::PLoop mmCIFLoop, mmdb::cpstr monID,
                         Container<Atom> & atoms );

      void write_mem ( PMemIO memIO, int version );
      bool read_mem  ( PMemIO memIO, int version, bool * Ok = NULL );

      void copy ( PChiCenter chicenter, mmdb::ivector anmatch=NULL );

    protected:
      ChiCenterID chicenter_id;    //!< chiral center id
      int         atom_center;     //!< chiral center atom [0,1,...]
      int         atom_1,atom_2;   //!< ordinal numbers of atoms that
      int         atom_3;          //!< form the center [0,1,...]
      int         volume_sign;     //!< chirality

  };

}  // namespace ccp4srs

#endif // CCP4SRS_CHICENTER_H
