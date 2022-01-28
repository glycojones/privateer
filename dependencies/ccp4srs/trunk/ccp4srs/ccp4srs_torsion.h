//  $Id: ccp4srs_torsion.h $
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
//  **** Module  :  ccp4srs_torsion  <interface>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Torsion  - torsion angle restraint class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#ifndef CCP4SRS_TORSION_H
#define CCP4SRS_TORSION_H

#include "ccp4srs_atom.h"

namespace ccp4srs  {

  DefineClass(Torsion);

  typedef char TorsionID[20];

  class Torsion {

    public:
      Torsion();
      virtual ~Torsion();

      inline mmdb::cpstr    id    () { return torsion_id;     }
      inline int            atom1 () { return atom_1;         }
      inline int            atom2 () { return atom_2;         }
      inline int            atom3 () { return atom_3;         }
      inline int            atom4 () { return atom_4;         }
      inline int            period() { return torsion_period; }
      inline mmdb::realtype value () { return torsion;        }
      inline mmdb::realtype esd   () { return torsion_esd;    }

      inline void set_period ( int period) { torsion_period = period; }

      static void makeCIFTags ( mmdb::mmcif::PLoop mmCIFLoop );
      int  readFromCIF ( mmdb::mmcif::PLoop mmCIFLoop, int angleNo,
                         Container<Atom> & atoms );
      void writeToCIF  ( mmdb::mmcif::PLoop mmCIFLoop, mmdb::cpstr monID,
                         Container<Atom> & atoms );

      void write_mem ( PMemIO memIO, int version );
      bool read_mem  ( PMemIO memIO, int version, bool * Ok = NULL );

      void copy ( PTorsion trs, mmdb::ivector anmatch=NULL );

    protected:
      TorsionID   torsion_id;     //!< torsion angle id
      int         atom_1,atom_2;  //!< ordinal numbers of atoms
      int         atom_3,atom_4;  ///  that form the angle [0,1,...]
      int         torsion_period; //!< torsion period
      mmdb::realtype torsion;     //!< torsion angle in degrees
      mmdb::realtype torsion_esd; //!< torsion angle esd in degrees;
                                  /// set to 0.0 if not provided

  };

}  // namespace ccp4srs

#endif // CCP4SRS_TORSION_H
