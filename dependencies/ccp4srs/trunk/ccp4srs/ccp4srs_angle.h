//  $Id: ccp4srs_angle.h $
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
//  **** Module  :  ccp4srs_angle  <interface>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Angle  - angle restraint class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2014
//
//  =================================================================
//

#ifndef CCP4SRS_ANGLE_H
#define CCP4SRS_ANGLE_H

#include "ccp4srs_atom.h"

namespace ccp4srs  {

  DefineClass(Angle);

  class Angle  {

    public:
      Angle();
      virtual ~Angle();

      inline int            atom1() { return atom_1;    }
      inline int            atom2() { return atom_2;    }
      inline int            atom3() { return atom_3;    }
      inline mmdb::realtype value() { return angle;     }
      inline mmdb::realtype esd  () { return angle_esd; }

      static void makeCIFTags ( mmdb::mmcif::PLoop mmCIFLoop );
      int  readFromCIF ( mmdb::mmcif::PLoop mmCIFLoop, int angleNo,
                         Container<Atom> & atoms );
      void writeToCIF  ( mmdb::mmcif::PLoop mmCIFLoop, mmdb::cpstr monID,
                         Container<Atom> & atoms );

      void write_mem ( PMemIO memIO, int version );
      bool read_mem  ( PMemIO memIO, int version,
                          bool * Ok = NULL );

      void copy ( PAngle ang, mmdb::ivector anmatch=NULL );

    protected:
      int  atom_1,atom_2,atom_3; // ordinal numbers of of atoms
                                 //   that form the angle [0,1,...]
      mmdb::realtype angle;      // angle in degrees
      mmdb::realtype angle_esd;  // angle esd in degrees; set to 0.0
                                 // if not provided

  };

}  // namespace ccp4srs

#endif // CCP4SRS_ANGLE_H
