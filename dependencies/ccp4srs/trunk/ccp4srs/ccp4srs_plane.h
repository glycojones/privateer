//  $Id: ccp4srs_plane.h $
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
//  **** Module  :  ccp4srs_plane  <interface>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Plane  - plain group of atoms
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#ifndef CCP4SRS_PLANE_H
#define CCP4SRS_PLANE_H

#include "ccp4srs_atom.h"

namespace ccp4srs  {

  DefineClass(Plane);

  typedef char PlaneID[20];

  class Plane  {

    public:

      enum { Ok=0,EmptyLoop=-11111,NoMorePlanes=-11112 };

      Plane();
      virtual ~Plane();

      inline mmdb::cpstr    id   () { return plane_id;    }
      inline int            size () { return n_atoms;     }
      inline mmdb::ivector  atoms() { return plane_atoms; }
      inline mmdb::rvector  esds () { return dist_esd;    }
      inline int            atom ( int atomNo ) { return plane_atoms[atomNo]; }
      inline mmdb::realtype esd  ( int atomNo ) { return dist_esd   [atomNo]; }

      static void makeCIFTags ( mmdb::mmcif::PLoop mmCIFLoop );
      int  readFromCIF ( mmdb::mmcif::PLoop mmCIFLoop, int planeNo,
                         Container<Atom> & atoms );
      void writeToCIF  ( mmdb::mmcif::PLoop mmCIFLoop, mmdb::cpstr monID,
                         Container<Atom> & atoms );

      void write_mem ( PMemIO memIO, int version );
      bool read_mem  ( PMemIO memIO, int version, bool * Ok = NULL );

      void copy ( PPlane plane, mmdb::ivector anmatch=NULL );

    protected:
      PlaneID       plane_id;     //!< plane id
      int           n_atoms;      //!< plane size
      mmdb::ivector plane_atoms;  //!< vector of plain's atom names
      mmdb::rvector dist_esd;     //!< vector of distance esds from the plane

      void  empty();

  };

}  // namespace ccp4srs

#endif // CCP4SRS_PLANE_H
