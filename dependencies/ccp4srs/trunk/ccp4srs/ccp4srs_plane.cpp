//  $Id: ccp4srs_plane.cpp $
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
//  **** Module  :  ccp4srs_plane  <implementation>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Plane  - plain group of atoms
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#include <string.h>

#include "ccp4srs_plane.h"
#include "ccp4srs_defs.h"

namespace ccp4srs  {

  Plane::Plane()  {
    plane_id[0] = char(0);
    n_atoms     = 0;
    plane_atoms = NULL;
    dist_esd    = NULL;
  }

  Plane::~Plane()  {
    empty();
  }

  void Plane::empty()  {
    mmdb::FreeVectorMemory ( plane_atoms,0 );
    mmdb::FreeVectorMemory ( dist_esd   ,0 );
    n_atoms = 0;
  }


  void Plane::makeCIFTags ( mmdb::mmcif::PLoop mmCIFLoop )  {
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_COMP_ID  );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_PLANE_ID );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_ATOM_ID  );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_DIST_ESD );
  }

  int Plane::readFromCIF ( mmdb::mmcif::PLoop mmCIFLoop,
                          int planeNo, Container<Atom> & atoms )  {
  PlaneID planeId;
  int     i,j,k,n,i1,loop_len;
  int     rc;

    empty();

    rc = 0;
    loop_len = mmCIFLoop->GetLoopLength();
    if (loop_len<=0)
      return Plane::EmptyLoop;

    i1 = -1;
    n  = 0;
    for (i=0;(i<loop_len) && (n<=planeNo);i++)  {
      mmCIFLoop->CopyString ( plane_id,sizeof(PlaneID)-1,
                              MMCIF_TAG_PLANE_ID,i,rc );
      k = -1;
      for (j=i-1;(j>=0) && (k<0);j--)  {
        mmCIFLoop->CopyString ( planeId,sizeof(PlaneID)-1,
                                MMCIF_TAG_PLANE_ID,j,rc );
        if (!strcmp(plane_id,planeId))
          k = j;
      }
      if (k<0)  {
        n++;
        i1 = i;
      }
    }

    if (n<=planeNo)
      return Plane::NoMorePlanes;

    mmdb::GetVectorMemory ( plane_atoms,loop_len-i1,0 );
    mmdb::GetVectorMemory ( dist_esd   ,loop_len-i1,0 );

    for (i=i1;i<loop_len;i++)  {
      mmCIFLoop->CopyString ( planeId,sizeof(PlaneID)-1,
                              MMCIF_TAG_PLANE_ID,i,rc );
      if (!strcmp(planeId,plane_id))  {
        get_atom_from_cif ( plane_atoms[n_atoms],mmCIFLoop,
                            MMCIF_TAG_ATOM_ID,i,atoms,rc );
        mmCIFLoop->CopyReal ( dist_esd[n_atoms],MMCIF_TAG_DIST_ESD,i,rc );
        n_atoms++;
      }
    }

    return rc;

  }

  /*
  loop_
  _chem_comp_plane_atom.comp_id
  _chem_comp_plane_atom.plane_id
  _chem_comp_plane_atom.atom_id
  _chem_comp_plane_atom.dist_esd
   ATP      plan-1    N9        0.020
  */

  void Plane::writeToCIF ( mmdb::mmcif::PLoop mmCIFLoop,
                           mmdb::cpstr monID,
                           Container<Atom> & atoms )  {
  int i;
    for (i=0;i<n_atoms;i++)  {
      mmCIFLoop->AddString ( monID,false     );
      mmCIFLoop->AddString ( plane_id,false  );
      add_atom_to_cif      ( mmCIFLoop,plane_atoms[i],atoms );
      mmCIFLoop->AddReal   ( dist_esd[i],"%10.3f" );
    }
  }


  void Plane::copy ( PPlane plane, mmdb::ivector anmatch )  {
  int i;

    empty();

    strcpy ( plane_id,plane->id() );
    n_atoms = plane->size();
    if (n_atoms>0)  {
      mmdb::GetVectorMemory ( plane_atoms,n_atoms,0 );
      mmdb::GetVectorMemory ( dist_esd   ,n_atoms,0 );
      if (anmatch)  {
        for (i=0;i<n_atoms;i++)  {
          plane_atoms[i] = anmatch[plane->atom(i)];
          dist_esd   [i] = plane->esd  ( i );
        }
      } else  {
        for (i=0;i<n_atoms;i++)  {
          plane_atoms[i] = plane->atom ( i );
          dist_esd   [i] = plane->esd  ( i );
        }
      }
    }

  }

  void  Plane::write_mem ( PMemIO memIO, int version )  {
  UNUSED_ARGUMENT(version);
  int  i;
    memIO->put_line   ( plane_id );
    memIO->put_ishort ( n_atoms  );
    for (i=0;i<n_atoms;i++)  {
      memIO->put_integer ( plane_atoms[i] );
      memIO->put_float   ( dist_esd   [i] );
    }
  }

  bool Plane::read_mem ( PMemIO memIO, int version, bool * Ok )  {
  UNUSED_ARGUMENT(version);
  int     i;
  bool success;
    empty();
    if (Ok)  success = *Ok;
       else  success = true;
    memIO->get_line   ( plane_id,&success );
    memIO->get_ishort ( n_atoms ,&success );
    if (n_atoms>0)  {
      mmdb::GetVectorMemory ( plane_atoms,n_atoms,0 );
      mmdb::GetVectorMemory ( dist_esd   ,n_atoms,0 );
      for (i=0;i<n_atoms;i++)  {
        memIO->get_integer ( plane_atoms[i],&success );
        memIO->get_float   ( dist_esd   [i],&success );
      }
    }
    if (Ok)  *Ok = success;
    return success;
  }

}  // namespace ccp4srs
