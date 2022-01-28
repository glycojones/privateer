//  $Id: ccp4srs_atom.h $
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
//  **** Module  :  ccp4srs_atom  <interface>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Atom  - atom description class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2014
//
//  =================================================================
//

#ifndef CCP4SRS_ATOM_H
#define CCP4SRS_ATOM_H

#include "ccp4srs_container.h"
#include "memio_.h"

#include "mmdb2/mmdb_atom.h"
#include "mmdb2/mmdb_mmcif_.h"
#include "mmdb2/mmdb_defs.h"

namespace ccp4srs  {

  DefineClass(Atom);

  class Atom  {

    public:
      Atom();
      virtual ~Atom();

      // function for searches in containers
      inline mmdb::cpstr      id         () { return atomName;     }

      inline mmdb::cpstr      name        () { return atomName;     }
      inline mmdb::cpstr      rcsb_name   () { return rcsbAtomName; }
      inline mmdb::cpstr      old_name    () { return oldAtomName;  }
      inline mmdb::cpstr      element     () { return chemElement;  }
      inline mmdb::cpstr      energy_type () { return energyType;   }
      mmdb::realtype x();
      mmdb::realtype y();
      mmdb::realtype z();
      inline mmdb::realtype   x_rcsb_cartn() { return rcsb_cartn_x; }
      inline mmdb::realtype   y_rcsb_cartn() { return rcsb_cartn_y; }
      inline mmdb::realtype   z_rcsb_cartn() { return rcsb_cartn_z; }
      inline mmdb::realtype   x_rcsb_ideal() { return rcsb_ideal_x; }
      inline mmdb::realtype   y_rcsb_ideal() { return rcsb_ideal_y; }
      inline mmdb::realtype   z_rcsb_ideal() { return rcsb_ideal_z; }
      inline mmdb::realtype   x_ccp4_mlib () { return ccp4_mlib_x;  }
      inline mmdb::realtype   y_ccp4_mlib () { return ccp4_mlib_y;  }
      inline mmdb::realtype   z_ccp4_mlib () { return ccp4_mlib_z;  }
      inline mmdb::realtype   charge      () { return atomCharge;   }
      inline mmdb::realtype   weight      () { return atomWeight;   }
      inline mmdb::realtype   vdw_radius  () { return vdwRadius;    }
      inline mmdb::realtype   vdwh_radius () { return vdwhRadius;   }
      inline mmdb::realtype   ion_radius  () { return ionRadius;    }
      inline int              valency     () { return atomValency;  }
      inline int              sp          () { return atomSP;       }

      inline char rcsb_chirality() const { return property[0];  }
      inline char ccp4_chirality() const { return property[1];  }
      char chirality() const; // ccp4's; if not set then rcsb's
      inline char leaving         () const { return property[2];  }
      inline bool isLeaving       () const { return (property[2]=='Y'); }
      inline char hb_type         () const { return property[3];  }
      inline bool ccp4_coordinates() const { return (property[4]=='Y'); }


      mmdb::cpstr name_pdb      ( mmdb::pstr aname );  // with significant spacing
      mmdb::cpstr rcsb_name_pdb ( mmdb::pstr aname );  // with significant spacing
      mmdb::cpstr old_name_pdb  ( mmdb::pstr aname );  // with significant spacing
      mmdb::cpstr element_pdb   ( mmdb::pstr elem  );  // with significant spacing

      void set_energy_type ( mmdb::cpstr etype );
      void set_old_name    ( mmdb::cpstr name  );
      inline void set_rcsb_chirality ( char c  )  { property[0] = c; }
      inline void set_chirality      ( char c  )  { property[1] = c; }
      inline void set_hb_type        ( char h  )  { property[3] = h; }
      void set_ccp4_coordinates      ( bool on );

      inline void set_vdw_radius  ( mmdb::realtype r )  { vdwRadius   = r; }
      inline void set_vdwh_radius ( mmdb::realtype r )  { vdwhRadius  = r; }
      inline void set_ion_radius  ( mmdb::realtype r )  { ionRadius   = r; }
      inline void set_valency     ( int v      )  { atomValency = v; }
      inline void set_sp          ( int s      )  { atomSP      = s; }

      void  makeAtom ( mmdb::RPAtom a );
      mmdb::PAtom makeAtom();

      static void makeCIFTags ( mmdb::mmcif::PLoop mmCIFLoop );
      int  readFromCIF_ccp4 ( mmdb::mmcif::PLoop mmCIFLoop, int atomNo );
      int  readFromCIF_rcsb ( mmdb::mmcif::PLoop mmCIFLoop, int atomNo );
      int  readFromCIF      ( mmdb::mmcif::PStruct mmCIFStruct );
      void writeToCIF       ( mmdb::mmcif::PLoop mmCIFLoop, mmdb::cpstr monID );

      void copy  ( PAtom atom );
      void merge ( PAtom atom );

      void write_mem ( PMemIO memIO, int version );
      bool read_mem  ( PMemIO memIO, int version, bool * Ok = NULL );

    protected:
      mmdb::AtomName    atomName;     //!< atom name
      mmdb::AtomName    rcsbAtomName; //!< alternative atom name
      mmdb::AtomName    oldAtomName;  //!< old atom name
      mmdb::Element     chemElement;  //!< chemical element name
      mmdb::EnergyType  energyType;   //!< energy type
      mmdb::realtype    ccp4_mlib_x ; //!< ccp4 monomer library
      mmdb::realtype    ccp4_mlib_y ; ///   cartesian
      mmdb::realtype    ccp4_mlib_z ; ///      coordinates
      mmdb::realtype    rcsb_cartn_x; //!< pdb
      mmdb::realtype    rcsb_cartn_y; ///   cartesian coordinates
      mmdb::realtype    rcsb_cartn_z; ///      coordinates
      mmdb::realtype    rcsb_ideal_x; //!< rcsb
      mmdb::realtype    rcsb_ideal_y; ///   idealised
      mmdb::realtype    rcsb_ideal_z; ///      coordinates
      mmdb::realtype    atomCharge;   //!< partial charge
      mmdb::realtype    atomWeight;   //!< atomic weight (in a.u.)
      mmdb::realtype    vdwRadius;    //!< Van-Der-Waals radius
      mmdb::realtype    vdwhRadius;   //!< Van-Der-Waals radius with hydrogen
      mmdb::realtype    ionRadius;    //!< ion radius
      int               atomValency;  //!< valency
      int               atomSP;       //!< sp-hybridization
      char              property[5];  //!< property vector:
                                /// [0]: RCSB chirality : 'R', 'S' or 'N'
                                /// [1]: MLib chirality : 'R', 'S' or 'N'
                                /// [2]: leaving atom: 'Y' or 'N'
                                /// [3]: hydrogen bond type:
                                ///                 'D' donor
                                ///                 'A' acceptor
                                ///                 'B' both
                                ///                 'H' hydrogen candidate
                                ///                 'N' neither
                                /// [4]: ccp4 coordinates:
                                ///    'Y'  - x,y,z return ccp4_mlib_*
                                ///    'N'  - x,y,z return rcsb_cartn_*

      void Init();

  };

  extern void get_atom_from_cif ( int & index, mmdb::mmcif::PLoop loop,
                                  mmdb::cpstr Tag, int row,
                                  Container<Atom> & atoms,
                                  int & rc );

  extern void get_atom_from_cif ( int & index, mmdb::mmcif::PStruct mmCIFStruct,
                                  mmdb::cpstr Tag,
                                  Container<Atom> & atoms,
                                  int & rc );

  extern void add_atom_to_cif  ( mmdb::mmcif::PLoop loop, int atom_no,
                                 Container<Atom> & atoms );

} // namespace ccp4srs


#endif // CCP4SRS_ATOM_H
