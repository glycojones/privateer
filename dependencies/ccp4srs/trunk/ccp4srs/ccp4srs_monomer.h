//  $Id: ccp4srs_monomer.h $
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
//  **** Module  :  ccp4srs_monomer  <interface>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Monomer  - monomer structure
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2014
//
//  =================================================================
//

#ifndef CCP4SRS_MONOMER_H
#define CCP4SRS_MONOMER_H

#include "mmdb2/mmdb_mmcif_.h"
#include "mmdb2/mmdb_math_graph.h"

#include "ccp4srs_entry.h"
#include "ccp4srs_bond.h"
#include "ccp4srs_angle.h"
#include "ccp4srs_torsion.h"
#include "ccp4srs_chicenter.h"
#include "ccp4srs_plane.h"
#include "ccp4srs_tree.h"
#include "ccp4srs_smiles.h"
#include "ccp4srs_types.h"

namespace ccp4srs  {

  DefineClass(Monomer);

  class Monomer : public Entry  {

    public:
      Monomer();
      ~Monomer();

      virtual int type()  { return Entry::Monomer; }

      void   reset();  // deletes all data

      inline mmdb::cpstr ID   ()  { return monID;    }
      inline mmdb::cpstr oldID()  { return oldMonID; }
      void   setOldID   ( mmdb::cpstr oldID );

      mmdb::cpstr  chem_name   ();
      mmdb::cpstr  chem_type   ();
      mmdb::cpstr  chem_formula();

      inline int n_atoms     () { return atoms     .numberOf(); }
      inline int n_bonds     () { return bonds     .numberOf(); }
      inline int n_angles    () { return angles    .numberOf(); }
      inline int n_torsions  () { return torsions  .numberOf(); }
      inline int n_chicenters() { return chicenters.numberOf(); }
      inline int n_planes    () { return planes    .numberOf(); }

      inline PAtom  atom  ( int n ) { return atoms   .at(n); }
      inline PBond  bond  ( int n ) { return bonds   .at(n); }
      inline PAngle angle ( int n ) { return angles  .at(n); }
      inline PPlane plane ( int n ) { return planes  .at(n); }
      inline PTree  tree  ()        { return &struct_tree;   }
      inline PTorsion torsion ( int n )
                                           { return torsions.at(n); }
      inline PChiCenter chicenter ( int n )
                                         { return chicenters.at(n); }

      inline void add ( PAtom  atom    ) { atoms .add(atom ); }
      inline void add ( PBond  bond    ) { bonds .add(bond ); }
      inline void add ( PAngle angle   ) { angles.add(angle); }
      inline void add ( PPlane plane   ) { planes.add(plane); }
      inline void add ( PTree  stree   ) { struct_tree.copy(stree); }
      inline void add ( PTorsion torsion )
                                            { torsions.add(torsion); }
      inline void add ( PChiCenter chicenter )
                                        { chicenters.add(chicenter); }

      inline Container<Atom>    *get_atoms   () { return &atoms;      }
      inline Container<Bond>    *get_bonds   () { return &bonds;      }
      inline Container<Angle>   *get_angles  () { return &angles;     }
      inline Container<Torsion> *get_torsions() { return &torsions;   }
      inline Container<ChiCenter> *get_chicenters() { return &chicenters; }
      inline Container<Plane>   *get_planes  () { return &planes;     }
      inline Tree *get_tree() { return &struct_tree; }

      inline Smiles *getACDLabs()  { return &ACDLabs; }
      inline Smiles *getCACTVS ()  { return &CACTVS;  }
      inline Smiles *getOpenEye()  { return &OpenEye; }
      inline Smiles *getInChI  ()  { return &InChI;   }

      int atom_no   ( mmdb::cpstr aname );
      int get_bound ( int  atomNo, mmdb::cpstr elem, int serNo );
      int get_bound ( int  atomNo, int n );
      int get_bound ( mmdb::cpstr aname, int n );

      CCP4SRS_RC addHydrogens ( mmdb::PResidue R );

      //  getAtomNameMatch(..) returns anmatch[i], i=0..nAtoms-1, equal
      // to j such that name(Atom[i])==name(A[j]). Note that atom names
      // are similarly aligned and space-padded in MMDB, but not in SRS.
      // However, the proper alignment is done on fly. If ith atom in
      // the structue is not found in A, anmatch[i] is set -1.
      //   If array A contains atoms in different alternative
      // conformations, the the value of altLoc is interpreted as
      // follows:
      //    NULL  - the highest occupancy atom will be taken
      //            if all occupancies are equal then atom with
      //            first altLoc taken
      //    other - atoms with given altLoc are taken. If such
      //            altLoc is not found, the function does as if
      //            NULL value for altLoc is given.
      //   A clean PDB file is anticipated, so that atoms with
      // alternative conformations are grouped together.
      //   It is Ok to have NULL pointers in A.
      void getAtomNameMatch ( mmdb::PPAtom A, int nat, mmdb::pstr altLoc,
                              mmdb::ivector anmatch );

      /// \brief Matches atoms in given container with atoms in the
      ///   monomer by atom names.
      ///  \param[in] A pointer to atom container
      ///  \param[out] anmatch allocated integer vector
      ///         [0..A->numberOf()-1], where the match will be returned:
      /// \arg \b anmatch[i]>=0 : atom A->at(i) matches anmathch[i]-th
      ///         atom in the monomer
      /// \arg \b anmatch[i]<0  : atom A->at(i) does not have a name match
      ///         in the monomer
      void getAtomNameMatch ( Container<Atom> *A, mmdb::ivector anmatch );

      mmdb::math::PGraph getGraph ( int *retCode );

      void getLeavingAtoms ( int & nLeavingAtoms,
                             mmdb::ivector & leavingAtom,
                             mmdb::ivector & bondedAtom );

      // Return code:
      //     CCP4SRS_Ok    success
      int readFromCIF ( mmdb::mmcif::PData mmCIFData,
                        mmdb::pstr error_desc = NULL );

      mmdb::mmcif::PData makeCIF();  // return to be deleted by application

      // service functions
      void merge ( PMonomer monomer );
      void check_ccp4_coordinates();

    protected:
      mmdb::ResName        monID;
      mmdb::ResName        oldMonID;
      mmdb::pstr           monName;
      mmdb::pstr           monType;
      mmdb::pstr           monFormula;
      Container<Atom>      atoms;
      Container<Bond>      bonds;
      Container<Angle>     angles;
      Container<Torsion>   torsions;
      Container<ChiCenter> chicenters;
      Container<Plane>     planes;
      Tree                 struct_tree;
      Smiles               ACDLabs;
      Smiles               CACTVS;
      Smiles               OpenEye;
      Smiles               InChI;

      void  Init();

      void  readMonIDFromCIF ( mmdb::mmcif::PLoop   loop );
      void  readMonIDFromCIF ( mmdb::mmcif::PStruct mmCIFStruct );
      int   readFromCIF_ccp4 ( mmdb::mmcif::PData   mmCIFData,
                               mmdb::pstr error_desc );
      int   readFromCIF_rcsb ( mmdb::mmcif::PData   mmCIFData,
                               mmdb::pstr error_desc );

      virtual void write_mem ( PMemIO memIO, int version );
      virtual bool read_mem  ( PMemIO memIO, int version,
                               bool * Ok = NULL );

  };

}  // namespace ccp4srs

#endif // CCP4SRS_MONOMER_H
