//  $Id: ccp4srs_bond.cpp $
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
//  **** Module  :  ccp4srs_bond  <implementation>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Bond  - atom description class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#include <string.h>

#include "ccp4srs_bond.h"
#include "ccp4srs_types.h"
#include "ccp4srs_defs.h"

namespace ccp4srs  {

  Bond::Bond()  {
    atom_1          = -1;
    atom_2          = -1;
    bond_order      = Bond::noOrder;
  //  bond_type       = Bond::noType;
    bond_length     = 0.0;
    bond_length_esd = 0.0;
  }

  Bond::~Bond()  {
  }


  void  get_bond_type_from_cif_ccp4 ( int & type, mmdb::mmcif::PLoop loop,
                                      mmdb::cpstr Tag, int row, int & rc )  {
  mmdb::pstr p;
     p = loop->GetString ( Tag,row,rc );
     if (!p)                         type = Bond::noOrder;
     else if (!strcmp(p,"single"))   type = Bond::Single;
     else if (!strcmp(p,"aromatic")) type = Bond::Aromatic;
     else if (!strcmp(p,"aromat"))   type = Bond::Aromatic;
     else if (!strcmp(p,"double"))   type = Bond::Double;
     else if (!strcmp(p,"triple"))   type = Bond::Triple;
     else if (!strcmp(p,"deloc"))    type = Bond::Deloc;
     else if (!strcmp(p,"coval"))    type = Bond::Covalent;
     else if (!strcmp(p,"metal"))    type = Bond::Metal;
     else  {
       type = Bond::noOrder;
       rc   = CCP4SRS_WrongBondType;
     }
  }

  void  get_bond_type_from_cif_rcsb ( int & type,
                                      mmdb::mmcif::PLoop loop,
                                      mmdb::cpstr Tag,
                                      int row, int & rc )  {
  mmdb::pstr p;
     p = loop->GetString ( Tag,row,rc );
     if (!p)                       type = Bond::noOrder;
     else if (!strcmp(p,"SING"))   type = Bond::Single;
     else if (!strcmp(p,"AROM"))   type = Bond::Aromatic;
     else if (!strcmp(p,"DOUB"))   type = Bond::Double;
     else if (!strcmp(p,"TRIP"))   type = Bond::Triple;
     else  {
       type = Bond::noOrder;
       rc   = CCP4SRS_WrongBondType;
     }
  }

  void  get_bond_type_from_cif ( int & type,
                                 mmdb::mmcif::PStruct mmCIFStruct,
                                 mmdb::cpstr Tag, int & rc )  {
  mmdb::pstr p;
     p = mmCIFStruct->GetString ( Tag,rc );
     if (!p)                       type = Bond::noOrder;
     else if (!strcmp(p,"SING"))   type = Bond::Single;
     else if (!strcmp(p,"AROM"))   type = Bond::Aromatic;
     else if (!strcmp(p,"DOUB"))   type = Bond::Double;
     else if (!strcmp(p,"TRIP"))   type = Bond::Triple;
     else  {
       type = Bond::noOrder;
       rc   = CCP4SRS_WrongBondType;
     }
  }

  void add_bond_type_to_cif ( mmdb::mmcif::PLoop loop, int type )  {
    switch (type)  {
      default:
      case Bond::noOrder:  loop->AddNoData ( mmdb::mmcif::CIF_NODATA_DOT );
                                break;
      case Bond::Single:   loop->AddString ( "single"  ,false );
                                break;
      case Bond::Aromatic: loop->AddString ( "aromatic",false );
                                break;
      case Bond::Double:   loop->AddString ( "double"  ,false );
                                break;
      case Bond::Triple:   loop->AddString ( "triple"  ,false );
                                break;
      case Bond::Deloc:    loop->AddString ( "deloc"   ,false );
    }
  }


  int Bond::readFromCIF_ccp4 ( mmdb::mmcif::PLoop mmCIFLoop,
                               int bondNo,
                               Container<Atom> & atoms )  {
  int rc;

    rc = 0;

    get_atom_from_cif ( atom_1,mmCIFLoop,MMCIF_TAG_ATOM_ID_1,
                        bondNo,atoms,rc );
    get_atom_from_cif ( atom_2,mmCIFLoop,MMCIF_TAG_ATOM_ID_2,
                        bondNo,atoms,rc );
    get_bond_type_from_cif_ccp4 ( bond_order,mmCIFLoop,MMCIF_TAG_TYPE,
                        bondNo,rc );

    if (mmCIFLoop->GetTagNo(MMCIF_TAG_VALUE_DIST)>=0)  {
      mmCIFLoop->CopyReal ( bond_length    ,MMCIF_TAG_VALUE_DIST    ,
                            bondNo,rc );
      mmCIFLoop->CopyReal ( bond_length_esd,MMCIF_TAG_VALUE_DIST_ESD,
                            bondNo,rc );
    } else  {
      bond_length     = 0.0;
      bond_length_esd = 0.0;
    }

    return rc;

  }


  int Bond::readFromCIF_rcsb ( mmdb::mmcif::PLoop mmCIFLoop,
                               int bondNo,
                               Container<Atom> & atoms )  {
  mmdb::pstr p;
  int        rc;

    rc = 0;

    get_atom_from_cif ( atom_1,mmCIFLoop,MMCIF_TAG_ATOM_ID_1,
                        bondNo,atoms,rc );
    get_atom_from_cif ( atom_2,mmCIFLoop,MMCIF_TAG_ATOM_ID_2,
                        bondNo,atoms,rc );
    get_bond_type_from_cif_rcsb ( bond_order,mmCIFLoop,
                        MMCIF_TAG_VALUE_ORDER,bondNo,rc );

    if ((bond_order==Bond::Single) ||
        (bond_order==Bond::Double))  {
      p = mmCIFLoop->GetString ( MMCIF_TAG_PDBX_AROMATIC_FLAG,bondNo,rc );
      if (p)  {
        if (p[0]=='Y')
          bond_order = Bond::Aromatic;
      }
    }

    bond_length     = 0.0;
    bond_length_esd = 0.0;

    return rc;

  }

  int Bond::readFromCIF ( mmdb::mmcif::PStruct mmCIFStruct,
                                 Container<Atom> & atoms ) {
  mmdb::pstr p;
  int        rc;

    rc = 0;

    get_atom_from_cif ( atom_1,mmCIFStruct,MMCIF_TAG_ATOM_ID_1,
                        atoms,rc );
    get_atom_from_cif ( atom_2,mmCIFStruct,MMCIF_TAG_ATOM_ID_2,
                        atoms,rc );
    get_bond_type_from_cif ( bond_order,mmCIFStruct,
                        MMCIF_TAG_VALUE_ORDER,rc );

    if ((bond_order==Bond::Single) ||
        (bond_order==Bond::Double))  {
      p = mmCIFStruct->GetString ( MMCIF_TAG_PDBX_AROMATIC_FLAG,rc );
      if (p)  {
        if (p[0]=='Y')
          bond_order = Bond::Aromatic;
      }
    }

    bond_length     = 0.0;
    bond_length_esd = 0.0;

    return rc;

  }

  /*
  loop_
  _chem_comp_bond.comp_id
  _chem_comp_bond.atom_id_1
  _chem_comp_bond.atom_id_2
  _chem_comp_bond.type
  _chem_comp_bond.value_dist
  _chem_comp_bond.value_dist_esd
   ATP      O1G    PG        deloc       1.510    0.020
  */


  bool Bond::compare ( Container<Atom> & atoms, PBond bond,
                       Container<Atom> & bond_atoms )  {
  mmdb::AtomName a11,a12, a21,a22;

    strcpy ( a11,atoms.at(atom_1)->id() );
    strcpy ( a12,atoms.at(atom_2)->id() );

    strcpy ( a21,bond_atoms.at(bond->atom1())->id() );
    strcpy ( a22,bond_atoms.at(bond->atom2())->id() );

    if (!strcmp(a11,a21))
      return (strcmp(a12,a22)==0);
    else if (!strcmp(a11,a22))
      return (strcmp(a12,a21)==0);

    return false;

  }

  bool Bond::copy ( Container<Atom> & atoms,
                         PBond bond,
                         Container<Atom> & bond_atoms )  {
  mmdb::AtomName a1,a2;
  int      i;

    strcpy ( a1,bond_atoms.at(bond->atom1())->id() );
    strcpy ( a2,bond_atoms.at(bond->atom2())->id() );

    atom_1 = -1;
    atom_2 = -1;

    for (i=0;(i<atoms.numberOf());i++)  {
      if (!strcmp(atoms.at(i)->id(),a1))  atom_1 = i;
      if (!strcmp(atoms.at(i)->id(),a2))  atom_2 = i;
    }

    if ((atom_1>=0) && (atom_2>=0))  {
      bond_order      = bond->order ();
      bond_length     = bond->length();
      bond_length_esd = bond->length_esd();
      return true;
    }

    return false;

  }

  void Bond::copy ( PBond bond, mmdb::ivector anmatch )  {
    if (anmatch)  {
      atom_1 = anmatch[bond->atom1()];
      atom_2 = anmatch[bond->atom2()];
    } else  {
      atom_1 = bond->atom1();
      atom_2 = bond->atom2();
    }
    bond_order      = bond->order     ();
    bond_length     = bond->length    ();
    bond_length_esd = bond->length_esd();
  }


  void Bond::makeCIFTags ( mmdb::mmcif::PLoop mmCIFLoop )  {
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_COMP_ID        );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_ATOM_ID_1      );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_ATOM_ID_2      );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_TYPE           );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_VALUE_DIST     );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_VALUE_DIST_ESD );
  }

  void Bond::writeToCIF ( mmdb::mmcif::PLoop mmCIFLoop,
                              mmdb::cpstr monID,
                              Container<Atom> & atoms )  {
    mmCIFLoop->AddString ( monID,false              );
    add_atom_to_cif      ( mmCIFLoop,atom_1,atoms   );
    add_atom_to_cif      ( mmCIFLoop,atom_2,atoms   );
    add_bond_type_to_cif ( mmCIFLoop,bond_order     );
    mmCIFLoop->AddReal   ( bond_length    ,"%10.3f" );
    mmCIFLoop->AddReal   ( bond_length_esd,"%10.3f" );
  }

  void Bond::write_mem ( PMemIO memIO, int version )  {
  UNUSED_ARGUMENT(version);
    memIO->put_integer ( atom_1          );
    memIO->put_integer ( atom_2          );
    memIO->put_ibyte   ( bond_order      );
  //  memIO->put_ibyte   ( bond_type       );
    memIO->put_float   ( bond_length     );
    memIO->put_float   ( bond_length_esd );
  }

  bool Bond::read_mem ( PMemIO memIO, int version, bool * Ok )  {
  UNUSED_ARGUMENT(version);
  bool success;
    if (Ok)  success = *Ok;
       else  success = true;
    memIO->get_integer ( atom_1         ,&success );
    memIO->get_integer ( atom_2         ,&success );
    memIO->get_ibyte   ( bond_order     ,&success );
  //  memIO->get_ibyte   ( bond_type      ,&success );
    memIO->get_float   ( bond_length    ,&success );
    memIO->get_float   ( bond_length_esd,&success );
    if (Ok)  *Ok = success;
    return success;
  }

}  // namespace ccp4srs
