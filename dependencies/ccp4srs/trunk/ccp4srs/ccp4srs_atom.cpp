//  $Id: ccp4srs_atom.cpp $
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
//    26.06.14   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  ccp4srs_atom  <implementation>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Atom  - atom description class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2014
//
//  =================================================================
//

#include <string.h>

#include "mmdb2/mmdb_tables.h"
#include "ccp4srs_atom.h"
#include "ccp4srs_types.h"
#include "ccp4srs_defs.h"

namespace ccp4srs  {

  Atom::Atom()  {
    Init();
  }

  Atom::~Atom()  {}

  void Atom::Init()  {

    atomName    [0] = char(0);       // atom name
    rcsbAtomName[0] = char(0);       // alternative atom name
    oldAtomName [0] = char(0);       // old atom name
    chemElement [0] = char(0);       // chemical element name
    energyType  [0] = char(0);       // energy type
    rcsb_cartn_x    = -mmdb::MaxShortReal/2.0; // cartesian coordinates
    rcsb_cartn_y    = rcsb_cartn_x;
    rcsb_cartn_z    = rcsb_cartn_x;
    rcsb_ideal_x    = rcsb_cartn_x;  // idealized coordinates
    rcsb_ideal_y    = rcsb_cartn_x;
    rcsb_ideal_z    = rcsb_cartn_x;
    ccp4_mlib_x     = rcsb_cartn_x;  // monomer library coordinates
    ccp4_mlib_y     = rcsb_cartn_x;
    ccp4_mlib_z     = rcsb_cartn_x;
    atomCharge      = 0.0;           // partial charge
    atomWeight      = 0.0;           // atomic weight (in a.u.)
    vdwRadius       = 0.0;           // Van-Der-Waals radius
    vdwhRadius      = 0.0;           // Van-Der-Waals radius with hydrogen
    ionRadius       = 0.0;           // ion radius
    atomValency     = 0;             // valency
    atomSP          = 0;             // sp-hybridization
    property[0]     = '-';           // RCSB chirality : 'R', 'S' or 'N'
    property[1]     = '-';           // MLib chirality : 'R', 'S' or 'N'
    property[2]     = '-';           // leaving atom: 'Y' or 'N'
    property[3]     = '-';           // hydrogen bond type:
                                     //           'D' donor
                                     //           'A' acceptor
                                     //           'B' both
                                     //           'H' hydrogen candidate
                                     //           'N' neither
    property[4]     = '-';           // ccp4 coordinates:
                                     //    'Y'  - x,y,z return ccp4_mlib_*
                                     //    'N'  - x,y,z return rcsb_cartn_*

  }


  mmdb::cpstr Atom::name_pdb ( mmdb::pstr aname )  {
  int  l;

    l = strlen(atomName);
    if (l>=4)
      strcpy ( aname,atomName );
    else  {
      strcpy ( aname,"    " );
      if (mmdb::isMetal(chemElement))
            strncpy ( aname,atomName,l );
      else  strncpy ( &(aname[1]),atomName,l );
    }

    return aname;

  }

  mmdb::cpstr Atom::rcsb_name_pdb ( mmdb::pstr aname )  {
  int  l;

    l = strlen(rcsbAtomName);
    if (l>=4)
      return rcsbAtomName;

    strcpy ( aname,"    " );

    if (mmdb::isMetal(chemElement))
          strncpy ( aname,rcsbAtomName,l       );
    else  strncpy ( &(aname[1]),rcsbAtomName,l );

    return aname;

  }

  mmdb::cpstr Atom::old_name_pdb ( mmdb::pstr aname )  {
  int  l;

    l = strlen(oldAtomName);
    if (l>=4)
      return oldAtomName;

    strcpy ( aname,"    " );

    if (mmdb::isMetal(chemElement))
          strncpy ( aname,oldAtomName,l       );
    else  strncpy ( &(aname[1]),oldAtomName,l );

    return aname;

  }

  mmdb::cpstr Atom::element_pdb ( mmdb::pstr elem )  {

    if (!chemElement[0])  {
      elem[0] = ' ';
      elem[1] = ' ';
    } else if (!chemElement[1])  {
      elem[0] = ' ';
      elem[1] = chemElement[0];
    } else  {
      elem[0] = chemElement[0];
      elem[1] = chemElement[1];
    }
    elem[2] = char(0);

    return elem;

  }

  void Atom::set_energy_type ( mmdb::cpstr etype )  {
    strcpy ( energyType,etype );
  }

  void Atom::set_old_name ( mmdb::cpstr name )  {
    strcpy ( oldAtomName,name );
  }

  void Atom::set_ccp4_coordinates ( bool on )  {
    if (on)  property[4] = 'Y';
       else  property[4] = 'N';
  }

  mmdb::realtype Atom::x()  {
    if (property[4]=='Y')  return ccp4_mlib_x;
                     else  return rcsb_cartn_x;
  }

  mmdb::realtype Atom::y()  {
    if (property[4]=='Y')  return ccp4_mlib_y;
                     else  return rcsb_cartn_y;
  }

  mmdb::realtype Atom::z()  {
    if (property[4]=='Y')  return ccp4_mlib_z;
                     else  return rcsb_cartn_z;
  }


  char Atom::chirality() const  {
  // ccp4's; if not set then rcsb's
    if (property[1]!='-')  return property[1];
    return property[0];
  }

  void Atom::makeAtom ( mmdb::RPAtom a )  {
  mmdb::AtomName aname;
  mmdb::Element  elem;
    if (!a)  a = mmdb::newAtom();
    a->SetAtomName    ( 0,0,name_pdb(aname),"","",element_pdb(elem) );
    a->SetCharge      ( mmdb::RMax(0.0,atomCharge) );
    a->SetCoordinates ( x(),y(),z(),1.0,0.0  );
    strcpy ( a->energyType,energyType  );
    a->sigOcc = -1.0;  // signal that atom was added from SBase
  }

  mmdb::PAtom Atom::makeAtom()  {
  mmdb::PAtom   a;
  mmdb::AtomName aname;
  mmdb::Element  elem;
    a = mmdb::newAtom();
    a->SetAtomName    ( 0,0,name_pdb(aname),"","",element_pdb(elem) );
    a->SetCharge      ( mmdb::RMax(0.0,atomCharge) );
    a->SetCoordinates ( x(),y(),z(),1.0,0.0 );
    strcpy ( a->energyType,energyType );
    a->sigOcc = -1.0;  // signal that atom was added from SBase
    return a;
  }


  void Atom::makeCIFTags ( mmdb::mmcif::PLoop mmCIFLoop )  {
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_COMP_ID        );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_ATOM_ID        );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_TYPE_SYMBOL    );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_TYPE_ENERGY    );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_PARTIAL_CHARGE );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_X              );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_Y              );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_Z              );
  }

  int Atom::readFromCIF_ccp4 ( mmdb::mmcif::PLoop mmCIFLoop,
                                      int atomNo )  {
  int  rc;

    rc = 0;
    mmCIFLoop->CopyString ( atomName,sizeof(mmdb::AtomName)-1,
                            MMCIF_TAG_ATOM_ID,atomNo,rc );
    mmCIFLoop->CopyString ( chemElement,sizeof(mmdb::Element)-1,
                            MMCIF_TAG_TYPE_SYMBOL,atomNo,rc );
    mmCIFLoop->CopyString ( energyType,sizeof(mmdb::EnergyType)-1,
                            MMCIF_TAG_TYPE_ENERGY,atomNo,rc );

    mmCIFLoop->CopyReal   ( atomCharge,MMCIF_TAG_PARTIAL_CHARGE,
                            atomNo,rc );

    if (mmCIFLoop->GetTagNo(MMCIF_TAG_X)>=0)  {
      mmCIFLoop->CopyReal ( ccp4_mlib_x,MMCIF_TAG_X,atomNo,rc );
      mmCIFLoop->CopyReal ( ccp4_mlib_y,MMCIF_TAG_Y,atomNo,rc );
      mmCIFLoop->CopyReal ( ccp4_mlib_z,MMCIF_TAG_Z,atomNo,rc );
    } else if (mmCIFLoop->GetTagNo(MMCIF_TAG_MODEL_CARTN_X)>=0)  {
      mmCIFLoop->CopyReal ( ccp4_mlib_x,MMCIF_TAG_MODEL_CARTN_X,atomNo,rc );
      mmCIFLoop->CopyReal ( ccp4_mlib_y,MMCIF_TAG_MODEL_CARTN_Y,atomNo,rc );
      mmCIFLoop->CopyReal ( ccp4_mlib_z,MMCIF_TAG_MODEL_CARTN_Z,atomNo,rc );
    }

    return rc;

  }

  int Atom::readFromCIF_rcsb ( mmdb::mmcif::PLoop mmCIFLoop,
                                      int atomNo )  {
  mmdb::pstr p;
  int        rc;

    rc = 0;
    mmCIFLoop->CopyString  ( atomName,sizeof(mmdb::AtomName)-1,
                             MMCIF_TAG_ATOM_ID,atomNo,rc );
    mmCIFLoop->CopyString  ( rcsbAtomName,sizeof(mmdb::AtomName)-1,
                             MMCIF_TAG_ALT_ATOM_ID,atomNo,rc );
    mmCIFLoop->CopyString  ( chemElement,sizeof(mmdb::Element)-1,
                             MMCIF_TAG_TYPE_SYMBOL,atomNo,rc );

    if (mmCIFLoop->GetTagNo(MMCIF_TAG_MODEL_CARTN_X)>=0)  {
      mmCIFLoop->CopyReal ( rcsb_cartn_x,MMCIF_TAG_MODEL_CARTN_X,atomNo,rc );
      mmCIFLoop->CopyReal ( rcsb_cartn_y,MMCIF_TAG_MODEL_CARTN_Y,atomNo,rc );
      mmCIFLoop->CopyReal ( rcsb_cartn_z,MMCIF_TAG_MODEL_CARTN_Z,atomNo,rc );
    }
    if (mmCIFLoop->GetTagNo(MMCIF_TAG_PDBX_MODEL_CARTN_X)>=0)  {
      mmCIFLoop->CopyReal ( rcsb_ideal_x,MMCIF_TAG_PDBX_MODEL_CARTN_X,atomNo,rc );
      mmCIFLoop->CopyReal ( rcsb_ideal_y,MMCIF_TAG_PDBX_MODEL_CARTN_Y,atomNo,rc );
      mmCIFLoop->CopyReal ( rcsb_ideal_z,MMCIF_TAG_PDBX_MODEL_CARTN_Z,atomNo,rc );
    }

    p = mmCIFLoop->GetString ( MMCIF_TAG_PDBX_STEREO_CONFIG,atomNo,rc );
    if (!p)  property[0] = '-';
       else  property[0] = p[0];


    p = mmCIFLoop->GetString ( MMCIF_TAG_PDBX_LEAVING_ATOM_FLAG,atomNo,rc );
    if (!p)              property[2] = 'N';
    else if (p[0]=='Y')  property[2] = 'Y';
                   else  property[2] = 'N';

    return rc;

  }

  int Atom::readFromCIF ( mmdb::mmcif::PStruct mmCIFStruct )  {
  mmdb::pstr p;
  int        rc;

    rc = 0;
    p = mmCIFStruct->GetString ( MMCIF_TAG_ATOM_ID,rc );
    if (p) strcpy ( atomName,p );
    p = mmCIFStruct->GetString ( MMCIF_TAG_ALT_ATOM_ID,rc );
    if (p) strcpy ( rcsbAtomName,p );
    p = mmCIFStruct->GetString ( MMCIF_TAG_TYPE_SYMBOL,rc );
    if (p) strcpy ( chemElement,p );
    p = mmCIFStruct->GetString ( MMCIF_TAG_TYPE_ENERGY,rc );
    if (p) strcpy ( energyType,p );

    rc = mmCIFStruct->GetReal ( atomCharge,MMCIF_TAG_PARTIAL_CHARGE );

    if (mmCIFStruct->GetTagNo(MMCIF_TAG_X)>=0)  {
      rc = mmCIFStruct->GetReal ( ccp4_mlib_x,MMCIF_TAG_X );
      rc = mmCIFStruct->GetReal ( ccp4_mlib_y,MMCIF_TAG_Y );
      rc = mmCIFStruct->GetReal ( ccp4_mlib_z,MMCIF_TAG_Y );
    } else if (mmCIFStruct->GetTagNo(MMCIF_TAG_MODEL_CARTN_X)>=0)  {
      rc = mmCIFStruct->GetReal ( rcsb_cartn_x,MMCIF_TAG_MODEL_CARTN_X );
      rc = mmCIFStruct->GetReal ( rcsb_cartn_y,MMCIF_TAG_MODEL_CARTN_Y );
      rc = mmCIFStruct->GetReal ( rcsb_cartn_z,MMCIF_TAG_MODEL_CARTN_Y );
      rc = mmCIFStruct->GetReal ( rcsb_ideal_x,MMCIF_TAG_PDBX_MODEL_CARTN_X );
      rc = mmCIFStruct->GetReal ( rcsb_ideal_y,MMCIF_TAG_PDBX_MODEL_CARTN_Y );
      rc = mmCIFStruct->GetReal ( rcsb_ideal_z,MMCIF_TAG_PDBX_MODEL_CARTN_Y );
      p  = mmCIFStruct->GetString ( MMCIF_TAG_PDBX_LEAVING_ATOM_FLAG,rc );
      if (!p)              property[2] = 'N';
      else if (p[0]=='Y')  property[2] = 'Y';
                     else  property[2] = 'N';
    }

    return rc;

  }

  void Atom::writeToCIF  ( mmdb::mmcif::PLoop mmCIFLoop, mmdb::cpstr monID )  {

    if (ccp4_mlib_x>-mmdb::MaxShortReal/2.1)  {

      mmCIFLoop->AddString ( monID      ,false );

      mmCIFLoop->AddString ( atomName   ,false );
      mmCIFLoop->AddString ( chemElement,false );
      mmCIFLoop->AddString ( energyType ,false );

      mmCIFLoop->AddReal   ( atomCharge ,"%10.3f" );
      mmCIFLoop->AddReal   ( ccp4_mlib_x,"%10.3f" );
      mmCIFLoop->AddReal   ( ccp4_mlib_y,"%10.3f" );
      mmCIFLoop->AddReal   ( ccp4_mlib_z,"%10.3f" );

    }

  }

  /*
  data_comp_ATP
  #
  loop_
  _chem_comp_atom.comp_id
  _chem_comp_atom.atom_id
  _chem_comp_atom.type_symbol
  _chem_comp_atom.type_energy
  _chem_comp_atom.partial_charge
  _chem_comp_atom.x
  _chem_comp_atom.y
  _chem_comp_atom.z
   ATP           O2A    O    OP       -0.500      0.000    0.000    0.000
  */

  void Atom::copy ( PAtom atom )  {

    strcpy ( atomName    ,atom->name       () );
    strcpy ( rcsbAtomName,atom->rcsb_name  () );
    strcpy ( oldAtomName ,atom->old_name   () );
    strcpy ( chemElement ,atom->element    () );
    strcpy ( energyType  ,atom->energy_type() );

    rcsb_cartn_x = atom->x_rcsb_cartn();
    rcsb_cartn_y = atom->y_rcsb_cartn();
    rcsb_cartn_z = atom->z_rcsb_cartn();
    rcsb_ideal_x = atom->x_rcsb_ideal();
    rcsb_ideal_y = atom->y_rcsb_ideal();
    rcsb_ideal_z = atom->z_rcsb_ideal();
    ccp4_mlib_x  = atom->x_ccp4_mlib ();
    ccp4_mlib_y  = atom->y_ccp4_mlib ();
    ccp4_mlib_z  = atom->z_ccp4_mlib ();

    atomCharge  = atom->charge        ();
    atomWeight  = atom->weight        ();
    vdwRadius   = atom->vdw_radius    ();
    vdwhRadius  = atom->vdwh_radius   ();
    ionRadius   = atom->ion_radius    ();
    atomValency = atom->valency       ();
    atomSP      = atom->sp            ();
    property[0] = atom->rcsb_chirality();
    property[1] = atom->chirality     ();
    property[2] = atom->leaving       ();
    property[3] = atom->hb_type       ();

  }


  void Atom::merge ( PAtom atom )  {

    if (!atomName    [0])  strcpy ( atomName    ,atom->name       () );
    if (!rcsbAtomName[0])  strcpy ( rcsbAtomName,atom->rcsb_name  () );
    if (!oldAtomName [0])  strcpy ( oldAtomName ,atom->old_name   () );
    if (!chemElement [0])  strcpy ( chemElement ,atom->element    () );
    if (!energyType  [0])  strcpy ( energyType  ,atom->energy_type() );

    if (rcsb_cartn_x<=-mmdb::MaxShortReal/2.0)  {
      rcsb_cartn_x = atom->x_rcsb_cartn();
      rcsb_cartn_y = atom->y_rcsb_cartn();
      rcsb_cartn_z = atom->z_rcsb_cartn();
    }
    if (rcsb_ideal_x<=-mmdb::MaxShortReal/2.0)  {
      rcsb_ideal_x = atom->x_rcsb_ideal();
      rcsb_ideal_y = atom->y_rcsb_ideal();
      rcsb_ideal_z = atom->z_rcsb_ideal();
    }
    if (ccp4_mlib_x<=-mmdb::MaxShortReal/2.0)  {
      ccp4_mlib_x = atom->x_ccp4_mlib();
      ccp4_mlib_y = atom->y_ccp4_mlib();
      ccp4_mlib_z = atom->z_ccp4_mlib();
    }

    if (atomCharge==0.0)  atomCharge  = atom->charge        ();
    if (atomWeight==0.0)  atomWeight  = atom->weight        ();
    if (vdwRadius ==0.0)  vdwRadius   = atom->vdw_radius    ();
    if (vdwhRadius==0.0)  vdwhRadius  = atom->vdwh_radius   ();
    if (ionRadius ==0.0)  ionRadius   = atom->ion_radius    ();
    if (!atomValency)     atomValency = atom->valency       ();
    if (!atomSP)          atomSP      = atom->sp            ();
    if (property[0]=='-') property[0] = atom->rcsb_chirality();
    if (property[1]=='-') property[1] = atom->chirality     ();
    if (property[2]=='-') property[2] = atom->leaving       ();
    if (property[3]=='-') property[3] = atom->hb_type       ();
    if (property[4]=='-')
      set_ccp4_coordinates ( atom->ccp4_coordinates() );

  }


  void Atom::write_mem ( PMemIO memIO, int version )  {
  UNUSED_ARGUMENT(version);
    memIO->put_line     ( atomName     );
    memIO->put_line     ( rcsbAtomName );
    memIO->put_line     ( oldAtomName  );
    memIO->put_line     ( chemElement  );
    memIO->put_line     ( energyType   );
    memIO->put_float    ( rcsb_cartn_x );
    memIO->put_float    ( rcsb_cartn_y );
    memIO->put_float    ( rcsb_cartn_z );
    memIO->put_float    ( rcsb_ideal_x );
    memIO->put_float    ( rcsb_ideal_y );
    memIO->put_float    ( rcsb_ideal_z );
    memIO->put_float    ( ccp4_mlib_x  );
    memIO->put_float    ( ccp4_mlib_y  );
    memIO->put_float    ( ccp4_mlib_z  );
    memIO->put_float    ( atomCharge   );
    memIO->put_float    ( atomWeight   );
    memIO->put_float    ( vdwRadius    );
    memIO->put_float    ( vdwhRadius   );
    memIO->put_float    ( ionRadius    );
    memIO->put_integer  ( atomValency  );
    memIO->put_integer  ( atomSP       );
    memIO->write_buffer ( property,sizeof(property) );
  }

  bool Atom::read_mem  ( PMemIO memIO, int version, bool * Ok )  {
  UNUSED_ARGUMENT(version);
  bool success;
    if (Ok)  success = *Ok;
       else  success = true;
    memIO->get_line    ( atomName    ,&success );
    memIO->get_line    ( rcsbAtomName,&success );
    memIO->get_line    ( oldAtomName ,&success );
    memIO->get_line    ( chemElement ,&success );
    memIO->get_line    ( energyType  ,&success );
    memIO->get_float   ( rcsb_cartn_x,&success );
    memIO->get_float   ( rcsb_cartn_y,&success );
    memIO->get_float   ( rcsb_cartn_z,&success );
    memIO->get_float   ( rcsb_ideal_x,&success );
    memIO->get_float   ( rcsb_ideal_y,&success );
    memIO->get_float   ( rcsb_ideal_z,&success );
    memIO->get_float   ( ccp4_mlib_x ,&success );
    memIO->get_float   ( ccp4_mlib_y ,&success );
    memIO->get_float   ( ccp4_mlib_z ,&success );
    memIO->get_float   ( atomCharge  ,&success );
    memIO->get_float   ( atomWeight  ,&success );
    memIO->get_float   ( vdwRadius   ,&success );
    memIO->get_float   ( vdwhRadius  ,&success );
    memIO->get_float   ( ionRadius   ,&success );
    memIO->get_integer ( atomValency ,&success );
    memIO->get_integer ( atomSP      ,&success );
    if (version<2)  {
      memIO->read_buffer ( &(property[1]),3,&success );
      property[0] = property[1];
    } else
      memIO->read_buffer ( property,sizeof(property),&success );
    if (Ok) *Ok = success;
    return success;
  }


  void  get_atom_from_cif ( int & index, mmdb::mmcif::PLoop loop,
                            mmdb::cpstr Tag,
                            int row,Container<Atom> & atoms,
                            int & rc )  {
  mmdb::pstr p;
     p = loop->GetString ( Tag,row,rc );
     if (!p)                     index = -1;
     else if (!strcmp(p,"n/a"))  index = AtomID_NA;
                           else  index = atoms.index ( p );
  }

  void  get_atom_from_cif ( int & index, mmdb::mmcif::PStruct mmCIFStruct,
                            mmdb::cpstr Tag,
                            Container<Atom> & atoms,
                            int & rc )  {
  mmdb::pstr p;
     p = mmCIFStruct->GetString ( Tag,rc );
     if (!p)                     index = -1;
     else if (!strcmp(p,"n/a"))  index = AtomID_NA;
                           else  index = atoms.index ( p );
  }

  void add_atom_to_cif ( mmdb::mmcif::PLoop loop, int atom_no,
                         Container<Atom> & atoms )  {
    if ((atom_no>=0) && (atom_no<atoms.numberOf()))
         loop->AddString ( atoms.at(atom_no)->name(),false );
    else if (atom_no==AtomID_NA)
         loop->AddString ( "n/a",false );
    else loop->AddNoData ( mmdb::mmcif::CIF_NODATA_DOT );
  }

}  // namespace ccp4srs
