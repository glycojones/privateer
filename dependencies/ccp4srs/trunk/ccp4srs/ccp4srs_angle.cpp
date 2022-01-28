//  $Id: ccp4srs_angle.cpp $
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
//  **** Module  :  ccp4srs_angle  <implementation>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Angle  - angle restraint class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#include "ccp4srs_angle.h"
#include "ccp4srs_defs.h"

namespace ccp4srs  {

  Angle::Angle()  {
    atom_1    = -1;
    atom_2    = -1;
    atom_3    = -1;
    angle     = 0.0;
    angle_esd = 0.0;
  }

  Angle::~Angle()  {
  }


  int Angle::readFromCIF ( mmdb::mmcif::PLoop mmCIFLoop, int angleNo,
                           Container<Atom> & atoms )  {
  int rc;

    rc = 0;

    get_atom_from_cif ( atom_1,mmCIFLoop,MMCIF_TAG_ATOM_ID_1,
                        angleNo,atoms,rc );
    get_atom_from_cif ( atom_2,mmCIFLoop,MMCIF_TAG_ATOM_ID_2,
                        angleNo,atoms,rc );
    get_atom_from_cif ( atom_3,mmCIFLoop,MMCIF_TAG_ATOM_ID_3,
                        angleNo,atoms,rc );

    mmCIFLoop->CopyReal(angle    ,MMCIF_TAG_VALUE_ANGLE    ,angleNo,rc);
    mmCIFLoop->CopyReal(angle_esd,MMCIF_TAG_VALUE_ANGLE_ESD,angleNo,rc);

    return rc;

  }

  /*
  loop_
  _chem_comp_angle.comp_id
  _chem_comp_angle.atom_id_1
  _chem_comp_angle.atom_id_2
  _chem_comp_angle.atom_id_3
  _chem_comp_angle.value_angle
  _chem_comp_angle.value_angle_esd
   ATP      O2A    PA     O1A     119.900    3.000
  */


  void Angle::makeCIFTags ( mmdb::mmcif::PLoop mmCIFLoop )  {
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_COMP_ID         );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_ATOM_ID_1       );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_ATOM_ID_2       );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_ATOM_ID_3       );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_VALUE_ANGLE     );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_VALUE_ANGLE_ESD );
  }

  void Angle::writeToCIF ( mmdb::mmcif::PLoop mmCIFLoop, mmdb::cpstr monID,
                              Container<Atom> & atoms )  {
    mmCIFLoop->AddString ( monID,false );
    add_atom_to_cif ( mmCIFLoop,atom_1,atoms );
    add_atom_to_cif ( mmCIFLoop,atom_2,atoms );
    add_atom_to_cif ( mmCIFLoop,atom_3,atoms );
    mmCIFLoop->AddReal ( angle    ,"%10.3f" );
    mmCIFLoop->AddReal ( angle_esd,"%10.3f" );
  }

  void  Angle::copy ( PAngle ang, mmdb::ivector anmatch )  {
    if (anmatch)  {
      atom_1 = anmatch[ang->atom1()];
      atom_2 = anmatch[ang->atom2()];
      atom_3 = anmatch[ang->atom3()];
    } else  {
      atom_1 = ang->atom1();
      atom_2 = ang->atom2();
      atom_3 = ang->atom3();
    }
    angle     = ang->value();
    angle_esd = ang->esd  ();
  }


  void  Angle::write_mem ( PMemIO memIO, int version )  {
  UNUSED_ARGUMENT(version);
    memIO->put_integer ( atom_1    );
    memIO->put_integer ( atom_2    );
    memIO->put_integer ( atom_3    );
    memIO->put_float   ( angle     );
    memIO->put_float   ( angle_esd );
  }

  bool Angle::read_mem ( PMemIO memIO, int version, bool * Ok ) {
  UNUSED_ARGUMENT(version);
  bool success;
    if (Ok)  success = *Ok;
       else  success = true;
    memIO->get_integer ( atom_1   ,&success );
    memIO->get_integer ( atom_2   ,&success );
    memIO->get_integer ( atom_3   ,&success );
    memIO->get_float   ( angle    ,&success );
    memIO->get_float   ( angle_esd,&success );
    if (Ok)  *Ok = success;
    return success;
  }

}  // namespace ccp4srs
