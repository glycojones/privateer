//  $Id: ccp4srs_torsion.cpp $
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
//  **** Module  :  ccp4srs_torsion  <implementation>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Torsion  - torsion angle restraint class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#include <string.h>

#include "ccp4srs_torsion.h"
#include "ccp4srs_defs.h"

namespace ccp4srs  {

  Torsion::Torsion()  {
    torsion_id[0]  = char(0);
    atom_1         = -1;
    atom_2         = -1;
    atom_3         = -1;
    atom_4         = -1;
    torsion_period = 0;
    torsion        = 0.0;
    torsion_esd    = 0.0;
  }

  Torsion::~Torsion()  {
  }


  void Torsion::makeCIFTags ( mmdb::mmcif::PLoop mmCIFLoop )  {
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_COMP_ID         );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_ID              );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_ATOM_ID_1       );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_ATOM_ID_2       );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_ATOM_ID_3       );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_ATOM_ID_4       );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_VALUE_ANGLE     );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_VALUE_ANGLE_ESD );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_PERIOD          );
  }

  int Torsion::readFromCIF ( mmdb::mmcif::PLoop mmCIFLoop,
                             int torsionNo, Container<Atom> & atoms ) {
  int rc;

    rc = 0;

    mmCIFLoop->CopyString ( torsion_id,sizeof(TorsionID)-1,
                            MMCIF_TAG_ID,torsionNo,rc );

    get_atom_from_cif ( atom_1,mmCIFLoop,MMCIF_TAG_ATOM_ID_1,
                        torsionNo,atoms,rc );
    get_atom_from_cif ( atom_2,mmCIFLoop,MMCIF_TAG_ATOM_ID_2,
                        torsionNo,atoms,rc );
    get_atom_from_cif ( atom_3,mmCIFLoop,MMCIF_TAG_ATOM_ID_3,
                        torsionNo,atoms,rc );
    get_atom_from_cif ( atom_4,mmCIFLoop,MMCIF_TAG_ATOM_ID_4,
                        torsionNo,atoms,rc );

    mmCIFLoop->CopyReal ( torsion,MMCIF_TAG_VALUE_ANGLE,torsionNo,rc );
    mmCIFLoop->CopyReal ( torsion_esd,MMCIF_TAG_VALUE_ANGLE_ESD,
                                                        torsionNo,rc );

    mmCIFLoop->CopyInteger(torsion_period,MMCIF_TAG_PERIOD,torsionNo,rc);

    return rc;

  }

  /*
  loop_
  _chem_comp_tor.comp_id
  _chem_comp_tor.id
  _chem_comp_tor.atom_id_1
  _chem_comp_tor.atom_id_2
  _chem_comp_tor.atom_id_3
  _chem_comp_tor.atom_id_4
  _chem_comp_tor.value_angle
  _chem_comp_tor.value_angle_esd
  _chem_comp_tor.period
   ATP      var_1    O2A    PA     O3A    PB       -59.994   20.000   1
  */


  void Torsion::writeToCIF ( mmdb::mmcif::PLoop mmCIFLoop,
                             mmdb::cpstr monID,
                             Container<Atom> & atoms )  {

    mmCIFLoop->AddString ( monID,false       );
    mmCIFLoop->AddString ( torsion_id,false  );

    add_atom_to_cif ( mmCIFLoop,atom_1,atoms );
    add_atom_to_cif ( mmCIFLoop,atom_2,atoms );
    add_atom_to_cif ( mmCIFLoop,atom_3,atoms );
    add_atom_to_cif ( mmCIFLoop,atom_4,atoms );
    mmCIFLoop->AddReal    ( torsion    ,"%10.3f" );
    mmCIFLoop->AddReal    ( torsion_esd,"%10.3f" );
    mmCIFLoop->AddInteger ( torsion_period );

  }

  void  Torsion::copy ( PTorsion trs,
                               mmdb::ivector anmatch )  {
    strcpy ( torsion_id,trs->id() );
    if (anmatch)  {
      atom_1 = anmatch[trs->atom1()];
      atom_2 = anmatch[trs->atom2()];
      atom_3 = anmatch[trs->atom3()];
      atom_4 = anmatch[trs->atom4()];
    } else  {
      atom_1 = trs->atom1();
      atom_2 = trs->atom2();
      atom_3 = trs->atom3();
      atom_4 = trs->atom4();
    }
    torsion_period = trs->period();
    torsion        = trs->value ();
    torsion_esd    = trs->esd   ();
  }

  void  Torsion::write_mem ( PMemIO memIO, int version )  {
  UNUSED_ARGUMENT(version);
    memIO->put_line    ( torsion_id     );
    memIO->put_integer ( atom_1         );
    memIO->put_integer ( atom_2         );
    memIO->put_integer ( atom_3         );
    memIO->put_integer ( atom_4         );
    memIO->put_ibyte   ( torsion_period );
    memIO->put_float   ( torsion        );
    memIO->put_float   ( torsion_esd    );
  }

  bool Torsion::read_mem ( PMemIO memIO, int version, bool * Ok )  {
  UNUSED_ARGUMENT(version);
  bool success;
    if (Ok)  success = *Ok;
       else  success = true;
    memIO->get_line    ( torsion_id    ,&success );
    memIO->get_integer ( atom_1        ,&success );
    memIO->get_integer ( atom_2        ,&success );
    memIO->get_integer ( atom_3        ,&success );
    memIO->get_integer ( atom_4        ,&success );
    memIO->get_ibyte   ( torsion_period,&success );
    memIO->get_float   ( torsion       ,&success );
    memIO->get_float   ( torsion_esd   ,&success );
    if (Ok)  *Ok = success;
    return success;
  }

}  // namespace ccp4srs
