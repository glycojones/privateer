//  $Id: ccp4srs_chicenter.cpp $
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
//  **** Module  :  ccp4srs_chicenter  <implementation>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::ChiCenter - chiral center class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#include <string.h>

#include "ccp4srs_chicenter.h"
#include "ccp4srs_defs.h"

namespace ccp4srs  {

  ChiCenter::ChiCenter()  {
    chicenter_id[0] = char(0);  // chiral center id
    atom_center     = -1;       // chiral center atom
    atom_1          = -1;       // ordinal numbers of atoms
    atom_2          = -1;       //   that form
    atom_3          = -1;       //     the center
    volume_sign     = ChiCenter::noSign;
  }

  ChiCenter::~ChiCenter()  {
  }


  #define  VOLUME_SIGN_POSITIVE  "positiv"
  #define  VOLUME_SIGN_NEGATIVE  "negativ"
  #define  VOLUME_SIGN_BOTH      "both"

  char ChiCenter::get_chirality()  {
    switch (volume_sign)  {
      default:
      case noSign   : break;
      case Positive : return 'R';
      case Negative : return 'S';
      case Both     : return 'B';
    }
    return 'N';
  }


  int ChiCenter::readFromCIF ( mmdb::mmcif::PLoop mmCIFLoop,
                               int chiNo, Container<Atom> & atoms )  {
  char S[100];
  int  rc;

    rc = 0;

    mmCIFLoop->CopyString ( chicenter_id,sizeof(ChiCenterID)-1,
                            MMCIF_TAG_ID,chiNo,rc );

    get_atom_from_cif ( atom_center,mmCIFLoop,MMCIF_TAG_ATOM_ID_CENTRE,
                        chiNo,atoms,rc );
    get_atom_from_cif ( atom_1,mmCIFLoop,MMCIF_TAG_ATOM_ID_1,
                        chiNo,atoms,rc );
    get_atom_from_cif ( atom_2,mmCIFLoop,MMCIF_TAG_ATOM_ID_2,
                        chiNo,atoms,rc );
    get_atom_from_cif ( atom_3,mmCIFLoop,MMCIF_TAG_ATOM_ID_3,
                        chiNo,atoms,rc );

    S[0] = char(0);
    mmCIFLoop->CopyString ( S,sizeof(S)-1,MMCIF_TAG_VOLUME_SIGN,
                            chiNo,rc );
    mmdb::LowerCase ( S );
    if (!strcmp(S,VOLUME_SIGN_POSITIVE))
      volume_sign = ChiCenter::Positive;
    else if (!strcmp(S,VOLUME_SIGN_NEGATIVE))
      volume_sign = ChiCenter::Negative;
    else if (!strcmp(S,VOLUME_SIGN_BOTH))
      volume_sign = ChiCenter::Both;
    else
      volume_sign = ChiCenter::noSign;

    return rc;

  }



  /*
  loop_
  _chem_comp_chir.comp_id
  _chem_comp_chir.id
  _chem_comp_chir.atom_id_centre
  _chem_comp_chir.atom_id_1
  _chem_comp_chir.atom_id_2
  _chem_comp_chir.atom_id_3
  _chem_comp_chir.volume_sign
   ATP      chir_01  C4*    C5*    O4*    C3*       negativ
   ATP      chir_02  C3*    C4*    O3*    C2*       negativ
   ATP      chir_03  C2*    C3*    O2*    C1*       negativ
   ATP      chir_04  C1*    O4*    C2*    N9        positiv
  */

  void ChiCenter::makeCIFTags ( mmdb::mmcif::PLoop mmCIFLoop )  {
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_COMP_ID        );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_ID             );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_ATOM_ID_CENTRE );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_ATOM_ID_1      );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_ATOM_ID_2      );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_ATOM_ID_3      );
    mmCIFLoop->AddLoopTag ( MMCIF_TAG_VOLUME_SIGN    );
  }

  void ChiCenter::writeToCIF ( mmdb::mmcif::PLoop mmCIFLoop,
                               mmdb::cpstr monID,
                               Container<Atom> & atoms )  {

    mmCIFLoop->AddString ( monID,false                 );
    mmCIFLoop->AddString ( chicenter_id,false          );

    add_atom_to_cif      ( mmCIFLoop,atom_center,atoms );
    add_atom_to_cif      ( mmCIFLoop,atom_1,atoms      );
    add_atom_to_cif      ( mmCIFLoop,atom_2,atoms      );
    add_atom_to_cif      ( mmCIFLoop,atom_3,atoms      );

    switch (volume_sign)  {
      case ChiCenter::Positive :
        mmCIFLoop->AddString ( VOLUME_SIGN_POSITIVE,false ); break;
      case ChiCenter::Negative :
        mmCIFLoop->AddString ( VOLUME_SIGN_NEGATIVE,false ); break;
      case ChiCenter::Both :
        mmCIFLoop->AddString ( VOLUME_SIGN_BOTH,false );     break;
      default :
        mmCIFLoop->AddNoData ( mmdb::mmcif::CIF_NODATA_DOT );
    }

  }

  void  ChiCenter::copy ( PChiCenter chicenter,
                          mmdb::ivector anmatch )  {
    strcpy ( chicenter_id,chicenter->id() );
    if (anmatch)  {
      atom_center = anmatch[chicenter->center()];
      atom_1 = anmatch[chicenter->atom1()];
      atom_2 = anmatch[chicenter->atom2()];
      atom_3 = anmatch[chicenter->atom3()];
    } else  {
      atom_1 = chicenter->atom1();
      atom_2 = chicenter->atom2();
      atom_3 = chicenter->atom3();
    }
    volume_sign = chicenter->sign();
  }

  void  ChiCenter::write_mem ( PMemIO memIO, int version )  {
  UNUSED_ARGUMENT(version);
    memIO->put_line    ( chicenter_id );
    memIO->put_integer ( atom_center  );
    memIO->put_integer ( atom_1       );
    memIO->put_integer ( atom_2       );
    memIO->put_integer ( atom_3       );
    memIO->put_ibyte   ( volume_sign  );
  }

  bool ChiCenter::read_mem ( PMemIO memIO, int version, bool * Ok )  {
  UNUSED_ARGUMENT(version);
  bool success;
    if (Ok)  success = *Ok;
       else  success = true;
    memIO->get_line    ( chicenter_id,&success );
    memIO->get_integer ( atom_center ,&success );
    memIO->get_integer ( atom_1      ,&success );
    memIO->get_integer ( atom_2      ,&success );
    memIO->get_integer ( atom_3      ,&success );
    memIO->get_ibyte   ( volume_sign ,&success );
    if (Ok)  *Ok = success;
    return success;
  }

}  // namespace ccp4srs
