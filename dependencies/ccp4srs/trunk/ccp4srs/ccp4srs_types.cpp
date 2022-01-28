//  $Id: ccp4srs_types.cpp $
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
//  **** Module  :  ccp4srs_types  <implementation>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::AtomPair  - atom pair
//       ~~~~~~~~~  ccp4srs::DASelHnds - donor-acceptor selection
//                                       handlers
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#include "ccp4srs_types.h"

namespace ccp4srs  {

  const int CCP4SRS_file_version = 3;

  mmdb::cpstr srsIndexFile  = "index.srs";
  mmdb::cpstr srsGraphFile  = "graph.srs";
  mmdb::cpstr srsStructFile = "struct.srs";

  const int AtomID_NA = -10000;

  void DASelHnds::getNewHandles ( mmdb::PManager MMDB )  {
    selHndDonor    = MMDB->NewSelection();
    selHndAcceptor = MMDB->NewSelection();
    selHndHydrogen = MMDB->NewSelection();
    selKey = mmdb::SKEY_OR;
  }

  void DASelHnds::makeSelIndexes ( mmdb::PManager MMDB )  {
    MMDB->MakeSelIndex ( selHndDonor    );
    MMDB->MakeSelIndex ( selHndAcceptor );
    MMDB->MakeSelIndex ( selHndHydrogen );
  }

  void DASelHnds::deleteSelections ( mmdb::PManager MMDB )  {
    MMDB->DeleteSelection ( selHndDonor    );
    MMDB->DeleteSelection ( selHndAcceptor );
    MMDB->DeleteSelection ( selHndHydrogen );
  }

}  // namespace ccp4srs
