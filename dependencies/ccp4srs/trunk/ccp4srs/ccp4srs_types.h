//  $Id: ccp4srs_types.h $
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
//  **** Module  :  ccp4srs_types  <interface>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::AtomPair  - atom pair
//       ~~~~~~~~~  ccp4srs::DASelHnds - donor-acceptor selection
//                                       handlers
//
//  (C) E. Krissinel 2010-2014
//
//  =================================================================
//

#ifndef CCP4SRS_TYPES_H
#define CCP4SRS_TYPES_H

#include "mmdb2/mmdb_manager.h"

namespace ccp4srs  {

  enum CCP4SRS_VERSION  {
    MAJOR_VERSION = 1,  //!< major version
    MINOR_VERSION = 0,  //!< minor version
    MICRO_VERSION = 0   //!< micro version
  };

  extern const int CCP4SRS_file_version;

  enum CCP4SRS_RC  {
    CCP4SRS_noDonors         =   4,
    CCP4SRS_noAcceptors      =   3,
    CCP4SRS_noHBonds         =   2,
    CCP4SRS_Incomplete       =   1,
    CCP4SRS_Ok               =   0,
    CCP4SRS_WrongBondType    =  -1,
    CCP4SRS_IndexCorrupt     =  -2,
    CCP4SRS_FileNotFound     =  -3,
    CCP4SRS_EntryNotFound    =  -4,
    CCP4SRS_ReadErrors       =  -5,
    CCP4SRS_EmptyResidue     =  -6,
    CCP4SRS_NoAtomsFound     =  -7,
    CCP4SRS_NoBonds          =  -8,
    CCP4SRS_EmptyResSet      =  -9,
    CCP4SRS_Fail             = -10,
    CCP4SRS_SVD_Fail         = -11,
    CCP4SRS_noCoordHierarchy = -12,
    CCP4SRS_BrokenBonds      = -13
  };

  extern mmdb::cpstr srsIndexFile;
  extern mmdb::cpstr srsGraphFile;
  extern mmdb::cpstr srsStructFile;


  //  special atom ids
  extern const int AtomID_NA;

  // ==================================================================

  DefineStructure(AtomPair);

  struct AtomPair  {
    mmdb::PAtom a1,a2;
  };

  // ==================================================================

  //   SDASelHandles is optionally used in MakeBonds(..), when
  // the latter works for hydrogen bond calculations.
  DefineStructure(DASelHnds);

  struct DASelHnds  {
    int  selHndDonor;
    int  selHndAcceptor;
    int  selHndHydrogen;
    mmdb::SELECTION_KEY selKey;
    void getNewHandles    ( mmdb::PManager MMDB );
    void makeSelIndexes   ( mmdb::PManager MMDB );
    void deleteSelections ( mmdb::PManager MMDB );
  };

}  // namespace ccp4srs

#endif // CCP4SRS_TYPES_H
