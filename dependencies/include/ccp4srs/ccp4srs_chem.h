//  $Id: ccp4srs_chem.h $
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
//  **** Module  :  ccp4srs_chem  <interface>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Chem - SRS chemistry manager class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#ifndef CCP4SRS_CHEM_H
#define CCP4SRS_CHEM_H

#include "ccp4srs_base.h"
#include "ccp4srs_types.h"

namespace ccp4srs  {

  // ==================================================================

  DefineClass(Chem);

  class Chem : public Base  {

    public:
      Chem ();
      ~Chem();

      CCP4SRS_RC getEnergyTypes ( mmdb::PResidue  R,
                                  mmdb::io::PFile structFile );
      CCP4SRS_RC getEnergyTypes ( mmdb::PPResidue R, int nRes,
                                  mmdb::io::PFile structFile );
      CCP4SRS_RC getEnergyTypes ( mmdb::PChain  chain,
                                  mmdb::io::PFile structFile );
      CCP4SRS_RC getEnergyTypes ( mmdb::PModel  model,
                                  mmdb::io::PFile structFile );
      CCP4SRS_RC getEnergyTypes ( mmdb::PManager MMDB,
                                  mmdb::io::PFile structFile );

      CCP4SRS_RC addHydrogens   ( mmdb::PResidue    R,
                                  mmdb::io::PFile structFile );
      CCP4SRS_RC addHydrogens   ( mmdb::PChain  chain,
                                  mmdb::io::PFile structFile );
      CCP4SRS_RC addHydrogens   ( mmdb::PModel  model,
                                  mmdb::io::PFile structFile );
      CCP4SRS_RC addHydrogens   ( mmdb::PManager MMDB,
                                  mmdb::io::PFile structFile );


      //   MakeBonds(..) makes bonds between atoms in MMDB's residue R
      // from data found in SBase. Residue R must be associated with
      // coordinate hierarchy. Data is retrieved from SBase on the basis
      // of residue name only. In case of multiple conformations, if
      // altLoc:
      //    NULL  - the highest occupancy atom will be taken
      //            if all occupancies are equal then atom with
      //            first altLoc taken
      //    other - atoms with given altLoc are taken. If such
      //            altLoc is not found, the function does as if
      //            NULL value for altLoc is given.
      //   If selHandles is not NULL, the function also selects atoms
      // in the residue according to their hydrogen bond attributes.
      // This is a special option for hydrogen bond calculations
      //   If ignoreNegSigOcc is set True then the function will ignore
      // atoms with negative occupancy standard deviation. Such atoms
      // may be hydrogens added by CSBase0::AddHydrogen(..) function,
      // in general any atoms added by CSBAtom::MakeCAtom(..) function.
      // Added hydrogens may be ignored if MakeBonds is used in
      // CSbase::CalcHBonds(..) function.
      //   Return:
      //     SBASE_Ok             success
      //     SBASE_FileNotFound   non-initiated SBase
      //     SBASE_StructNotFound the residue's name is not found
      //                          in SBase
      //     SBASE_EmptyResidue   residue R does not contain atoms
      //     SBASE_NoAtomsFound   SBase entry does not contain atoms
      //     SBASE_BrokenBonds    some bonds could not be set up because
      //                          of missing atoms in R. This could be
      //                          a result of residue R named wrongly.
      CCP4SRS_RC makeBonds ( mmdb::PResidue  R,
                             mmdb::pstr      altLoc,
                             mmdb::io::PFile structFile,
                             PDASelHnds      selHandles,
                             bool            ignoreNegSigOcc );

      //   CalcHBonds(..) calculates hydrogen bonds and salt bridges
      // between atoms of residues given in Res1[0..nres1-1] and
      // Res2[0..nres2-1]. H-bonds are returned in HBond[0..nHBonds-1]
      // and salt bridges are returned in SBridge[0..nSBridges-1].
      //
      //   On input:
      //     Res1, Res2         must belong to the same MMDB
      //     HBond, SBridge     should be set NULL, otherwise the
      //                        function attempts to deallocate them
      //     nHBonds, nSBridges ignored
      //     structFile         may be a pointer to open SBase stucture
      //                        file in order to save on file open
      //                        operation. If structFile is set NULL,
      //                        the function will open the SBase
      //                        structure file by itself
      //     altLoc             specifies the alternative location for
      //                        donors and acceptors.
      //                          altLoc==NULL  the highest occupancy
      //                                        atoms are taken. If all
      //                                        occupancies are equal
      //                                        equal then atom with
      //                                        first altLoc taken
      //                          altLoc!=NULL  atoms with given altLoc
      //                                        are taken. If such
      //                                        altLoc is not found,
      //                                        the function acts as
      //                                        if NULL value for altLoc
      //                                        were given.
      //     ignoreNegSigOcc    if it is set True, then the function
      //                        ignores atoms with negative occupancy
      //                        standard deviation. Such atoms may be
      //                        hydrogens added by addHydrogen(..)
      //                        function, or, in general, any atoms added
      //                        by CCP4SRSAtom::MakeCAtom(..) function.
      //                        Such added hydrogens are not guaranteed to
      //                        be in correct place, therefore the
      //                        function may mistake on some hydrogen
      //                        bonds if they are not neglected.
      //
      //   On output:
      //     Allocated arrays HBond[0..nHBonds-1] and
      //     SBridge[0..nSBridges-1]. If no bonds/bridges were found,
      //     the corresponding array is not allocated and set NULL.
      //     Application is responsible for deallocation of the arrays,
      //     when not needed, using statements
      //       if (HBond)    delete[] HBond;
      //       if (SBridge)  delete[] SBridge;
      //     HBond[i].a1, SBridge[i].a1 always refer atom from Res1[],
      //     and HBond[i].a2, SBridge[i].a2 always refer atom from
      //     Res2[].
      //
      CCP4SRS_RC CalcHBonds ( mmdb::PPResidue Res1, int nres1,
                              mmdb::PPResidue Res2, int nres2,
                              RPAtomPair     HBond, int & nHBonds,
                              RPAtomPair   SBridge, int & nSBridges,
                              mmdb::io::PFile structFile=NULL,
                              mmdb::pstr      altLoc=NULL,
                              bool ignoreNegSigOcc=false );

    protected:
      // Geometry parameters for H-bond calculations
      mmdb::realtype minDAdist,maxSBdist,maxDAdist,maxHAdist2;
      mmdb::realtype maxDHAcos,maxHAAcos,maxDAAcos,maxDDAcos;

      void init_chem();

  };

}  // namespace ccp4srs

#endif // CCP4SRS_CHEM_H
