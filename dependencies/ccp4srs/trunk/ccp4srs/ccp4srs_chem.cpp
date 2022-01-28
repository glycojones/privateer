//  $Id: ccp4srs_chem.cpp $
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
//  **** Module  :  ccp4srs_chem  <implementation>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Chem - SRS chemistry manager class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#include <string.h>

#include "ccp4srs_chem.h"
#include "ccp4srs_defs.h"

namespace ccp4srs  {

  Chem::Chem() : Base()  {
    init_chem();
  }

  Chem::~Chem()  {}

  void Chem::init_chem()  {

    minDAdist  = 2.0;
    maxDAdist  = 3.9;  // max distance for H-bonds
    maxSBdist  = 4.0;  // max distance for salt bridges
    maxHAdist2 = 2.5*2.5;
    maxDHAcos  = 0.0;
    maxHAAcos  = 0.0;
    maxDAAcos  = 0.0;
    maxDDAcos  = 0.0;

  }


  CCP4SRS_RC Chem::getEnergyTypes ( mmdb::PResidue R,
                                    mmdb::io::PFile structFile ) {
  PMonomer       monomer;
  PAtom          srs_atom;
  mmdb::PPAtom   A;
  mmdb::AtomName aname;
  int            i,j,natoms,nat;
  CCP4SRS_RC     rc;

    R->GetAtomTable ( A,natoms );
    if (!A)  return CCP4SRS_EmptyResidue;

    for (i=0;i<natoms;i++)
      if (A[i])  A[i]->energyType[0] = char(0);

    monomer = getMonomer ( R->GetResName(),structFile );
    if (monomer)  {
      nat = monomer->n_atoms();
      if (nat>0)  {
        for (i=0;i<natoms;i++)
          if (A[i])  {
            for (j=0;j<nat;j++)  {
              srs_atom = monomer->atom(j);
              if (srs_atom)  {
                if (!strcmp(A[i]->name,srs_atom->name_pdb(aname)))  {
                  strcpy ( A[i]->energyType,srs_atom->energy_type() );
                  A[i]->charge = srs_atom->charge();
                }
              }
            }
          }
        rc = CCP4SRS_Ok;
      } else
        rc = CCP4SRS_NoAtomsFound;
      delete monomer;
    } else
      rc = CCP4SRS_EntryNotFound;

    return rc;

  }


  CCP4SRS_RC Chem::getEnergyTypes ( mmdb::PPResidue R, int nRes,
                                    mmdb::io::PFile structFile )  {
  mmdb::io::PFile sFile;
  PMonomer        monomer;
  PAtom           srs_atom;
  mmdb::PPAtom    A;
  mmdb::AtomName  aname;
  mmdb::ovector   B;
  int             i,j,k,n,natoms,nat;

    mmdb::GetVectorMemory ( B,nRes,0 );
    for (i=0;i<nRes;i++)
      B[i] = (!R[i]);

    if (structFile)  sFile = structFile;
               else  sFile = getStructFile();
    if (!sFile)  return CCP4SRS_FileNotFound;

    i = 0;
    while (i<nRes)  {
      while (i<nRes)
        if (B[i])  i++;
             else  break;
      if (i<nRes)  {
        monomer = getMonomer ( R[i]->GetResName(),sFile );
        j   = i;
        while (j<nRes)  {
          B[j] = true;
          R[j]->GetAtomTable ( A,natoms );
          if (A)  {
            for (k=0;k<natoms;k++)
              if (A[k])  A[k]->energyType[0] = char(0);
            if (monomer)  {
              nat = monomer->n_atoms();
              if (nat>0)  {
                for (k=0;k<natoms;k++)
                  if (A[k])  {
                    for (n=0;n<nat;n++)  {
                      srs_atom = monomer->atom(n);
                      if (srs_atom) {
                        if (!strcmp(A[k]->name,
                                    srs_atom->name_pdb(aname))) {
                          strcpy ( A[k]->energyType,
                                   srs_atom->energy_type() );
                          A[k]->charge = srs_atom->charge();
                        }
                      }
                    }
                  }
              }
            }
          }
          j++;
          while (j<nRes)
            if (!B[j])  {
              if (!strcmp(R[i]->name,R[j]->name))  break;
                                             else  j++;
            } else
              j++;
        }
        if (monomer)  {
          delete monomer;
          monomer = NULL;
        }
        i++;
      }
    }

    if (!structFile)  delete sFile;
    mmdb::FreeVectorMemory ( B,0 );

    return CCP4SRS_Ok;

  }


  CCP4SRS_RC Chem::getEnergyTypes ( mmdb::PChain chain,
                                    mmdb::io::PFile structFile )  {
  mmdb::PPResidue Res;
  int             nRes;
    chain->GetResidueTable ( Res,nRes );
    if (nRes>0)
          return getEnergyTypes ( Res,nRes,structFile );
    else  return CCP4SRS_EmptyResSet;
  }

  CCP4SRS_RC Chem::getEnergyTypes ( mmdb::PModel model,
                                    mmdb::io::PFile structFile )  {
  mmdb::PPResidue Res;
  mmdb::PPChain   chain;
  mmdb::PManager  MMDB;
  int             selHnd,i,nRes,nChains;
  CCP4SRS_RC      rc;

    rc   = CCP4SRS_Ok;
    MMDB = mmdb::PManager(model->GetCoordHierarchy());

    if (MMDB)  {

      selHnd = MMDB->NewSelection();
      MMDB->Select ( selHnd,mmdb::STYPE_RESIDUE,model->GetSerNum(),
                     "*",mmdb::ANY_RES,"*",mmdb::ANY_RES,"*",
                     "*","*","*","*",mmdb::SKEY_NEW );
      MMDB->GetSelIndex ( selHnd,Res,nRes );
      if (nRes>0)  rc = getEnergyTypes ( Res,nRes,structFile );
             else  rc = CCP4SRS_EmptyResSet;
      MMDB->DeleteSelection ( selHnd );

    } else  {

      model->GetChainTable ( chain,nChains );
      for (i=0;i<nChains;i++)
        if (chain[i])  {
          chain[i]->GetResidueTable ( Res,nRes );
          if (nRes>0) getEnergyTypes ( Res,nRes,structFile );
        }

    }

    return rc;

  }

  CCP4SRS_RC Chem::getEnergyTypes ( mmdb::PManager  MMDB,
                                    mmdb::io::PFile structFile ) {
  mmdb::PPResidue Res;
  int             selHnd,nRes;
  CCP4SRS_RC      rc;

    rc = CCP4SRS_Ok;

    selHnd = MMDB->NewSelection();
    MMDB->Select ( selHnd,mmdb::STYPE_RESIDUE,0,
                   "*",mmdb::ANY_RES,"*",mmdb::ANY_RES,"*",
                   "*","*","*","*",mmdb::SKEY_NEW );
    MMDB->GetSelIndex ( selHnd,Res,nRes );
    if (nRes>0)  rc = getEnergyTypes ( Res,nRes,structFile );
           else  rc = CCP4SRS_EmptyResSet;

    MMDB->DeleteSelection ( selHnd );

    return rc;

  }


  CCP4SRS_RC Chem::addHydrogens ( mmdb::PResidue R,
                                  mmdb::io::PFile structFile )  {
  //   Return:
  //     CCP4SRS_Ok             success
  //     CCP4SRS_EmptyResidue   residue R does not contain atoms
  //     CCP4SRS_NoAtomsFound   SBStructure does not contain atoms
  //     CCP4SRS_NoBonds        SBStructure does not contain bonds
  //     CCP4SRS_NoAtomsData    SBStructure is not complete
  //     CCP4SRS_NoSimilarity   too few coomon atom names in R and SBase
  //                          entry with the same structure name
  //     CCP4SRS_SuperpositionFailed  failed residue superposition
  // NOTE: the function does not rearranges existing atoms in the
  // residue, but places the hydrogens on top of them (leaving the
  // Ter pseudoatom, if found, on top of the list)
  mmdb::io::PFile sFile;
  PMonomer        monomer;
  CCP4SRS_RC      rc;

    if (structFile)  sFile = structFile;
               else  sFile = getStructFile();

    if (!sFile)  return CCP4SRS_FileNotFound;

    monomer = getMonomer ( R->GetResName(),sFile );
    if (!structFile)  delete sFile;

    if (!monomer)  return CCP4SRS_EntryNotFound;

    rc = monomer->addHydrogens ( R );

    delete monomer;

    return rc;

  }

  CCP4SRS_RC Chem::addHydrogens ( mmdb::PChain chain,
                                  mmdb::io::PFile structFile )  {
  mmdb::io::PFile sFile;
  mmdb::PPResidue Res;
  int             i,k,nRes;
  bool            B;
  CCP4SRS_RC      rc;

    if (structFile)  sFile = structFile;
               else  sFile = getStructFile();

    if (!sFile)  return CCP4SRS_FileNotFound;

    rc = CCP4SRS_Ok;
    B  = false;
    chain->GetResidueTable ( Res,nRes );
    for (i=0;i<nRes;i++)
      if (Res[i])  {
        k = addHydrogens ( Res[i],sFile );
        if (k==CCP4SRS_Ok)  B = true;
                      else  rc = CCP4SRS_Incomplete;
      }
    if (!B)  rc = CCP4SRS_Fail;

    if (!structFile)  delete sFile;

    return rc;

  }


  CCP4SRS_RC Chem::addHydrogens ( mmdb::PModel    model,
                                  mmdb::io::PFile structFile )  {
  mmdb::io::PFile sFile;
  mmdb::PPChain   chain;
  mmdb::PPResidue Res;
  int             i,j,k,nChains,nRes;
  bool            B;
  CCP4SRS_RC      rc;

    if (structFile)  sFile = structFile;
               else  sFile = getStructFile();

    if (!sFile)  return CCP4SRS_FileNotFound;

    rc = CCP4SRS_Ok;
    B  = false;
    model->GetChainTable ( chain,nChains );
    for (i=0;i<nChains;i++)
      if (chain[i])  {
        chain[i]->GetResidueTable ( Res,nRes );
        for (j=0;j<nRes;j++)
          if (Res[j])  {
            k = addHydrogens ( Res[j],sFile );
            if (k==CCP4SRS_Ok)  B = true;
                          else  rc = CCP4SRS_Incomplete;
          }
      }
    if (!B)  rc = CCP4SRS_Fail;

    if (!structFile)  delete sFile;

    return rc;

  }

  CCP4SRS_RC Chem::addHydrogens ( mmdb::PManager MMDB,
                                  mmdb::io::PFile structFile )  {
  mmdb::io::PFile sFile;
  mmdb::PPModel   model;
  mmdb::PPChain   chain;
  mmdb::PPResidue Res;
  int             i,j,n,k,nModels,nChains,nRes;
  bool            B;
  CCP4SRS_RC      rc;

    if (structFile)  sFile = structFile;
               else  sFile = getStructFile();

    if (!sFile)  return CCP4SRS_FileNotFound;

    rc = CCP4SRS_Ok;
    B  = false;
    MMDB->GetModelTable ( model,nModels );
    for (i=0;i<nModels;i++)
      if (model[i])  {
        model[i]->GetChainTable ( chain,nChains );
        for (j=0;j<nChains;j++)
          if (chain[j])  {
            chain[j]->GetResidueTable ( Res,nRes );
            for (n=0;n<nRes;n++)
              if (Res[n])  {
                k = addHydrogens ( Res[n],sFile );
                if (k==CCP4SRS_Ok)  B = true;
                              else  rc = CCP4SRS_Incomplete;
              }
          }
      }
    if (!B)  rc = CCP4SRS_Fail;

    if (!structFile)  delete sFile;

    return rc;

  }



  CCP4SRS_RC Chem::makeBonds ( mmdb::PResidue  R,
                               mmdb::pstr      altLoc,
                               mmdb::io::PFile structFile,
                               PDASelHnds      selHandles,
                               bool            ignoreNegSigOcc )  {
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
  // This is a special option for hydrogen bond calculations.
  //   If ignoreNegSigOcc is set true then the function will ignore
  // atoms with negative occupancy standard deviation. Such atoms
  // may be hydrogens added by CSBase0::AddHydrogens(..) function,
  // in general any atoms added by CSBAtom::MakeCAtom(..) function.
  // Added hydrogens may be ignored if MakeBonds is used in
  // CSbase::CalcHBonds(..) function.
  //   Return:
  //     SBASE_Ok             success
  //     SBASE_FileNotFound   non-initiated SBase
  //     SBASE_StructNotFound the residue's name is not found in SBase
  //     SBASE_EmptyResidue   residue R does not contain atoms
  //     SBASE_NoAtomsFound   SBase entry does not contain atoms
  //     SBASE_BrokenBonds    some bonds could not be set up because
  //                          of missing atoms in R. This could be
  //                          a result of residue R named wrongly.
  PMonomer       monomer;
  PBond          bondi;
  mmdb::PManager MMDB;
  mmdb::PPAtom   A;
  mmdb::ivector  anmatch;
  int            natoms, mon_natoms, i,i1,i2;
  mmdb::SELECTION_KEY sk;
  CCP4SRS_RC     rc;

    R->GetAtomTable ( A,natoms );
    if (!A)  return CCP4SRS_EmptyResidue;

    for (i=0;i<natoms;i++)
      if (A[i])  A[i]->FreeBonds();

    monomer = getMonomer ( R->GetResName(),structFile );
    if (!monomer)  return CCP4SRS_EntryNotFound;

    mon_natoms = monomer->n_atoms();
    if (mon_natoms<=0)  {
      delete monomer;
      return CCP4SRS_NoAtomsFound;
    }

    mmdb::GetVectorMemory ( anmatch,mon_natoms,0 );
    monomer->getAtomNameMatch ( A,natoms,altLoc,anmatch );
    if (ignoreNegSigOcc)
      for (i=0;i<mon_natoms;i++)  {
        i1 = anmatch[i];
        if (i1>=0)  {
          if (A[i1]->sigOcc<0.0)  anmatch[i] = -1;
        }
      }

    if (selHandles)  {
      MMDB = mmdb::PManager(R->GetCoordHierarchy());
      if (MMDB)  {
        sk = selHandles->selKey;
        for (i=0;i<mon_natoms;i++)  {
          i1 = anmatch[i];
          if (i1>=0)
            switch (monomer->atom(i)->hb_type())  {
              case 'D' :  MMDB->SelectAtom ( selHandles->selHndDonor,
                                             A[i1],sk,false );
                        break;
              case 'A' :  MMDB->SelectAtom ( selHandles->selHndAcceptor,
                                             A[i1],sk,false );
                        break;
              case 'B' :  MMDB->SelectAtom ( selHandles->selHndDonor,
                                             A[i1],sk,false );
                          MMDB->SelectAtom ( selHandles->selHndAcceptor,
                                             A[i1],sk,false );
                        break;
              case 'H' :  MMDB->SelectAtom ( selHandles->selHndHydrogen,
                                             A[i1],sk,false );
                        break;
              default  :
              case 'N' : ;
            }
        }
      }
    }

    rc = CCP4SRS_Ok;
    for (i=0;i<monomer->n_bonds();i++)  {
      bondi = monomer->bond(i);
      i1 = anmatch[bondi->atom1()];
      i2 = anmatch[bondi->atom2()];
      if ((i1>=0) && (i2>=0))  {
        A[i1]->AddBond ( A[i2],bondi->order(),2 );
        A[i2]->AddBond ( A[i1],bondi->order(),2 );
      } else
        rc = CCP4SRS_BrokenBonds;
    }

    mmdb::FreeVectorMemory ( anmatch,0 );
    delete monomer;

    return rc;

  }

  void _addAtomPair ( mmdb::PAtom a1, mmdb::PAtom a2,
                      PAtomPair & Pair, int & nPairs, int & nPAlloc )  {
  PAtomPair AP;
  int       i,i1,i2,j1,j2;
  bool      Found;

    Found = false;
    i1    = a1->GetIndex();
    i2    = a2->GetIndex();
    for (i=0;(i<nPairs) && (!Found);i++)  {
      j1 = Pair[i].a1->GetIndex();
      j2 = Pair[i].a2->GetIndex();
      Found = (((i1==j1) && (i2==j2)) ||
               ((i1==j2) && (i2==j1)));
    }

    if (!Found)  {
      if (nPairs>=nPAlloc)  {
        nPAlloc = nPairs+20;
        AP = new AtomPair[nPAlloc];
        for (i=0;i<nPairs;i++)  {
          AP[i].a1 = Pair[i].a1;
          AP[i].a2 = Pair[i].a2;
        }
        if (Pair)  delete[] Pair;
        Pair = AP;
      }
      Pair[nPairs].a1 = a1;
      Pair[nPairs].a2 = a2;
      nPairs++;
    }

  }

  CCP4SRS_RC Chem::CalcHBonds ( mmdb::PPResidue Res1, int nres1,
                                mmdb::PPResidue Res2, int nres2,
                                RPAtomPair     HBond, int & nHBonds,
                                RPAtomPair   SBridge, int & nSBridges,
                                mmdb::io::PFile structFile,
                                mmdb::pstr altLoc,
                                bool ignoreNegSigOcc )  {
  mmdb::io::PFile sFile;
  mmdb::PManager  MMDB;
  mmdb::PPAtom    Donor1,Acceptor1, Donor2,Acceptor2;
  mmdb::PAtom     D,A,H;
  mmdb::PContact  Contact;
  mmdb::PAtomBond ABond,DBond;
  DASelHnds       selHandles1,selHandles2;
  mmdb::pstr      resName;
  int             nDonors1,nAcceptors1, nDonors2,nAcceptors2, nContacts;
  int             i,j,k,nDBonds,nABonds,nHBAlloc,nSBAlloc;
  bool            isHBond,isSBridge;

    if (HBond)  {
      delete[] HBond;
      HBond = NULL;
    }
    nHBonds  = 0;
    nHBAlloc = 0;

    if (SBridge)  {
      delete[] SBridge;
      SBridge = NULL;
    }
    nSBridges = 0;
    nSBAlloc  = 0;

    //  1. Calculate bonds between atoms in given residues and
    //     select donors and acceptors

    i = 0;
    while ((i<nres1) && (!Res1[i]))  i++;
    if (i>=nres1)  return CCP4SRS_EmptyResSet;

    MMDB = mmdb::PManager(Res1[i]->GetCoordHierarchy());
    if (!MMDB)  return CCP4SRS_noCoordHierarchy;

    if (structFile)  sFile = structFile;
               else  sFile = getStructFile();

    if (!sFile)  return CCP4SRS_FileNotFound;

    selHandles1.getNewHandles ( MMDB );
    for (i=0;i<nres1;i++)
      if (Res1[i])
        makeBonds ( Res1[i],altLoc,sFile,&selHandles1,ignoreNegSigOcc );
    selHandles1.makeSelIndexes ( MMDB );

    selHandles2.getNewHandles ( MMDB );
    for (i=0;i<nres2;i++)
      if (Res2[i])
        makeBonds ( Res2[i],altLoc,sFile,&selHandles2,ignoreNegSigOcc );
    selHandles2.makeSelIndexes ( MMDB );

    if (!structFile)  delete sFile;

    // 2. Calculate contacts between donors and acceptors as
    //    potential hydrogen bond contacts

    MMDB->GetSelIndex ( selHandles1.selHndDonor,Donor1,nDonors1 );
    MMDB->GetSelIndex ( selHandles2.selHndDonor,Donor2,nDonors2 );
    if ((nDonors1<=0) && (nDonors2<=0))  {
      selHandles1.deleteSelections ( MMDB );
      selHandles2.deleteSelections ( MMDB );
      return CCP4SRS_noDonors;
    }

    MMDB->GetSelIndex(selHandles1.selHndAcceptor,Acceptor1,nAcceptors1);
    MMDB->GetSelIndex(selHandles2.selHndAcceptor,Acceptor2,nAcceptors2);
    if ((nAcceptors1<=0) && (nAcceptors2<=0))  {
      selHandles1.deleteSelections ( MMDB );
      selHandles2.deleteSelections ( MMDB );
      return CCP4SRS_noAcceptors;
    }

    if ((nDonors1*nAcceptors2<=0) && (nDonors2*nAcceptors1<=0))  {
      selHandles1.deleteSelections ( MMDB );
      selHandles2.deleteSelections ( MMDB );
      return CCP4SRS_noHBonds;
    }

    //   We now calculate contacts such that 1st contacting atom, either
    // acceptor or donor, belongs to 1st set of residues, and the second
    // one - to 2nd set of residues. Therefore we run SeekContacts on
    // two sets of donors and acceptors, identified by different group
    // id, merging the array of contacts for convenience.
    Contact   = NULL;
    nContacts = 0;
    MMDB->SeekContacts ( Donor1,nDonors1,Acceptor2,nAcceptors2,
                         minDAdist,mmdb::RMax(maxDAdist,maxSBdist),0,
                         Contact,nContacts,0,NULL,1,0 );
    MMDB->SeekContacts ( Acceptor1,nAcceptors1,Donor2,nDonors2,
                         minDAdist,mmdb::RMax(maxDAdist,maxSBdist),0,
                         Contact,nContacts,0,NULL,2,0 );
    if (nContacts<=0)  {
      selHandles1.deleteSelections ( MMDB );
      selHandles2.deleteSelections ( MMDB );
      return CCP4SRS_noHBonds;
    }


    // 3. Check all contacts for h-bond geometry

    // merge all hydrogens into one selection as it is used
    // for checking with only
    MMDB->Select ( selHandles1.selHndHydrogen,mmdb::STYPE_ATOM,
                   selHandles2.selHndHydrogen,mmdb::SKEY_OR );

    for (i=0;i<nContacts;i++)  {

      if (Contact[i].group<=1)  {
        D = Donor1   [Contact[i].id1];
        A = Acceptor2[Contact[i].id2];
      } else  {
        A = Acceptor1[Contact[i].id1];
        D = Donor2   [Contact[i].id2];
      }

      isSBridge = false;

      if ((Contact[i].dist<=maxSBdist)       &&
          (D->GetResidue()!=A->GetResidue()) &&
          (!strcmp(D->element," N")) && (!strcmp(A->element," O")))  {
        // Check for salt bridge, which may be formed only by N-O
        // pairs of aminoacid atoms at distances less then maxSBdist
        if (!strcmp(D->name," N  "))  {
          // mainchain nitrogen can form salt bridge only at N-terminus
          isSBridge = D->isNTerminus();
        } else  {
          // other nitrogens can form salt bridge only in LYS, ARG
          // and HIS
          resName   = D->GetResName();
          isSBridge = ((!strcmp(resName,"LYS")) ||
                       (!strcmp(resName,"ARG")) ||
                       (!strcmp(resName,"HIS")));
        }
        if (isSBridge)  {
          if ((!strcmp(A->name," O  ")) || (!strcmp(A->name," OXT")))  {
            // mainchain oxygens can form salt bridge only at C-terminus
            isSBridge = A->isCTerminus();
          } else  {
            // other oxygens can form salt bridge only in GLU and ASP
            resName   = A->GetResName();
            isSBridge = ((!strcmp(resName,"GLU")) ||
                         (!strcmp(resName,"ASP")));
          }
          if (isSBridge)  {
            if (Contact[i].group<=1)
                  _addAtomPair ( D,A,SBridge,nSBridges,nSBAlloc );
            else  _addAtomPair ( A,D,SBridge,nSBridges,nSBAlloc );
          }
        }
      }

      if ((Contact[i].dist<=maxDAdist) && (!isSBridge))  {
        // Check for H-bond if salt bridge was not identified
        D->GetBonds ( DBond,nDBonds );
        A->GetBonds ( ABond,nABonds );
        if (nABonds>0)  {
          // Check whether there are hydrogens bound to the donor,
          // and if they are then calculate h-bonds using them
          H = NULL;
          for (j=0;j<nDBonds;j++)
            if ((DBond[j].atom->occupancy>0.0)  &&
                 DBond[j].atom->isInSelection(
                                  selHandles1.selHndHydrogen))  {
              H = DBond[j].atom;
              if ((H->GetDist2(A)<maxHAdist2) &&
                  (H->GetCosine(D,A)<=maxDHAcos))  {
                // Check angles with all acceptor neighbours
                isHBond = true;
                for (k=0;(k<nABonds) && isHBond;k++)
                  isHBond = (A->GetCosine(H,ABond[k].atom)<=maxHAAcos);
                if (isHBond)  {
                  if (Contact[i].group<=1)
                        _addAtomPair ( H,A,HBond,nHBonds,nHBAlloc );
                  else  _addAtomPair ( A,H,HBond,nHBonds,nHBAlloc );
                }
              }
            }
          if ((!H) && (nDBonds>0))  {
            // There were no hydrogens bonded to donor, assume that
            // the structure is incomplete and check donor-acceptor
            // geometry for possible h-bonding.
            isHBond = true;
            for (j=0;(j<nDBonds) && isHBond;j++)
              isHBond = (D->GetCosine(DBond[j].atom,A)<=maxDDAcos);
            for (j=0;(j<nABonds) && isHBond;j++)
              isHBond = (A->GetCosine(D,ABond[j].atom)<=maxDAAcos);
            if (isHBond)  {
              if (Contact[i].group<=1)
                    _addAtomPair ( D,A,HBond,nHBonds,nHBAlloc );
              else  _addAtomPair ( A,D,HBond,nHBonds,nHBAlloc );
            }
          }
        }
      }

    }

    if (Contact)  delete[] Contact;

    selHandles1.deleteSelections ( MMDB );
    selHandles2.deleteSelections ( MMDB );

    return CCP4SRS_Ok;

  }

}  // namespace ccp4srs
