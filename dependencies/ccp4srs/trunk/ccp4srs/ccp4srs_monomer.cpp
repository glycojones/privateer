//  $Id: ccp4srs_monomer.cpp $
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
//  **** Module  :  ccp4srs_monomer  <implementation>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Monomer  - monomer structure
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2014
//
//  =================================================================
//

#include <string.h>

#include "ccp4srs_monomer.h"
#include "ccp4srs_defs.h"

#include "mmdb2/mmdb_math_linalg.h"
#include "mmdb2/mmdb_tables.h"

namespace ccp4srs  {

  Monomer::Monomer() : Entry()  {
    Init();
  }

  Monomer::~Monomer()  {
    reset();
  }

  void Monomer::Init()  {
    monID   [0] = char(0);
    oldMonID[0] = char(0);
    monName     = NULL;
    monType     = NULL;
    monFormula  = NULL;
  }

  void Monomer::reset()  {

    atoms      .empty();
    bonds      .empty();
    angles     .empty();
    torsions   .empty();
    chicenters .empty();
    planes     .empty();
    struct_tree.empty();

    if (monName)     delete[] monName;
    if (monType)     delete[] monType;
    if (monFormula)  delete[] monFormula;

    Init();

  }

  void Monomer::setOldID ( mmdb::cpstr oldID )  {
    strcpy ( oldMonID,oldID );
  }

  mmdb::cpstr  Monomer::chem_name()  {
    if (monName)  return monName;
    return "";
  }

  mmdb::cpstr  Monomer::chem_type()  {
    if (monType)  return monType;
    return "";
  }

  mmdb::cpstr  Monomer::chem_formula()  {
    if (monFormula)  return monFormula;
    return "";
  }

  void Monomer::getLeavingAtoms ( int & nLeavingAtoms,
                                         mmdb::ivector & leavingAtom,
                                         mmdb::ivector & bondedAtom )  {
  //  should be called after all atoms and bonds are set up
  int i,j;

    mmdb::FreeVectorMemory ( leavingAtom,0 );
    mmdb::FreeVectorMemory ( bondedAtom ,0 );

    nLeavingAtoms = 0;
    for (i=0;i<atoms.numberOf();i++)
      if (atom(i)->isLeaving())
        nLeavingAtoms++;

    if (nLeavingAtoms>0)  {
      mmdb::GetVectorMemory ( leavingAtom,nLeavingAtoms,0 );
      mmdb::GetVectorMemory ( bondedAtom ,nLeavingAtoms,0 );
      j = 0;
      for (i=0;i<atoms.numberOf();i++)
        if (atom(i)->isLeaving())  {
          leavingAtom[j] =  i;
          bondedAtom [j] = -1;
          j++;
        }
      for (i=0;i<nLeavingAtoms;i++)
        for (j=0;j<bonds.numberOf();j++)  {
          if ((leavingAtom[i]==bond(j)->atom1()) &&
              (strcmp(atom(bond(j)->atom2())->element(),"H")))
            bondedAtom[i] = bond(j)->atom2();
          else if ((leavingAtom[i]==bond(j)->atom2()) &&
              (strcmp(atom(bond(j)->atom1())->element(),"H")))
            bondedAtom[i] = bond(j)->atom1();
        }
    }

  }


  int  superpose_atoms ( mmdb::mat44 & T, PPAtom A1,
                  mmdb::PPAtom A2, int nAtoms,
                  mmdb::rmatrix & A, mmdb::rmatrix & U, mmdb::rmatrix & V,
                  mmdb::rvector & W, mmdb::rvector & RV1 )  {
  //   Given two sets of atoms, A1 and A2, superpose_atoms(...)
  // calculates the rotational-translational matrix T such that
  // |T*A1 - A2| is minimal in least-square terms. The transfomation
  // superposes exactly the atoms A1[0] and A2[0].
  mmdb::realtype det,B;
  mmdb::vect3    vc1,vc2;
  int            i,j,k;

    //  1.  Calculate the correlation matrix. The rotation will be
    //      done around

    for (i=1;i<=3;i++)
      for (j=1;j<=3;j++)
        A[i][j] = 0.0;

    for (k=1;k<nAtoms;k++)  {
      vc1[0] = A1[k]->x() - A1[0]->x();
      vc1[1] = A1[k]->y() - A1[0]->y();
      vc1[2] = A1[k]->z() - A1[0]->z();
      vc2[0] = A2[k]->x - A2[0]->x;
      vc2[1] = A2[k]->y - A2[0]->y;
      vc2[2] = A2[k]->z - A2[0]->z;
      for (i=1;i<=3;i++)
        for (j=1;j<=3;j++)
          A[i][j] += vc1[j-1]*vc2[i-1];
    }

    //  2. Calculate transformation matrix (to be applied to A1)

    det = A[1][1]*A[2][2]*A[3][3] +
          A[1][2]*A[2][3]*A[3][1] +
          A[2][1]*A[3][2]*A[1][3] -
          A[1][3]*A[2][2]*A[3][1] -
          A[1][1]*A[2][3]*A[3][2] -
          A[3][3]*A[1][2]*A[2][1];

    //  2.1 SV-decompose the correlation matrix

    mmdb::math::SVD ( 3,3,3,A,U,V,W,RV1,true,true,i );

    if (i!=0)  return CCP4SRS_SVD_Fail;

    //  2.2 Check for parasite inversion and fix it if found

    if (det<=0.0)  {
      k = 0;
      B = mmdb::MaxReal;
      for (j=1;j<=3;j++)
        if (W[j]<B)  {
          B = W[j];
          k = j;
        }
      for (j=1;j<=3;j++)
        V[j][k] = -V[j][k];
    }

    //  2.3 Calculate rotational part of T

    for (j=1;j<=3;j++)
      for (k=1;k<=3;k++)  {
        B = 0.0;
        for (i=1;i<=3;i++)
          B += U[j][i]*V[k][i];
        T[j-1][k-1] = B;
      }

    //  2.4 Add translational part to T

    T[0][3] = A2[0]->x - T[0][0]*A1[0]->x() - T[0][1]*A1[0]->y() -
                         T[0][2]*A1[0]->z();
    T[1][3] = A2[0]->y - T[1][0]*A1[0]->x() - T[1][1]*A1[0]->y() -
                         T[1][2]*A1[0]->z();
    T[2][3] = A2[0]->z - T[2][0]*A1[0]->x() - T[2][1]*A1[0]->y() -
                         T[2][2]*A1[0]->z();

    return CCP4SRS_Ok;

  }

  int Monomer::atom_no ( mmdb::cpstr aname )  {
  // returns number of atom with given name [0...], or -1 if atom
  // is not found
  int i,n;

    n = -1;
    for (i=0;(i<atoms.numberOf()) && (n<0);i++)
      if (!strcmp(atoms.at(i)->name(),aname))
        n = i;

    return n;

  }

  int Monomer::get_bound ( int atomNo, mmdb::cpstr elem, int serNo ) {
  // returns number of atom with chemical element 'elem', bound to
  // atom number 'atomNo' [0,1,..]. serNo==0 returns 1st atom found,
  // serNo==1 returns second atom and so on. If atom not found, -1 is
  // returned.
  int i,j,n, atm1,atm2;

    n = -1;
    j = 0;
    for (i=0;(i<bonds.numberOf()) && (n<0);i++)  {
      atm2 = bonds.at(i)->atom2();
      if (atm2==atomNo)  {
        atm1 = atm2;
        atm2 = bonds.at(i)->atom1();
      } else
        atm1 = bonds.at(i)->atom1();
      if (atm1==atomNo)  {
        if (!strcmp(atoms.at(atm2)->element(),elem))  {
          if (j==serNo)  n = atm2;
                   else  j++;
        }
      }
    }

    return n;

  }

  int Monomer::get_bound ( int atomNo, int n ) {
  // returns nth [0,1,..] atom bound to atom with name 'aname'
  int  i,j,k,atm1,atm2;

    k = -1;
    j = 0;
    for (i=0;(i<bonds.numberOf()) && (k<0);i++)  {
      atm1 = bonds.at(i)->atom1();
      atm2 = bonds.at(i)->atom2();
      if (atm1==atomNo)  {
        if (j==n) k = atm2;
             else j++;
      } else if (atm2==atomNo)  {
        if (j==n) k = atm1;
             else j++;
      }
    }

    return k;

  }

  int Monomer::get_bound ( mmdb::cpstr aname, int n ) {
  // returns nth [0,1,..] atom bound to atom with name 'aname'
  int  i,j,k,atm1,atm2;

    k = -1;
    j = 0;
    for (i=0;(i<bonds.numberOf()) && (k<0);i++)  {
      atm1 = bonds.at(i)->atom1();
      atm2 = bonds.at(i)->atom2();
      if (!strcmp(atoms.at(atm1)->name(),aname))  {
        if (j==n) k = atm2;
             else j++;
      } else if (!strcmp(atoms.at(atm2)->name(),aname))  {
        if (j==n) k = atm1;
             else j++;
      }
    }

    return k;

  }


  CCP4SRS_RC Monomer::addHydrogens ( mmdb::PResidue R )  {
  //
  //   Return:
  //     _Ok             success
  //     _EmptyResidue   residue R does not contain atoms
  //     _NoAtomsFound   monomer does not contain atoms
  //     _NoBonds        monomer does not contain bonds
  //     _NoAtomsData    monomer is not complete
  //     _NoSimilarity   too few coomon atom names in R and SRS
  //                            entry with the same structure name
  //     _SuperpositionFailed  failed residue superposition
  //
  // NOTE1: the function does not rearranges existing atoms in the
  // residue, but places the hydrogens on top of them (leaving the
  // Ter pseudoatom, if found, on top of the list).
  //
  // NOTE2: in case of alternative locations, the first one in the
  // residue is chosen.
  //
  mmdb::PPAtom   H;
  PAtom          sa1[10];
  PAtom          atomi;
  mmdb::PAtom    sa2[10];
  mmdb::AtomName aname;
  mmdb::rmatrix  D,U,V;
  mmdb::rvector  W,RV1;
  mmdb::imatrix  c;
  mmdb::ivector  A,nb;
  mmdb::mat44    T;
  int            natoms,nbonds, i,j,k,m,n,mm,nTer,nH;

    //  1.  Make simple checks first

    natoms = atoms.numberOf();
    nbonds = bonds.numberOf();

    if (natoms<=2)          return CCP4SRS_NoAtomsFound;
    if (nbonds<=0)          return CCP4SRS_NoBonds;
    if (R->nAtoms<=2)       return CCP4SRS_EmptyResidue;
    if (natoms==R->nAtoms)  return CCP4SRS_Ok;


    //  2.  Map existing atoms from the residue onto a local array

    mmdb::GetVectorMemory ( A,natoms,0 );
    nTer = 20000;
    for (i=0;i<natoms;i++)  {
      A[i] = -1;  // signal "no atom"
      for (j=0;(j<R->nAtoms) && (A[i]<0);j++)
        if (R->atom[j])  {
          if (R->atom[j]->Ter)
            nTer = j;
          else if (!strcmp(atoms.at(i)->name_pdb(aname),R->atom[j]->name))
            A[i] = j;  // here is the place to check for altlocs
        }
    }


    //  3.  Make bond matrix

    mmdb::GetMatrixMemory ( c ,natoms,10,0,0 );
    mmdb::GetVectorMemory ( nb,natoms,0 );

    for (i=0;i<natoms;i++)
      nb[i] = 0;  // number of bonds at ith atom

    for (i=0;i<nbonds;i++)  {
      j = bonds.at(i)->atom1();
      k = bonds.at(i)->atom2();
      c[j][nb[j]] = k;
      c[k][nb[k]] = j;
      nb[j]++;
      nb[k]++;
    }


    //  4.  Loop over all hydrogens. Locate core atoms bonded to
    //      hydrogen in SBStructure and superpose them with the
    //      corresponding atoms in the residue. Using the superposition
    //      matrix, add hydrogens to the residue.

    mmdb::GetMatrixMemory ( D  ,3,3,1,1 );
    mmdb::GetMatrixMemory ( U  ,3,3,1,1 );
    mmdb::GetMatrixMemory ( V  ,3,3,1,1 );
    mmdb::GetVectorMemory ( W  ,3,1 );
    mmdb::GetVectorMemory ( RV1,3,1 );

    H  = new mmdb::PAtom[natoms];
    nH = 0;

    for (i=0;i<natoms;i++)  {
      atomi = atoms.at(i);
      if ((!strcmp(atomi->element(),"H")) && (A[i]<0) && (nb[i]>0) &&
          (atomi->x()>-mmdb::MaxShortReal/2.0)) {
        // ith atom is a hydrogen which is not found in the residue.
        // Find 3+ core atoms that are most closely bonded to this one.
        n = c[i][0];  // core atom bonded to the hydrogen
        if (A[n]>=0)  {
          sa1[0] = atoms.at(n);
          sa2[0] = R->atom[A[n]];
          mm = 1;
          m  = 1;
          while ((m<3) && (mm<2))  {
            k = n;
            for (j=0;(j<nb[k]) && (m<10);j++)  {
              n = c[k][j];
              if (A[n]>=0)  {
                sa1[m] = atoms.at(n);
                sa2[m] = R->atom[A[n]];
                m++;
              }
            }
            mm++;
          }
          if (m>=3)  {
            // superpose atoms and add the hydrogen to the residue
            k = superpose_atoms ( T,sa1,sa2,m,D,U,V,W,RV1 );
            if (k==CCP4SRS_Ok)  {
              H[nH] = atomi->makeAtom();
              H[nH]->Transform ( T );
              nH++;
            }
          }
        }
      }
    }

    //  5.  Put hydrogens into the residue

    for (i=0;i<nH;i++)  {
      R->InsertAtom ( H[i],nTer );
      nTer++;
    }


    //  6.  Release memory and return

    if (H)  delete[] H;
    mmdb::FreeVectorMemory ( A,0 );

    mmdb::FreeVectorMemory ( RV1,1 );
    mmdb::FreeVectorMemory ( W  ,1 );
    mmdb::FreeMatrixMemory ( V  ,3,1,1 );
    mmdb::FreeMatrixMemory ( U  ,3,1,1 );
    mmdb::FreeMatrixMemory ( D  ,3,1,1 );

    mmdb::FreeVectorMemory ( nb,0 );
    mmdb::FreeMatrixMemory ( c ,natoms,0,0 );

    return CCP4SRS_Ok;

  }


  void Monomer::getAtomNameMatch ( mmdb::PPAtom  A,
                                   int           nat,
                                   mmdb::pstr    altLoc,
                                   mmdb::ivector anmatch )  {
  mmdb::AtomName anamei;
  int            i,j,k, natoms;
  bool           done;

    natoms = atoms.numberOf();
    for (i=0;i<natoms;i++)  {
      k    = -1;
      j    = 0;
      done = false;
      while ((j<nat) && (!done))  {
        if (A[j])  {
          atoms.at(i)->name_pdb ( anamei );
          if ((!A[j]->Ter) && (!strcmp(A[j]->name,anamei)))  {
            k = j;  // atom found
            // check now for altLocs
            j++;
            while ((j<nat) && (!done))  {
              if (A[j])  {
                done = (A[j]->Ter || (strcmp(A[j]->name,anamei)));
                if (!done)  {
                  if (A[j]->occupancy>A[k]->occupancy)  k = j;
                  if (altLoc)  {
                    if (!strcmp(A[j]->altLoc,altLoc))  {
                      k    = j;
                      done = true;
                    }
                  }
                }
              }
              j++;
            }
          }
        }
        j++;
      }
      anmatch[i] = k;
    }

  }

  void Monomer::getAtomNameMatch ( Container<Atom> *A,
                                   mmdb::ivector anmatch )  {
  int i,j,n,nat;

    n   = A->numberOf();
    nat = atoms.numberOf();
    for (i=0;i<n;i++)  {
      anmatch[i] = -1;
      for (j=0;(j<nat) && (anmatch[i]<0);j++)
        if (!strcmp(A->at(i)->id(),atoms.at(j)->id()))
          anmatch[i] = j;
    }

  }


  void Monomer::readMonIDFromCIF ( mmdb::mmcif::PLoop loop )  {
  mmdb::pstr p;
  int        rc;
    if (loop && (!monID[0]))  {
      p  = loop->GetString ( MMCIF_TAG_COMP_ID,0,rc );
      if (p) strcpy ( monID,p );
    }
  }

  void Monomer::readMonIDFromCIF ( mmdb::mmcif::PStruct mmCIFStruct )  {
  mmdb::pstr p;
  int        rc;
    if (mmCIFStruct && (!monID[0]))  {
      p  = mmCIFStruct->GetString ( MMCIF_TAG_COMP_ID,rc );
      if (p) strcpy ( monID,p );
    }
  }

  // Return code:
  //     _Ok    success
  int Monomer::readFromCIF ( mmdb::mmcif::PData mmCIFData,
                                    mmdb::pstr error_desc )  {
    reset();
    if (mmCIFData->GetStructure(MMCIF_STRUCT_CHEM_COMP))
          return readFromCIF_rcsb ( mmCIFData,error_desc );
    else  return readFromCIF_ccp4 ( mmCIFData,error_desc );
  }

  int Monomer::readFromCIF_rcsb ( mmdb::mmcif::PData mmCIFData,
                                         mmdb::pstr error_desc )  {
  mmdb::mmcif::PStruct mmCIFStruct;
  mmdb::mmcif::PLoop   mmCIFLoop;
  PAtom                atom;
  PBond                bond;
  int                  rc,loop_len,i;

    rc = 0;

    mmCIFStruct = mmCIFData->GetStructure ( MMCIF_STRUCT_CHEM_COMP );
    if (mmCIFStruct)  {
      mmdb::CreateCopy ( monName,mmCIFStruct->GetString(MMCIF_TAG_Name,rc) );
      mmdb::CreateCopy ( monType,mmCIFStruct->GetString(MMCIF_TAG_Type,rc) );
      mmdb::CreateCopy ( monFormula,mmCIFStruct->GetString(MMCIF_TAG_Formula,
                                                     rc) );
    }

    // take atoms
    mmCIFLoop = mmCIFData->GetLoop ( MMCIF_LOOP_CHEM_COMP_ATOM );
    if (!mmCIFLoop)
      mmCIFStruct = mmCIFData->GetStructure ( MMCIF_LOOP_CHEM_COMP_ATOM );

    if (mmCIFLoop)  {

      readMonIDFromCIF ( mmCIFLoop );

      loop_len = mmCIFLoop->GetLoopLength();
      rc       = 0;
      for (i=0;(i<loop_len) && (!rc);i++)  {
        atom = new Atom();
        rc   = atom->readFromCIF_rcsb ( mmCIFLoop,i );
        if (!rc)  atoms.add ( atom );
        else  {
          if (error_desc)
            sprintf ( error_desc,"Error %i in row %i of %s loop",
                      rc,i+1,MMCIF_LOOP_CHEM_COMP_ATOM );
          delete atom;
        }
      }

    } else if (mmCIFStruct)  {

      readMonIDFromCIF ( mmCIFStruct );

      atom = new Atom();
      rc = atom->readFromCIF ( mmCIFStruct );
      if (!rc)  atoms.add ( atom );
      else  {
        if (error_desc)
          sprintf ( error_desc,"Error %i in %s structure",
                    rc,MMCIF_LOOP_CHEM_COMP_ATOM );
        delete atom;
      }

    }

    if (rc)  {
      reset();
      return rc;
    }


    // take bonds
    mmCIFLoop =  mmCIFData->GetLoop ( MMCIF_LOOP_CHEM_COMP_BOND );
    if (!mmCIFLoop)
      mmCIFStruct = mmCIFData->GetStructure ( MMCIF_LOOP_CHEM_COMP_BOND );

    if (mmCIFLoop)  {

      readMonIDFromCIF ( mmCIFLoop );

      loop_len = mmCIFLoop->GetLoopLength();
      rc       = 0;
      for (i=0;(i<loop_len) && (!rc);i++)  {
        bond = new Bond();
        rc = bond->readFromCIF_rcsb ( mmCIFLoop,i,atoms );
        if (!rc)  bonds.add ( bond );
        else  {
          if (error_desc)
            sprintf ( error_desc,"Error %i in row %i of %s loop",
                      rc,i+1,MMCIF_LOOP_CHEM_COMP_BOND );
          delete bond;
        }
      }

    } else if (mmCIFStruct)  {

      readMonIDFromCIF ( mmCIFStruct );

      bond = new Bond();
      rc = bond->readFromCIF ( mmCIFStruct,atoms );
      if (!rc)  bonds.add ( bond );
      else  {
        if (error_desc)
          sprintf ( error_desc,"Error %i in %s structure",
                    rc,MMCIF_LOOP_CHEM_COMP_BOND );
        delete bond;
      }

    }

    // take bonds
    mmCIFLoop =  mmCIFData->GetLoop ( MMCIF_LOOP_CHEM_COMP_DESCRIPTOR );
    if (mmCIFLoop)  {
      ACDLabs.readFromCIF_rcsb ( mmCIFLoop,"ACDLabs" );
      CACTVS .readFromCIF_rcsb ( mmCIFLoop,"CACTVS"  );
      OpenEye.readFromCIF_rcsb ( mmCIFLoop,"OpenEye OEToolkits" );
      InChI  .readFromCIF_rcsb ( mmCIFLoop,"InChI"   );
    } else  {
      ACDLabs.freeSmiles();
      CACTVS .freeSmiles();
      OpenEye.freeSmiles();
      InChI  .freeSmiles();
    }

/*
loop_
_pdbx_chem_comp_descriptor.comp_id
_pdbx_chem_comp_descriptor.type
_pdbx_chem_comp_descriptor.program
_pdbx_chem_comp_descriptor.program_version
_pdbx_chem_comp_descriptor.descriptor
000 SMILES           ACDLabs              12.01 "O=C(O)OC"
000 SMILES_CANONICAL CACTVS               3.370 "COC(O)=O"
000 SMILES           CACTVS               3.370 "COC(O)=O"
000 SMILES_CANONICAL "OpenEye OEToolkits" 1.7.0 "COC(=O)O"
000 SMILES           "OpenEye OEToolkits" 1.7.0 "COC(=O)O"
000 InChI            InChI                1.03  "InChI=1S/C2H4O3/c1-5-2(3)4/h1H3,(H,3,4)"
000 InChIKey         InChI                1.03  CXHHBNMLPJOKQD-UHFFFAOYSA-N
*/

    if (rc)  {
      reset();
      return rc;
    }

    return CCP4SRS_Ok;

  }

  int Monomer::readFromCIF_ccp4 ( mmdb::mmcif::PData mmCIFData,
                                  mmdb::pstr error_desc )  {
  mmdb::mmcif::PLoop loop;
  PAtom       atom;
  PBond       bond;
  PAngle      angle;
  PChiCenter  chicenter;
  PTorsion    torsion;
  PPlane      plane;
  int                loop_len, i,n, rc;

    // take atoms
    loop =  mmCIFData->GetLoop ( MMCIF_LOOP_CHEM_COMP_ATOM );
    readMonIDFromCIF ( loop );

    if (loop)  {

      loop_len = loop->GetLoopLength();
      rc       = 0;
      for (i=0;(i<loop_len) && (!rc);i++)  {
        atom = new Atom();
        rc = atom->readFromCIF_ccp4 ( loop,i );
        if (!rc)  atoms.add ( atom );
        else  {
          if (error_desc)
            sprintf ( error_desc,"Error %i in row %i of %s loop",
                      rc,i+1,MMCIF_LOOP_CHEM_COMP_ATOM );
          delete atom;
        }
      }

      if (rc)  {
        reset();
        return rc;
      }

    }

    // take tree records
    loop =  mmCIFData->GetLoop ( MMCIF_LOOP_CHEM_COMP_TREE );
    readMonIDFromCIF ( loop );
    rc = struct_tree.readFromCIF ( loop,atoms );

    if (rc)  {
      reset();
      return rc;
    }

    // take bonds
    loop =  mmCIFData->GetLoop ( MMCIF_LOOP_CHEM_COMP_BOND );
    readMonIDFromCIF ( loop );

    if (loop)  {

      loop_len = loop->GetLoopLength();
      rc       = 0;
      for (i=0;(i<loop_len) && (!rc);i++)  {
        bond = new Bond();
        rc = bond->readFromCIF_ccp4 ( loop,i,atoms );
        if (!rc)  bonds.add ( bond );
        else  {
          if (error_desc)
            sprintf ( error_desc,"Error %i in row %i of %s loop",
                      rc,i+1,MMCIF_LOOP_CHEM_COMP_BOND );
          delete bond;
        }
        /*
        if ((!rc) || (rc==_WrongBondType))  {
          bonds.add ( bond );
          rc = 0;
        } else  {
          if (error_desc)
            sprintf ( error_desc,"Error %i in row %i of %s loop",
                      rc,i+1,MMCIF_LOOP_CHEM_COMP_BOND );
          delete bond;
        }
        */
      }

      if (rc)  {
        reset();
        return rc;
      }

    }

    // take angles
    loop =  mmCIFData->GetLoop ( MMCIF_LOOP_CHEM_COMP_ANGLE );
    readMonIDFromCIF ( loop );

    if (loop)  {

      loop_len = loop->GetLoopLength();
      rc       = 0;
      for (i=0;(i<loop_len) && (!rc);i++)  {
        angle = new Angle();
        rc = angle->readFromCIF ( loop,i,atoms );
        if (!rc)  angles.add ( angle );
        else  {
          if (error_desc)
            sprintf ( error_desc,"Error %i in row %i of %s loop",
                      rc,i+1,MMCIF_LOOP_CHEM_COMP_ANGLE );
          delete angle;
        }
      }

      if (rc)  {
        reset();
        return rc;
      }

    }

    // take chiral centers
    loop =  mmCIFData->GetLoop ( MMCIF_LOOP_CHEM_COMP_CHIR );
    readMonIDFromCIF ( loop );

    if (loop)  {

      // reset atom chirality signs
      for (i=0;i<atoms.numberOf();i++)
        atoms.at(i)->set_chirality ( 'N' );

      loop_len = loop->GetLoopLength();
      rc       = 0;
      for (i=0;(i<loop_len) && (!rc);i++)  {
        chicenter = new ChiCenter();
        rc = chicenter->readFromCIF ( loop,i,atoms );
        if (!rc)  {
          chicenters.add ( chicenter );
          atoms.at(chicenter->center())->set_chirality (
                                            chicenter->get_chirality() );
        } else  {
          if (error_desc)
            sprintf ( error_desc,"Error %i in row %i of %s loop",
                      rc,i+1,MMCIF_LOOP_CHEM_COMP_CHIR );
          delete chicenter;
        }
      }

      if (rc)  {
        reset();
        return rc;
      }

    }


    // take torsions
    loop =  mmCIFData->GetLoop ( MMCIF_LOOP_CHEM_COMP_TOR );
    readMonIDFromCIF ( loop );

    if (loop)  {

      loop_len = loop->GetLoopLength();
      rc       = 0;
      for (i=0;(i<loop_len) && (!rc);i++)  {
        torsion = new Torsion();
        rc = torsion->readFromCIF ( loop,i,atoms );
        if (!rc)  torsions.add ( torsion );
        else  {
          if (error_desc)
            sprintf ( error_desc,"Error %i in row %i of %s loop",
                      rc,i+1,MMCIF_LOOP_CHEM_COMP_TOR );
          delete torsion;
        }
      }

      if (rc)  {
        reset();
        return rc;
      }

    }


    // take planes
    loop =  mmCIFData->GetLoop ( MMCIF_LOOP_CHEM_COMP_PLANE_ATOM );
    readMonIDFromCIF ( loop );

    if (loop)  {

      n  = 0;
      rc = Plane::Ok;
      while (rc==Plane::Ok)  {
        plane = new Plane();
        rc = plane->readFromCIF ( loop,n,atoms );
        if (rc==Plane::Ok)  {
          planes.add ( plane );
          n++;
        } else  {
          if (error_desc && (rc!=Plane::NoMorePlanes))
            sprintf ( error_desc,"Error %i in plane %s of %s loop",
                      rc,plane->id(),MMCIF_LOOP_CHEM_COMP_PLANE_ATOM );
          delete plane;
        }
      }

      if (rc!=Plane::NoMorePlanes)  {
        reset();
        return rc;
      }

    }

    return CCP4SRS_Ok;

  }


  mmdb::mmcif::PData Monomer::makeCIF()  {
  mmdb::mmcif::PData data;
  mmdb::mmcif::PLoop loop;
  char        S[100];
  int         i;

    strcpy ( S,MMCIF_COMP );
    strcat ( S,monID );
    data = new mmdb::mmcif::Data ( S );

    data->AddLoop ( MMCIF_LOOP_CHEM_COMP_ATOM,loop );
    Atom::makeCIFTags ( loop );
    for (i=0;i<atoms.numberOf();i++)
      atoms.at(i)->writeToCIF ( loop,monID );

    struct_tree.writeToCIF ( data,monID,atoms );

    data->AddLoop ( MMCIF_LOOP_CHEM_COMP_BOND,loop );
    Bond::makeCIFTags ( loop );
    for (i=0;i<bonds.numberOf();i++)
      bonds.at(i)->writeToCIF ( loop,monID,atoms );

    data->AddLoop ( MMCIF_LOOP_CHEM_COMP_ANGLE,loop );
    Angle::makeCIFTags ( loop );
    for (i=0;i<angles.numberOf();i++)
      angles.at(i)->writeToCIF ( loop,monID,atoms );

    data->AddLoop ( MMCIF_LOOP_CHEM_COMP_CHIR,loop );
    ChiCenter::makeCIFTags ( loop );
    for (i=0;i<chicenters.numberOf();i++)
      chicenters.at(i)->writeToCIF ( loop,monID,atoms );

    data->AddLoop ( MMCIF_LOOP_CHEM_COMP_TOR,loop );
    Torsion::makeCIFTags ( loop );
    for (i=0;i<torsions.numberOf();i++)
      torsions.at(i)->writeToCIF ( loop,monID,atoms );

    data->AddLoop ( MMCIF_LOOP_CHEM_COMP_PLANE_ATOM,loop );
    Plane::makeCIFTags ( loop );
    for (i=0;i<planes.numberOf();i++)
      planes.at(i)->writeToCIF ( loop,monID,atoms );

    return data;

  }

  void Monomer::write_mem ( PMemIO memIO, int version )  {
    memIO->put_line   ( monID      );
    memIO->put_line   ( oldMonID   );
    memIO->put_string ( monName    );
    memIO->put_string ( monType    );
    memIO->put_string ( monFormula );
    atoms      .write_mem ( memIO,version );
    bonds      .write_mem ( memIO,version );
    angles     .write_mem ( memIO,version );
    torsions   .write_mem ( memIO,version );
    chicenters .write_mem ( memIO,version );
    planes     .write_mem ( memIO,version );
    struct_tree.write_mem ( memIO,version );
    ACDLabs    .write_mem ( memIO,version );
    CACTVS     .write_mem ( memIO,version );
    OpenEye    .write_mem ( memIO,version );
    InChI      .write_mem ( memIO,version );
  }

  bool Monomer::read_mem ( PMemIO memIO, int version, bool * Ok )  {
  bool success;
    if (Ok)  success = *Ok;
       else  success = true;
    memIO->get_line   ( monID     ,&success );
    memIO->get_line   ( oldMonID  ,&success );
    memIO->get_string ( monName   ,&success );
    memIO->get_string ( monType   ,&success );
    memIO->get_string ( monFormula,&success );
    atoms      .read_mem ( memIO,version,&success );
    bonds      .read_mem ( memIO,version,&success );
    angles     .read_mem ( memIO,version,&success );
    torsions   .read_mem ( memIO,version,&success );
    chicenters .read_mem ( memIO,version,&success );
    planes     .read_mem ( memIO,version,&success );
    struct_tree.read_mem ( memIO,version,&success );
    ACDLabs    .read_mem ( memIO,version,&success );
    CACTVS     .read_mem ( memIO,version,&success );
    OpenEye    .read_mem ( memIO,version,&success );
    InChI      .read_mem ( memIO,version,&success );

    if (Ok)  *Ok = success;
    return success;
  }


  int makeElementType ( int ElType, char chirality )  {
    if (chirality=='S')  return ElType | mmdb::math::CHIRAL_LEFT;
    if (chirality=='R')  return ElType | mmdb::math::CHIRAL_RIGHT;
    return ElType;
  }


  mmdb::math::PGraph Monomer::getGraph ( int *retCode )  {
  // Generate chemical graph
  mmdb::math::PGraph  graph;
  mmdb::math::PVertex vertex;
  mmdb::math::PEdge   edge;
  PAtom        a;
  PBond        b;
  int                 i,n,rc,vtype,vx1,vx2;

    graph = new mmdb::math::Graph();
    graph->SetName ( ID() );
    rc = 0;

    n = n_atoms();
    for (i=0;(i<n) && (!rc);i++)  {

      a = atom(i);

      vtype = mmdb::getElementNo ( a->element() );
      if (vtype!=mmdb::ELEMENT_UNKNOWN)  {

        if (a->leaving()=='Y')
          vtype |= mmdb::math::ATOM_LEAVING;

        vertex = new mmdb::math::Vertex (
                     makeElementType(vtype,a->chirality()) );
        vertex->SetName   ( a->name() );
        vertex->SetUserID ( i );
        graph->AddVertex  ( vertex );

      } else
        rc = 10000 + i;

    }

    n = n_bonds();
    for (i=0;(i<n) && (!rc);i++)  {

      b   = bond(i);
      vx1 = b->atom1();
      vx2 = b->atom2();

      if (vx1<0)  rc = 20000 + i;
      if (vx2<0)  rc = 30000 + i;

      vx1++;
      vx2++;

      if (!rc)  {
        edge = new mmdb::math::Edge ( vx1,vx2,b->order() );
        graph->AddEdge ( edge );
      }

    }

    if (rc)  {
      delete graph;
      graph = NULL;
    }

    if (retCode)  *retCode = rc;

    return graph;

  }

  // service function, assume merging rcsb monomer with monlib one
  void Monomer::merge ( PMonomer monomer )  {
  PAtom            tatom;
  PBond            mbond,tbond;
  PAngle           angle;
  PPlane           plane;
  PTorsion         torsion;
  PChiCenter       chicenter;
  Container<Atom> *matoms;
  mmdb::AtomName   aname;
  mmdb::ivector    anmatch;
  int              i,j,k;

    if (!monID[0])  strcpy ( monID,monomer->ID() );
    else if (strcmp(monID,monomer->ID()))
      return;  // can't merge dufferent monomers

    if (!oldMonID[0])  strcpy ( oldMonID,monomer->oldID() );

    if (!chem_name()[0])  mmdb::CreateCopy ( monName,monomer->chem_name() );
    if (!chem_type()[0])  mmdb::CreateCopy ( monType,monomer->chem_type() );
    if (!chem_formula()[0])
                    mmdb::CreateCopy ( monFormula,monomer->chem_formula() );
    if (ACDLabs.isEmpty())  ACDLabs.copy ( monomer->getACDLabs() );
    if (CACTVS .isEmpty())  CACTVS .copy ( monomer->getCACTVS () );
    if (OpenEye.isEmpty())  OpenEye.copy ( monomer->getOpenEye() );
    if (InChI  .isEmpty())  InChI  .copy ( monomer->getInChI  () );

    for (i=0;i<monomer->n_atoms();i++)  {
      k = -1;
      strcpy ( aname,monomer->atom(i)->id() );
      for (j=0;(j<atoms.numberOf()) && (k<0);j++)
        if (!strcmp(atoms.at(j)->id(),aname))
          k = j;
      if (k<0)  {
        tatom = new Atom();
        tatom->copy ( monomer->atom(i) );
        add ( tatom );
      } else
        atoms.at(k)->merge ( monomer->atom(i) );
    }

    mmdb::GetVectorMemory ( anmatch,monomer->n_atoms(),0 );
    getAtomNameMatch ( monomer->get_atoms(),anmatch );

    matoms = monomer->get_atoms();
    for (i=0;i<monomer->n_bonds();i++)  {
      k     = -1;
      mbond = monomer->bond(i);
      for (j=0;(j<bonds.numberOf()) && (k<0);j++)
        if (mbond->compare(*matoms,bonds.at(j),atoms))
          k = j;
      if (k<0)  {
        tbond = new Bond();
        if (tbond->copy(atoms,mbond,*matoms))
              bonds.add ( tbond );
        else  delete tbond;
      } else if (bonds.at(k)->length()<=1.0e-10)
        bonds.at(k)->copy ( mbond,anmatch );
    }

    if (angles.numberOf()<=0)  {
      for (i=0;i<monomer->n_angles();i++)  {
        angle = new Angle();
        angle->copy ( monomer->angle(i),anmatch );
        angles.add ( angle );
      }
    }

    if (torsions.numberOf()<=0)  {
      for (i=0;i<monomer->n_torsions();i++)  {
        torsion = new Torsion();
        torsion->copy ( monomer->torsion(i),anmatch );
        torsions.add ( torsion );
      }
    }

    if (chicenters.numberOf()<=0)  {
      for (i=0;i<monomer->n_chicenters();i++)  {
        chicenter = new ChiCenter();
        chicenter->copy ( monomer->chicenter(i),anmatch );
        chicenters.add ( chicenter );
      }
    }

    if (planes.numberOf()<=0)  {
      for (i=0;i<monomer->n_planes();i++)  {
        plane = new Plane();
        plane->copy ( monomer->plane(i),anmatch );
        planes.add ( plane );
      }
    }

    if (struct_tree.size()<=0)
      struct_tree.copy ( monomer->tree(),anmatch );

    mmdb::FreeVectorMemory ( anmatch,0 );

  }

  void Monomer::check_ccp4_coordinates()  {
  int  i;
  bool ccp4_coors = false;

    for (i=0;(i<atoms.numberOf()) && (!ccp4_coors);i++)
      ccp4_coors =
          (atoms.at(i)->x_ccp4_mlib()>-mmdb::MaxShortReal/2.0) &&
          (atoms.at(i)->y_ccp4_mlib()>-mmdb::MaxShortReal/2.0) &&
          (atoms.at(i)->z_ccp4_mlib()>-mmdb::MaxShortReal/2.0);

    for (i=0;i<atoms.numberOf();i++)
      atoms.at(i)->set_ccp4_coordinates ( ccp4_coors );

  }


}  // namespace
