
/* ------------------------------------------------------------------

   17.04.12   <--  Date of Last Modification.

   --- Module      :  srsgen.cpp

   --- Purpose     :

   This program compiles graph representations of chemical compounds
   from CCP4 Monomer Library cifs. The output is arranged in 3 binary
   files:

   index.srs   the index. The file contains a sequence of
               CCP4SRSIndex structure (cf. ccp4srs_index.h)

   graph.srs   the graphs: arrays of vertices and edges.
               For each entry mentioned in index.srs, this file
               contains a CGraph structure (cf. mmdb_graph.h).
               Offset of this structure for a given compound, as
               well as values of nAtoms and nBonds are found in
               the CCP4SRSIndex, which is read from index.srs

   struct.srs  descriptions of chemical structures. For each
               entry mentioned in index.srs, this file contains
               a CCP4SRSMonomer structure (cf. ccp4srs_monomer.h).
               Offset of this structure for a given compound, as
               well as values of nAtoms and nBonds are found in
               the CCP4SRSIndex, which is read from index.srs


   --- Invokation     :

  srsgen rcsblib_dir monlib_dir [old_mon_lib] [-d srs_dir]


  rcsblib_dir is dirctory with RCSB's chemical compound library

  monlib_dir is CCP4's monomer library directory. In the original
  CCP4 distribution, it is found in $CCP4/lib/data/monomers

  srs_dir is an optional parameter specifying where SRS should be
  placed. If SRS directory is not specified, the current one will
  be used.


  E. Krissinel 2011


  ------------------------------------------------------------------- */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(__DEC) || defined(_WIN32)
#define USE_DIRENT_H 1
#endif

#if USE_DIRENT_H
 #include <dirent.h>
#else
 #include <sys/dir.h>
#endif


#include "mmdb2/mmdb_tables.h"
#include "mmdb2/mmdb_mmcif_.h"
#include "mmdb2/mmdb_math_graph.h"


#ifdef  _WIN32
//  for DOS/WINDOWS machines:
#define  _dir_sep    "\\"
#define  _dir_sep_c  '\\'
#else
//  for UNIX machines:
#define  _dir_sep    "/"
#define  _dir_sep_c  '/'
#endif

#include "ccp4srs/ccp4srs_monomer.h"
#include "ccp4srs/ccp4srs_index.h"
#include "ccp4srs/ccp4srs_defs.h"

/*
void  MakePDBName ( pstr pdb_name, pstr name,
                    int pdb_rcsb_align )  {
int  i,k;
  k = 0;
  while (k<pdb_rcsb_align)
    pdb_name[k++] = ' ';
  i = 0;
  while (name[i] && (k<4))
    pdb_name[k++] = name[i++];
  while (k<4)
    pdb_name[k++] = ' ';
  pdb_name[k] = '\0';
}
*/

/*
int GetBondOrder ( pstr border )  {
int etype;
  switch (border[0])  {
    case 's' : case 'S' : etype = BOND_SINGLE;   break;
    case 'd' : case 'D' : etype = BOND_DOUBLE;   break;
    case 'a' : case 'A' : etype = BOND_AROMATIC; break;
    case 't' : case 'T' : etype = BOND_TRIPLE;   break;
    default  : etype = -1;
  }
  return etype;
}
*/

//  --------  Calculating the Energy-Type Derivable Data  --------

/*
PCSBAngle GetBondAngle ( PCMMCIFLoop Loop,
                         pstr  at1, pstr  at2, pstr  at3,
                         int atom1, int atom2, int atom3 )  {
PCSBAngle angle;
pstr      et1,et2,et3;
int       i,l, a1,a3, found, rc;

  angle = NULL;
  l = Loop->GetLoopLength();
  found = 0;
  for (i=0;(i<l) && (found<2);i++)  {
    et2 = Loop->GetString ( "atom_type_2",i,rc );
    if (et2)  {
      if (!strcmp(et2,at2))  {
        et1 = Loop->GetString ( "atom_type_1",i,rc );
        et3 = Loop->GetString ( "atom_type_3",i,rc );
        if (et1)  {
          if (!strcmp(et1,at1))      a1 = 1;
          else if (!strcmp(et1,at3)) a1 = 3;
          else a1 = 0;
        } else
          a1 = 2;
        if (a1==0)  a3 = 0;
        else if (et3)  {
          if ((a1==1) && (!strcmp(et3,at3)))      a3 = 3;
          else if ((a1==3) && (!strcmp(et3,at1))) a3 = 1;
          else if (a1==2)  {
            if (!strcmp(et3,at1))      a3 = 1;
            else if (!strcmp(et3,at3)) a3 = 3;
            else a3 = 2;
          } else
            a3 = 0;
        } else
          a3 = 2;
        if (a3>0)  {
          if (!angle)  angle = new CSBAngle();
          angle->atom1 = atom1;
          angle->atom2 = atom2;
          angle->atom3 = atom3;
          if (Loop->GetReal(angle->angle,"value",i,False))
              angle->angle = 0.0;
          if (angle->angle>0.0)  {
            if ((a1!=2) && (a3!=2))  found = 2;
                               else  found = 1;
          }
        }
      }
    }
  }

  if (!found)
    printf ( " cannot find angle for energy types %s, %s and %s\n",
             at1,at2,at3 );

  return angle;

}



void FillETLibData ( PCSBStructure SBS, RCMMCIFData ETLib )  {
PCMMCIFLoop  Loop;
PCSBAtom     Atom;
PCSBBond     Bond,Bond2;
PCSBAngle    angle;
realtype     blength,blength_esd;
pstr         at1,at2, bt, et1,et2;
pstr         at21,at22;
int          i,j,l,a1,a2,order,btype, rc,found,nb,nb_esd;
int          a21,a22;

  // 1. Set atom properties depending on energy type

  Loop = ETLib.GetLoop ( "_lib_atom" );
  if (!Loop)  {
    printf ( " energy lib error:  _lib_atom not found\n" );
    return;
  }

  l = Loop->GetLoopLength();
  for (i=0;i<SBS->nAtoms;i++)  {
    Atom = SBS->Atom[i];
    at1  = Atom->energyType;
    a1   = -1;
    if (at1[0])  {
      for (j=0;(j<l) && (a1<0);j++)  {
        et1 = Loop->GetString ( "type",j,rc );
        if ((!rc) && et1)  {
          if (!strcmp(et1,at1))  a1 = j;
        }
      }
      if (a1>=0)  {
        if (Loop->GetReal(Atom->vdw_radius,"vdw_radius",a1,False))
          Atom->vdw_radius = getVdWaalsRadius  ( Atom->element );
        if (Loop->GetReal(Atom->vdwh_radius,"vdwh_radius",a1,False))
          Atom->vdw_radius = getVdWaalsRadius  ( Atom->element );
        if (Loop->GetReal(Atom->ion_radius,"ion_radius",a1,False))
          Atom->ion_radius  = getCovalentRadius ( Atom->element );
        if (Loop->GetInteger(Atom->valency,"valency",a1,False))
          Atom->valency = 0;
        bt = Loop->GetString ( "hb_type",a1,rc );
        if ((!rc) && bt)  {
          if (strpbrk(bt,"ADHBN"))  Atom->hb_type = bt[0];
                              else  Atom->hb_type = 'N';
// temporary?? terrible fix:
if ((!strcmp(SBS->compoundID,"DC")) &&
    (!strcmp(Atom->pdb_name," N3 ")))
  Atom->hb_type = 'A';

        }
      }
    }
    if (a1<0)  {
      // take whatever is possible from MMDB tables
      Atom->vdw_radius  = getVdWaalsRadius  ( Atom->element );
      Atom->vdwh_radius = Atom->vdw_radius + getVdWaalsRadius  ( "H" );
      Atom->ion_radius  = getCovalentRadius ( Atom->element );
      Atom->valency     = 0;    // valency is unknown
      Atom->hb_type     = 'N';  // no hydrogen bonding
    }
  }


  // 2. Fill in bond lengths

  Loop = ETLib.GetLoop ( "_lib_bond" );
  if (!Loop)  {
    printf ( " energy lib error:  _lib_bond not found\n" );
    return;
  }

  l = Loop->GetLoopLength();
  for (i=0;i<SBS->nBonds;i++)  {
    Bond = SBS->Bond[i];
    Bond->length     = 0.0;
    Bond->length_esd = 0.0;
    if (Bond->atom1>Bond->atom2)
      ISwap ( Bond->atom1,Bond->atom2 );
    at1 = SBS->Atom[Bond->atom1-1]->energyType;
    at2 = SBS->Atom[Bond->atom2-1]->energyType;
    if ((at1[0]) && (at2[0]))  {
      order       = Bond->order;
      blength     = 0.0;
      blength_esd = 0.0;
      nb          = 0;
      nb_esd      = 0;
      found       = 0;
      for (j=0;(j<l) && (found<2);j++)  {
        et1 = Loop->GetString ( "atom_type_1",j,rc );
        if (!rc)  {
          if (et1)  {
            if (!strcmp(et1,at1))      a1 = 1; // 1st atom identified
            else if (!strcmp(et1,at2)) a1 = 2; // 2nd atom identified
                                  else a1 = 0; // unsuitable pair in
                                               // ETLib
          } else
            a1 = 3;  // mmCIF '.' meaning any energy type for 1st atom
        } else
          a1 = 0;
        if (a1>0)  {
          et2 = Loop->GetString ( "atom_type_2",j,rc );
          if (!rc)  {
            if (et2)  {
              if ((a1==1) && (!strcmp(et2,at2)))      a2 = 2;
              else if ((a1==2) && (!strcmp(et2,at1))) a2 = 1;
              else if (a1==3)  {
                // first atom was identified as that of any energy type
                if (!strcmp(et2,at1))      a2 = 1;
                else if (!strcmp(et2,at2)) a2 = 2;
                                      else a2 = 0; // unsuitable pair
                                                   // in ETLib
              } else
                 a2 = 0;
            } else
              a2 = 3; // mmCIF '.' meaning any energy type for 2nd atom
            if (a2>0)  {
              bt = Loop->GetString ( "type",j,rc );
              if (!rc)  {
                if (bt)  {
                  btype = GetBondOrder ( bt );
                  if (btype!=order)  btype = -1;
                } else
                  btype = MaxInt;  // mmCIF '.' meaning any bond order
              } else
                btype = -MaxInt;
              if (btype>=0)  {
                if (!Loop->GetReal(blength,"length",j,False))  {
                  if ((a1<3) || (a2<3) || (nb==0))  {
                    nb++;                     // count
                    Bond->length += blength;  // and accumulate
                  }
                } else
                  blength = 0.0;  // unsuccessful retrieval
                if (!Loop->GetReal(blength_esd,"value_esd",j,False))  {
                  if ((a1<3) || (a2<3) || (nb_esd==0))  {
                    nb_esd++;                        // count
                    Bond->length_esd += blength_esd; // and accumulate
                  }
                } else  // unsuccessful retrieval
                  blength_esd = 0.0;
                if ((blength>0.0) && (btype<MaxInt) &&
                    (a2<3) && (a1<3))  {
                  found = 2;
                  Bond->length     = blength;
                  Bond->length_esd = blength_esd;
                } else
                  found = 1;
              }
            }
          }
        }
      }
      if (found<2)  {
        // average over incomplete energy type pairs
        if (nb>0)     Bond->length     /= nb;
        if (nb_esd>0) Bond->length_esd /= nb_esd;
      }
      if (!found)
        printf ( " can't find bond data for energy types %s and %s\n",
                 at1,at2 );
    }
  }


  // 3. Make angles

  Loop = ETLib.GetLoop ( "_lib_angle" );
  if (!Loop)  {
    printf ( " energy lib error:  _lib_angle not found\n" );
    return;
  }

  for (i=0;i<SBS->nBonds;i++)  {
    Bond = SBS->Bond[i];
    a1   = Bond->atom1;
    a2   = Bond->atom2;
    at1  = SBS->Atom[a1-1]->energyType;
    at2  = SBS->Atom[a2-1]->energyType;
    if ((at1[0]) && (at2[0]))
      for (j=i+1;j<SBS->nBonds;j++)  {
        Bond2 = SBS->Bond[j];
        a21   = Bond2->atom1;
        a22   = Bond2->atom2;
        at21  = SBS->Atom[a21-1]->energyType;
        at22  = SBS->Atom[a22-1]->energyType;
        if ((at21[0]) && (at22[0]))  {
          if (a21==a1)
            angle = GetBondAngle ( Loop, at22,at1,at2, a22,a1,a2 );
          else if (a22==a1)
            angle = GetBondAngle ( Loop, at21,at1,at2, a21,a1,a2 );
          else if (a21==a2)
            angle = GetBondAngle ( Loop, at1,at2,at22, a1,a2,a22 );
          else if (a22==a2)
            angle = GetBondAngle ( Loop, at1,at2,at21, a1,a2,a21 );
          else angle = NULL;
          if (angle)
            SBS->AddAngle ( angle );
        }
      }
  }

}


void FillMonLibData ( PCSBStructure SBS, RCMMCIFFile MonLib )  {
PCMMCIFData  mmCIFData;
PCMMCIFLoop  Loop;
char         mmCIFDataName[100];
char         pdb_name[100];
pstr         atname;
int          i,j,k,l;

  strcpy ( mmCIFDataName,"comp_" );
  strcat ( mmCIFDataName,SBS->compoundID );
  mmCIFData = MonLib.GetCIFData ( mmCIFDataName );

  if (!mmCIFData)  return;

  printf ( "; monomer data found " );

  Loop = mmCIFData->GetLoop ( "_chem_comp_atom" );
  if (!Loop)  {
    printf ( " monomer lib error:  %s:_chem_comp_atom not found\n",
             SBS->compoundID );
    return;
  }

  l = Loop->GetLoopLength();
  for (i=0;i<SBS->nAtoms;i++)  {
    strcpy_css ( pdb_name,SBS->Atom[i]->pdb_name );
    if ((!strcmp(SBS->Atom[i]->element," H" )) &&
        (pdb_name[0]>='0') && (pdb_name[0]<='9'))  {
      j = strlen(pdb_name);
      pdb_name[j]   = pdb_name[0];
      pdb_name[0]   = ' ';
      pdb_name[j+1] = char(0);
      DelSpaces ( pdb_name );
    }
    atname = NULL;
    for (j=0;(j<l) && (!atname);j++)  {
      atname = Loop->GetString ( "atom_id",j,k );
      if (atname)  {
        if (!strcmp(atname,pdb_name))  {
          if (Loop->GetReal(SBS->Atom[i]->ccp4_charge,
              "partial_charge",j,False))
            SBS->Atom[i]->ccp4_charge = 0.0;  // unsuccessful retrieval
        } else
          atname = NULL;
      }
    }
  }

}



pstr getString ( pstr S, PCMMCIFLoop mmCIFLoop,
                 PCMMCIFStruct mmCIFStruct, cpstr tag,
                 int row, cpstr defS )  {
pstr p  = NULL;
int  rc = 0;
  if (mmCIFLoop)         p = mmCIFLoop->GetString   ( tag,row,rc );
  else if (mmCIFStruct)  p = mmCIFStruct->GetString ( tag,rc );
  if (p)  strcpy ( S,p    );
    else  strcpy ( S,defS );
  return p;
}

int getInteger ( int & v, PCMMCIFLoop mmCIFLoop,
                 PCMMCIFStruct mmCIFStruct, cpstr tag,
                 int row, int defV )  {
int rc=1;
  if (mmCIFLoop)         rc = mmCIFLoop->GetInteger   ( v,tag,row );
  else if (mmCIFStruct)  rc = mmCIFStruct->GetInteger ( v,tag );
  if (rc)  v = defV;
  return rc;
}

int getReal ( realtype & v, PCMMCIFLoop mmCIFLoop,
              PCMMCIFStruct mmCIFStruct, cpstr tag,
              int row, realtype defV )  {
int rc=1;
  if (mmCIFLoop)         rc = mmCIFLoop->GetReal   ( v,tag,row );
  else if (mmCIFStruct)  rc = mmCIFStruct->GetReal ( v,tag );
  if (rc)  v = defV;
  return rc;
}



DefineClass(CSortSBIndex)

class CSortSBIndex : public CQuickSort  {
  public :
    CSortSBIndex() : CQuickSort() {}
    int  Compare ( int i, int j );
    void Swap    ( int i, int j );
    void Sort    ( PPCSBIndex index, int nCompounds );
};

int CSortSBIndex::Compare ( int i, int j )  {
  return strcmp ( ((PPCSBIndex)data)[i]->compoundID,
                  ((PPCSBIndex)data)[j]->compoundID );
}

void CSortSBIndex::Swap ( int i, int j )  {
PCSBIndex I;
  I = ((PPCSBIndex)data)[i];
  ((PPCSBIndex)data)[i] = ((PPCSBIndex)data)[j];
  ((PPCSBIndex)data)[j] = I;
}


void CSortSBIndex::Sort ( PPCSBIndex index, int nCompounds )  {
  CQuickSort::Sort ( &(index[0]),nCompounds );
}
*/

/*
int MakeSBase ( int procKey, psvector compFileName, int nNames,
                pstr monLib, pstr energyLib, pstr sbaseDir )  {
CSortSBIndex  sortIndex;
CFile         graphFile;
CFile         structFile;
CFile         indexFile;
CFile         compFile;
CMMCIFData    ETLib;
CMMCIFData    compound;
CMMCIFFile    MonomerLib;  // used only for fetching atom charges
PCMMCIFStruct mmCIFStruct;
PCMMCIFLoop   mmCIFLoop;
PPCSBIndex    Index;
PPCSBIndex    Kndex;
PCSBIndex     SBIndex;
PCGraph       Graph;
PCVertex      Vertex;
PCEdge        Edge;
PCSBStructure SBStructure;
PCSBAtom      SBAtom;
PCSBBond      SBBond;
pstr          p;
char          S[1000];
char          chirality[10];
char          scb_bond_order[20];
AtomName      pdb_name;
AtomName      pdb_name_old;
AtomName      sca_name;
AtomName      sca1_name;
AtomName      sca2_name;
ResName       cif_id;
Element       element;
realtype      sca_charge,coor_x,coor_y,coor_z;
int           pdb_rcsb_align;
int           lcount,nCompounds,nCAlloc, n,i,natoms,nbonds;
int           vtype, etype,vx1,vx2, nXTs;
int           k,rc;
Boolean       takeit,rel,obs;

  Index       = NULL;
  SBIndex     = NULL;
  Graph       = NULL;
  Vertex      = NULL;
  Edge        = NULL;
  SBStructure = NULL;
  SBAtom      = NULL;
  SBBond      = NULL;

  rel = procKey & PKEY_Rel;
  obs = procKey & PKEY_Obs;

  ETLib.SetFlag ( CIFFL_PrintWarnings );
  k = ETLib.ReadMMCIFData ( energyLib );
  if (k)  {
    printf ( " *** cannot read energy type library:\n"
             " %s\n",GetCIFMessage(S,k) );
    return -1;
  }

  MonomerLib.SetPrintWarnings ( True );
  k = MonomerLib.ReadMMCIFFile ( monLib );
  if (k)  {
    printf ( " *** cannot read monomer library:\n"
             " %s\n",GetCIFMessage(S,k) );
    return -1;
  }

  printf ( " ... total %i compounds in monomer library\n",
           MonomerLib.GetNofData() );

  strcpy ( S,sbaseDir    );
  strcat ( S,sbGraphFile );
  graphFile.assign ( S,False,True );
  if (!graphFile.rewrite())  {
    printf ( " *** can't open sbase file '%s' file.\n",S );
    return -1;
  }

  strcpy ( S,sbaseDir     );
  strcat ( S,sbStructFile );
  structFile.assign ( S,False,True );
  if (!structFile.rewrite())  {
    printf ( " *** can't open sbase file '%s' file.\n",S );
    return -1;
  }

  strcpy ( S,sbaseDir    );
  strcat ( S,sbIndexFile );
  indexFile.assign ( S,False,True );
  if (!indexFile.rewrite())  {
    printf ( " *** can't open sbase file '%s' file.\n",S );
    return -1;
  }

  SBIndex     = new CSBIndex();
  Graph       = new CGraph();
  SBStructure = new CSBStructure();

  nCompounds = 0;
  nCAlloc    = 0;

  for (n=0;n<nNames;n++)  {

    compFile.assign ( compFileName[n],True,False );
    if (!compFile.reset(True))  {
      printf ( " *** cannot open components file %s\n",
               compFileName[n] );
      return -1;
    } else  {

      printf ( " --- try file %i/%i '%s'\n",
               n+1,nNames,compFileName[n] );

      lcount     = 0;
      rc         = 0;

      while (!compFile.FileEnd())  {

        S[0] = char(0);
        rc = compound.ReadMMCIFData ( compFile,S,lcount );

        if (rc!=CIFRC_Ok)  {  // error or warning
          if ((rc<0) && (!compFile.FileEnd()))  { // error
            printf ( " *** error reading file %s:\n"
                   "     %s\n",compFileName[n],GetCIFMessage(S,rc) );
            exit(rc);
          } else if (rc>0)  { // warning
            printf ( " ... warning on reading file %s:\n"
                   "     %s\n",compFileName[n],GetCIFMessage(S,rc) );
          }
        } else  {

          SBIndex->fGraphPos  = graphFile .Position();
          SBIndex->fStructPos = structFile.Position();

          SBStructure->Reset();
          Graph      ->Reset();
          nXTs = 0;

          mmCIFStruct = compound.GetStructure ( "_chem_comp" );
          strcpy ( cif_id,compound.GetDataName() );

          if (mmCIFStruct)  {
            rc = 0;
            p = mmCIFStruct->GetString ( "id",rc );
            if (p)  strcpy ( cif_id,p );
            SBStructure->PutName    (
               mmCIFStruct->GetString ( "name",rc ) );
            SBStructure->PutFormula (
               mmCIFStruct->GetString ( "formula",rc ) );
            SBStructure->PutCharge  (
               mmCIFStruct->GetString ( "pdbx_formal_charge",rc ) );
            SBStructure->PutSynonym  (
               mmCIFStruct->GetString ( "pdbx_synonyms",rc ) );
            takeit = True;
            if (rel || obs)  {
              p = mmCIFStruct->GetString ( "pdbx_release_status",rc );
              if (rel && strcmp(p,"REL"))  {
                printf ( " +++ [%4i] compound '%s' not released\n",
                         nCompounds+1,cif_id );
                takeit = False;
              }
              if (obs && ((!strcmp(p,"OBS")) ||
                          (!strcmp(p,"REF_ONLY"))))  {
                printf ( " +++ [%4i] compound '%s' is obsolete\n",
                         nCompounds+1,cif_id );
                takeit = False;
              }
            }
          } else  {
            printf ( " *** [%4i] compound '%s' is not annotated\n",
                     nCompounds+1,cif_id );
            takeit = False;
          }

          //  1. Fetch and store atom data

          mmCIFStruct = NULL;
          mmCIFLoop   = compound.GetLoop ( "_chem_comp_atom" );
          if (!mmCIFLoop)  {
            mmCIFStruct = compound.GetStructure ( "_chem_comp_atom" );
            if (!mmCIFStruct)  {
              printf ( " *** no atom data in entry '%s'\n",cif_id );
              takeit = false;
            }
          }

          if (takeit && (mmCIFLoop || mmCIFStruct))  {

            printf ( " ... [%4i] compound '%s'",
                     nCompounds+1,compound.GetDataName() );

            SBStructure->xyz_source = 'R';

            if (mmCIFLoop)  natoms = mmCIFLoop->GetLoopLength();
                      else  natoms = 1;
            for (i=0;i<natoms;i++)  {

              if (getString(element,mmCIFLoop,mmCIFStruct,
                            "type_symbol",i,""))  {
                if (element[0])  {
                  if (!element[1])  {
                    element[1] = element[0];
                    element[0] = ' ';
                    element[2] = char(0);
                  }
                  vtype = getElementNo ( element );
                  if (getString(S,mmCIFLoop,mmCIFStruct,
                                "pdbx_leaving_atom_flag",i,"N"))  {
                    if (S[0]=='Y')  vtype |= ATOM_LEAVING;
                  }
                } else
                  vtype = -11;
              } else
                vtype = -11;

              if (vtype>=0)  {

                getString ( chirality,mmCIFLoop,mmCIFStruct,
                            "pdbx_stereo_config",i,"N" );
                UpperCase ( chirality );

                getString ( pdb_name_old,mmCIFLoop,mmCIFStruct,
                            "alt_atom_id",i,"" );
                getString ( sca_name,mmCIFLoop,mmCIFStruct,
                            "atom_id",i,"" );
                getString ( pdb_name,mmCIFLoop,mmCIFStruct,
                            "pdbx_component_atom_id",i,"" );

                getInteger ( pdb_rcsb_align,mmCIFLoop,mmCIFStruct,
                             "pdbx_align",i,2-strlen(element) );
                getReal    ( sca_charge,mmCIFLoop,mmCIFStruct,
                             "charge",i,0.0 );

                // 1.1  Generate vertex of structural graph
                Vertex = new CVertex ( MakeElementType(vtype,
                                             chirality[0],True) );
                Graph->AddVertex ( Vertex );

                // 1.2  Generate SB atom
                SBAtom = new CSBAtom();

                strcpy ( SBAtom->sca_name,sca_name );

                MakePDBName ( SBAtom->pdb_name,pdb_name,
                              pdb_rcsb_align );

                strcpy      ( SBAtom->old_pdb_name,pdb_name_old );
                if (!strcmp(&(SBAtom->pdb_name[2]),"XT"))  nXTs++;
                if (!element[1])  {
                  SBAtom->element[0] = ' ';
                  strcpy ( &(SBAtom->element[1]),element );
                } else
                  strcpy ( SBAtom->element ,element  );
                SBAtom->chirality = chirality[0];
                if (vtype & ATOM_LEAVING)
                  SBAtom->leaving = 'Y';
                SBAtom->sca_charge     = sca_charge;
                SBAtom->partial_charge = sca_charge;

                getReal ( coor_x,mmCIFLoop,mmCIFStruct,
                          "pdbx_model_Cartn_x_ideal",i,0.0 );
                getReal ( coor_y,mmCIFLoop,mmCIFStruct,
                          "pdbx_model_Cartn_y_ideal",i,0.0 );
                getReal ( coor_z,mmCIFLoop,mmCIFStruct,
                          "pdbx_model_Cartn_z_ideal",i,0.0 );

                SBAtom->x     = coor_x;
                SBAtom->y     = coor_y;
                SBAtom->z     = coor_z;
                SBAtom->x_esd = 0.0;
                SBAtom->y_esd = 0.0;
                SBAtom->z_esd = 0.0;

                SBStructure->AddAtom ( SBAtom );

              } else
                printf ( "\n"
                     "    *** unknown element %s in _chem_comp_atom"
                     " in entry %s.\n",element,cif_id );
            }

            //  2. Fetch and store bond data

            mmCIFStruct = NULL;
            mmCIFLoop   = compound.GetLoop ( "_chem_comp_bond" );
            if (!mmCIFLoop)  {
              mmCIFStruct = compound.GetStructure ( "_chem_comp_bond" );
              if (!mmCIFStruct)  {
                if (natoms>1)
                  printf ( "\n    *** no bond loop in entry '%s'\n",
                         cif_id );
//              takeit = false;
              }
            }

            if (takeit && (mmCIFLoop || mmCIFStruct))  {

              if (mmCIFLoop)  nbonds = mmCIFLoop->GetLoopLength();
                        else  nbonds = 1;

              for (i=0;i<nbonds;i++)  {

                getString ( sca1_name,mmCIFLoop,mmCIFStruct,
                            "atom_id_1",i,"" );
                getString ( sca2_name,mmCIFLoop,mmCIFStruct,
                            "atom_id_2",i,"" );
                getString ( scb_bond_order,mmCIFLoop,mmCIFStruct,
                            "value_order",i,"SING" );

                //  2.1 Find ordinal numbers of bonded atoms

                rc  = 0;
                vx1 = SBStructure->GetAtomNo_nss ( sca1_name );
                vx2 = SBStructure->GetAtomNo_nss ( sca2_name );
                if (vx1<=0)  {
                  printf ( " ... wrong bond vertex %s in "
                           "_chem_comp_bond at cif_id=%s.\n",
                           sca1_name,cif_id );
                  rc = 1;
                }
                if (vx2<=0)  {
                  printf ( " ... wrong bond vertex %s in "
                           "_chem_comp_bond at cif_id=%s.\n",
                           sca2_name,cif_id );
                  rc = 1;
                }

                //  2.2 Encode the bond order

                etype = GetBondOrder ( scb_bond_order );
                if (etype<0)  {
                  printf ( " ... wrong bond type '%s' in "
                           "_chem_comp_bond at cif_id=%s.\n",
                           scb_bond_order,cif_id );
                  rc = 1;
                }

                if (!rc)  {
                  //  2.3  Add edge to the graph
                  Edge = new CEdge ( vx1,vx2,etype );
                  Graph->AddEdge ( Edge );
                  //  2.4  Add bond to the structure
                  SBBond = new CSBBond();
                  SBBond->SetBond ( vx1,vx2,etype );
                  SBStructure->AddBond ( SBBond );
                }

              }

            }

            SBStructure->MakeLeavingAtoms ();
            SBStructure->RemoveEnergyTypes();

          }

          //  5. Store data in binary files

          if (!takeit)  {

            //  4.1  Unsuccessful retrieval -- backspace files
            graphFile .seek ( SBIndex->fGraphPos  );
            structFile.seek ( SBIndex->fStructPos );

          } else  {

            //  4.2  Store graph and structure
            Graph->SetName ( cif_id                  );
            StreamWrite    ( graphFile  ,Graph       );

            strcpy         ( SBStructure->compoundID,cif_id );
            FillETLibData  ( SBStructure,ETLib       );
            FillMonLibData ( SBStructure,MonomerLib  );
            StreamWrite    ( structFile ,SBStructure );

            if (nCompounds>=nCAlloc)  {
              nCAlloc = nCompounds + 5000;
              Kndex = new PCSBIndex[nCAlloc];
              for (i=0;i<nCompounds;i++)
                Kndex[i] = Index[i];
              for (i=nCompounds;i<nCAlloc;i++)
                Kndex[i] = NULL;
              if (Index)  delete[] Index;
              Index = Kndex;
            }
            SBIndex->MakeCompositions ( SBStructure );
            strcpy ( SBIndex->compoundID,cif_id );
            SBIndex->nXTs = nXTs;
            Index[nCompounds++] = SBIndex;
            SBIndex = new CSBIndex();

            printf ( " -- done\n" );

          }

        }

      }

    }


  }

  printf ( "\n\n -------- total %i compounds.\n",nCompounds );

  compFile  .shut();
  graphFile .shut();
  structFile.shut();

  sortIndex.Sort ( Index,nCompounds );

  for (i=0;i<nCompounds;i++)
    if (Index[i])  {
      StreamWrite ( indexFile,Index[i] );
      delete Index[i];
    }
  indexFile .shut();

  if (Index)        delete[] Index;
  if (Graph)        delete   Graph;
  if (SBIndex)      delete   SBIndex;
  if (SBStructure)  delete   SBStructure;

  return rc;

}
*/

void printInstructions ( mmdb::pstr argv0 )  {

  printf (
    "\n"
    " CCP4SRS generator\n"
    " -----------------\n"
    "\n"
    " USAGE:\n"
    "\n"
    "%s rcsblib_dir [monlib_dir [old_monlib_dir]] [-d srs_dir]\n"
    "\n"
    " where 'rcsblib_dir' is directory with RCSB's chemical\n"
    " compound library directory;\n"
    "\n"
    " 'monlib_dir' is CCP4 optional monomer library's directory. In\n"
    " the original CCP4 distribution, it is found in\n"
    " $CCP4/lib/data/monomers;\n"
    "\n"
    " 'old_monlib_dir' is optional CCP4 monomer library's directory.\n"
    " Old entries are graph-matched to new ones, and name maps are\n"
    " stored in SRS for renaming utilities;\n"
    "\n"
    " 'srs_dir' is an optional parameter specifying where SRS\n"
    " should be placed. If SRS directory is not specified,\n"
    " the current one will be used.\n"
    ,argv0
   );

}


DefineClass(CSortSRSIndex)

class CSortSRSIndex : public mmdb::QuickSort  {
  public :
    CSortSRSIndex() : mmdb::QuickSort() {}
    int  Compare ( int i, int j );
    void Swap    ( int i, int j );
    void Sort    ( ccp4srs::PPIndex    index,
                   ccp4srs::PPMonomer  monomer,
                   mmdb::math::PPGraph graph,
                   int nEntries );
  protected:
    ccp4srs::PPMonomer  m_array;
    mmdb::math::PPGraph g_array;
};

int CSortSRSIndex::Compare ( int i, int j )  {
  return strcmp ( ((ccp4srs::PPIndex)data)[i]->entryID,
                  ((ccp4srs::PPIndex)data)[j]->entryID );
}

void CSortSRSIndex::Swap ( int i, int j )  {
ccp4srs::PIndex    I;
ccp4srs::PMonomer  M;
mmdb::math::PGraph G;

  I = ((ccp4srs::PPIndex)data)[i];
  ((ccp4srs::PPIndex)data)[i] = ((ccp4srs::PPIndex)data)[j];
  ((ccp4srs::PPIndex)data)[j] = I;

  M = m_array[i];
  m_array[i] = m_array[j];
  m_array[j] = M;

  G = g_array[i];
  g_array[i] = g_array[j];
  g_array[j] = G;

}


void CSortSRSIndex::Sort ( ccp4srs::PPIndex    index,
                           ccp4srs::PPMonomer  monomer,
                           mmdb::math::PPGraph graph,
                           int nEntries )  {
  m_array = monomer;
  g_array = graph;
  mmdb::QuickSort::Sort ( &(index[0]),nEntries );
}


class CMakeSRS  {

  public:
    CMakeSRS ();
    ~CMakeSRS();

    void FreeMemory();

    int  openLogs ( mmdb::cpstr dir );
    void prints   ( mmdb::cpstr S, int key );
    void closeLogs();

    int readEnergyLib      ( mmdb::cpstr monLib  );
    int scanRCSBLibrary    ( mmdb::cpstr rcsbLib );
    int scanMonomerLibrary ( mmdb::cpstr monLib, bool match_only );

    void sortIndex         ();
    int  writeSRS          ( mmdb::cpstr srsDir );

  protected:
    mmdb::mmcif::Data      energyLib;
    mmdb::math::GraphMatch graphMatch;
    mmdb::io::File         fout;
    mmdb::io::File         ferr;
    ccp4srs::PPIndex       srsIndex;
    ccp4srs::PPMonomer     srsMonomer;
    mmdb::math::PPGraph    srsGraph;
    ccp4srs::MemIO         memIO;
    int                    srsVersion;
    int                    nEntries,nEntries0;
    bool                   matchOnly;

    // error counters
    int              n_rcsb_read_err;
    int              n_ml_read_err;
    int              n_graph_err;
    int              n_elib_err;
    int              n_oldml_umbiguous;
    int              n_oldml_no_structure;
    int              n_oldml_not_found;
    int              n_oldml_no_matches;

    int  scanMonLibDirectory ( mmdb::cpstr dirPath  );
    int  processMonLibFile   ( mmdb::cpstr filePath );
    int  processRCSBFile     ( mmdb::cpstr filePath );
    void FillETLibData       ( ccp4srs::PMonomer monomer );
    int  matchAtomNames      ( ccp4srs::PMonomer mon,
                               ccp4srs::PMonomer mon_old );
    int  processMonomer      ( ccp4srs::PMonomer monomer );
    int  getMonomerNo        ( mmdb::cpstr monID );

  private:
    int index_alloc;

};



CMakeSRS::CMakeSRS()  {

  srsVersion  = ccp4srs::CCP4SRS_file_version;
  srsIndex    = NULL;
  srsMonomer  = NULL;
  srsGraph    = NULL;
  nEntries    = 0;
  index_alloc = 0;

  n_rcsb_read_err      = 0;
  n_ml_read_err        = 0;
  n_graph_err          = 0;
  n_elib_err           = 0;
//  n_oldml_read_err     = 0;
  n_oldml_umbiguous    = 0;
  n_oldml_no_structure = 0;
  n_oldml_not_found    = 0;
  n_oldml_no_matches   = 0;

}

CMakeSRS::~CMakeSRS()  {
  FreeMemory();
}

void CMakeSRS::FreeMemory()  {
int i;

  if (srsIndex)  {
    for (i=0;i<index_alloc;i++)
      if (srsIndex[i])  delete srsIndex[i];
    delete[] srsIndex;
  }

  if (srsMonomer)  {
    for (i=0;i<index_alloc;i++)
      if (srsMonomer[i])  delete srsMonomer[i];
    delete[] srsMonomer;
  }

  if (srsGraph)  {
    for (i=0;i<index_alloc;i++)
      if (srsGraph[i])  delete srsGraph[i];
    delete[] srsGraph;
  }

  nEntries    = 0;
  index_alloc = 0;

}

int CMakeSRS::openLogs ( mmdb::cpstr dir )  {
mmdb::pstr fname;
int        rc;

  fname = NULL;
  rc    = 0;

  mmdb::CreateCopCat ( fname,dir,"srsgen_out.txt" );
  fout.assign ( fname,true );
  if (!fout.rewrite())  rc = -1;

  mmdb::CreateCopCat ( fname,dir,"srsgen_err.txt" );
  ferr.assign ( fname,true );
  if (!ferr.rewrite())  rc = -2;

  delete[] fname;

  n_rcsb_read_err      = 0;
  n_ml_read_err        = 0;
  n_graph_err          = 0;
  n_elib_err           = 0;
//  n_oldml_read_err     = 0;
  n_oldml_umbiguous    = 0;
  n_oldml_no_structure = 0;
  n_oldml_not_found    = 0;
  n_oldml_no_matches   = 0;

  return rc;

}

void CMakeSRS::closeLogs()  {
char S[1000];

  sprintf ( S,"\n ---------------------------------------------\n"
              " RCSB read errors           : %i\n"
              " Monomer Library read errors: %i\n"
              " Graph erorrs               : %i\n"
              " Energy Library read errors : %i\n"
              " Umbiguous graoh matches    : %i\n"
              " Entries with no structure  : %i\n"
              " Not found old ML entries   : %i\n"
              " Not matched old ML entries : %i\n",
              n_rcsb_read_err  , n_ml_read_err       ,
              n_graph_err      , n_elib_err          ,
              n_oldml_umbiguous, n_oldml_no_structure,
              n_oldml_not_found, n_oldml_no_matches );

  prints ( S,3 );

  fout.shut();
  ferr.shut();

}

void CMakeSRS::prints ( mmdb::cpstr S, int key )  {
  printf ( "%s",S );
  if (key & 0x00000001)  fout.Write ( S );
  if (key & 0x00000002)  ferr.Write ( S );
}

int CMakeSRS::getMonomerNo ( mmdb::cpstr monID )  {
mmdb::ResName id;
int           l1,l2,l,k;

  mmdb::strcpy_css ( id,monID );

  k = -1;

  if (nEntries0<=3)  {
    for (l=0;(l<nEntries0) && (k<0);l++)
      if (!strcmp(srsIndex[l]->entryID,id))
        k = l;
    if (k>=0)  return k;
    return ccp4srs::CCP4SRS_EntryNotFound;
  }

  l  = 0;
  l1 = 0;
  l2 = nEntries0-1;
  while (l1<l2-1)  {
    l = (l1+l2)/2;
    k = strcmp ( srsIndex[l]->entryID,id );
    if (k<0)       l1 = l;
    else if (k>0)  l2 = l;
    else {
      l1 = l;
      l2 = l;
    }
  }

  if (k==0)  return l;
  else if (l==l1)  {
    if (!strcmp(srsIndex[l2]->entryID,id))  return l2;
  } else if (l==l2)  {
    if (!strcmp(srsIndex[l1]->entryID,id))  return l1;
  }

  return ccp4srs::CCP4SRS_EntryNotFound;

}


void CMakeSRS::FillETLibData ( ccp4srs::PMonomer monomer )  {
mmdb::mmcif::PLoop  Loop;
ccp4srs::PAtom      Atom;
mmdb::realtype      vdwRadius,vdwhRadius,ionRadius;
mmdb::cpstr         at1,et1,bt;
int                 valency,sp;
int                 i,j,l,a1, rc;

  // 1. Set atom properties depending on energy type

  Loop = energyLib.GetLoop ( "_lib_atom" );
  if (!Loop)  {
    prints ( " *** energy lib error:  _lib_atom not found\n",3 );
    return;
  }

  l = Loop->GetLoopLength();
  for (i=0;i<monomer->n_atoms();i++)  {
    Atom = monomer->atom(i);
    at1  = Atom->energy_type();
    a1   = -1;
    if (at1[0])  {
      for (j=0;(j<l) && (a1<0);j++)  {
        et1 = Loop->GetString ( "type",j,rc );
        if ((!rc) && et1)  {
          if (!strcmp(et1,at1))  a1 = j;
        }
      }
      if (a1>=0)  {
        vdwRadius  = Atom->vdw_radius ();
        vdwhRadius = Atom->vdwh_radius();
        ionRadius  = Atom->ion_radius ();
        valency    = Atom->valency    ();
        sp         = Atom->sp         ();
        if (vdwRadius<=0.01)  {
          if (Loop->GetReal(vdwRadius,"vdw_radius",a1,false))
            vdwRadius = mmdb::getVdWaalsRadius  ( Atom->element() );
          Atom->set_vdw_radius ( vdwRadius );
        }
        if (vdwhRadius<=0.01)  {
          if (Loop->GetReal(vdwhRadius,"vdwh_radius",a1,false))
            vdwhRadius = vdwRadius + mmdb::getVdWaalsRadius  ( "H" );
          Atom->set_vdwh_radius ( vdwhRadius );
        }
        if (ionRadius<=0.01)  {
          if (Loop->GetReal(ionRadius,"ion_radius",a1,false))
            ionRadius = mmdb::getCovalentRadius ( Atom->element() );
          Atom->set_ion_radius ( ionRadius );
        }
        if (valency<=0)  {
          Loop->GetInteger ( valency,"valency",a1,false );
          Atom->set_valency ( valency );
        }
        if (sp<=0)  {
          Loop->GetInteger ( sp,"sp",a1,false );
          if (sp<0)  sp = 0;
          Atom->set_sp ( sp );
        }
        bt = Loop->GetString ( "hb_type",a1,rc );
        if ((!rc) && bt)  {
          if (strpbrk(bt,"ADHBN"))  Atom->set_hb_type (  bt[0] );
                              else  Atom->set_hb_type ( 'N' );
// temporary?? terrible fix:
//if ((!strcmp(SBS->compoundID,"DC")) &&
//    (!strcmp(Atom->pdb_name," N3 ")))
//  Atom->hb_type = 'A';
        }
      }
    }
    if (a1<0)  {
      // take whatever is possible from MMDB tables
      Atom->set_vdw_radius  ( mmdb::getVdWaalsRadius(Atom->element()) );
      Atom->set_vdwh_radius ( Atom->vdw_radius() +
                              mmdb::getVdWaalsRadius("H") );
      Atom->set_ion_radius  ( mmdb::getCovalentRadius(Atom->element()) );
      Atom->set_valency ( 0 );    // valency is unknown
      Atom->set_sp      ( 0 );    // sp is unknown
      Atom->set_hb_type ( 'N' );  // no hydrogen bonding
    }
  }

  /*

  // 2. Fill in bond lengths

  Loop = ETLib.GetLoop ( "_lib_bond" );
  if (!Loop)  {
    printf ( " *** energy lib error:  _lib_bond not found\n" );
    return;
  }

  l = Loop->GetLoopLength();
  for (i=0;i<SBS->nBonds;i++)  {
    Bond = SBS->Bond[i];
    Bond->length     = 0.0;
    Bond->length_esd = 0.0;
    if (Bond->atom1>Bond->atom2)
      ISwap ( Bond->atom1,Bond->atom2 );
    at1 = SBS->Atom[Bond->atom1-1]->energyType;
    at2 = SBS->Atom[Bond->atom2-1]->energyType;
    if ((at1[0]) && (at2[0]))  {
      order       = Bond->order;
      blength     = 0.0;
      blength_esd = 0.0;
      nb          = 0;
      nb_esd      = 0;
      found       = 0;
      for (j=0;(j<l) && (found<2);j++)  {
        et1 = Loop->GetString ( "atom_type_1",j,rc );
        if (!rc)  {
          if (et1)  {
            if (!strcmp(et1,at1))      a1 = 1; // 1st atom identified
            else if (!strcmp(et1,at2)) a1 = 2; // 2nd atom identified
                                  else a1 = 0; // unsuitable pair in
                                               // ETLib
          } else
            a1 = 3;  // mmCIF '.' means any energy type for 1st atom
        } else
          a1 = 0;
        if (a1>0)  {
          et2 = Loop->GetString ( "atom_type_2",j,rc );
          if (!rc)  {
            if (et2)  {
              if ((a1==1) && (!strcmp(et2,at2)))      a2 = 2;
              else if ((a1==2) && (!strcmp(et2,at1))) a2 = 1;
              else if (a1==3)  {
                // first atom was identified as that of any energy type
                if (!strcmp(et2,at1))      a2 = 1;
                else if (!strcmp(et2,at2)) a2 = 2;
                                      else a2 = 0; // unsuitable pair
                                                   // in ETLib
              } else
                 a2 = 0;
            } else
              a2 = 3; // mmCIF '.' meaning any energy type for 2nd atom
            if (a2>0)  {
              bt = Loop->GetString ( "type",j,rc );
              if (!rc)  {
                if (bt)  {
                  btype = GetBondOrder ( bt );
                  if (btype!=order)  btype = -1;
                } else
                  btype = MaxInt;  // mmCIF '.' meaning any bond order
              } else
                btype = -MaxInt;
              if (btype>=0)  {
                if (!Loop->GetReal(blength,"length",j,False))  {
                  if ((a1<3) || (a2<3) || (nb==0))  {
                    nb++;                     // count
                    Bond->length += blength;  // and accumulate
                  }
                } else
                  blength = 0.0;  // unsuccessful retrieval
                if (!Loop->GetReal(blength_esd,"value_esd",j,False))  {
                  if ((a1<3) || (a2<3) || (nb_esd==0))  {
                    nb_esd++;                        // count
                    Bond->length_esd += blength_esd; // and accumulate
                  }
                } else  // unsuccessful retrieval
                  blength_esd = 0.0;
                if ((blength>0.0) && (btype<MaxInt) &&
                    (a2<3) && (a1<3))  {
                  found = 2;
                  Bond->length     = blength;
                  Bond->length_esd = blength_esd;
                } else
                  found = 1;
              }
            }
          }
        }
      }
      if (found<2)  {
        // average over incomplete energy type pairs
        if (nb>0)     Bond->length     /= nb;
        if (nb_esd>0) Bond->length_esd /= nb_esd;
      }
      if (!found)
        printf ( " can't find bond data for energy types %s and %s\n",
                 at1,at2 );
    }
  }


  // 3. Make angles

  Loop = ETLib.GetLoop ( "_lib_angle" );
  if (!Loop)  {
    printf ( " energy lib error:  _lib_angle not found\n" );
    return;
  }

  for (i=0;i<SBS->nBonds;i++)  {
    Bond = SBS->Bond[i];
    a1   = Bond->atom1;
    a2   = Bond->atom2;
    at1  = SBS->Atom[a1-1]->energyType;
    at2  = SBS->Atom[a2-1]->energyType;
    if ((at1[0]) && (at2[0]))
      for (j=i+1;j<SBS->nBonds;j++)  {
        Bond2 = SBS->Bond[j];
        a21   = Bond2->atom1;
        a22   = Bond2->atom2;
        at21  = SBS->Atom[a21-1]->energyType;
        at22  = SBS->Atom[a22-1]->energyType;
        if ((at21[0]) && (at22[0]))  {
          if (a21==a1)
            angle = GetBondAngle ( Loop, at22,at1,at2, a22,a1,a2 );
          else if (a22==a1)
            angle = GetBondAngle ( Loop, at21,at1,at2, a21,a1,a2 );
          else if (a21==a2)
            angle = GetBondAngle ( Loop, at1,at2,at22, a1,a2,a22 );
          else if (a22==a2)
            angle = GetBondAngle ( Loop, at1,at2,at21, a1,a2,a21 );
          else angle = NULL;
          if (angle)
            SBS->AddAngle ( angle );
        }
      }
  }

  */

}


int CMakeSRS::processMonomer ( ccp4srs::PMonomer monomer )  {
ccp4srs::PIndex     index;
ccp4srs::PPIndex    kndex;
ccp4srs::PPMonomer  marray;
mmdb::math::PPGraph garray;
mmdb::math::PGraph  graph;
char                S[500];
int                 rc,i;

  rc    = 0;

  index = new ccp4srs::Index();

  strcpy ( index->entryID,monomer->ID() );
  index->nAtoms = monomer->n_atoms();
  index->nBonds = monomer->n_bonds();
  index->MakeCompositions ( monomer );

  graph = monomer->getGraph ( &rc );

  sprintf ( S," ... [%4i] compound '%s'",nEntries+1,monomer->ID() );
  prints  ( S,1 );

  if (rc>=30000)  {
    sprintf ( S,"\n +++ wrong bond vertex #2, bond no. %i in entry %s.\n",
                rc-30000,monomer->ID() );
    prints  ( S,3 );
    n_graph_err++;
    rc = 2;
  } else if (rc>=20000)  {
    sprintf ( S,"\n +++ wrong bond vertex #1, bond no. %i in entry %s.\n",
                rc-20000,monomer->ID() );
    n_graph_err++;
    prints  ( S,3 );
    rc = 2;
  } else if (rc>=10000)  {
    rc -= 10000;
    sprintf ( S,"\n"
             "    *** unknown element %s in _chem_comp_atom"
             " in entry %s.\n",monomer->atom(rc)->element(),
             monomer->ID() );
    prints  ( S,3 );
    n_graph_err++;
    rc = 1;
  }

  //  Store SRS data

  if (nEntries>=index_alloc)  {
    index_alloc = nEntries + 5000;
    kndex  = new ccp4srs::PIndex   [index_alloc];
    marray = new ccp4srs::PMonomer [index_alloc];
    garray = new mmdb::math::PGraph[index_alloc];
    for (i=0;i<nEntries;i++)  {
      kndex [i] = srsIndex  [i];
      marray[i] = srsMonomer[i];
      garray[i] = srsGraph  [i];
    }
    for (i=nEntries;i<index_alloc;i++)  {
      kndex [i] = NULL;
      marray[i] = NULL;
      garray[i] = NULL;
    }
    if (srsIndex)   delete[] srsIndex;
    if (srsMonomer) delete[] srsMonomer;
    if (srsGraph)   delete[] srsGraph;
    srsIndex   = kndex;
    srsMonomer = marray;
    srsGraph   = garray;
  }
  srsIndex  [nEntries] = index;
  srsMonomer[nEntries] = monomer;
  srsGraph  [nEntries] = graph;
  nEntries++;
  prints ( " -- added\n",1 );

  return rc;

}


int CMakeSRS::matchAtomNames ( ccp4srs::PMonomer mon,
                               ccp4srs::PMonomer mon_old )  {
mmdb::math::PGraph  graph;
mmdb::math::PGraph  graph_old;
mmdb::AtomName      aname1,aname2;
mmdb::ivector       FV1,FV2;
mmdb::ivector       amap;
mmdb::ivector       hcnt;
mmdb::realtype      p1,p2;
char                S[500];
int                 htype,nv,i,n1,n2, atm_old,atm_new,atm_h;

  if (mon->n_bonds()<=0)  {
    sprintf ( S," +++++ compound '%s' does not have a structure\n",
                mon->ID() );
    prints  ( S,3 );
    n_oldml_no_structure++;
    return -1;
  }

  if (mon_old->n_bonds()<=0)  {
    sprintf ( S," +++++ compound '%s' from old ML does not have "
                "a structure\n",mon->ID() );
    prints  ( S,3 );
    n_oldml_no_structure++;
    return -1;
  }

  graph = mon->getGraph(&i);
  if (!graph)  {
    sprintf ( S," +++++ cannot build graph for compound '%s'\n",
                mon->ID() );
    prints  ( S,3 );
    n_graph_err++;
    return -1;
  }

  graph_old = mon_old->getGraph(&i);
  if (!graph_old)  {
    sprintf ( S," +++++ cannot build graph for compound '%s' "
                "from old ML\n",mon_old->ID() );
    prints  ( S,3 );
    n_graph_err++;
    return -2;
  }

  htype = mmdb::getElementNo(mmdb::pstr("H"));

  graph->ExcludeType ( htype );
  graph->MakeVertexIDs();
  graph->MakeSymmetryRelief ( false );
  graph->Build ( false );

  graph_old->ExcludeType ( htype );
  graph_old->MakeVertexIDs();
  graph_old->MakeSymmetryRelief ( false );
  graph_old->Build ( false );

  graphMatch.SetMaxNofMatches ( 2,true );
  graphMatch.MatchGraphs ( graph,graph_old,
                           mmdb::IMin(graph->GetNofVertices(),
                                      graph_old->GetNofVertices()),
                           true,mmdb::math::EXTTYPE_Ignore );

  if (graphMatch.GetNofMatches()<=0)  {

    sprintf ( S," +++++ no matches for compound '%s' from old ML "
                "(%i:%i)\n",
                mon_old->ID(),graph->GetNofVertices(),
                graph_old->GetNofVertices() );
    n_oldml_no_matches++;
    prints  ( S,3 );
    nv = -1;

  } else  {

    mon->setOldID ( mon_old->ID() );

    if (graphMatch.GetNofMatches()>1)  {
      sprintf ( S," ----- umbiguous matching for compound '%s' "
                  "from old ML\n",mon_old->ID() );
      prints  ( S,3 );
      n_oldml_umbiguous++;
    }

    FV1 = NULL;
    FV2 = NULL;
    graphMatch.GetMatch ( 0,FV1,FV2,nv,p1,p2 );

    n1 = mmdb::IMax(mon_old->n_atoms(),mon->n_atoms());
    mmdb::GetVectorMemory ( amap,n1,0 );
    mmdb::GetVectorMemory ( hcnt,n1,0 );
    for (i=0;i<n1;i++)  {
      amap[i] = -1;
      hcnt[i] = 0;
    }

    for (i=1;i<=nv;i++)  {
      n1 = graph->GetVertex(FV1[i])->GetUserID();
      n2 = graph_old->GetVertex(FV2[i])->GetUserID();
      sprintf ( S,"        map '%4s':%i -> '%4s':%i\n",
                  mon_old->atom(n2)->name_pdb(aname2),n2,
                  mon->atom(n1)->name_pdb(aname1),n1 );
      prints  ( S,1 );
      mon->atom(n1)->set_old_name ( mon_old->atom(n2)->name() );
      amap[n2] = n1;
    }

    // map hydrogens

    for (i=0;i<mon_old->n_atoms();i++)
      if (!strcmp(mon_old->atom(i)->element(),"H"))  {
        atm_old = mon_old->get_bound ( i,0 );
        if (atm_old>=0)  {
          atm_new = amap[atm_old];
          if (atm_new>=0)  {
            atm_h = mon->get_bound ( atm_new,"H",hcnt[atm_new] );
            if (atm_h>=0)  {
              hcnt[atm_new]++;
              sprintf ( S,"        map '%4s':%i -> '%4s':%i\n",
                          mon_old->atom(i)->name_pdb(aname2),i,
                          mon->atom(atm_h)->name_pdb(aname1),atm_h );
              prints  ( S,1 );
              mon->atom(atm_h)->set_old_name (
                                         mon_old->atom(i)->name() );
            } else  {
              sprintf ( S," +++ inconsistency #3 for atom '%s':%i:%i\n",
                          mon_old->atom(i)->name(),atm_new,
                          hcnt[atm_new] );
              prints  ( S,3 );
            }
          } else  {
            sprintf ( S," +++ inconsistency #2 for atom '%s':%i\n",
                        mon_old->atom(i)->name(),atm_old );
            prints  ( S,3 );
          }
        } else  {
          sprintf ( S," +++ inconsistency #1 for atom '%s'\n",
                      mon_old->atom(i)->name() );
          prints  ( S,3 );
        }
      }

    mmdb::FreeVectorMemory ( hcnt,0 );
    mmdb::FreeVectorMemory ( amap,0 );

  }

  delete graph;
  delete graph_old;

  if (nv<0)  return -3;
  return 0;

}


int CMakeSRS::processMonLibFile ( mmdb::cpstr filePath )  {
// processing an entry from Monomer Library
mmdb::mmcif::File  mmCIFFile;
mmdb::mmcif::PData mmCIFData;
ccp4srs::PMonomer  monomer;
ccp4srs::PAtom     atom;
char               ErrLog[1000];
mmdb::ResName      rname;
mmdb::pstr         p;
char               S[500];
int                rc, nData,i,j, nMon;

  rc = mmCIFFile.ReadMMCIFFile ( filePath,mmdb::io::GZM_CHECK );

  if (rc)  {
    sprintf ( S," +++ file %s read with return code %i\n",
                filePath,rc );
    prints  ( S,3 );
    n_ml_read_err++;
    return 1;
  }

  nData = mmCIFFile.GetNofData();

  for (i=0;i<nData;i++)  {

    mmCIFData = mmCIFFile.GetCIFData ( i );
    if (strcmp(mmCIFData->GetDataName(),"comp_list"))  {

      monomer = new ccp4srs::Monomer();
      rc      = monomer->readFromCIF ( mmCIFData,ErrLog );
      if (rc)  {
        sprintf ( S," +++ error reading monomer %s:\n"
                    "     %s\n",filePath,ErrLog );
        prints  ( S,3 );
        delete monomer;
        n_ml_read_err++;
        monomer = NULL;
      } else
        FillETLibData ( monomer );

      if (monomer)  {

        strcpy ( rname,mmCIFData->GetDataName() );
        p = mmdb::FirstOccurence ( rname,'_' );
        if (!p)  p  = rname;
           else  p++;

        nMon = getMonomerNo ( p );
        if (nMon==ccp4srs::CCP4SRS_EntryNotFound)  {

          if (matchOnly)  {
            sprintf ( S," +++ [%4i] compound '%s' from old ML is not "
                        "found in new ML\n",
                        nMon,monomer->ID() );
            prints  ( S,3 );
            n_oldml_not_found++;
            delete monomer;
          } else  {
            processMonomer ( monomer );
//            nEntries0 = nEntries;
//            sortIndex();
          }

        } else if (!matchOnly)  {

          srsMonomer[nMon]->merge ( monomer );

          // a hack
          for (j=0;j<srsMonomer[nMon]->n_atoms();j++)  {
            atom = srsMonomer[nMon]->atom(j);
            if (!strcmp(atom->name(),"OXT"))  {
              atom->set_energy_type ( "O" );
              atom->set_chirality   ( 'N' );
              atom->set_hb_type     ( 'A' );
              atom->set_valency     ( 2   );
              atom->set_vdw_radius  ( 1.520 );
              atom->set_vdwh_radius ( 1.520 );
              atom->set_ion_radius  ( 1.280 );
            }
            if (!strcmp(atom->element(),"H"))  {
              if (atom->chirality()=='-')  atom->set_chirality   ( 'N' );
              if (!atom->energy_type()[0]) atom->set_energy_type ( "H" );
              if (atom->hb_type  ()=='-')  atom->set_hb_type     ( 'H' );
              atom->set_valency ( 1 );
              if (atom->vdw_radius()==0.0)
                atom->set_vdw_radius  ( 1.200 );
              if (atom->vdwh_radius()==0.0)
                atom->set_vdwh_radius ( 2.400 );
              if (atom->ion_radius()==0.0)
                atom->set_ion_radius  ( 0.320 );
            }
          }

          sprintf ( S," ::: [%4i] compound '%s' merged\n",
                      nMon,monomer->ID() );
          prints  ( S,1 );

          delete monomer;

        } else  {
          // just match atom names
          sprintf ( S," ::: [%4i] compound '%s' is being matched\n",
                      nMon,monomer->ID() );
          prints  ( S,1 );
          matchAtomNames ( srsMonomer[nMon],monomer );
          delete monomer;
        }

      }

    }

  }

  return 0;

}


int CMakeSRS::scanMonLibDirectory ( mmdb::cpstr dirPath )  {
DIR    *dir;
#if USE_DIRENT_H
dirent *dp;
#else
direct *dp;
#endif
mmdb::pstr p;
mmdb::pstr filePath;
int        rc;

  filePath = NULL;
  rc       = 0;

  dir = opendir(dirPath);
  if (!dir)
    return -111112;

  dp = readdir ( dir );
  while (dp)  {
    p = mmdb::LastOccurence ( dp->d_name,'.' );
    if (p)  {
      if (!strcmp(p,".cif"))  {
        mmdb::CreateCopCat ( filePath,dirPath,dp->d_name );
        rc = processMonLibFile ( filePath );
      }
    }
    dp = readdir ( dir );
  }

  closedir ( dir );

  return rc;

}

int CMakeSRS::processRCSBFile ( mmdb::cpstr filePath )  {
mmdb::io::File    compFile;
mmdb::mmcif::Data mmCIFData;
ccp4srs::PMonomer monomer;
char              S[1000];
char              L[1000];
char              ErrLog[1000];
int               lcount,rc;

  compFile.assign ( filePath,true,false );
  if (!compFile.reset(true))  {
    sprintf ( S," *** cannot open components file %s\n",filePath );
    prints  ( S,3 );
    n_rcsb_read_err++;
    return -1;
  } else  {

//    printf ( " --- try file '%s'\n",filePath );

    lcount     = 0;
    rc         = 0;
    S[0] = char(0);

    while (!compFile.FileEnd())  {

      rc = mmCIFData.ReadMMCIFData ( compFile,S,lcount );

      if (rc!=mmdb::mmcif::CIFRC_Ok)  {

        if (rc<0)  {
          sprintf ( S," *** error reading file (%i):\n"
                      "     %s\n",rc,mmdb::mmcif::GetCIFMessage(L,rc) );
          prints  ( S,3 );
          n_rcsb_read_err++;
        } else  {
          sprintf ( S," ... warning on reading file:\n"
                      "     %s\n",mmdb::mmcif::GetCIFMessage(L,rc) );
          prints  ( S,3 );
        }

      } else  {

        monomer = new ccp4srs::Monomer();
        rc = monomer->readFromCIF ( &mmCIFData,ErrLog );
        if (rc)  {
          sprintf ( S," +++ error reading monomer %s:\n"
                      "     %s\n",filePath,ErrLog );
          prints  ( S,3 );
          n_rcsb_read_err++;
          delete monomer;
        } else
          processMonomer ( monomer );

      }

    }

  }

  return rc;

}

int CMakeSRS::readEnergyLib ( mmdb::cpstr monLib )  {
mmdb::pstr S;
char       L[1000];
int        rc;

  S  = NULL;
  rc = 0;

  S = NULL;
  mmdb::CreateCopCat ( S,monLib,"ener_lib.cif" );
  energyLib.SetFlag  ( mmdb::mmcif::CIFFL_PrintWarnings );
  rc = energyLib.ReadMMCIFData ( S );
  if (rc)  {
    sprintf ( L," *** cannot read energy type library:\n"
                " file '%s':\n",S );
    prints  ( L,3 );
    delete[] S;
    S = new char[2000];
    sprintf ( L," %s\n",mmdb::mmcif::GetCIFMessage(S,rc) );
    prints  ( L,3 );
    delete S;
    n_elib_err++;
    return -111111;
  }

  sprintf ( L," ... energy library : %s\n",S );
  prints  ( L,1 );
  sprintf ( L,"                      %i energy type definitions\n",
              energyLib.GetLoopLength("_lib_atom") );
  prints  ( L,1 );

  return rc;

}

int CMakeSRS::scanRCSBLibrary ( mmdb::cpstr rcsbLib )  {
DIR    *dir;
#if USE_DIRENT_H
dirent *dp;
#else
direct *dp;
#endif
mmdb::pstr S;
char       L[1000];
int        rc;

  S  = NULL;
  rc = 0;

  dir = opendir(rcsbLib);
  if (!dir)  {
    mmdb::CreateCopy ( S,rcsbLib );
    if (rcsbLib[strlen(S)-1]==_dir_sep_c)
      S[strlen(S)-1] = char(0);
    rc = processRCSBFile ( S );
    nEntries0 = nEntries;
    if (rc==-1)  {
      sprintf ( L,"\n *** RCSB library '%s' does not exist\n",rcsbLib );
      prints  ( L,3 );
      n_rcsb_read_err++;
      return -111112;
    } else
      return rc;
  }

  dp = readdir ( dir );
  while (dp)  {
    if (dp->d_name[0]!='.')  {
      mmdb::CreateCopCat ( S,rcsbLib,dp->d_name );
      rc = processRCSBFile ( S );
    }
    dp = readdir ( dir );
  }

  nEntries0 = nEntries;

  closedir ( dir );

  return rc;

}


int CMakeSRS::scanMonomerLibrary ( mmdb::cpstr monLib,
                                   bool match_only )  {
DIR    *dir;
#if USE_DIRENT_H
dirent *dp;
#else
direct *dp;
#endif
mmdb::pstr S;
char       L[1000];
int        rc;

  matchOnly = match_only;

  S  = NULL;
  rc = 0;

  dir = opendir(monLib);
  if (!dir)  {
    sprintf ( L,"\n *** monomer library '%s' does not exist\n",monLib );
    prints  ( L,3 );
    n_ml_read_err++;
    return -111112;
  }

  dp = readdir ( dir );
  while (dp)  {
    if ((dp->d_name[0]!='.') && (!dp->d_name[1]))  {
      mmdb::CreateCopCat ( S,monLib,dp->d_name,_dir_sep );
      rc = scanMonLibDirectory ( S );
    }
    dp = readdir ( dir );
  }

  closedir ( dir );

  nEntries0 = nEntries;

  return rc;

}

void  CMakeSRS::sortIndex()  {
CSortSRSIndex sortSRSIndex;
  prints ( "\n ... sorting ...\n",1 );
  sortSRSIndex.Sort ( srsIndex,srsMonomer,srsGraph,nEntries );
}

int  CMakeSRS::writeSRS ( mmdb::cpstr srsDir )  {
mmdb::io::File graphFile;
mmdb::io::File structFile;
mmdb::io::File indexFile;
mmdb::pstr     S;
mmdb::pstr     buffer;
char           L[1000];
int            i,n,k, buffer_length;

  S = NULL;

  mmdb::CreateCopCat ( S,srsDir,ccp4srs::srsGraphFile );
  graphFile.assign ( S,false,true );
  if (!graphFile.rewrite())  {
    sprintf ( L," *** can't open srs file '%s' file.\n",S );
    prints  ( L,3 );
    return -111113;
  }

  mmdb::CreateCopCat ( S,srsDir,ccp4srs::srsStructFile );
  structFile.assign ( S,false,true );
  if (!structFile.rewrite())  {
    sprintf ( L," *** can't open srs file '%s' file.\n",S );
    prints  ( L,3 );
    return -111114;
  }

  mmdb::CreateCopCat ( S,srsDir,ccp4srs::srsIndexFile );
  indexFile.assign ( S,false,true );
  if (!indexFile.rewrite())  {
    sprintf ( L," *** can't open srs file '%s' file.\n",S );
    prints  ( L,3 );
    return -111115;
  }

  prints ( "\n ... writing SRS files: ",1 );
  k = mmdb::IMin(50,nEntries);
  for (i=0;i<k;i++)
    prints ( ".",1 );
  printf ( "\r ... writing SRS files: " );
  fflush ( stdout );

  n = nEntries/k;
  k = 0;


  memIO.setCompressionLevel ( 1 );

  indexFile.WriteInt ( &srsVersion );
  indexFile.WriteInt ( &nEntries   );

  for (i=0;i<nEntries;i++)  {

    srsIndex[i]->fGraphPos  = graphFile .Position();
    srsIndex[i]->fStructPos = structFile.Position();

    buffer = (char *)memIO.get_buffer ( 50000 );
    buffer_length = 0;
    if (srsGraph[i])
      srsGraph[i]->mem_write ( buffer,buffer_length );
    memIO.set_buffer_length ( buffer_length );
    memIO.write ( graphFile );

//    StreamWrite ( graphFile,srsGraph[i] );

    memIO.reset();
    srsMonomer[i]->check_ccp4_coordinates();
    srsMonomer[i]->write ( structFile,srsVersion,&memIO );

    memIO.reset();
    srsIndex[i]->write ( indexFile,srsVersion,&memIO );

    k++;
    if (k>=n)  {
      printf ( ">" );
      fflush ( stdout );
      k = 0;
    }

  }

  graphFile .shut();
  structFile.shut();
  indexFile .shut();

  prints ( "\n --- done.\n",1 );

  return 0;

}


int main ( int argc, char ** argv, char ** env )  {
UNUSED_ARGUMENT(env);
CMakeSRS makeSRS;
mmdb::pstr rcsbLib;
mmdb::pstr monLib;
mmdb::pstr oldMonLib;
mmdb::pstr srsDir;
char       S[1000];
int        argNo,rc;

  if (argc<2)  {
    printInstructions ( argv[0] );
    return 1;
  }

  mmdb::InitMatType();
  makeSRS.openLogs ( NULL );

  rcsbLib   = NULL;
  monLib    = NULL;
  oldMonLib = NULL;
  srsDir    = NULL;

  rc    = 0;
  argNo = 1;

  mmdb::CreateCopy ( rcsbLib,argv[argNo++] );
  if (rcsbLib[strlen(rcsbLib)-1]!=_dir_sep_c)
    mmdb::CreateConcat ( rcsbLib,_dir_sep );

  if (argNo<argc)  {
    mmdb::CreateCopy ( monLib,argv[argNo++] );
    if (monLib[strlen(monLib)-1]!=_dir_sep_c)
      mmdb::CreateConcat ( monLib,_dir_sep );
  }

  if (argNo<argc)  {
    if (!strcmp(argv[argNo],"-d"))  {
      argNo++;
      mmdb::CreateCopy ( srsDir,argv[argNo] );
    } else {
      mmdb::CreateCopy ( oldMonLib,argv[argNo++] );
      if (oldMonLib[strlen(oldMonLib)-1]!=_dir_sep_c)
        mmdb::CreateConcat ( oldMonLib,_dir_sep );
      if (argNo<argc)  {
        if (!strcmp(argv[argNo++],"-d"))
          mmdb::CreateCopy ( srsDir,argv[argNo] );
      }
    }
  }

  if (srsDir)  {
    if (srsDir[strlen(srsDir)-1]!=_dir_sep_c)
      mmdb::CreateConcat ( srsDir,_dir_sep );
  } else
    mmdb::CreateCopy ( srsDir,"" );

  sprintf ( S," ... SRS Generator version 1.0.1\n\n"
              " ... RCSB library            : %s\n",rcsbLib );
  makeSRS.prints ( S,1 );
  if (monLib) sprintf ( S,
              " ... CCP4 monomer library    : %s\n",monLib );
  makeSRS.prints ( S,1 );
  if (oldMonLib)  {
    sprintf ( S," ... old CCP4 monomer library: %s\n",oldMonLib );
    makeSRS.prints ( S,1 );
  }
  rc = makeSRS.readEnergyLib ( monLib );
  if (rc)  return rc;

  sprintf ( S," ... SRS directory           : %s\n"
              "\n ... machine properties:\n"
              "      max_short_real = %10.3f\n"
              "      sizeof(short)  = %lu\n"
              "      sizeof(int)    = %lu\n"
              "      sizeof(long)   = %lu\n\n",
              srsDir,mmdb::MaxShortReal,
              sizeof(short),sizeof(int),sizeof(long)  );
  makeSRS.prints ( S,1 );

  rc = makeSRS.scanRCSBLibrary ( rcsbLib );
  makeSRS.sortIndex();

  if (monLib)  {
    rc = makeSRS.scanMonomerLibrary ( monLib,false );
    makeSRS.sortIndex();
    if (oldMonLib)  {
      rc = makeSRS.scanMonomerLibrary ( oldMonLib,true );
      makeSRS.sortIndex();
    }
  }

  makeSRS.writeSRS ( srsDir );
  makeSRS.closeLogs();

  if (rcsbLib)   delete[] rcsbLib;
  if (monLib)    delete[] monLib;
  if (oldMonLib) delete[] oldMonLib;
  if (srsDir)    delete[] srsDir;

  return rc;

}

