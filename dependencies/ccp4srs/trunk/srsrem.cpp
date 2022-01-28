//
// =================================================================
//
//  CCP4 Storage, Retrieval and Searches System (SRS)
//
//  File  srsrem.cpp
//
//  SRS-based PDB remediator.
//
//  Eugene Krissinel (2011)
//
// =================================================================
//

#include <string.h>

#include "ccp4srs/ccp4srs_manager.h"
#include "ccp4srs/ccp4srs_defs.h"

#include "mmdb2/mmdb_tables.h"


void printInstructions ( mmdb::cpstr argv0 )  {

  printf (
    "\n"
    " CCP4SRS remediator\n"
    " ------------------\n"
    "\n"
    " USAGE:\n"
    "\n"
    "%s [-f srs_dir] pdbin.pdb pdbout.pdb\n"
    "\n"
    "\n", argv0
   );

}


mmdb::PManager getPDB ( mmdb::cpstr fileName )  {
mmdb::PManager mmdb;
char           S[1000];
int            lcount;
mmdb::ERROR_CODE rc;

  //  1.  Create an instance of MMDB
  mmdb = new mmdb::Manager();

  //  2.  Read coordinate file.
  //  2.1 Set all necessary read flags -- check with the top of
  //      file  mmdb_file.h  as needed
  mmdb->SetFlag ( mmdb::MMDBF_PrintCIFWarnings |
                  mmdb::MMDBF_FixSpaceGroup    |
                  mmdb::MMDBF_IgnoreDuplSeqNum |
                  mmdb::MMDBF_IgnoreHash       |
                  mmdb::MMDBF_IgnoreNonCoorPDBErrors );

  //  3. Read file
  rc = mmdb->ReadCoorFile ( fileName );

  //  4. Check for possible errors:
  if (rc) {
    //  An error was encountered. MMDB provides an error messenger
    //  function for easy error message printing.
    printf ( " ***** error #%i on reading file %s:\n\n %s\n\n",
             rc,fileName,
             mmdb::GetErrorDescription(rc) );
    //  Location of the error may be identified as precise as line
    //  number and the line itself (PDB only. Errors in mmCIF are
    //  located by category/item name). This information is now
    //  retrieved from MMDB input buffer:
    mmdb->GetInputBuffer ( S,lcount );
    if (lcount>=0)
      printf ( "       offending line #%i:\n%s\n\n",lcount,S );
    else if (lcount==-1)
      printf ( "       offending cif item: %s\n\n",S );
    //  dispose instance of MMDB and quit:
    delete mmdb;
    return NULL;
  } else  {
    //  MMDB allows to identify the type of file that has been just
    //  read:
    switch (mmdb->GetFileType())  {
      case mmdb::MMDB_FILE_PDB    : printf ( " PDB"         );  break;
      case mmdb::MMDB_FILE_CIF    : printf ( " mmCIF"       );  break;
      case mmdb::MMDB_FILE_Binary : printf ( " MMDB binary" );  break;
      default : printf ( " Unknown (report as a bug!)" );
    }
    printf ( " ... file %s has been read in.\n",fileName );
  }

  return mmdb;

}

int translateHydrogens ( ccp4srs::PMonomer monomer,
                         mmdb::PPAtom atom, mmdb::AtomName aname[],
                         mmdb::ivector elem, mmdb::ivector amap,
                         int nAtoms )  {
mmdb::AtomName a;
mmdb::realtype dx,dy,dz,d0;
mmdb::ivector  amapm;
int            i,j,k, mon_atom, rc;
bool           done;

  // calculate nearest neighbours for hydrogens

  for (i=0;i<nAtoms;i++)
    if (elem[i]==1)  {  // check only hydrogens
      d0 = mmdb::MaxReal;
      k  = -1;
      for (j=0;j<nAtoms;j++)
        if (elem[j]>1)  {
          dx = atom[j]->x - atom[i]->x;
          dy = atom[j]->y - atom[i]->y;
          dz = atom[j]->z - atom[i]->z;
          dx = dx*dx + dy*dy + dz*dz;
          if (dx<d0)  {
            d0 = dx;
            k  = j;
          }
        }
      if (k<0)  {
        printf ( " +++ monomer contains only hydrogens, translation"
                 " is not possible\n" );
        return -10;
      }
      amap[i] = k;  // nearest neighbour
    } else
      amap[i] = 0;  // will be used as counter

  // make all possible name-based renamings

  amapm = NULL;
  mmdb::GetVectorMemory ( amapm,monomer->n_atoms(),0 );
  for (j=0;j<monomer->n_atoms();j++)
    amapm[j] = 0;

  done = true;
  for (i=0;i<nAtoms;i++)
    if (elem[i]==1)  {  // check only hydrogens
      k = -1;
      for (j=0;(j<monomer->n_atoms()) && (k<0);j++)
        if (!strcmp(aname[i],monomer->atom(j)->name()))
          k = j;
        else if (!strcmp(aname[i],monomer->atom(j)->rcsb_name()))
          k = j;
        else if (!strcmp(aname[i],monomer->atom(j)->old_name()))
          k = j;
      if (k>=0)  {
        atom [i]->SetAtomName ( monomer->atom(k)->name_pdb(a) );
        strcpy ( aname[i],monomer->atom(k)->name() );
        amap [i] = -1; // mapped
        amapm[k] = -1; // used
      } else  {
        done = false;
      }
    }

  rc = 1;

  if (!done)  {

    for (i=0;i<nAtoms;i++)
      if ((elem[i]==1) && (amap[i]>=0))  {
        mmdb::strcpy_css ( a,atom[amap[i]]->name );
        mon_atom = monomer->atom_no ( a );
        if (mon_atom>=0)  {
          k = 0;
          do {
            j = monomer->get_bound ( mon_atom,"H",k++ );
            if (j>=0)  {
              if (amapm[j]>=0)  {
                atom [i]->SetAtomName ( monomer->atom(j)->name_pdb(a) );
                strcpy ( aname[i],monomer->atom(j)->name() );
                amap [i] = -1; // mapped
                amapm[j] = -1; // used
                j        = -2; // proper terminate
              }
            }
          } while (j>=0);
          if (j==-1)  {
            printf ( " +++ atom '%s' was not matched\n",aname[i] );
            rc = 2;
          }
        } else  {
          printf ( " +++ atom '%s' was not identified\n",
                   atom[amap[i]]->name );
          rc = 3;
        }
      }

  }

  mmdb::FreeVectorMemory ( amapm,0 );

  return rc;

}

int processMonomer ( ccp4srs::PMonomer monomer,
                     mmdb::PPAtom atom, mmdb::AtomName aname[],
                     mmdb::ivector elem, mmdb::ivector amap,
                     int nAtoms )  {
mmdb::AtomName a;
int            i,j,k;
bool           areHydrogens;

  areHydrogens = false;

  // named as in current Monomer Library?
  k = 1;
  for (i=0;(i<nAtoms) && (k>=0);i++)
    if (elem[i]>1)  {  // check only non-hydrogens
      k = -1;
      for (j=0;(j<monomer->n_atoms()) && (k<0);j++)
        if (!strcmp(aname[i],monomer->atom(j)->name()))
          k = j;
    } else
      areHydrogens = true;

  if (k>=0)  {
    if (areHydrogens)
      return translateHydrogens ( monomer,atom,aname,
                                  elem,amap,nAtoms );
    return 0;  // no translation needed
  }

  // named as in PDB 3?
  k = 1;
  for (i=0;(i<nAtoms) && (k>=0);i++)
    if (elem[i]>1)  {
      k = -1;
      for (j=0;(j<monomer->n_atoms()) && (k<0);j++)
        if (!strcmp(aname[i],monomer->atom(j)->rcsb_name()))  {
          k       = j;
          amap[i] = k;
         }
    }

  // named as in old Monomer Library?
  if (k<0)  {
    k = 1;
    for (i=0;(i<nAtoms) && (k>=0);i++)
      if (elem[i]>1)  {
        k = -1;
        for (j=0;(j<monomer->n_atoms()) && (k<0);j++)
          if (!strcmp(aname[i],monomer->atom(j)->old_name()))  {
            k       = j;
            amap[i] = k;
          }
      }

  }

  if (k<0)
    return -1;  // equivalence with old names not found

  for (i=0;i<nAtoms;i++)
    if (elem[i]>1)
      atom[i]->SetAtomName ( monomer->atom(amap[i])->name_pdb(a) );

  if (areHydrogens)
    return translateHydrogens ( monomer,atom,aname,
                                elem,amap,nAtoms );

  return 1;  // atoms renamed

}



mmdb::cpstr special_name[] = {
    "Ar", "A",
    "Cr", "C",
    "Tr", "T",
    "Gr", "G",
    "Ur", "U",
    "Ad", "DA",
    "Cd", "DC",
    "Td", "DT",
    "Gd", "DG",
    "Ud", "DU",

    "TRY", "TRP",
    "DAL", "ALA-D",
    "DVA", "VAL-D",
    "DPN", "PHE-D",
    "DPR", "PRO",
    "DTR", "TRP-D",
    "DTH", "THR",
    "DLE", "LEU-D",
    "+A",  "DA",
    "+C",  "DC",
    "+G",  "DG",
    "+T",  "DT",
    "+U",  "DU",
    "XLS", "XYL",
    "SUL", "SO4",
    "SO1", "SO4",
    "SO2", "SO4",
    "PHO", "PO4",
    "SPS", "P",
    "IPS", "P",
    "ZN1", "ZN",
    "ZN2", "ZN",
    "GMP", "3GP",
    "WAT", "HOH",
    "H2O", "HOH",
    "OH2", "HOH",
    "DOD", "HOH",
    "BOX", "BEZ",
    "AMT", "ATA",
    "GSH", "GTT",
    NULL
};


void translateResidue ( mmdb::PResidue res, ccp4srs::PManager srs,
                        mmdb::io::PFile structFile )  {
ccp4srs::PMonomer monomer;
mmdb::PPAtom      atom;
mmdb::AtomName   *aname;
mmdb::ivector     elem,amap;
char              S[1000];
mmdb::cpstr       p;
int               nAtoms,nAtoms0,nNonHAtoms0, i,rc;

  atom   = NULL;
  elem   = NULL;
  nAtoms = 0;
  res->GetAtomTable1 ( atom,nAtoms );
  aname       = new mmdb::AtomName[nAtoms];
  mmdb::GetVectorMemory ( elem,nAtoms,0 );
  mmdb::GetVectorMemory ( amap,nAtoms,0 );

  nAtoms0     = 0; // number of atoms in main location
  nNonHAtoms0 = 0;
  for (i=0;i<nAtoms;i++)  {
    mmdb::strcpy_css ( aname[i],atom[i]->name );
    elem[i] = mmdb::getElementNo ( atom[i]->element );
    if ((!atom[i]->altLoc[0]) || (atom[i]->altLoc[0]==' '))  {
      nAtoms0++;
      if (elem[i]>1)
        nNonHAtoms0++;
    }
  }

  rc = -1;

  // check for special cases first

  p = NULL;
  for (i=0;special_name[i] && (!p);i+=2)
    if (!strcmp(res->GetResName(),special_name[i]))
      p = special_name[i+1];

  if (p)  {
    monomer = srs->getMonomer ( p,structFile );
    if (!monomer)  {
      printf ( " +++ residue %s was not found in SRS\n",
               res->GetResidueID(S) );
      rc = -1;
    } else  {
      rc = processMonomer  ( monomer,atom,aname,elem,amap,nAtoms );
      if (rc>=0)  {
        printf ( " +++ residue %s was matched and renamed to %s\n",
                 res->GetResidueID(S),monomer->ID() );
        res->SetResName ( monomer->ID() );
      }
    }
  }

  // try own name

  if (rc<0)  {
    monomer = srs->getMonomer ( res->GetResName(),structFile );
    if (!monomer)  {
      printf ( " +++ residue %s was not found in SRS\n",
               res->GetResName() );
      rc = -1;
    } else  {
      rc = processMonomer ( monomer,atom,aname,elem,amap,nAtoms );
      if (rc==1)
        printf ( " ... residue %s was successfully translated\n",
                 res->GetResidueID(S) );
      else if (rc>1)
        printf ( " *** residue %s was translated, hydrogens need "
                 "to be checked manually\n",
                 res->GetResidueID(S) );
    }
  }

  /*
  // do full srs search
  for (i=0;(i<srs->n_entries()) && (rc<0);i++)  {
    index = srs->getIndex ( i );
    if ((nAtoms0<=index->nAtoms) &&
        (nNonHAtoms0<=index->nNonHAtoms))  {
      monomer = srs->getMonomer ( i,structFile );
      rc      = processMonomer  ( monomer,atom,aname,elem,amap,nAtoms );
      if (rc>=0)  {
        printf ( " +++ residue %s was matched and renamed to %s, "
                 "check all LINK records.\n",
                 res->GetResidueID(S),monomer->ID() );

        /// change resname everywhere else in PDB?

        res->SetResName ( monomer->ID() );
      }
    }
  }
  */

  if (rc<0)
    printf ( " *** residue %s could not be identified\n",
             res->GetResidueID(S) );

  mmdb::FreeVectorMemory ( amap,0 );
  mmdb::FreeVectorMemory ( elem,0 );

  delete[] aname;
  if (atom)  delete[] atom;

}


int main ( int argc, char ** argv, char ** env )  {
UNUSED_ARGUMENT(env);
ccp4srs::PManager srs;
mmdb::PManager    mmdb;
mmdb::PPResidue   res;
mmdb::io::PFile   structFile;
char              S[1000];
int               selHnd,nResidues;
int               argNo,rc, i;

  if (argc<3)  {
    printInstructions ( argv[0] );
    return 1;
  }

  //  1.  Make routine initializations, which must always be done
  //      before working with MMDB
  mmdb::InitMatType();

  // $CLIB/ccp4srs/srsdata

  //  2.  Load SRS index
  if (!strcmp(argv[1],"-f"))  {
    if (argc<4)  {
      printInstructions ( argv[0] );
      return 1;
    }
    strcpy ( S,argv[2] );
    argNo = 3;
  } else if (argc<3)  {
    printInstructions ( argv[0] );
    return 1;
  } else  {
    strcpy ( S,"" );
    argNo = 1;
  }

  srs = new ccp4srs::Manager();
  rc = srs->loadIndex ( S );
  if (rc!=ccp4srs::CCP4SRS_Ok)  {
    printf ( " *** *.srs files not found at '%s'.\n",S );
    delete srs;
    return 1;
  }

  // 3. Read file in
  mmdb = getPDB ( argv[argNo++] );
  if (!mmdb)
    return 1;  // unsuccessful read of file

  // 4. Select all residues
  selHnd = mmdb->NewSelection();
  mmdb->Select ( selHnd,mmdb::STYPE_RESIDUE,"*",mmdb::SKEY_NEW );

  // 5. Obtain selected residues
  mmdb->GetSelIndex ( selHnd,res,nResidues );

  printf ( " ... total %i residues selected\n",nResidues );

  structFile = srs->getStructFile();
  if (!structFile)  {
    printf ( " *** unable to open SRS structure file (error)\n" );
    return 2;
  }

  for (i=0;i<nResidues;i++)
    translateResidue ( res[i],srs,structFile );

  mmdb->WritePDBASCII ( argv[argNo] );

  delete structFile;
  delete mmdb;
  delete srs;

  return 0;

}



