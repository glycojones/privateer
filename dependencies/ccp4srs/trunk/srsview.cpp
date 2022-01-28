//
// =================================================================
//
//  CCP4 Storage, Retrieval and Searches System (SRS)
//
//  File  srsview.cpp
//
//  SRS contents viewer.
//
//  Eugene Krissinel (2013)
//
// =================================================================
//

#include <string.h>

#include "ccp4srs/ccp4srs_manager.h"
#include "ccp4srs/ccp4srs_defs.h"



void printInstructions ( mmdb::cpstr argv0 )  {

  printf (
    "\n"
    " CCP4SRS viewer\n"
    " --------------\n"
    "\n"
    " USAGE:\n"
    "\n"
    "%s [-f srs_dir] -l\n"
    "\n"
    "lists all compounds found in SRS.\n"
    "\n"
    "%s [-f srs_dir] -c entry_id\n"
    "\n"
    "lists properties of compound 'entry_id'.\n"
    "\n"
    "%s [-f srs_dir] -e entry_id\n"
    "\n"
    "prints mmCIF extract for compound 'entry_id'.\n"
    "\n"
    "%s [-f srs_dir] -g entry_id\n"
    "\n"
    "lists properties of chemical graph for 'entry_id'.\n"
    "\n"
    "\n", argv0,argv0,argv0,argv0
   );

}

void printSmiles ( mmdb::cpstr title, ccp4srs::PSmiles smiles )  {
mmdb::cpstr  ver,plain,canonical;
  if (!smiles->isEmpty())  {
    ver       = smiles->getVersion();
    plain     = smiles->getSmiles();
    canonical = smiles->getCanonical();
    printf ( " Smiles from %s version %s\n",title,ver );
    if (plain)
      printf ( "      plain: %s\n",plain );
    if (canonical)
      printf ( "  canonical: %s\n",canonical );
  }
}

void printInChiI ( ccp4srs::PSmiles smiles )  {
mmdb::cpstr  ver,plain,canonical;
  if (!smiles->isEmpty())  {
    ver       = smiles->getVersion();
    plain     = smiles->getInChI();
    canonical = smiles->getInChIKey();
    printf ( " InChI version %s\n",ver );
    if (plain)
      printf ( "      InChI: %s\n",plain );
    if (canonical)
      printf ( "  InChI Key: %s\n",canonical );
  }
}

void printListOfEntries ( ccp4srs::PManager srs )  {
ccp4srs::PMonomer  Monomer;
int                l;

  printf ( "\n"
           "   ##   Code  EntryID   Natoms\n\n" );
  l = srs->n_entries();
  for (int i=0;i<l;i++)  {
    Monomer = srs->getMonomer ( i,NULL );
    if (Monomer)  {
      printf ( " %5i   %3s  %3s   %3i\n",i+1,
               Monomer->ID(),srs->getIndex(i)->entryID,
               Monomer->n_atoms() );
      delete Monomer;
      Monomer = NULL;
    } else
      printf ( " %5i   --- can't be retrieved (report as a bug).\n",
               i+1 );
  }

}

void printMonomerData ( ccp4srs::PMonomer Monomer )  {
ccp4srs::PAtom   Atom;
ccp4srs::PBond   Bond;
ccp4srs::PAngle  Angle;
mmdb::AtomName   aname;
mmdb::Element    elem;
mmdb::realtype   total_charge;
char             S[1000];
int              i;

  printf ( "\n Compound: %s\n"
           " ~~~~~~~~~~~~~\n\n",Monomer->ID() );

  printf ( "\n"
"Atom| ML |Chem| Energy |Chi| Lv|    X     |    Y     |    Z     \n"
" No |Name|elem|  Type  |   |   |          |          |          \n"
"====|====|====|========|===|===|==========|==========|==========\n"
         );
  for (i=0;i<Monomer->n_atoms();i++)  {
    Atom = Monomer->atom(i);
    printf ( " %3i|%4s| %s |%8s| %1c | %1c ",i+1,
             Atom->name_pdb   ( aname ),
             Atom->element_pdb( elem  ),
             Atom->energy_type(),
             Atom->chirality  (),
             Atom->leaving    () );
    if (Atom->x()>-mmdb::MaxShortReal/2.1)
          printf ( "|%10.3f|%10.3f|%10.3f\n",
                   Atom->x(),Atom->y(),Atom->z() );
    else  printf ( "|   ****   |   ****   |   ****   \n" );
  }


  printf ( "\n"
" X - RCSB | Y - RCSB | Z - RCSB | X - RCSB | Y - RCSB | Z - RCSB |    X     |    Y     |    Z     \n"
"  CARTN   |   CARTN  |   CARTN  |  IDEAL   |   IDEAL  |   IDEAL  |   CCP4   |   CCP4   |   CCP4   \n"
"==========|==========|==========|==========|==========|==========|==========|==========|==========\n"
         );
  for (i=0;i<Monomer->n_atoms();i++)  {
    Atom = Monomer->atom(i);
    if (Atom->x_rcsb_cartn()>-mmdb::MaxShortReal/2.1)
          printf ( "%10.3f|%10.3f|%10.3f",
                   Atom->x_rcsb_cartn(),Atom->y_rcsb_cartn(),
                   Atom->z_rcsb_cartn() );
    else  printf ( "   ****   |   ****   |   ****   " );
    if (Atom->x_rcsb_ideal()>-mmdb::MaxShortReal/2.1)
          printf ( "|%10.3f|%10.3f|%10.3f",
                   Atom->x_rcsb_ideal(),Atom->y_rcsb_ideal(),
                   Atom->z_rcsb_ideal() );
    else  printf ( "|   ****   |   ****   |   ****   " );
    if (Atom->x_ccp4_mlib()>-mmdb::MaxShortReal/2.1)
          printf ( "|%10.3f|%10.3f|%10.3f\n",
                   Atom->x_ccp4_mlib(),Atom->y_ccp4_mlib(),
                   Atom->z_ccp4_mlib() );
    else  printf ( "|   ****   |   ****   |   ****   \n" );
  }


  printf ( "\n"
"Atom| ML | HB |Valency| Energy |   VdW    |  VdW+H   |   Ion\n"
" No |Name|type|       |  Type  |  radius  |  radius  |  radius\n"
"====|====|====|=======|========|==========|==========|============\n"
         );
  for (i=0;i<Monomer->n_atoms();i++)  {
    Atom = Monomer->atom(i);
    printf ( " %3i|%4s|  %1c |  %3i  |%8s",i+1,
             Atom->name_pdb   ( aname ),
             Atom->hb_type    (),
             Atom->valency    (),
             Atom->energy_type() );
    printf ( "|%10.3f|%10.3f|%10.3f\n",Atom->vdw_radius(),
             Atom->vdwh_radius(),Atom->ion_radius() );
  }

  total_charge = 0.0;
  printf ( "\n"
"Atom| ML |RCSB| Old| Chi|RCSB|  Charge  | SP |\n"
" No |Name|Name|Name|    | Chi|          |    |\n"
"====|====|====|====|====|====|==========|====|\n" );
  for (i=0;i<Monomer->n_atoms();i++)  {
    Atom = Monomer->atom(i);
    printf ( " %3i|%4s",i+1,Atom->name_pdb(aname) );
    printf ( "|%4s",Atom->rcsb_name_pdb(aname) );
    printf ( "|%4s|  %1c |  %1c |%10.3f| %2i |\n",
             Atom->old_name_pdb(aname),
             Atom->chirality(),Atom->rcsb_chirality(),
             Atom->charge(),Atom->sp() );
    total_charge += Atom->charge();
  }
  printf (
"---------|------------------------------|----|\n"
"  Total: |                    %10.3f|    |\n",total_charge );

  printf ( "\n"
"Bond|Atom|Atom|Bond |  Length  | Length ESD\n"
" No | 1  | 2  |order|          |\n"
"====|=========|=====|==========|=============\n" );
  for (i=0;i<Monomer->n_bonds();i++)  {
    Bond = Monomer->bond(i);
    switch (Bond->order())  {
      case ccp4srs::Bond::noOrder : strcpy ( S,"nord" );  break;
      case ccp4srs::Bond::Single  : strcpy ( S,"sing" );  break;
      case ccp4srs::Bond::Aromatic: strcpy ( S,"arom" );  break;
      case ccp4srs::Bond::Double  : strcpy ( S,"doub" );  break;
      case ccp4srs::Bond::Triple  : strcpy ( S,"trip" );  break;
      case ccp4srs::Bond::Deloc   : strcpy ( S,"delc" );  break;
      case ccp4srs::Bond::Covalent: strcpy ( S,"covl" );  break;
      case ccp4srs::Bond::Metal   : strcpy ( S,"metl" );  break;
      default                     : strcpy ( S,"unkn" );
    }
    printf ( "%4i| %3i  %3i|%5s|%10.3f|%10.3f\n",
             i+1,Bond->atom1()+1,Bond->atom2()+1,S,
             Bond->length(),Bond->length_esd() );
  }

  if (Monomer->n_angles()>0)  {
    printf ( "\n"
"Angle|Atom|Atom|Atom|  Angle   |  Angle ESD\n"
"  No | 1  | 2  | 3  |          |\n"
"=====|==============|==========|=============\n" );
    for (i=0;i<Monomer->n_angles();i++)  {
      Angle = Monomer->angle(i);
      printf ( " %4i| %2i   %2i   %2i |%10.3f|%10.3f\n",
               i+1,Angle->atom1()+1,Angle->atom2()+1,
               Angle->atom3()+1,
               Angle->value(),Angle->esd() );
    }
  }


  printSmiles ( "ACDLabs",Monomer->getACDLabs() );
  printSmiles ( "CACTVS" ,Monomer->getCACTVS () );
  printSmiles ( "OpenEye",Monomer->getOpenEye() );
  printInChiI ( Monomer->getInChI() );

  /*
  if (Monomer->nSBS->nLeavingAtoms>0)  {
    printf ( "\n Leaving atoms:\n" );
    for (i=0;i<SBS->nLeavingAtoms;i++)
      printf ( " %2i.  %3i -> %3i\n",i+1,SBS->leavingAtom[i],
               SBS->bondedAtom[i] );
  }
  */
}

void extractMonomerData ( ccp4srs::PMonomer Monomer )  {
mmdb::mmcif::PData mmCIFData = Monomer->makeCIF();
mmdb::io::File     file;

  if (mmCIFData)  {
    file.assign ( "stdout",true );
    if (file.rewrite())  {
      mmCIFData->WriteMMCIF ( file );
      file.shut();
    } else
      printf ( " *** cannot write into standard output.\n" );
    delete mmCIFData;
  }

}

void printGraphData ( mmdb::math::PGraph Graph )  {

  printf ( "\n"
    " Number of vertices: %i\n"
    " Number of edges;    %i\n",
    Graph->GetNofVertices(),
    Graph->GetNofEdges   () );

}

int main ( int argc, char ** argv, char ** env )  {
UNUSED_ARGUMENT(env);
ccp4srs::PManager  srs;
ccp4srs::PMonomer  Monomer;
mmdb::math::PGraph Graph;
char               S[1000];
int                rc,i, argNo;

  if (argc<2)  {
    printInstructions ( argv[0] );
    return 1;
  }

  //  1.  Make routine initializations, which must always be done
  //      before working with MMDB
  mmdb::InitMatType();

  //  2.  Load SRS index
  if (!strcmp(argv[1],"-f"))  {
    if (argc<4)  {
      printInstructions ( argv[0] );
      return 1;
    }
    strcpy ( S,argv[2] );
    argNo = 3;
  } else if (argc<2)  {
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

  Monomer = NULL;
  Graph   = NULL;


  if (!strcmp(argv[argNo],"-l"))  {

    printListOfEntries ( srs );

  } else if (!strcmp(argv[argNo],"-c"))  {

    Monomer = srs->getMonomer ( argv[argNo+1],NULL );
    if (!Monomer)
          printf ( " --- compound %s not found in SRS.\n",argv[argNo+1]);
    else  printMonomerData ( Monomer );

  } else if (!strcmp(argv[argNo],"-e"))  {

    Monomer = srs->getMonomer ( argv[argNo+1],NULL );
    if (!Monomer)
          printf ( " --- compound %s not found in SRS.\n",argv[argNo+1]);
    else  extractMonomerData ( Monomer );

  } else if (!strcmp(argv[argNo],"-g"))  {

    Graph = NULL;  // must!
    i     = srs->getGraph ( argv[argNo+1],Graph,0 );

    if (i==ccp4srs::CCP4SRS_EntryNotFound)
      printf ( " --- compound %s not found in SRS.\n",argv[argNo+1]);
    else if (i==ccp4srs::CCP4SRS_FileNotFound)
      printf ( " --- graph descriptions are not found in SRS.\n" );
    else
      printGraphData ( Graph );

  }

  if (Monomer)  delete Monomer;
  if (Graph)    delete Graph;
  if (srs)      delete srs;

  return 0;

}



