//  $Id: ccp4srs_base.cpp $
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
//  **** Module  :  ccp4srs_base  <implementation>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Base - SRS base manager class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2014
//
//  =================================================================
//

#include <string.h>
#include <stdlib.h>

#include "ccp4srs_base.h"
#include "ccp4srs_defs.h"
#include "mmdb2/mmdb_tables.h"

namespace ccp4srs  {

  Base::Base()  {
    init_base();
  }

  Base::~Base()  {
    freeMemory();
  }

  void Base::init_base()  {

    index       = NULL;
    srsDir      = NULL;
    srsVersion  = CCP4SRS_file_version;
    nEntries    = 0;
    index_alloc = 0;

  }

  void Base::freeMemory()  {
  int i;

    if (index)  {
      for (i=0;i<index_alloc;i++)
        if (index[i])  delete index[i];
      delete[] index;
      index = NULL;
    }
    nEntries    = 0;
    index_alloc = 0;

    if (srsDir)  delete[] srsDir;
    srsDir = NULL;

  }

  int Base::loadIndex ( mmdb::cpstr srsPath )  {
  mmdb::io::File f;
  MemIO          memIO;
  mmdb::pstr     S;
  int            i;

    freeMemory();

    if (srsPath)  {
      mmdb::CreateCopy ( srsDir,srsPath );
      if (srsPath[0])  {
        if (srsPath[strlen(srsPath)-1]!=mmdb::io::_dir_sep_c)
          mmdb::CreateConcat ( srsDir,mmdb::io::_dir_sep );
      }
    }

    S = NULL;
    f.assign ( getPath(S,srsIndexFile),false,true );
    if (S)  delete[] S;

    if (f.reset(true))  {

      f.ReadInt ( &srsVersion  );
      f.ReadInt ( &index_alloc );

      if (index_alloc>0)  {

        index = new PIndex[index_alloc];
        for (i=0;(i<index_alloc) && (!f.FileEnd()) && f.Success();i++)  {
          memIO.reset();
          index[nEntries] = new Index();
          if (index[i]->read(f,srsVersion,&memIO)!=CCP4SRS_Ok)
                delete index[nEntries];
          else  nEntries++;
        }

        for (i=nEntries;i<index_alloc;i++)
          index[i] = NULL;

      }

      f.shut();

      if (nEntries==index_alloc)  return CCP4SRS_Ok;
      return CCP4SRS_IndexCorrupt;

    } else
      return CCP4SRS_FileNotFound;

  }

  PIndex Base::getIndex ( mmdb::cpstr entryID )  {
  int n = getEntryNo ( entryID );
    if ((n>=0) && (n<nEntries))  return index[n];
    return NULL;
  }

  PIndex Base::getIndex ( int n )  {
    if ((n>=0) && (n<nEntries))  return index[n];
    return NULL;
  }


  int  Base::loadStructure ( mmdb::cpstr entryID )  {
  UNUSED_ARGUMENT(entryID);
    return CCP4SRS_Ok;
  }

  int  Base::unloadStructure ( mmdb::cpstr entryID )  {
  UNUSED_ARGUMENT(entryID);
    return CCP4SRS_Ok;
  }


  mmdb::pstr Base::getPath ( mmdb::pstr & S, mmdb::cpstr FName )  {
    if (srsDir)  mmdb::CreateCopCat ( S,srsDir,FName );
           else  mmdb::CreateCopy   ( S,FName );
    return S;
  }


  int Base::getEntryNo ( mmdb::cpstr entryID )  {
  mmdb::ResName id;
  int           l1,l2,l,k;

    mmdb::strcpy_css ( id,entryID );

    k = -1;

    if (nEntries<=3)  {
      for (l=0;(l<nEntries) && (k<0);l++)
        if (!strcmp(index[l]->entryID,id))
          k = l;
      if (k>=0)  return k;
      return CCP4SRS_EntryNotFound;
    }

    l  = 0;
    l1 = 0;
    l2 = nEntries-1;
    while (l1<l2-1)  {
      l = (l1+l2)/2;
      k = strcmp ( index[l]->entryID,id );
      if (k<0)       l1 = l;
      else if (k>0)  l2 = l;
      else {
        l1 = l;
        l2 = l;
      }
    }

    if (k==0)  return l;
    else if (l==l1)  {
      if (!strcmp(index[l2]->entryID,id))  return l2;
    } else if (l==l2)  {
      if (!strcmp(index[l1]->entryID,id))  return l1;
    }

    return CCP4SRS_EntryNotFound;

  }

  mmdb::io::PFile Base::getStructFile()  {
  mmdb::io::PFile structFile;
  mmdb::pstr      S;
    structFile = new mmdb::io::File();
    S = NULL;
    structFile->assign ( getPath(S,srsStructFile),false,true );
    if (S)  delete[] S;
    if (!structFile->reset(true))  {
      delete structFile;
      structFile = NULL;
    }
    return structFile;
  }

  mmdb::io::PFile Base::getGraphFile()  {
  mmdb::io::PFile graphFile;
  mmdb::pstr      S;
    graphFile = new mmdb::io::File();
    S = NULL;
    graphFile->assign ( getPath(S,srsGraphFile),false,true );
    if (S)  delete[] S;
    if (!graphFile->reset(true))  {
      delete graphFile;
      graphFile = NULL;
    }
    return graphFile;
  }


  PMonomer Base::getMonomer ( mmdb::cpstr entryID )  {
  //   getMonomer returns pointer to the monomer structure
  // identified by 3-letter entryID. If such structure is not
  // found, the function returns NULL.
  //   The function returns a pointer to a private copy of the
  // structure. Modifying it will not change data in the structural
  // database. The application is responsible for deallocating
  // the structure after use (simply use delete).
  //   See description of CCP4SRSMonomer for the explanation of
  // its fields.
  mmdb::io::PFile structFile;
  PMonomer        srsMonomer;
  int             entryNo;

    srsMonomer = NULL;
    entryNo    = getEntryNo ( entryID );
    if (entryNo!=CCP4SRS_EntryNotFound)  {
      structFile = getStructFile();
      if (structFile)  {
        structFile->seek ( index[entryNo]->fStructPos );
        srsMonomer = new Monomer();
        memIO.reset();
        if (srsMonomer->read(*structFile,srsVersion,&memIO)!=CCP4SRS_Ok) {
          delete srsMonomer;
          srsMonomer = NULL;
        }
        structFile->shut ();
        delete structFile;
      }
    }

    return srsMonomer;

  }

  PMonomer Base::getMonomer ( int entryNo, mmdb::io::PFile structFile ) {
  mmdb::io::PFile sFile;
  PMonomer        srsMonomer;

    srsMonomer = NULL;
    if ((0<=entryNo) && (entryNo<nEntries))  {
      if (!structFile)  sFile = getStructFile();
                  else  sFile = structFile;
      if (sFile)  {
        sFile->seek ( index[entryNo]->fStructPos );
        srsMonomer = new Monomer();
        memIO.reset();
        if (srsMonomer->read(*sFile,srsVersion,&memIO)!=CCP4SRS_Ok) {
          delete srsMonomer;
          srsMonomer = NULL;
        }
        if (!structFile)  delete sFile;
      }
    }

    return srsMonomer;

  }

  PMonomer Base::getMonomer ( mmdb::cpstr entryID,
                              mmdb::io::PFile structFile )  {
  //   Another form of getMonomer(..) uses an open structure
  // file, which allows to save on opening/closing file if
  // multiple access to SRS structures is required.
  int entryNo;

    entryNo    = getEntryNo ( entryID );
    if (entryNo!=CCP4SRS_EntryNotFound)
          return getMonomer ( entryNo,structFile );
    else  return NULL;

  }

  int  Base::getGraph ( mmdb::io::PFile graphFile,
                        mmdb::math::RPGraph G )  {
  mmdb::pstr buffer;
  int        buffer_length,rc,n;

    rc = CCP4SRS_Ok;

    memIO.reset();
    if (memIO.read(*graphFile)==MemIO::Ok)  {
      memIO.get_buffer ( &buffer,&buffer_length );
      if (buffer_length<=0)  {
        if (G)  delete G;
        G = NULL;
      } else  {
        if (!G)  G = new mmdb::math::Graph();
        n = 0;
        G->mem_read ( buffer,n );
        if (n!=buffer_length)  {
          delete G;
          G  = NULL;
          rc = CCP4SRS_ReadErrors;
        }
      }
    } else  {
      rc = CCP4SRS_ReadErrors;
      if (G)  delete G;
      G = NULL;
    }

    return rc;

  }

  int  Base::getGraph ( mmdb::io::PFile graphFile, int entryNo,
                               mmdb::math::RPGraph G, int Hflag )  {
  //   getGraph(..) retrieves data for chemical structure number
  // entryNo (as in Index) from graph file graphFile, then allocates
  // and builds the corresponding graph, which is returned in G.
  //   If Hflag is set to 1, all hydrogens are removed from the graph.
  //   If Hflag is set to 2, element types of atoms, to which hydrogens
  // are bonded, are modified with flag HYDROGEN_BOND and moved to
  // the end.
  //   Returns 0 in case of success.
  int rc,htype;

    if ((entryNo<0) || (entryNo>=nEntries))
      return CCP4SRS_EntryNotFound;

    rc = CCP4SRS_Ok;

  /*
    if (Index[structNo]->loadPos>=0)  {

      if (!G)  G = new CGraph();
      G->Copy ( ldGraph[Index[structNo]->loadPos] );

    } else  {
  */
    graphFile->seek ( index[entryNo]->fGraphPos );
    graphFile->SetSuccess();
    rc = getGraph ( graphFile,G );
  /*
    if (memIO.read(*graphFile)==CMemIO::Ok)  {
  //  StreamRead ( *graphFile,G );
      pstr buffer;
      int  buffer_length,n;
      memIO.get_buffer ( &buffer,&buffer_length );
      if (buffer_length<=0)  {
        if (G)  delete G;
        G = NULL;
      } else  {
        if (!G)  G = new CGraph();
        n = 0;
        G->mem_read ( buffer,n );
        if (n!=buffer_length)  {
          delete G;
          G  = NULL;
          rc = CCP4SRS_ReadErrors;
        }
      }
    } else  {
  //  if (!graphFile->Success())  {
      rc = CCP4SRS_ReadErrors;
      if (G)  delete G;
      G = NULL;
    }
  */

  /*
    }
  */

    if (G)  {
      G->MakeVertexIDs();
      if (Hflag>=1)  {
        htype = mmdb::getElementNo(mmdb::pstr("H"));
        if (Hflag==2)  G->HideType    ( htype );
                 else  G->ExcludeType ( htype );
      }
      G->Build ( false );
    }

    return rc;

  }

  int  Base::getGraph ( mmdb::io::PFile graphFile,
                               mmdb::math::RPGraph G,
                               int Hflag )  {
  //   getGraph(..) retrieves data for chemical structure, which is
  // next in the graph file, then allocates and builds the corresponding
  // graph, which is returned in G.
  //   If Hflag is set to 1, all hydrogens are removed from the graph.
  //   If Hflag is set to 2, element types of atoms, to which hydrogens
  // are bonded, are modified with flag HYDROGEN_BOND and moved to
  // the end.
  //   Returns 0 in case of success.
  int rc,htype;

  //  rc = 0;

    graphFile->SetSuccess();
    rc = getGraph ( graphFile,G );
  //  StreamRead ( *graphFile,G );

    if (!graphFile->Success())  {
      rc = CCP4SRS_ReadErrors;
      if (G)  delete G;
      G = NULL;
    }

    if (G)  {
      G->MakeVertexIDs();
      if (Hflag>=1)  {
        htype = mmdb::getElementNo(mmdb::pstr("H"));
        if (Hflag==2)  G->HideType    ( htype );
                 else  G->ExcludeType ( htype );
      }
      G->Build ( false );
    }

    return rc;

  }

  int  Base::getGraph ( int entryNo, mmdb::math::RPGraph G,
                               int Hflag ) {
  mmdb::io::PFile graphFile;
  int             rc,htype;

    if ((0<=entryNo) && (entryNo<nEntries))  {
      /*
      if (Index[structNo]->loadPos>=0)  {
        if (!G) G = new CGraph();
        G->Copy ( ldGraph[Index[structNo]->loadPos] );
      } else  {
      */
      graphFile = getGraphFile();
      if (graphFile)  {
        graphFile->seek ( index[entryNo]->fGraphPos );
        rc = getGraph ( graphFile,G );
  //      StreamRead  ( *graphFile,G );
        graphFile->shut ();
        delete graphFile;
        if (rc!=CCP4SRS_Ok)
          return rc;
      } else  {
        if (G)  delete G;
        G = NULL;
        return CCP4SRS_FileNotFound;
      }
      /*
      }
      */
    } else  {
      if (G)  delete G;
      G = NULL;
      return CCP4SRS_EntryNotFound;
    }

    if (G)  {
      G->MakeVertexIDs();
      if (Hflag>=1)  {
        htype = mmdb::getElementNo(mmdb::pstr("H"));
        if (Hflag==2)  G->HideType    ( htype );
                 else  G->ExcludeType ( htype );
      }
      G->Build ( false );
    }

    return CCP4SRS_Ok;

  }

  int  Base::getGraph ( mmdb::cpstr entryID, mmdb::math::RPGraph G, int Hflag ) {
  mmdb::io::PFile graphFile;
  int    rc,entryNo,htype;

    entryNo = getEntryNo ( entryID );

    if (entryNo!=CCP4SRS_EntryNotFound)  {
      /*
      if (Index[structNo]->loadPos>=0)  {
        if (!G )  G = new CGraph();
        G->Copy ( ldGraph[Index[structNo]->loadPos] );
      } else  {
      */
      graphFile = getGraphFile();
      if (graphFile)  {
        graphFile->seek ( index[entryNo]->fGraphPos );
        rc = getGraph ( graphFile,G );
  //      StreamRead      ( *graphFile,G );
        graphFile->shut ();
        delete graphFile;
        if (rc!=CCP4SRS_Ok)
          entryNo = CCP4SRS_FileNotFound;
      } else  {
        entryNo = CCP4SRS_FileNotFound;
        if (G)  delete G;
        G = NULL;
      }
      /*
      }
      */
    } else  {
      if (G)  delete G;
      G = NULL;
    }

    if (G)  {
      G->MakeVertexIDs();
      if (Hflag>=1)  {
        htype = mmdb::getElementNo(mmdb::pstr("H"));
        if (Hflag==2)  G->HideType    ( htype );
                 else  G->ExcludeType ( htype );
      }
      G->Build ( false );
    }

    return  entryNo;

    /*
    if (structNo!=SBASE_StructNotFound)  {
      return  GetGraph ( structNo,G,Hflag );
    } else  {
      if (G)  delete G;
      G = NULL;
      return  structNo;
    }
    */

  }


  int Base::getNofAtoms ( mmdb::cpstr entryID )  {
  // Returns the number of atoms in the references monomer. 'entryID' is
  // the 3-letter cif code. Non-negative return delivers the number of
  // atoms, otherwise CCP4SRS_EntryNotFound indicates that the requested
  // structure was not found in the database.
  int entryNo;
    entryNo = getEntryNo ( entryID );
    if (entryNo>=0)
      return index[entryNo]->nAtoms;
    return entryNo;
  }

  int Base::getNofAtoms ( int entryNo )  {
  // Returns the number of atoms in the references monomer. 'entryNo' is
  // the entry number (0...n_entries()-1). Non-negative return delivers
  // the number of atoms, otherwise CCP4SRS_EntryNotFound indicates that
  // the requested structure was not found in the database
    if ((0<=entryNo) && (entryNo<nEntries))
      return index[entryNo]->nAtoms;
    return CCP4SRS_EntryNotFound;
  }


  int Base::getNofAtoms ( int entryNo,  int & nNonHAtoms,
                                 int & nHAtoms )  {
  mmdb::pstr  p1,p2;

    nNonHAtoms = 0;
    nHAtoms    = 0;

    if ((entryNo<0) || (entryNo>=nEntries))
      return CCP4SRS_EntryNotFound;

    nNonHAtoms = index[entryNo]->nAtoms;
    nHAtoms    = 0;
    if (index[entryNo]->Comp1)  {
      p1 = strstr ( index[entryNo]->Comp1,"H(" );
      if (p1)  {
        p1 += 2;
        p2  = p1;
        while ((*p2) && (*p2!=')'))  p2++;
        if (*p2==')')  {
          *p2 = char(0);
          nHAtoms = mmdb::mround(strtod(p1,NULL));
          *p2 = ')';
        }
      }
    }

    nNonHAtoms -= nHAtoms;

    return CCP4SRS_Ok;

  }


  int Base::getAtNames ( int  entryNo, mmdb::PAtomName AtName,
                                int & nAtoms, int & nH )  {
  mmdb::io::PFile structFile;
  int             rc;
    nAtoms = 0;
    nH     = 0;
    if ((entryNo<0) || (entryNo>=nEntries))
      return CCP4SRS_EntryNotFound;
    structFile = getStructFile();
    if (structFile)  {
      rc = getAtNames ( structFile,entryNo,AtName,nAtoms,nH );
      structFile->shut();
      delete structFile;
    } else
      return CCP4SRS_FileNotFound;
    return rc;
  }


  int Base::getAtNames ( mmdb::io::PFile structFile, int entryNo,
                         mmdb::PAtomName AtName,
                         int & nAtoms, int & nH )  {
  PMonomer       monomer;
  PAtom          atom;
  mmdb::AtomName aname;
  int            i,j;
  mmdb::pstr     p1,p2;
  bool           removeHydrogens;

    if ((entryNo<0) || (entryNo>=nEntries))
      return CCP4SRS_EntryNotFound;

    removeHydrogens = (nH==-1);

    monomer = getMonomer ( entryNo,structFile );
    if (!monomer)
      return CCP4SRS_ReadErrors;

    nAtoms = index[entryNo]->nAtoms;
    nH     = 0;
    if (index[entryNo]->Comp1)  {
      p1 = strstr ( index[entryNo]->Comp1,"H(" );
      if (p1)  {
        p1 += 2;
        p2  = p1;
        while ((*p2) && (*p2!=')'))  p2++;
        if (*p2==')')  {
          *p2 = char(0);
          nH  = mmdb::mround(strtod(p1,NULL));
          *p2 = ')';
        }
      }
    }

    if (removeHydrogens)  {
      j = 0;
      for (i=0;i<nAtoms;i++)  {
        atom = monomer->atom(i);
        if (atom)  {
          if ((atom->element()[0]!='H') || (atom->element()[1]))
            strcpy ( AtName[j++],atom->name_pdb(aname) );
        }
      }
    } else  {
      for (i=0;i<nAtoms;i++)  {
        atom = monomer->atom(i);
        if (atom)
          strcpy ( AtName[i],atom->name_pdb(aname) );
      }
    }

    delete monomer;

    return CCP4SRS_Ok;

  }


  int _makeChirInd ( char chirality )  {
    if (chirality=='S')  return -1;
    if (chirality=='R')  return +1;
    return 0;
  }

  int _findName ( mmdb::PAtomName Nams, mmdb::cpstr N, int len )  {
  mmdb::AtomName Nam;
  int            i;
    i = 0;
    while (i<len)  {
      strcpy ( Nam,Nams[i] );
      if (!strcmp(N,Nam))  break;
      i++;
    }
    if (i>=len)  i = -1;
    return i;
  }

  int Base::getAtoms ( mmdb::cpstr     entryID,
                              int           & nNonHAtoms,
                              mmdb::PAtomName NonHAtName,
                              int           & nHAtoms,
                              mmdb::PAtomName HAtName,
                              mmdb::ivector   Hconnect,
                              mmdb::ivector   Elem,
                              mmdb::ivector   Chiral )  {
  PMonomer       monomer;
  PAtom          atom,atom2;
  PBond          bond;
  mmdb::AtomName aname;
  int            i,j,entryNo,rc;

    nNonHAtoms = 0;
    nHAtoms    = 0;

    entryNo = getEntryNo ( entryID );
    if (entryNo<0)  return CCP4SRS_EntryNotFound;

    monomer = getMonomer ( entryNo,NULL );
    if (!monomer)
      return CCP4SRS_ReadErrors;

    rc = CCP4SRS_Ok;

    for (i=0;i<index[entryNo]->nAtoms;i++)  {
      atom = monomer->atom(i);
      if (atom)  {
        Elem[i] = mmdb::getElementNo ( atom->element() );
        if ((atom->element()[0]=='H') && (!atom->element()[1]))  {
          strcpy ( HAtName[nHAtoms],atom->name_pdb(aname) );
          nHAtoms++;
        } else  {
          strcpy ( NonHAtName[nNonHAtoms],atom->name_pdb(aname) );
          nNonHAtoms++;
        }
        Chiral[i] = _makeChirInd ( atom->chirality() );
      }
    }

    if (nHAtoms>0)  {
      for (j=0;j<nHAtoms;j++)
        Hconnect[j] = -1;
      for (i=0;i<index[entryNo]->nBonds;i++)  {
        bond = monomer->bond(i);
        if (bond)  {
          atom  = monomer->atom(bond->atom1());
          atom2 = monomer->atom(bond->atom2());
          if (atom && atom2)  {
            j = _findName ( HAtName,atom->name_pdb(aname),nHAtoms );
            if (j>=0)
              Hconnect[j] = _findName ( NonHAtName,
                                        atom2->name_pdb(aname),
                                        nNonHAtoms );
            else  {
              j = _findName ( HAtName,atom2->name_pdb(aname),nHAtoms );
              if (j>=0)
                Hconnect[j] = _findName ( NonHAtName,
                                          atom->name_pdb(aname),
                                          nNonHAtoms );
            }
          }
        }
      }
      j = 0;
      while ((j<nHAtoms) && (Hconnect[j]>=0)) j++;
      if (j<nHAtoms)
        rc = CCP4SRS_BrokenBonds;
    }

    delete monomer;

    return rc;

  }


  int  Base::getBonds ( mmdb::cpstr   entryID,
                               mmdb::ivector nBonds,
                               mmdb::imatrix bondPair,
                               int         & nAtoms,
                               int           maxNAtoms,
                               int           maxNBonds )  {
  mmdb::math::PGraph G;
  mmdb::math::PEdge  edge;
  int                entryNo,i, a1,a2,a;

    for (i=0;i<maxNAtoms;i++)
      nBonds[i] = 0;

    nAtoms = 0;

    entryNo = getEntryNo ( entryID );
    if (entryNo<0)  return CCP4SRS_EntryNotFound;

    nAtoms = index[entryNo]->nAtoms;
    if (nAtoms<=0)  return CCP4SRS_Ok;

    G = NULL;
    i = getGraph ( entryNo,G,0 );
    if (i!=CCP4SRS_Ok)  return i;
    if (!G)             return CCP4SRS_ReadErrors;

    for (i=1;i<=G->GetNofEdges();i++) {
      edge = G->GetEdge(i);
      if (edge)  {
        a1 = edge->GetVertex1() - 1;
        a2 = edge->GetVertex2() - 1;
        if (a1>a2)  {
          a = a1;  a1 = a2;  a2 = a;
        }
        if (nBonds[a1]<maxNBonds)  {
          bondPair[a1][nBonds[a1]] = a2;
          nBonds[a1]++;
        }
      }
    }

    delete G;

    return CCP4SRS_Ok;

  }


  int  Base::getHetInfo ( mmdb::cpstr       entryID,
                                 mmdb::pstr        Formula,
                                 mmdb::pstr        Hname,
                                 mmdb::pstr        Hsynonym,
                                 mmdb::pstr        Hcharge,
                                 mmdb::PAtomName & ClinkAtom,
                                 mmdb::PElement  & ClinkEle,
                                 mmdb::PAtomName & SlinkAtom,
                                 mmdb::PElement  & SlinkEle,
                                 int             & nLeavingAtoms )  {
  PMonomer       monomer;
  PAtom          atom;
  mmdb::AtomName aname;
  mmdb::Element  elem;
  mmdb::ivector  leavingAtom;
  mmdb::ivector  bondedAtom;
  int            i,entryNo;

    Formula [0] = char(0);
    Hname   [0] = char(0);
    Hsynonym[0] = char(0);
    Hcharge [0] = char(0);
    ClinkAtom   = NULL;
    ClinkEle    = NULL;
    SlinkAtom   = NULL;
    SlinkEle    = NULL;
    nLeavingAtoms = 0;

    entryNo = getEntryNo ( entryID );
    if (entryNo<0)  return CCP4SRS_EntryNotFound;

    monomer = getMonomer ( entryNo,NULL );
    if (!monomer)
      return CCP4SRS_ReadErrors;

    if (monomer->chem_formula()) strcpy ( Formula ,monomer->chem_formula() );
    if (monomer->chem_formula()) strcpy ( Hsynonym,monomer->chem_formula() );
    if (monomer->chem_type   ()) strcpy ( Hsynonym,monomer->chem_type   () );
    if (monomer->chem_name   ()) strcpy ( Hname   ,monomer->chem_name   () );
  //  if (monomer->Charge)  strcpy ( Hcharge ,monomer->Charge  );

    leavingAtom = NULL;
    bondedAtom  = NULL;
    monomer->getLeavingAtoms ( nLeavingAtoms,leavingAtom,bondedAtom );

    if (nLeavingAtoms>0)  {
      SlinkAtom = new mmdb::AtomName[nLeavingAtoms];
      SlinkEle  = new mmdb::Element [nLeavingAtoms];
      ClinkAtom = new mmdb::AtomName[nLeavingAtoms];
      ClinkEle  = new mmdb::Element [nLeavingAtoms];
      for (i=0;i<nLeavingAtoms;i++)  {
        atom = monomer->atom ( leavingAtom[i] );
        strcpy (SlinkAtom[i],atom->name_pdb(aname) );
        strcpy (SlinkEle [i],atom->element_pdb(elem) );
        if (bondedAtom[i]>=0)  {
          atom = monomer->atom ( bondedAtom[i] );
          strcpy(ClinkAtom[i],atom->name_pdb(aname) );
          strcpy(ClinkEle [i],atom->element_pdb(elem) );
        } else  {
          ClinkAtom[i][0] = char(0);
          ClinkEle [i][0] = char(0);
        }
      }
    }

    mmdb::FreeVectorMemory ( leavingAtom,0 );
    mmdb::FreeVectorMemory ( bondedAtom ,0 );
    delete monomer;

    return CCP4SRS_Ok;

  }

}  // namespace ccp4srs
