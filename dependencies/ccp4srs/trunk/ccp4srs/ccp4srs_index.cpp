//  $Id: ccp4srs_index.cpp $
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
//  **** Module  :  ccp4srs_index  <implementation>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Index - SRS data index class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2014
//
//  =================================================================
//

#include <string.h>

#include "mmdb2/mmdb_tables.h"
#include "ccp4srs_index.h"

namespace ccp4srs  {

  Index::Index()  {
    IndexInit();
  }

  Index::~Index() {
    if (Comp1)  delete[] Comp1;
    if (Comp2)  delete[] Comp2;
  }

  void Index::IndexInit ()  {
    entryID[0] = char(0);
    nAtoms     = 0;
    nNonHAtoms = 0;
    nBonds     = 0;
    fGraphPos  = -1;
    fStructPos = -1;
    Comp1      = NULL;
    Comp2      = NULL;
  }

  int Index::MakeCompositions ( PMonomer mon )  {
  //    MakeCompositions(..) makes the compositions strings for the
  //  given structure.
  //    A composition string consists of records E(N), where E stands
  //  for chemical element name, and N - for the number of atoms of this
  //  element in the structure. The records E(N) follow in alphabetical
  //  order of E without spaces and contain no spaces, records with N=0
  //  are excluded.
  //    Comp2 differs of Comp1 only if there are leaving atoms and
  //  represents the composition with leaving atoms taken into account.
  //  If there is no leaving atoms, Comp2==Comp1.
  //    The function returns the number of leaving atoms in the.
  //  structure.
  char          Cmp1[1000];
  char          Cmp2[1000];
  char          N[50];
  mmdb::ivector elem;
  int           i,k,x,l,nl;
  short         n0,n;

    nAtoms = mon->n_atoms();
    nBonds = mon->n_bonds();

    strcpy ( Cmp1,"" );
    strcpy ( Cmp2,"" );

    mmdb::GetVectorMemory ( elem,nAtoms,0 );
    nNonHAtoms = 0;
    for (i=0;i<nAtoms;i++)  {
      elem[i] = mmdb::getElementNo ( mon->atom(i)->element() );
      if (elem[i]>1)  nNonHAtoms++;
    }

    n0 = -1;
    n  = 0;
    nl = 0;
    for (i=0;(i<nAtoms) && (n<mmdb::nElementNames);i++)  {
      n = 10000;
      for (k=0;k<nAtoms;k++)
        if ((elem[k]>n0) && (elem[k]<n))  n = elem[k];
      if (n<mmdb::nElementNames)   {
        x = 0;
        l = 0;
        for (k=0;k<nAtoms;k++)
          if (elem[k]==n)  {
            x++;
            if (mon->atom(k)->leaving()=='Y')  l++;
          }
        nl += l;
        sprintf ( N,"%s(%i)",mmdb::ElementName[n-1],x );
        strcat  ( Cmp1,N );
        sprintf ( N,"%s(%i)",mmdb::ElementName[n-1],x-l );
        strcat  ( Cmp2,N );
        n0 = n;
      }
    }

    mmdb::FreeVectorMemory ( elem,0 );

    mmdb::CreateCopy ( Comp1,Cmp1 );
    mmdb::CreateCopy ( Comp2,Cmp2 );

    return nl;

  }



  int Index::write ( mmdb::io::RFile f, int version, PMemIO memIO ) {
  PMemIO mem_io;
  int    rc;

    if (memIO)  mem_io = memIO;
          else  mem_io = new MemIO();

    write_mem ( mem_io,version );
    rc = mem_io->write ( f );

    if (!memIO)
      delete mem_io;

    return rc;

  }

  int Index::read ( mmdb::io::RFile f, int version, PMemIO memIO )  {
  PMemIO mem_io;
  int    rc;

    if (memIO)  mem_io = memIO;
          else  mem_io = new MemIO();

    rc = mem_io->read ( f );
    if (rc==MemIO::Ok)  {
      if (!read_mem(mem_io,version))
        rc = MemIO::MemReadError;
    }

    if (!memIO)
      delete mem_io;

    return rc;

  }

  void Index::write_mem ( PMemIO memIO, int version )  {
  UNUSED_ARGUMENT(version);

    memIO->put_line    ( entryID    );
    memIO->put_integer ( nAtoms     );
    memIO->put_integer ( nNonHAtoms );
    memIO->put_integer ( nBonds     );
    memIO->put_long    ( fGraphPos  );
    memIO->put_long    ( fStructPos );
    memIO->put_string  ( Comp1      );
    memIO->put_string  ( Comp2      );

  }

  bool Index::read_mem ( PMemIO memIO, int version, bool * Ok )  {
  UNUSED_ARGUMENT(version);
  bool success;

    if (Ok)  success = *Ok;
       else  success = true;

    memIO->get_line    ( entryID   ,&success );
    memIO->get_integer ( nAtoms    ,&success );
    memIO->get_integer ( nNonHAtoms,&success );
    memIO->get_integer ( nBonds    ,&success );
    memIO->get_long    ( fGraphPos ,&success );
    memIO->get_long    ( fStructPos,&success );
    memIO->get_string  ( Comp1     ,&success );
    memIO->get_string  ( Comp2     ,&success );

    if (Ok) *Ok = success;
    return success;

  }

}  // namespace ccp4srs
