//  $Id: ccp4srs_tree.cpp $
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
//  **** Module  :  ccp4srs_tree  <implementation>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Tree  - structure tree
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#include <string.h>

#include "ccp4srs_tree.h"
#include "ccp4srs_defs.h"

namespace ccp4srs  {

  Tree::Tree()  {
    n_sections    = 0;
    atom_tree     = NULL;  // tree atom
    atom_backward = NULL;  // previous atom
    atom_forward  = NULL;  // next atom
    connect_type  = NULL;  // type of connection
  }

  Tree::~Tree()  {
    empty();
  }

  void Tree::empty()  {
    mmdb::FreeVectorMemory ( atom_tree    ,0 );
    mmdb::FreeVectorMemory ( atom_backward,0 );
    mmdb::FreeVectorMemory ( atom_forward ,0 );
    n_sections = 0;
  }


  int Tree::readFromCIF ( mmdb::mmcif::PLoop mmCIFLoop,
                          Container<Atom> & atoms ) {
  mmdb::pstr p;
  int        rc,i;

    rc = 0;

    empty();

    if (mmCIFLoop)  {

      n_sections = mmCIFLoop->GetLoopLength();

      if (n_sections>0)  {

        mmdb::GetVectorMemory ( atom_tree,    n_sections,0 );
        mmdb::GetVectorMemory ( atom_backward,n_sections,0 );
        mmdb::GetVectorMemory ( atom_forward, n_sections,0 );
        mmdb::GetVectorMemory ( connect_type, n_sections,0 );
        for (i=0;(i<n_sections) && (!rc);i++)  {
          get_atom_from_cif ( atom_tree[i],mmCIFLoop,
                              MMCIF_TAG_ATOM_ID,i,atoms,rc );
          get_atom_from_cif ( atom_backward[i],mmCIFLoop,
                              MMCIF_TAG_ATOM_BACK,i,atoms,rc );
          get_atom_from_cif ( atom_forward[i],mmCIFLoop,
                              MMCIF_TAG_ATOM_FORWARD,i,atoms,rc );
          p = mmCIFLoop->GetString ( MMCIF_TAG_CONNECT_TYPE,i,rc );
          if (!p)   connect_type[i] = Tree::None;
          else if (!strcmp(p,"START"))
                    connect_type[i] = Tree::Start;
          else if (!strcmp(p,"END"))
                    connect_type[i] = Tree::End;
          else if (!strcmp(p,"ADD"))
                    connect_type[i] = Tree::Add;
          else      connect_type[i] = Tree::None;
        }

      }

    }

    return rc;

  }

  void  Tree::writeToCIF ( mmdb::mmcif::PData mmCIFData,
                           mmdb::cpstr monID,
                           Container<Atom> & atoms )  {
  mmdb::mmcif::PLoop loop;
  int                i;

    if (n_sections>0)  {

      mmCIFData->AddLoop ( MMCIF_LOOP_CHEM_COMP_TREE,loop );

      loop->AddLoopTag ( MMCIF_TAG_COMP_ID      );
      loop->AddLoopTag ( MMCIF_TAG_ATOM_ID      );
      loop->AddLoopTag ( MMCIF_TAG_ATOM_BACK    );
      loop->AddLoopTag ( MMCIF_TAG_ATOM_FORWARD );
      loop->AddLoopTag ( MMCIF_TAG_CONNECT_TYPE );

      for (i=0;i<n_sections;i++)  {
        loop->AddString ( monID,false );
        add_atom_to_cif ( loop, atom_tree    [i], atoms );
        add_atom_to_cif ( loop, atom_backward[i], atoms );
        add_atom_to_cif ( loop, atom_forward [i], atoms );
        switch (connect_type[i])  {
          default :
          case Tree::None  : loop->AddNoData ( mmdb::mmcif::CIF_NODATA_DOT );
                                    break;
          case Tree::Start : loop->AddString ( "START",false );
                                    break;
          case Tree::End   : loop->AddString ( "END",false );
                                    break;
          case Tree::Add   : loop->AddString ( "ADD",false );
        }
      }

    }

  }

  /*
  loop_
  _chem_comp_tree.comp_id
  _chem_comp_tree.atom_id
  _chem_comp_tree.atom_back
  _chem_comp_tree.atom_forward
  _chem_comp_tree.connect_type
   ATP      O2A    n/a    PA     START

  */


  void  Tree::write_mem ( PMemIO memIO, int version )  {
  UNUSED_ARGUMENT(version);
  int i;
    memIO->put_integer ( n_sections );
    for (i=0;i<n_sections;i++)  {
      memIO->put_integer ( atom_tree    [i] );
      memIO->put_integer ( atom_backward[i] );
      memIO->put_integer ( atom_forward [i] );
      memIO->put_integer ( connect_type [i] );
    }
  }

  void  Tree::copy ( PTree tree, mmdb::ivector anmatch )  {
  int i;
    empty();
    n_sections = tree->size();
    if (n_sections>0)  {
      mmdb::GetVectorMemory ( atom_tree,    n_sections,0 );
      mmdb::GetVectorMemory ( atom_backward,n_sections,0 );
      mmdb::GetVectorMemory ( atom_forward, n_sections,0 );
      mmdb::GetVectorMemory ( connect_type, n_sections,0 );
      if (anmatch)  {
        for (i=0;i<n_sections;i++)  {
          if (tree->atom(i)>=0)
            atom_tree    [i] = anmatch [ tree->atom    (i) ];
          if (tree->backward(i)>=0)
            atom_backward[i] = anmatch [ tree->backward(i) ];
          if (tree->forward(i)>=0)
            atom_forward [i] = anmatch [ tree->forward (i) ];
          connect_type [i] = tree->connect (i);
        }
      } else  {
        for (i=0;i<n_sections;i++)  {
          atom_tree    [i] = tree->atom    (i);
          atom_backward[i] = tree->backward(i);
          atom_forward [i] = tree->forward (i);
          connect_type [i] = tree->connect (i);
        }
      }
    }
  }


  bool  Tree::read_mem ( PMemIO memIO, int version, bool * Ok )  {
  UNUSED_ARGUMENT(version);
  int  i;
  bool success;
    empty();
    if (Ok)  success = *Ok;
       else  success = true;
    memIO->get_integer ( n_sections,&success );
    if (n_sections>0)  {
      mmdb::GetVectorMemory ( atom_tree,    n_sections,0 );
      mmdb::GetVectorMemory ( atom_backward,n_sections,0 );
      mmdb::GetVectorMemory ( atom_forward, n_sections,0 );
      mmdb::GetVectorMemory ( connect_type, n_sections,0 );
      for (i=0;i<n_sections;i++)  {
        memIO->get_integer ( atom_tree    [i],&success );
        memIO->get_integer ( atom_backward[i],&success );
        memIO->get_integer ( atom_forward [i],&success );
        memIO->get_integer ( connect_type [i],&success );
      }
    }
    if (Ok)  *Ok = success;
    return success;
  }

}  // namespace ccp4srs
