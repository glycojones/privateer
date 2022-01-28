//  $Id: ccp4srs_tree.h $
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
//  **** Module  :  ccp4srs_tree  <interface>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Tree  - structure tree
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2013
//
//  =================================================================
//

#ifndef CCP4SRS_TREE_H
#define CCP4SRS_TREE_H

#include "ccp4srs_atom.h"

namespace ccp4srs  {

  DefineClass(Tree);

  class Tree  {

    // atom connectivity types
    enum  {None,Start,End,Add};

    public:
      Tree();
      virtual ~Tree();

      void  empty();

      inline int size() { return n_sections;    }
      inline int atom     ( int secNo ) { return atom_tree    [secNo]; }
      inline int backward ( int secNo ) { return atom_backward[secNo]; }
      inline int forward  ( int secNo ) { return atom_forward [secNo]; }
      inline int connect  ( int secNo ) { return connect_type [secNo]; }

      int   readFromCIF   ( mmdb::mmcif::PLoop mmCIFLoop,
                            Container<Atom> & atoms );
      void  writeToCIF    ( mmdb::mmcif::PData mmCIFData,
                            mmdb::cpstr monID,
                            Container<Atom> & atoms );

      void  write_mem   ( PMemIO memIO, int version );
      bool  read_mem    ( PMemIO memIO, int version,
                          bool * Ok = NULL );

      void  copy ( PTree tree, mmdb::ivector anmatch=NULL );

    protected:
      int            n_sections;     //!< number of sections
      mmdb::ivector  atom_tree;      //!< tree atom name
      mmdb::ivector  atom_backward;  //!< previous atom name
      mmdb::ivector  atom_forward;   //!< next atom name
      mmdb::ivector  connect_type;   //!< type of connection

  };

}  // namespace ccp4srs

#endif // CCP4SRS_TREE_H
