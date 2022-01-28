//  $Id: ccp4srs_base.h $
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
//  **** Module  :  ccp4srs_base  <interface>
//       ~~~~~~~~~
//  **** Classes :  ccp4srs::Base - SRS base manager class
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2010-2014
//
//  =================================================================
//

#ifndef CCP4SRS_BASE_H
#define CCP4SRS_BASE_H

#include "ccp4srs_index.h"
#include "mmdb2/mmdb_math_graph.h"

namespace ccp4srs  {

  // ==================================================================

  DefineClass(Base);

  class Base  {

    public:

      Base ();  //!< Constructor
      ~Base();  //!< Destructor

      //   LoadIndex() loads index of the structural database. 'srsPath'
      // must point on the directory containing the database files.
      // The index must be loaded once before retrieving any
      // information from the database.
      //   LoadIndex() may return either CCP4SRS_Ok, CCP4SRS_CorruptIndex
      // or CCP4SRS_FileNotFound.
      int loadIndex ( mmdb::cpstr srsPath );

      //   LoadStructure(..) reads structure from *.srs files and
      // stores it in RAM for faster access. There is no special
      // functions to access loaded structures, all requests to
      // *.srs files and RAM-storage are dispatched automatically.
      int  loadStructure   ( mmdb::cpstr entryID );

      PIndex getIndex ( mmdb::cpstr entryID );
      PIndex getIndex ( int n         );

      //   UnloadStructure(..) deletes strtucture from RAM and releases
      // its memory. The structure is then accessible in the normal
      // way from *.srs files, which is slower.
      int  unloadStructure ( mmdb::cpstr entryID );

      //   GetPath() returns full path to a file with file name FName
      // in the database directory.  Length of S should suffice for
      // accomodating the path. The function returns S.
      //   GetPath() will work only after loading the database index.
      mmdb::pstr getPath ( mmdb::pstr & S, mmdb::cpstr FName );

      //   GetStructFile() creates and opens the database structure
      // file and returns its pointer. In the case of errors returns
      // NULL. Application is responsible for deleting this class.
      mmdb::io::PFile getStructFile();

      //   GetGraphFile() creates and opens the database graph
      // file and returns its pointer. In the case of errors returns
      // NULL. Application is responsible for deleting this class.
      mmdb::io::PFile getGraphFile();

      //   n_entries() returns number of entries in the database index.
      inline int  n_entries() { return nEntries; }

      //   getEntryNo() returns position of the structure with
      // (3-letter) name 'name' as found in the database index.
      //   Non-negative return means success, otherwise
      // CCP4SRS_EntryNotFound indicates that the requested structure
      // was not found in the database index.
      int  getEntryNo ( mmdb::cpstr entryID );

      //   getMonomer returns pointer to the monomer structure
      // identified by 3-letter entryID. If such structure is not
      // found, the function returns NULL.
      //   The function returns a pointer to a private copy of the
      // structure. Modifying it will not change data in the structural
      // database. The application is responsible for deallocating
      // the structure after use (simply use delete).
      //   See description of CCP4SRSMonomer for the explanation of
      // its fields.
      PMonomer getMonomer ( mmdb::cpstr entryID );

      //   Another form of getMonomer(..) uses an open structure
      // file, which allows to save on opening/closing file if
      // multiple access to SRS structures is required.
      PMonomer getMonomer ( int entryNo,  // 0...nEntries-1
                            mmdb::io::PFile structFile );
      PMonomer getMonomer ( mmdb::cpstr entryID,
                            mmdb::io::PFile structFile );


      /// \brief Returns the number of atoms in the references monomer
      /// \param[inp] entryID 3-letter cif code
      /// \return Non-negative return delivers the number of atoms,
      /// otherwise CCP4SRS_EntryNotFound indicates that the requested
      /// structure was not found in the database
      int getNofAtoms ( mmdb::cpstr entryID );

      /// \brief Returns the number of atoms in the references monomer
      /// \param[in] entryNo entry number (0...n_entries()-1)
      /// \return Non-negative return delivers the number of atoms,
      /// otherwise CCP4SRS_EntryNotFound indicates that the requested
      /// structure was not found in the database
      int getNofAtoms ( int entryNo );

      /// \brief Returns number of non-hydrogen and hydrogen atoms
      ///        for referenced structure.
      /// \param[in] entryNo entry number (0...n_entries()-1)
      /// \param[out] nNonHAtoms number of non-hydrogen atoms
      /// \param[out] nHAtoms number of hydrogen atoms
      /// \return CCP4SRS_Ok or CCP4SRS_EntryNotFound.
      int  getNofAtoms ( int entryNo, int & nNonHAtoms, int & nHAtoms );

      /// \brief Returns atom names, total number of atoms and number of
      /// hydrogens (nH) for referenced monomer.
      /// \param[in] entryNo entry number (0...n_entries()-1)
      /// \param[out] AtName array of atom names, prealocated with
      /// sufficient length
      /// \param[out] nAtoms total number of atoms
      /// \param[out] nH number of hydrogens
      /// \return CCP4SRS_Ok, CCP4SRS_FileNotFound, CCP4SRS_EntryNotFound
      /// and CCP4SRS_ReadErrors.
      int  getAtNames ( int  entryNo, mmdb::PAtomName AtName,
                        int & nAtoms, int & nH );

      /// \brief Returns atom names, total number of atoms and number of
      /// hydrogens (nH) for referenced monomer.
      /// \param[in] structFile pre-open SRS structure file
      /// \param[in] entryNo entry number (0...n_entries()-1)
      /// \param[out] AtName array of atom names, prealocated with
      /// sufficient length
      /// \param[out] nAtoms total number of atoms
      /// \param[out] nH number of hydrogens
      /// \return CCP4SRS_Ok, CCP4SRS_FileNotFound, CCP4SRS_EntryNotFound
      /// and CCP4SRS_ReadErrors.
      int  getAtNames  ( mmdb::io::PFile structFile, int entryNo,
                         mmdb::PAtomName  AtName, int & nAtoms, int & nH  );

      /// \brief Retrieves various atom properties for referenced
      ///        structure.
      /// \param[in] entryID monomer's cif id
      /// \param[out] nNonHAtoms number of non-hydrogen atoms
      /// \param[out] nNonHAtomName non-hydrogen atom names
      /// \param[out] nHAtoms number of hydrogen atoms
      /// \param[out] nNonHAtomName hydrogen atom names
      /// \param[out] Hconnect connectivity to non-hydrogen atoms:
      ///             hydrogen HAtName[i] is connected to non-hydrogen
      ///             NonHAtom[Hconnect[i]], if Hconnect[i]>=0.
      /// \param[out] Elem element IDs
      /// \param[out] Chiral chirality IDs
      /// \return CCP4SRS_Ok, CCP4SRS_FileNotFound, CCP4SRS_EntryNotFound,
      /// CCP4SRS_ReadErrors and CCP4SRS_BrokenBonds.
      int  getAtoms ( mmdb::cpstr entryID,
                      int & nNonHAtoms, mmdb::PAtomName NonHAtName,
                      int & nHAtoms,    mmdb::PAtomName HAtName,
                      mmdb::ivector Hconnect, mmdb::ivector Elem,
                      mmdb::ivector Chiral );

      /// \brief Retrieves bonds for referenced structure.
      /// \param[in] entryID monomer's cif id
      /// \param[out] nBonds   number of bonds for ith atom
      /// \param[out] bondPair bonded atoms: bondPair[i][j] (0<=i<nAtoms,
      ///             0<=j<=nBonds[i]) gives the number of atom connected
      ///             to ith atom; only pairs i<j are returned
      /// \param[in] maxNAtoms length of nBonds[] and bondPairs[]
      /// \param[in] maxNBonds length of bondPairs[][]
      /// \return CCP4SRS_Ok, CCP4SRS_FileNotFound, CCP4SRS_EntryNotFound
      /// and CCP4SRS_ReadErrors.
      int  getBonds ( mmdb::cpstr   entryID,
                      mmdb::ivector nBonds, mmdb::imatrix bondPair,
                      int &   nAtoms, int maxNAtoms,
                      int     maxNBonds );

      int  getHetInfo ( mmdb::cpstr       entryID,
                        mmdb::pstr        Formula,
                        mmdb::pstr        Hname,
                        mmdb::pstr        Hsynonym,
                        mmdb::pstr        Hcharge,
                        mmdb::PAtomName & ClinkAtom,  // will
                        mmdb::PElement  & ClinkEle,   //   be
                        mmdb::PAtomName & SlinkAtom,  //     allocated
                        mmdb::PElement  & SlinkEle,   //       or NULL
                        int             & nLeavingAtoms );

      //   getGraph(..) retrieves data for chemical structure number
      // entryNo (as in Index) from graph file graphFile, then
      // allocates and builds the corresponding graph, which is
      // returned in G.
      //   If Hflag is set >= 1, all hydrogens are removed from
      // the graph. If Hflag is set to 2, element types of atoms,
      // to which hydrogens are bonded, are modified with flag
      // HYDROGEN_BOND.
      //   Returns CCP4SRS_Ok in case of success. Other return code are
      // CCP4SRS_WrongIndex and CCP4SRS_ReadError.
      int  getGraph ( mmdb::io::PFile graphFile, mmdb::math::RPGraph G );
      int  getGraph ( mmdb::io::PFile graphFile, int entryNo,
                      mmdb::math::RPGraph G, int Hflag );
      int  getGraph ( mmdb::io::PFile graphFile, mmdb::math::RPGraph G,
                      int Hflag );
      int  getGraph ( int    entryNo, mmdb::math::RPGraph G, int Hflag );
      int  getGraph ( mmdb::cpstr entryID, mmdb::math::RPGraph G,
                      int Hflag );


    protected:
      MemIO      memIO;
      PPIndex    index;
      mmdb::pstr srsDir;
      int        srsVersion;
      int        nEntries;

      void init_base ();
      void freeMemory();

    private:
      int index_alloc;

  };

} // namespace ccp4srs


#endif // CCP4SRS_BASE_H
