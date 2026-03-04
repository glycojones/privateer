/*! \file minimol_io_gemmi.h
  Header file for atomic model io types GEMMI*/

// C Copyright (C) 2000-2006 Kevin Cowtan and University of York
// L
// L  This library is free software and is distributed under the terms
// L  and conditions of version 2.1 of the GNU Lesser General Public
// L  Licence (LGPL) with the following additional clause:
// L
// L     `You may also combine or link a "work that uses the Library" to
// L     produce a work containing portions of the Library, and distribute
// L     that work under terms of your choice, provided that you give
// L     prominent notice with each copy of the work that the specified
// L     version of the Library is used in it, and that you include or
// L     provide public access to the complete corresponding
// L     machine-readable source code for the Library including whatever
// L     changes were used in the work. (i.e. If you make changes to the
// L     Library you must distribute those, but you do not need to
// L     distribute source or object code to those portions of the work
// L     not covered by this licence.)'
// L
// L  Note that this clause grants an additional right and does not impose
// L  any additional restriction, and so does not affect compatibility
// L  with the GNU General Public Licence (GPL). If you wish to negotiate
// L  other terms, please contact the maintainer.
// L
// L  You can redistribute it and/or modify the library under the terms of
// L  the GNU Lesser General Public License as published by the Free Software
// L  Foundation; either version 2.1 of the License, or (at your option) any
// L  later version.
// L
// L  This library is distributed in the hope that it will be useful, but
// L  WITHOUT ANY WARRANTY; without even the implied warranty of
// L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// L  Lesser General Public License for more details.
// L
// L  You should have received a copy of the CCP4 licence and/or GNU
// L  Lesser General Public License along with this library; if not, write
// L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
// L  The GNU Lesser General Public can also be obtained by writing to the
// L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// L  MA 02111-1307 USA

#ifndef CLIPPER_MINIMOL_IO_GEMMI
#define CLIPPER_MINIMOL_IO_GEMMI

#include <gemmi/mmread.hpp>  // coor_format_from_ext, read_pdb,

#include "../gemmi/clipper_gemmi_model.h"
#include "minimol.h"
// make_structure, read_mmjson, read
#include <gemmi/align.hpp>     // assign_label_seq_id
#include <gemmi/assembly.hpp>  // shorten_chain_name
#include <gemmi/fstream.hpp>   // Ofstream file
#include <gemmi/select.hpp>    // Selection
#include <gemmi/to_json.hpp>   // JsonWriter

namespace clipper {

//! Write options for writing PDB/CIF file
/*! A few simplified quick options for writing PDB/CIF files.
    Options such as Minimal and AllAuth will always overwrite
    what is set in GEMMIFile::set_mmcif_output_groups
    if they contradicts. */
struct QuickWriteOptions {
  bool Minimal = false;             //!< disable many records (HEADER, TITLE, ...) for PDB and mmcif
  bool PdbShortTer = false;         //!< write PDB TER records without numbers
  bool PdbLinkR = false;            //!< use non-standard Refmac LINKR record instead of LINK
  bool PdbShortChainNames = false;  //!< shorten chain names for chain name size > 2
  bool CifAllAuth = false;          //!< include _atom_site.auth_atom_id and auth_comp_id
  bool UpdateCifDoc = false;        //!< update existing mmcif blocks from read in document
};

//! GEMMI file object for MiniMol i/o
/*! This object is an i/o object for MiniMol, representing an
  interface between MiniMol and PDB or CIF files. The model
  is read using GEMMI implemented methods and then imported to
  MiniMol object and vice versa. Works fine for gemmi>=v0.6.4.*/
class GEMMIfile
{
 public:
  enum TYPE { Default = -1, PDB, CIF, Mmjson };

  //! null constructor
  GEMMIfile() {}
  //! read structure from file
  void read_file( const String& file,
                  gemmi::PdbReadOptions pdb_read_options = gemmi::PdbReadOptions(),
                  bool force_label = false );
  //! write structure to file
  void write_file( const String& file, TYPE type = Default,
                   QuickWriteOptions write_options = QuickWriteOptions() );
  //! import MiniMol from Gemmi Structure.
  void import_minimol( MiniMol& minimol, const int model_num = 1 );
  //! export MiniMol to Gemmi Structure.
  void export_minimol( MiniMol& minimol );
  //! set structure if not read from file
  void set_gemmi_structure( const GemmiStructure& st ) { structure_ = st; };
  //! get structure held in GEMMIFile
  GemmiStructure get_gemmi_structure() const { return structure_; };
  //! set mmcif output groups, default: GemmiMmcifOutputGroups groups(true).
  /*! For more details, see MmcifOutputGroups in gemmi's to_mmcif.hpp. */
  void set_mmcif_output_groups( gemmi::MmcifOutputGroups groups ) { groups_ = groups; };
  //! set pdb write options
  /*! For more details, see PdbWriteOptions in gemmi's to_pdb.hpp */
  void set_pdb_write_options( gemmi::PdbWriteOptions opts ) { pdb_write_opts_ = opts; };
  //! set spacegroup to structure
  void set_spacegroup( const Spacegroup& sg );
  //! set cell to structure
  void set_cell( const Cell& cell_in );
  //! get spacegroup in clipper::Spacegroup format
  Spacegroup spacegroup() const { return structure_.spacegroup(); };
  //! get cell in clipper::Cell format
  Cell cell() const { return structure_.get_cell(); };
  //! set cif output style.
  /*! For more details, see WriteOptions in gemmi's to_cif.hpp */
  void set_cif_style( gemmi::cif::WriteOptions cif_write_style ) {
    cif_write_style_ = cif_write_style;
  };

 private:
  GemmiStructure structure_;
  gemmi::cif::Document st_doc_;
  bool entities_set = false;
  gemmi::cif::WriteOptions cif_write_style_ = gemmi::cif::WriteOptions();
  gemmi::MmcifOutputGroups groups_ = gemmi::MmcifOutputGroups( true );
  gemmi::PdbWriteOptions pdb_write_opts_ = gemmi::PdbWriteOptions();
};

}  // namespace clipper

#endif
