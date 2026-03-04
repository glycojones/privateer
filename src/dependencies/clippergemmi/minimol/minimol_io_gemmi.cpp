/* minimol_io.cpp: atomic model io types GEMMI*/
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

#include "minimol_io_gemmi.h"

#include <gemmi/chemcomp.hpp>
#include <gemmi/resinfo.hpp>

extern "C" {
#include <string.h>
}

namespace clipper {

// GEMMIfile
// local types for referencing between GEMMI and MiniMol
struct GAtom {
  gemmi::Atom* db;
  const MAtom* mm;
};
struct GMono {
  gemmi::Residue* db;
  const MMonomer* mm;
  std::vector<GAtom> data;
};
struct GPoly {
  gemmi::Chain* db;
  const MPolymer* mm;
  std::vector<GMono> data;
};
struct GModl {
  gemmi::Model* db;
  const MModel* mm;
  std::vector<GPoly> data;
};

/*! This reads file of either PDB, mmCIF or mmjson format.
  Input file will be determined from extension. i.e. ".pdb", ".ent", ".cif", ".mmcif", ".json"
  For more details on PdbReadOptions, see gemmi's model.hpp
  \param file The filename (or pathname) of the file to read.
  \param pdb_read_opts PdbReadOptions struct
  \parblock
  max_line_length=0, split_chain_on_ter=false, skip_remarks=false
  \endparblock
  \param force_label If true, label_seq in Mmccif will be assigned even if full sequence is not
  known*/
void GEMMIfile::read_file( const String& file, gemmi::PdbReadOptions pdb_read_opts, bool force_label )
{
  gemmi::CoorFormat format = ::gemmi::coor_format_from_ext( file.trim().tail() );
  switch ( format ) {
    case gemmi::CoorFormat::Pdb: {
      structure_ = ::gemmi::read_pdb( ::gemmi::BasicInput( file ), pdb_read_opts );
      // setup entities for mmcif format
      setup_entities( structure_ );
      entities_set = true;
      // force_label=true will assign label_seq even if full sequence is not known (assumes no gaps)
      // refer GEMMI documentation
      assign_label_seq_id( structure_, force_label );
      break;
    }
    case gemmi::CoorFormat::Mmcif: {
      st_doc_.clear();
      structure_ =
          ::gemmi::make_structure( ::gemmi::cif::read( ::gemmi::BasicInput( file ) ), &st_doc_ );
      break;
    }
    case gemmi::CoorFormat::Mmjson: {
      st_doc_.clear();
      structure_ = ::gemmi::make_structure(
          ::gemmi::cif::read_mmjson( ::gemmi::BasicInput( file ) ), &st_doc_ );
      break;
    }
    case gemmi::CoorFormat::ChemComp: {
      st_doc_ = ::gemmi::cif::read( ::gemmi::BasicInput( file ) );
      structure_ = ::gemmi::make_structure_from_chemcomp_doc( st_doc_ );
    }
    case gemmi::CoorFormat::Unknown:
    case gemmi::CoorFormat::Detect: {
      String msg = "GEMMIfile: Error reading file, unknown format: " + file + "\n";
      Message::message( Message_fatal( msg ) );
    }
  }
}

/*! The default output file type will be the same as the file read, otherwise PDB.
  QuickWriteOptions.Minimal = true will always take precedence even when
  GemmiMmcifOutputGroups is set through GEMMIfile::set_mmcif_output_groups or
  PdbWriteOptions is set through GEMMIfile::set_pdb_write_options.
  \param file The filename (or pathname) of the file to be written.
  \param type 0=PDB, 1=CIF, 2=Mmjson, default=same as input file or PDB
  \param write_options
  \parblock
  QuickWriteOptions Minimal=false, PdbShortTER=false, PdbLinkR=false, PdbShortChainNames = false,
  CifAllAuth=false, UpdateCifDoc = false;
  \endparblock */
void GEMMIfile::write_file( const String& file, TYPE type, QuickWriteOptions write_options )
{
  // file types in Gemmi enum
  // { Unknown, Detect, Pdb, Mmcif, Mmjson, ChemComp }
  const TYPE types[4] = { PDB, CIF, Mmjson };
  auto intype = structure_.input_format;
  if ( type == Default ) {
    switch ( intype ) {
      case gemmi::CoorFormat::Mmcif:
        type = types[1];
        break;
      case gemmi::CoorFormat::Mmjson:
        type = types[2];
        break;
      case gemmi::CoorFormat::Pdb:
      default:
        type = types[0];
    }
  }
  gemmi::cif::Document doc = st_doc_;
  ::gemmi::Ofstream os( file, &std::cout );
  switch ( type ) {
    case CIF:  // write to mmcif/mmjson refer Gemmi's convert program
    case Mmjson: {
      if ( doc.blocks.size() == 0 ) {
        doc.clear();
        doc = ::gemmi::make_mmcif_document( structure_ );
      }
      if ( write_options.Minimal ) {
        doc.clear();
        doc.blocks.resize( 1 );
        ::gemmi::add_minimal_mmcif_data( structure_, doc.blocks[0] );
        gemmi::MmcifOutputGroups groups( false );
        groups.atoms = true;
        groups.group_pdb = true;
        ::gemmi::update_mmcif_block( structure_, doc.blocks[0], groups );
      } else {
        if ( !write_options.UpdateCifDoc ) {  // UpdateCifDoc = false so update only atoms
          // this will overwrite set_mmcif_output_groups
          gemmi::MmcifOutputGroups groups( false );
          groups.atoms = true;
          groups.group_pdb = true;
          groups.auth_all = write_options.CifAllAuth;
          ::gemmi::update_mmcif_block( structure_, doc.blocks[0], groups );
        } else {  // cif groups in document can be updated individually if set to true
          // i.e. groups.symmetry = true
          // use GEMMIfile::set_mmcif_output_groups to set groups
          // refer ::gemmi::MmcifOutputGroups or gemmi's to_mmcif.hpp
          ::gemmi::update_mmcif_block( structure_, doc.blocks[0], groups_ );
        }
      }
      if ( type == CIF ) {
        ::gemmi::cif::write_cif_to_stream( os.ref(), doc, cif_write_style_ );
      } else {  // mmjson
        ::gemmi::cif::write_json_to_stream( os.ref(), doc,
                                            ::gemmi::cif::JsonWriteOptions::mmjson() );
      }
    } break;
    case PDB:
    default: {
      ::gemmi::shorten_ccd_codes( structure_ );
      if ( write_options.Minimal ) pdb_write_opts_ = gemmi::PdbWriteOptions::minimal();
      if ( write_options.PdbShortTer ) pdb_write_opts_.numbered_ter = false;
      if ( write_options.PdbLinkR ) pdb_write_opts_.use_linkr = true;
      if ( write_options.PdbShortChainNames ) shorten_chain_names( structure_ );
      try {
        if ( pdb_write_opts_.ter_records && !entities_set ) ::gemmi::setup_entities( structure_ );
        ::gemmi::write_pdb( structure_, os.ref(), pdb_write_opts_ );
      } catch ( const std::exception& err ) {
        String msg = "GEMMIfile: write file error: " + file + "\n>> ";
        msg += err.what();
        Message::message( Message_fatal( msg ) );
      }
    } break;
  }
}

/*! Import data from the GEMMI hierarchy into a MiniMol object. All the
  atoms, residues and chains are imported unless a selection handle is
  supplied, in which case only the selected atoms and their containers
  are imported.

  Each data is tagged with an optional Property<String> (see
  clipper::PropertyManager), called "CID", which describes the
  position in the GEMMI hierarchy from which it was taken. This may
  (optionally) be used later to put a modified element back into the
  hierarchy.

  Any atoms which have alternate conformation codes all also be given
  an "AltConf" Property<String>.

  For more details, see clipper::PropertyManager::exists_property(),
  clipper::PropertyManager::get_property(). (e.g. label = "CID")
  \param minimol The minimol to import.
  \param model_num (optional) Model selection handle when number of models > 1 in file, starts
  with 1. */
void GEMMIfile::import_minimol( MiniMol& minimol, const int model_num )
{
  // clear the model
  minimol = MiniMol( structure_.spacegroup(), structure_.get_cell() );
  // import objects
  String txt;
  MModel& mol = minimol.model();
  // get name of first model or model number from sel_hnd
  gemmi::Model* p_mod = nullptr;
  if ( model_num == 0 )  // just in case
    p_mod = structure_.find_model( 1 );
  else
    p_mod = structure_.find_model( model_num );
  mol.set_property( "StrucName", Property<String>( structure_.name ) );
  gemmi::CRA cra;

  if ( p_mod != NULL ) {
    for ( int c = 0; c < p_mod->chains.size(); c++ ) {
      gemmi::Chain* p_chn = &p_mod->chains.at( c );
      if ( p_chn != NULL ) {
        // import the chain
        cra.chain = p_chn;
        MPolymer pol;
        for ( gemmi::Residue res : p_chn->whole() ) {
          // import residue
          cra.residue = &res;
          MMonomer mon;
          for ( gemmi::Atom a : res.children() ) {
            // import atoms
            cra.atom = &a;
            MAtom atm( Atom::null() );
            String altLoc = a.has_altloc() ? String( &a.altloc, 1 ) : "";
            atm.set_name( a.padded_name(), altLoc, true );
            atm.set_element( a.element.uname() );
            atm.set_coord_orth( Coord_orth( a.pos.at( 0 ), a.pos.at( 1 ), a.pos.at( 2 ) ) );
            atm.set_occupancy( a.occ );
            atm.set_u_iso( Util::b2u( a.b_iso ) );
            if ( a.aniso.u11 != 0 )
              atm.set_u_aniso_orth( U_aniso_orth( a.aniso.u11, a.aniso.u22, a.aniso.u33,
                                                  a.aniso.u12, a.aniso.u13, a.aniso.u23 ) );
            txt = structure_.GetID_str( p_mod->num, cra, a.what() );  // Get atom ID
            atm.set_property( "CID", Property<String>( txt ) );
            if ( a.has_altloc() ) atm.set_property( "AltConf", Property<String>( altLoc ) );
            if ( a.charge ) atm.set_property( "Charge", Property<signed char>( a.charge ) );
            mon.insert( atm );  // store the atom
          }
          cra.atom = nullptr;
          String inscode = res.seqid.has_icode() ? String( &res.seqid.icode, 1 ) : "";
          mon.set_seqnum( res.seqid.num.value, inscode );  // String(&res.seqid.icode, 1));
          mon.set_type( res.name );
          txt = structure_.GetID_str( p_mod->num, cra, res.what() );  // Get residue ID
          mon.set_property( "CID", Property<String>( String( txt ) ) );
          pol.insert( mon );  // store the residue
        }
        pol.set_id( p_chn->name );
        cra.residue = nullptr;
        txt = structure_.GetID_str( p_mod->num, cra, p_chn->what() );  // Get chain ID
        pol.set_property( "CID", Property<String>( String( txt ) ) );
        mol.insert( pol );  // store the chain
      }
    }
    cra.chain = nullptr;
    txt = structure_.GetID_str( p_mod->num, cra, p_mod->what() );  // Get model ID
    mol.set_property( "CID", Property<String>( String( txt ) ) );
  }
}

/*! Export data to the GEMMI hierachy from a MiniMol object. All the
  atoms, residues and chains are exported.

  If any MiniMol object has a "CID" Property<String> (see
  clipper::PropertyManager), then the information from that object
  will be used to update the corresponding object in the MMDB
  hierarchy, if it exists. If there is no such entry in the MMDB
  hierarchy, or if no "CID" Property<String> exists, then a new object
  will be created in the MMDB hierarchy.
  \param minimol The MiniMol to be exported into.*/
void GEMMIfile::export_minimol( MiniMol& minimol )
{
  // export spacegroup/cell
  if ( !minimol.spacegroup().is_null() ) structure_.set_spacegroup( minimol.spacegroup() );
  if ( !minimol.cell().is_null() ) structure_.set_cell( minimol.cell() );

  // create structure for relationships between MiniMol and GEMMI
  GModl gmod;
  clipper::String cid;

  // fill structure
  gmod.mm = &( minimol.model() );
  gmod.db = NULL;
  gmod.data.resize( gmod.mm->size() );
  for ( int p = 0; p < gmod.data.size(); p++ ) {  // loop over chains
    GPoly& gpol = gmod.data[p];
    gpol.mm = &( ( *gmod.mm )[p] );
    gpol.db = NULL;
    gpol.data.resize( gpol.mm->size() );
    for ( int r = 0; r < gpol.data.size(); r++ ) {  // loop over residues
      GMono& gmon = gpol.data[r];
      gmon.mm = &( ( *gpol.mm )[r] );
      gmon.data.resize( gmon.mm->size() );
      for ( int a = 0; a < gmon.data.size(); a++ ) {  // loop over atoms
        GAtom& gatm = gmon.data[a];
        gatm.mm = &( ( *gmon.mm )[a] );
        gatm.db = NULL;
      }
    }
  }

  // make the GEMMI reference by CID if present
  if ( gmod.mm->exists_property( "CID" ) ) {
    cid = dynamic_cast<const Property<String>&>( gmod.mm->get_property( "CID" ) )
              .value();  // + "/A/0/A";
    std::pair<gemmi::Model*, gemmi::CRA> PM_CRA =
        ::gemmi::Selection( cid.trim() ).first( structure_ );
    if ( PM_CRA.first != nullptr ) gmod.db = PM_CRA.first;

    if ( gmod.db != NULL ) {
      for ( int p = 0; p < gmod.data.size(); p++ ) {
        GPoly& gpol = gmod.data[p];
        if ( gpol.mm->exists_property( "CID" ) ) {
          cid = dynamic_cast<const Property<String>&>( gpol.mm->get_property( "CID" ) )
                    .value();  // + "/0/A";
          PM_CRA = ::gemmi::Selection( cid.trim() ).first( structure_ );
          gpol.db = PM_CRA.second.chain;
        }
        if ( gpol.db != NULL ) {
          for ( int r = 0; r < gpol.data.size(); r++ ) {
            GMono& gmon = gpol.data[r];
            if ( gmon.mm->exists_property( "CID" ) ) {
              cid = dynamic_cast<const Property<String>&>( gmon.mm->get_property( "CID" ) )
                        .value();  // + "/A";
              PM_CRA = ::gemmi::Selection( cid.trim() ).first( structure_ );
              gmon.db = PM_CRA.second.residue;
            }
            if ( gmon.db != NULL ) {
              for ( int a = 0; a < gmon.data.size(); a++ ) {
                GAtom& gatm = gmon.data[a];
                if ( gatm.mm->exists_property( "CID" ) ) {
                  cid = dynamic_cast<const Property<String>&>( gatm.mm->get_property( "CID" ) )
                            .value();
                  PM_CRA = ::gemmi::Selection( cid.trim() ).first( structure_ );
                  gatm.db = PM_CRA.second.atom;
                }
              }
            }
          }
        }
      }
    }
  }
  // create GEMMI object for anything that is missing
  // fill information in Gemmi from MiniMol
  if ( gmod.db == NULL ) {
    structure_.models.emplace_back( 1 );
    gmod.db = &structure_.models.back();
  }
  for ( int p = 0; p < gmod.data.size(); p++ ) {  // chains
    GPoly& gpol = gmod.data[p];
    if ( gpol.db == NULL ) {  // nullptr
      String chn_name = gpol.mm->id().trim();
      gmod.db->chains.emplace_back( chn_name );
      gpol.db = &gmod.db->chains.back();
    } else
      gpol.db->name = gpol.mm->id().trim();
    for ( int r = 0; r < gpol.data.size(); r++ ) {  // residues
      GMono& gmon = gpol.data[r];
      if ( gmon.db == NULL )  // nullptr
      {
        gemmi::ResidueId rid;
        gpol.db->residues.emplace_back( rid );
        gmon.db = &gpol.db->residues.back();
      }
      gmon.db->seqid.num = gmon.mm->seqnum();
      gmon.db->name = gmon.mm->type();  //.substr(0, 19);
      int pos = gmon.mm->id().find( ":" );
      if ( pos != String::npos ) {
        char inscode = gmon.mm->id()[pos + 1];
        if ( inscode != '\r' && inscode != '\n' ) gmon.db->seqid.icode = inscode;
      }
      if ( !gmon.db->het_flag ) {
        bool std_res = ::gemmi::find_tabulated_residue( gmon.db->name ).is_standard();
        gmon.db->het_flag = std_res ? 'A' : 'H';
      }
      for ( int a = 0; a < gmon.data.size(); a++ ) {  // loop through atoms
        GAtom& gatm = gmon.data[a];
        if ( gatm.db == NULL )  // nullptr
        {
          gemmi::Atom atm1;
          gmon.db->atoms.emplace_back( atm1 );
          gatm.db = &gmon.db->atoms.back();
        }

        if ( !gatm.mm->coord_orth().is_null() ) {  // set atom coord
          gatm.db->pos.x = gatm.mm->coord_orth().x();
          gatm.db->pos.y = gatm.mm->coord_orth().y();
          gatm.db->pos.z = gatm.mm->coord_orth().z();
        }
        if ( !Util::is_nan( gatm.mm->occupancy() ) )  // set atom occupancy
          gatm.db->occ = gatm.mm->occupancy();
        if ( !Util::is_nan( gatm.mm->u_iso() ) )  // set atom u_iso
          gatm.db->b_iso = Util::u2b( gatm.mm->u_iso() );
        if ( !gatm.mm->u_aniso_orth().is_null() ) {  // set atom u_aniso
          gatm.db->aniso.u11 = gatm.mm->u_aniso_orth().mat00();
          gatm.db->aniso.u22 = gatm.mm->u_aniso_orth().mat11();
          gatm.db->aniso.u33 = gatm.mm->u_aniso_orth().mat22();
          gatm.db->aniso.u12 = gatm.mm->u_aniso_orth().mat01();
          gatm.db->aniso.u13 = gatm.mm->u_aniso_orth().mat02();
          gatm.db->aniso.u23 = gatm.mm->u_aniso_orth().mat12();
        }
        if ( gatm.mm->name() != "" )  // set atom id/name for Gemmi
          gatm.db->name = gatm.mm->name().trim();
        if ( gatm.mm->element() != "" )  // set atom element
          gatm.db->element = ::gemmi::Element( gatm.mm->element().c_str() );
        if ( gatm.mm->exists_property( "Charge" ) ) {  // set atom charge if any
          signed char a =
              dynamic_cast<const Property<signed char>&>( gatm.mm->get_property( "Charge" ) )
                  .value();
          gatm.db->charge = a;
        }
        if ( gatm.mm->exists_property( "AltConf" ) ) {  // set alt conf code
          String a =
              dynamic_cast<const Property<String>&>( gatm.mm->get_property( "AltConf" ) ).value();
          if ( a != "" ) gatm.db->altloc = ( a[0] == ' ' ) ? '\0' : a[0];
        }
      }
    }
  }
}

/*! Set spacegroup to structure held in GEMMIfile object.
  \param sg The spacegroup to set */
void GEMMIfile::set_spacegroup( const Spacegroup& sg ) { structure_.set_spacegroup( sg ); }

/*! Set unit cell to structure held in GEMMIfile object.
  \param cell The cell to set */
void GEMMIfile::set_cell( const Cell& cell_in ) { structure_.set_cell( cell_in ); }

}  // namespace clipper
