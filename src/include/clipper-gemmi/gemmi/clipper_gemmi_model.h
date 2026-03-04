/*! \file gemmi/clipper_gemmi_model.h
    \author Soon Wen Hoh
    Header file for clipper-gemmi model class type conversion
*/

#ifndef CLIPPER_GEMMI_MODEL
#define CLIPPER_GEMMI_MODEL

#include <gemmi/model.hpp>  // for Gemmi's model data structures
#include <gemmi/to_cif.hpp>
#include <gemmi/to_mmcif.hpp>
#include <gemmi/to_pdb.hpp>

#include "clipper_gemmi.h"

namespace clipper {

//! GEMMI Atom object wrapper
/*! This class is a trivial derivation of the corresponding GEMMI,
  providing access in terms of Clipper types. Thus, when you need
  such access, simply cast the GEMMI object reference to this type
  to access the additional functions. */
class GemmiAtom : public gemmi::Atom
{
 public:
  //! null constructor
  GemmiAtom() {}
  //! constructor: from GEMMI atom
  GemmiAtom( const gemmi::Atom& a ) : gemmi::Atom( a ) {}

  // standard atom properties
  String id() const;                               //!< Atom id, e.g. CA, CB, CH3
  String element() const;                          //!< Atom element, e.g. C, H, Zn2+
  Coord_orth coord_orth() const;                   //!< Atom coordinate (orthogonal Angstroms)
  ftype occupancy() const;                         //!< Atom occupancy (0...1)
  ftype u_iso() const;                             //!< Atom isotropic U
  U_aniso_orth u_aniso_orth() const;               //!< Atom anisotropic U (orthogonal As)
  void set_id( const String& n );                  //!< set id
  void set_element( const String& n );             //!< set element
  void set_coord_orth( const Coord_orth& v );      //!< set coordinate
  void set_occupancy( const ftype& v );            //!< set occupancy
  void set_u_iso( const ftype& v );                //!< set iso U
  void set_u_aniso_orth( const U_aniso_orth& v );  //!< set aniso U
  // other atom properties
  String altconf() const;               //!< get atom alternate conformation code
  void set_altconf( const String& n );  //!< set altconf
  int serial_num() const;               //!< get atom serial number
  String charge() const;                //!< get atom charge
};

//! GEMMI residue object wrapper
/*! This class is a trivial derivation of the corresponding GEMMI,
    providing access in terms of Clipper types. Thus, when you need
    such access, simply cast the GEMMI object reference to this type
    to access the additional functions. */
class GemmiResidue : public gemmi::Residue
{
 public:
  //! null constructor
  GemmiResidue() {}
  //! constructor: from GEMMI residue
  GemmiResidue( const gemmi::Residue& r ) : gemmi::Residue( r ) {}

  // residue properties
  String type() const;                  //!< Residue type
  int seqnum() const;                   //!< Sequence number
  String inscode() const;               //!< Insertion code
  void set_type( const String& n );     //!< Set residue type
  void set_seqnum( const int& n );      //!< Set sequence number
  void set_inscode( const String& n );  //!< Set insertion code
};

//! GEMMI chain object wrapper
/*! This class is a trivial derivation of the corresponding GEMMI,
  providing access in terms of Clipper types. Thus, when you need
  such access, simply cast the GEMMI object reference to this type
  to access the additional functions. */
class GemmiChain : public gemmi::Chain
{
 public:
  //! null constructor
  // GemmiChain() {}
  //! constructor: from GEMMI chain
  GemmiChain( const gemmi::Chain& c ) : gemmi::Chain( c ) {}

  // chain properties
  String id() const;               //!< Chain id
  void set_id( const String& n );  //!< Set chain id
};

//! GEMMI model object wrapper
/*! This class is a trivial derivation of the corresponding GEMMI,
  providing access in terms of Clipper types. Thus, when you need
  such access, simply cast the GEMMI object reference to this type
  to access the additional functions. */
class GemmiModel : public gemmi::Model {
 public:
  //! null constructor
  // GemmiModel() {}
  //! constructor: from GEMMI model
  GemmiModel( const gemmi::Model& m ) : gemmi::Model( m ) {}

  // model properties
  String id() const;               //!< Model num as String
  void set_id( const String& n );  //!< Set Model num from String
  void set_id( const int& n );     //!< Set Model num
};

//! Gemmi Structure object wrapper
/*! This class is a trivial derivation of the corresponding Gemmi,
  providing access in terms of Clipper types. Thus, when you need
  such access, simply cast the Gemmi object reference to this type
  to access additional functions. */
class GemmiStructure : public gemmi::Structure
{
 public:
  //! null constructor
  GemmiStructure() {}
  //! copy constructor: from GEMMI Structure
  GemmiStructure( const gemmi::Structure& s ) : gemmi::Structure( s ) {}

  //! get IDs for model, chain, residue or atom
  static String GetID_str( const String& model_num, const gemmi::CRA& cra, const String& entity );
  //! get IDs for model, chain, residue or atom
  static String GetID_str( const int& model_num, const gemmi::CRA& cra, const String& entity );
  //! get spacegroup and return in clipper::Spacegroup format
  Spacegroup spacegroup() const;
  //! get cell and return in clipper::Cell format
  Cell get_cell() const;
  //! set spacegroup to Gemmi Structure
  void set_spacegroup( const Spacegroup& spacegroup );
  //! set cell to Gemmi Structure
  void set_cell( const Cell& cell );
  //! print experimental method
  String get_exptlmethod() const;
  //! get transformation/origx in clipper::RTop_orth format: ORIGXn; _database_PDB_matrix.origx*
  RTop_orth get_transform() const;
  //! set transformation/origx to Gemmi Structure
  void set_transform( const RTop_orth& rtop );
  //! get resolution in clipper::Resolution format
  Resolution get_resolution() const;
  //! set resolution to Gemmi Structure
  void set_resolution( const Resolution& resol );
};

//! GEMMI atom list class
/*! This class is used to convert gemmi CraProxy to Clipper Atom_list
 It is a trivial derivation of clipper::Atom_list, and may be used
 wherever an Atom_list is required. */
class GemmiAtom_list : public Atom_list
{
 public:
  //! constructor: from CraProxy
  GemmiAtom_list( gemmi::CraProxy cra );
};

}  // namespace clipper

#endif
