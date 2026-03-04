/* clipper_gemmi_model.cpp: clipper-gemmi model class type conversion */

#include "clipper_gemmi_model.h"

extern "C" {
#include <string.h>
}

namespace clipper {
// GEMMI wrapper types
// Atom
/*! \return Atom id as String*/
String GemmiAtom::id() const
{ 
  return String( name );
}

/*! \return Element name as String*/
String GemmiAtom::element() const
{ 
  return String( gemmi::Atom::element.name() );
}

/*! \return Orthogonal coordinates in clipper Coord_orth.
  Returns null values if NAN.*/
Coord_orth GemmiAtom::coord_orth() const
{
  if ( !pos.has_nan() )
    return Coord_orth( pos.x, pos.y, pos.z );
  else
    return Coord_orth( Coord_orth::null() );
}

ftype GemmiAtom::occupancy() const
{ 
  return occ;
}

ftype GemmiAtom::u_iso() const
{ 
  return Util::b2u( b_iso );
}

/*! \return Anisotropic U values in clipper U_aniso_orth format.
  Returns null values if u11+u22+u33=0. */
U_aniso_orth GemmiAtom::u_aniso_orth() const
{
  if ( this->aniso.nonzero() )
    return U_aniso_orth( aniso.u11, aniso.u22, aniso.u33, aniso.u12, aniso.u13, aniso.u23 );
  else
    return U_aniso_orth( U_aniso_orth::null() );
}

void GemmiAtom::set_id( const String& n )
{ 
  name = ( char* )n.c_str();
}

void GemmiAtom::set_element( const String& n )
{
  gemmi::Atom::element = ::gemmi::Element( ( char* )n.c_str() );
}

void GemmiAtom::set_coord_orth( const Coord_orth& v )
{
  if ( !v.is_null() ) {
    pos.x = v.x();
    pos.y = v.y();
    pos.z = v.z();
  }
}

void GemmiAtom::set_occupancy( const ftype& v )
{
  if ( !Util::is_nan( v ) ) {
    occ = float( v );
  }
}

void GemmiAtom::set_u_iso( const ftype& v )
{
  if ( !Util::is_nan( v ) ) {
    b_iso = Util::u2b( v );
  }
}

void GemmiAtom::set_u_aniso_orth( const U_aniso_orth& v )
{
  if ( !v.is_null() ) {
    aniso.u11 = v( 0, 0 );
    aniso.u22 = v( 1, 1 );
    aniso.u33 = v( 2, 2 );
    aniso.u12 = v( 0, 1 );
    aniso.u13 = v( 0, 2 );
    aniso.u23 = v( 1, 2 );
  }
}

/*! \return The atom alternate conformation code as String. */
String GemmiAtom::altconf() const
{ 
  return has_altloc() ? String( &altloc, 1 ) : "";
}

void GemmiAtom::set_altconf( const String& n )
{ 
  altloc = ( n[0] == ' ' ) ? '\0' : n[0];
}

/*! \return The atom serial number. */
int GemmiAtom::serial_num() const
{ 
  return serial;
}

/*! \return The atomic charge as clipper String with sign; sign+digit. */
String GemmiAtom::charge() const
{
  signed char c = gemmi::Atom::charge;
  char digit = c ? c > 0 ? '0' + c : '0' - c : ' ';
  char sign = c ? c > 0 ? '+' : '-' : ' ';
  return String( &sign, 1 ) + String( &digit, 1 );
}

// Residue
String GemmiResidue::type() const
{ 
  return String( name );
}

int GemmiResidue::seqnum() const
{ 
  return int( seqid.num );
}

String GemmiResidue::inscode() const
{ 
  return seqid.has_icode() ? String( &seqid.icode, 1 ) : "";
}

void GemmiResidue::set_type( const String& n )
{ 
  name = ( char* )n.c_str();
}

void GemmiResidue::set_seqnum( const int& n )
{ 
  seqid.num = n;
}

void GemmiResidue::set_inscode( const String& n )
{
  if ( n[0] != '\r' && n[0] != '\n' ) seqid.icode = n[0];
}

// Chain
String GemmiChain::id() const
{ 
  return String( name );
}

void GemmiChain::set_id( const String& n )
{ 
  name = ( char* )n.c_str();
}

// Model
String GemmiModel::id() const
{ 
  return String( num, 0 );
}

void GemmiModel::set_id( const String& n )
{
  num = n.i();  // clipper String convert to int
}

void GemmiModel::set_id( const int& n )
{ 
  num = n;
}

// Structure
/*! Return a String ID for atom, residue, chain or model.
  similar to gemmi's - model: /1, chn: /1/C, res: /1/C/LYS.ins,
  atom: /1/C/LYS.ins/NZ[N]:altloc

  Currently used in GEMMIFile::import_minimol to set "CID" labelled property
  \param model_num Model number
  \param cra chain, residue, atom - CRA struct object
  \param entity entity - Model, Chain, Residue, or Atom */
String GemmiStructure::GetID_str( const String& model_num, const gemmi::CRA& cra, const String& entity )
{
  String r = "/";  // clipper String start with "/"
  if ( !model_num.empty() ) {
    r += model_num;
  }
  if ( entity == "Model" ) return r;  // return just /modelID
  r += '/';
  if ( !cra.chain->name.empty() ) {
    r += cra.chain->name;
  }
  if ( entity == "Chain" ) return r;  // return /modelID/chnID
  r += '/';
  r += String( int( cra.residue->seqid.num ), 0 );
  if ( !cra.residue->seqid.has_icode() ) {
    r += '(';
    r += cra.residue->name;
    r += ')';
  } else  // /seqnum(restype).Ins causes syntax error in GEMMI selection
  {
    r += '.';
    r += cra.residue->seqid.icode;
  }
  if ( entity == "Residue" ) return r;  // return /modelID/chnID/seqnum.Ins
  r += '/';
  r += cra.atom->name;
  r += '[';
  r += cra.atom->element.uname();
  r += ']';
  if ( cra.atom->has_altloc() ) {
    r += ':';
    r += cra.atom->altloc;
  }
  return r;  // return /modelID/chnID/seqnum.Ins/atomName[elem]:altloc
}

/*! Return a String ID for atom, residue, chain or model.
  similar to gemmi's - model: /1, chn: /1/C, res: /1/C/LYS.ins,
  atom: /1/C/LYS.ins/NZ[N]:altloc

  Currently used in GEMMIFile::import_minimol to set "CID" labelled property
  \param model_num Model number
  \param cra chain, residue, atom - CRA struct object
  \param entity entity - Model, Chain, Residue, or Atom */
String GemmiStructure::GetID_str( const int& model_num, const gemmi::CRA& cra, const String& entity )
{
  return GetID_str( String( model_num, 0 ), cra, entity );
}

/*! Set spacegroup
  \param spacegroup Clipper spacegroup. */
void GemmiStructure::set_spacegroup( const Spacegroup& spacegroup )
{
  spacegroup_hm = ( char* )spacegroup.symbol_hm().c_str();
}

/*! Set cell
  \param cell_in Clipper cell. */
void GemmiStructure::set_cell( const Cell& cell_in )
{ 
  this->cell = GEMMI::cell( cell_in );
}

/*! \return Get spacegroup from GEMMI structure and return in
  clipper::Spacegroup format. */
Spacegroup GemmiStructure::spacegroup() const
{
  if ( String( spacegroup_hm ) != "" ) {
    const ::gemmi::SpaceGroup* sg = this->find_spacegroup();
    return GEMMI::spacegroup( *sg );
  } else
    return Spacegroup::null();
}

/*! \return Get cell from GEMMI structure and return in
  clipper::Cell format. */
Cell GemmiStructure::get_cell() const
{ 
  return GEMMI::cell( this->cell );
}

/*! \return Experimental method (e.g. X-ray, EM, etc.). */
String GemmiStructure::get_exptlmethod() const
{ 
  return get_info( "_exptl.method" );
}

/*! \return Get transformation/origx from GEMMI Structure
  and return in clipper::RTop_orth format*/
RTop_orth GemmiStructure::get_transform() const
{ 
  return GEMMI::transform( origx );
}

/*! Set transformation/origx to GEMMI Structure
  \param rtop Rotational and translation operator to be set. */
void GemmiStructure::set_transform( const RTop_orth& rtop )
{ 
  origx = GEMMI::transform( rtop );
}

/*! \return Get resolution from GEMMI Structure. */
Resolution GemmiStructure::get_resolution() const
{ 
  return Resolution( resolution );
}

/*! Set resolution to GEMMI Structure.
\param resol The resolution limit in Angstrom. */
void GemmiStructure::set_resolution( const Resolution& resol )
{ 
  resolution = resol.limit();
}

/*! Conversion of gemmi CraProxy to Atom_list
  \param cra gemmi's CraProxy. */
GemmiAtom_list::GemmiAtom_list( gemmi::CraProxy cra )
{
  for ( auto a : cra ) push_back( Atom( *( ( const GemmiAtom* )a.atom ) ) );
}
}  // namespace clipper
