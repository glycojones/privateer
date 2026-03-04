/* clipper_gemmi.cpp: Gemmi wrapper class file. */

#include "clipper_gemmi.h"
#include <map>

namespace clipper {
/* Internal: Get dataset name and a list of column names from clipper path name. */
const std::vector<String> extract_column_names( const String& assign, const int& f_size ) {
  std::vector<String> col_names( f_size, "MNF" );
  // interpret list name in terms of columns
  if ( assign.find( "[" ) == String::npos ) {
    // name is a single path
    col_names[0] = assign.split( "/" )[2];
  } else {
    // name is a list of mtz columns: extract column names from list
    String pref = assign.substr( 0, assign.find( "[" ) );
    String post =
        assign.substr( assign.find( "[" ) + 1, assign.find( "]" ) - assign.find( "[" ) - 1 );
    std::vector<String> list = post.split( ", " );
    for ( int i = 0; i < list.size(); i++ ) col_names[i] = list[i];
  }
  return col_names;
}

/*! Convert Clipper Spacegroup to Gemmi SpaceGroup through symmetry operations
  search from the spacegroup given.
  \param spgr Clipper::Spacegroup object. */
const gemmi::SpaceGroup* GEMMI::spacegroup( const Spacegroup& spgr ) {
  gemmi::GroupOps symops = gemmi::symops_from_hall( spgr.symbol_hall().c_str() );
  return gemmi::find_spacegroup_by_ops( symops );
}

/*! Convert Gemmi SpaceGroup to Clipper Spacegroup through symmetry operations
  of the spacegroup given.
  \param spgr gemmi::SpaceGroup object. */
Spacegroup GEMMI::spacegroup( const gemmi::SpaceGroup& spgr ) {
  String ops;
  for ( auto op : spgr.operations() ) ops += op.triplet() + ";";
  return Spacegroup( Spgr_descr( ops, Spgr_descr::Symops ) );
}

/*! Return Gemmi Transform from Clipper RTop_orth.
  \param rtop Rotation-translation operator. */
gemmi::Transform GEMMI::transform( const RTop_orth& rtop ) {
  gemmi::Transform grtop;  // identity mat33, zero vec3
  if ( !rtop.rot().is_null() ) {
    for ( int i = 0; i < 3; i++ ) {
      grtop.mat[i][0] = rtop.rot()( i, 0 );
      grtop.mat[i][1] = rtop.rot()( i, 1 );
      grtop.mat[i][2] = rtop.rot()( i, 2 );
    }
  }
  if ( !rtop.trn().is_null() ) {
    for ( int i = 0; i < 3; i++ ) {
      grtop.vec.at( i ) = rtop.trn()[i];
    }
  }
  return grtop;
}

/*! Return Clipper RTop_orth from gemmi Transfrom.
  \param grtop gemmi Transform object. */
RTop_orth GEMMI::transform( const gemmi::Transform& grtop ) {
  Mat33<> m( grtop.mat[0][0], grtop.mat[0][1], grtop.mat[0][2], grtop.mat[1][0], grtop.mat[1][1],
             grtop.mat[1][2], grtop.mat[2][0], grtop.mat[2][1], grtop.mat[2][2] );
  Vec3<> v( grtop.vec.at( 0 ), grtop.vec.at( 1 ), grtop.vec.at( 2 ) );
  return RTop_orth( m, v );
}

/*! Return Gemmi Miller from Clipper HKL
  \param hkl HKL object. */
Miller GEMMI::Hkl( const HKL& hkl ) { return Miller{ hkl.h(), hkl.k(), hkl.l() }; }

/*! Return Clipper HKL from Gemmi Miller
  \param hkl gemmi Miller object; std::array<int,3>. */
HKL GEMMI::Hkl( const Miller& hkl ) { return HKL( hkl[0], hkl[1], hkl[2] ); }

/*! Return gemmi UnitCell from clipper Cell.
  \param cell Cell object. */
gemmi::UnitCell GEMMI::cell( const Cell& cell ) {
  return gemmi::UnitCell( cell.a(), cell.b(), cell.c(), cell.alpha_deg(), cell.beta_deg(),
                          cell.gamma_deg() );
}

/*! Return clipper Cell from gemmi UnitCell.
  \param cell gemmi UnitCell object. */
Cell GEMMI::cell( const gemmi::UnitCell& cell ) {
  return Cell( Cell_descr( cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma ) );
}

/*! Return initialised HKL_info using data from Gemmi::Mtz which
  then can be used to construct HKL_data. By default, reflections
  list is taken frmo Mtz obj.
  \param mtzobj Gemmi Mtz class where the data is store.
  \param tolerance Tolerance limit.
  \param generate If true, a reflection list will be generated for an ASU. */
HKL_info GEMMI::as_HKL_info( const gemmi::Mtz& mtzobj, double tolerance, const bool& generate ) {
  if ( mtzobj.max_1_d2 == 0 ) Message::message( Message_fatal( "GEMMI: max_d_star_sq is zero." ) );
  if ( generate )
    return as_HKL_info( mtzobj.cell, *mtzobj.spacegroup, mtzobj.resolution_high(), tolerance,
                        generate );
  else {
    std::vector<HKL> hkllist;
    hkllist.reserve( mtzobj.nreflections );
    for ( size_t i = 0; i < mtzobj.data.size(); i += mtzobj.columns.size() )
      hkllist.push_back( Hkl( mtzobj.get_hkl( i ) ) );

    HKL_info result = as_HKL_info( mtzobj.cell, *mtzobj.spacegroup, mtzobj.resolution_high(),
                                   tolerance, generate );
    result.add_hkl_list( hkllist );
    return result;
  }
}

/*! Return initialised HKL_info using Gemmi unit cell, spacegroup and resolution limit
  which then can be used to construct HKL_data.
  \param unit_cell The unit cell
  \param sg The spacegroup
  \param d_min The resolution limit
  \param tolerance Tolerance limit.
  \param generate If true, a reflection list will be generated for an ASU. */
HKL_info GEMMI::as_HKL_info( const gemmi::UnitCell& unit_cell, const gemmi::SpaceGroup& sg,
                             double d_min, double tolerance, const bool& generate ) {
  return HKL_info( spacegroup( sg ), cell( unit_cell ), Resolution( d_min - tolerance ), generate );
}

/*! Return initialised HKL_info using Gemmi unit cell, spacegroup, and miller indices
  which then can be used to construct HKL_data.
  \param unit_cell The unit cell
  \param sg The spacegroup
  \param miller Vector of Miller indices
  \param tolerance Tolerance limit. */
HKL_info GEMMI::as_HKL_info( const gemmi::UnitCell& unit_cell, const gemmi::SpaceGroup& sg,
                             const std::vector<gemmi::Miller>& miller_indices, double tolerance ) {
  std::vector<HKL> hkl_list;
  hkl_list.reserve( miller_indices.size() );
  double max_d_star_sq = 0;
  for ( std::size_t i = 0; i < miller_indices.size(); i++ ) {
    max_d_star_sq = std::max( max_d_star_sq, unit_cell.calculate_1_d2( miller_indices[i] ) );
    hkl_list.push_back( Hkl( miller_indices[i] ) );
  }
  if ( max_d_star_sq == 0 ) Message::message( Message_fatal( "GEMMI: max_d_star_sq is zero." ) );
  HKL_info result = as_HKL_info( unit_cell, sg, 1 / std::sqrt( max_d_star_sq ), tolerance );
  result.add_hkl_list( hkl_list );
  return result;
}

/*! Import data from gemmi's Mtz object into HKL_data object,
  given the column paths. Preferred general way to import HKL_data from gemmi Mtz.

  An MTZ column type must be present in the MTZ_type_registry for the
  HKL_data type element name concerned.

  \param cdata The HKL_data object into which data is to be imported.
  \param mtzobj Gemmi Mtz class where the data is store.
  \param mtzpath The MTZ column names, as a path. See \ref MTZpaths for details. */
void GEMMI::import_hkl_data( HKL_data_base& cdata, const gemmi::Mtz& mtzobj,
                             const String mtzpath ) {
  // check if HKL_data is initialised
  if ( cdata.is_null() ) Message::message( Message_fatal( "GEMMI: HKL_data is not initialised." ) );
  String colpath = mtzpath;
  // add exported data columns to local list
  int ncols = cdata.data_size();
  size_t cols[ncols];
  std::vector<String> col_names = extract_column_names( colpath, ncols );
  std::vector<String> dat_names = cdata.data_names().split( " " );
  std::vector<ftype> scls( ncols, 1.0 );
  String dataset_name = colpath.split( "/" )[1];
  // get the column indices
  for ( int c = 0; c < ncols; c++ ) {
    if ( col_names[c] != "MNF" && col_names[c] != "NAN" && col_names[c] != "mnf" &&
         col_names[c] != "nan" ) {
      if ( dataset_name != "*" )
        cols[c] = mtzobj
                      .column_with_label( col_names[c].c_str(),
                                          mtzobj.dataset_with_name( dataset_name.c_str() ) )
                      ->idx;
      else
        cols[c] = mtzobj.column_with_label( col_names[c].c_str() )->idx;
    }
    if ( dat_names[c] == "phi" ) scls[c] = Util::rad2d( 1.0 );
  }
  // loop through nreflections,data
  xtype values[ncols];
  for ( size_t i = 0; i < mtzobj.data.size(); i += mtzobj.columns.size() ) {
    // import data
    for ( int c = 0; c < ncols; c++ ) {
      // set to Nan unless readable and present
      Util::set_null( values[c] );
      if ( !Util::is_nan( xtype( mtzobj.data[i + cols[c]] ) ) )
        values[c] = xtype( mtzobj.data[i + cols[c]] / scls[c] );
    }
    cdata.data_import( Hkl( mtzobj.get_hkl( i ) ), values );
  }
}

/*! Convert vectors of F and sigF values as HKL_data
  \param hkl_info The HKL_info object
  \param miller_indices std::vector containing Gemmi's Miller indices
  \param data std::vector containing F values
  \param sigmas std::vector containing sigF values. */
HKL_data<data64::F_sigF> GEMMI::as_HKL_data_fsigf( HKL_info& hkl_info,
                                                   const std::vector<gemmi::Miller>& miller_indices,
                                                   const std::vector<double>& data,
                                                   const std::vector<double>& sigmas ) {
  if ( hkl_info.is_null() )
    Message::message( Message_fatal( "GEMMI: HKL_info is not initialised." ) );
  if ( data.size() != miller_indices.size() )
    Message::message(
        Message_fatal( "GEMMI: Vectors for data and miller indices have different lengths." ) );
  if ( sigmas.size() != miller_indices.size() )
    Message::message(
        Message_fatal( "GEMMI: Vectors for sigmas and miller indices have different lengths." ) );
  HKL_data<data64::F_sigF> hkl_data( hkl_info );
  for ( std::size_t i = 0; i < miller_indices.size(); i++ ) {
    data64::F_sigF datum;
    datum.f() = data[i];
    datum.sigf() = sigmas[i];
    if ( !hkl_data.set_data( Hkl( miller_indices[i] ), datum ) )
      Message::message(
          Message_fatal( "GEMMI: Unable to set data for " + Hkl( miller_indices[i] ).format() ) );
  }
  return hkl_data;
}

/*! Convert complex structure factors, gemmi::AsuData<complex<double>> to HKL_data.
    HKL_info has to be initialiased with reflections list.
  \param hkl_info The HKL_info object
  \param data gemmi's AsuData containing complex numbers for amplitude and phase. */
HKL_data<data64::F_phi> GEMMI::as_HKL_data_fphi(
    HKL_info& hkl_info, const gemmi::AsuData<std::complex<double> >& data ) {
  if ( hkl_info.is_null() )
    Message::message( Message_fatal( "GEMMI: HKL_info is not initialised." ) );
  HKL_data<data64::F_phi> hkl_data( hkl_info );
  for ( std::size_t i = 0; i < data.size(); i++ ) {
    if ( !hkl_data.set_data( Hkl( data.get_hkl( i ) ), data64::F_phi( data.v[i].value ) ) )
      Message::message(
          Message_fatal( "GEMMI: Unable to set data for " + Hkl( data.get_hkl( i ) ).format() ) );
  }
  return hkl_data;
}

/*! Convert vectors of amplitude and phases to HKL_data.
  \param hkl_info The HKL_info object
  \param miller_indicues Vector of reflections
  \param data_f Vector of amplitude values
  \param data_phi Vector of amplitude values.*/
HKL_data<data64::F_phi> GEMMI::as_HKL_data_fphi( HKL_info& hkl_info,
                                                 const std::vector<gemmi::Miller>& miller_indices,
                                                 const std::vector<double>& data_f,
                                                 const std::vector<double>& data_phi ) {
  if ( hkl_info.is_null() )
    Message::message( Message_fatal( "GEMMI: HKL_info is not initialised." ) );
  if ( ( data_f.size() != miller_indices.size() ) || ( data_phi.size() != miller_indices.size() ) )
    Message::message(
        Message_fatal( "GEMMI: Vectors for data and miller indices have different lengths." ) );
  HKL_data<data64::F_phi> hkl_data( hkl_info );
  for ( std::size_t i = 0; i < data_f.size(); i++ ) {
    if ( !hkl_data.set_data( Hkl( miller_indices[i] ), data64::F_phi( data_f[i], data_phi[i] ) ) )
      Message::message(
          Message_fatal( "GEMMI: Unable to set data for " + Hkl( miller_indices[i] ).format() ) );
  }
  return hkl_data;
}

/*! Extract complex structure factors from HKL_data into std::vector.
    \param hkl_data HKL_data containing structure factors
    \param miller_indices Vector of reflections. */
std::vector<std::complex<double> > GEMMI::extract_complex(
    const HKL_data<data64::F_phi>& hkl_data, const std::vector<Miller>& miller_indices ) {
  std::vector<std::complex<double> > result( miller_indices.size() );
  data64::F_phi datum;
  for ( std::size_t i = 0; i < miller_indices.size(); i++ ) {
    if ( !hkl_data.get_data( Hkl( miller_indices[i] ), datum ) )
      Message::message(
          Message_fatal( "CLIPPER: Unable to get data for " + Hkl( miller_indices[i] ).format() ) );
    if ( datum.missing() )
      Message::message(
          Message_fatal( "CLIPPER: Missing F_phi data for " + Hkl( miller_indices[i] ).format() ) );
    result.push_back( datum );
  }
  return result;
}

/*! Converts Hendrickson-Lattman coefficients in std::vector<const gemmi::Mtz::Column*> returned
  from gemmi method Mtz::columns_with_type('A') to clipper's HKL_data<data64::ABCD>.
  HKL_info has to be initialiased with reflections list.
  \param hkl_info HKL_info object
  \param miller_indices Vector of reflections
  \param data Vector of Mtz::Column*. */
HKL_data<data64::ABCD> GEMMI::as_HKL_data( HKL_info& hkl_info,
                                           const std::vector<gemmi::Miller>& miller_indices,
                                           const std::vector<const gemmi::Mtz::Column*>& data ) {
  if ( hkl_info.is_null() )
    Message::message( Message_fatal( "GEMMI: HKL_info is not initialised." ) );
  if ( data.size() != 4 )
    Message::message(
        Message_fatal( "GEMMI: Size for vector<const gemmi::Mtz::Column*> data not 4." ) );
  if ( !std::all_of( data.cbegin(), data.cend(), [&]( const gemmi::Mtz::Column* c ) {
         return c->size() == miller_indices.size();
       } ) )
    Message::message(
        Message_fatal( "GEMMI: Vectors for data and miller indices have different lengths." ) );

  HKL_data<data64::ABCD> hkl_data( hkl_info );
  data64::ABCD abcd;
  for ( std::size_t i = 0; i < miller_indices.size(); i++ ) {
    double datum[4] = { data[0]->at( i ), data[1]->at( i ), data[2]->at( i ), data[3]->at( i ) };
    abcd.data_import( datum );
    if ( !hkl_data.set_data( Hkl( miller_indices[i] ), abcd ) )
      Message::message(
          Message_fatal( "GEMMI: Unable to set data for " + Hkl( miller_indices[i] ).format() ) );
  }
  return hkl_data;
}

/*! Extract Hendrickson-Lattman coefficients from HKL_data into std::vector<std::array<double,4>>.
  \param hkl_data HKL_data containing Hendrickson-Lattman coefficients
  \param miller_indices Vector of reflections. */
std::vector<std::array<double, 4> > GEMMI::extract_hendrickson_lattman(
    const HKL_data<data64::ABCD>& hkl_data, const std::vector<gemmi::Miller>& miller_indices ) {
  std::vector<std::array<double, 4> > result( miller_indices.size() );
  data64::ABCD datum;
  for ( std::size_t i = 0; i < miller_indices.size(); i++ ) {
    if ( !hkl_data.get_data( Hkl( miller_indices[i] ), datum ) )
      Message::message(
          Message_fatal( "CLIPPER: Failed to retrieve Hendrickson-Lattman coefficients for " +
                         Hkl( miller_indices[i] ).format() ) );
    if ( datum.missing() )
      Message::message( Message_fatal( "CLIPPER: Missing Hendrickson-Lattman coefficients for " +
                                       Hkl( miller_indices[i] ).format() ) );
    std::array<double, 4> hl_coef = { datum.a(), datum.b(), datum.c(), datum.d() };
    result.push_back( hl_coef );
  }
  return result;
}

/*! Extracting centroid phases from HKL_data into std::vector<double>.
  \param hkl_data HKL_data containing centroid phases
  \param miller_indices Vector of reflections. */
std::vector<double> GEMMI::extract_centroid_phases(
    const HKL_data<data64::Phi_fom>& hkl_data, const std::vector<gemmi::Miller>& miller_indices ) {
  std::vector<double> result( miller_indices.size() );
  data64::Phi_fom datum;
  for ( std::size_t i = 0; i < miller_indices.size(); i++ ) {
    if ( !hkl_data.get_data( Hkl( miller_indices[i] ), datum ) )
      Message::message( Message_fatal( "CLIPPER: Failed to retrieve Phi for " +
                                       Hkl( miller_indices[i] ).format() ) );
    if ( datum.missing() )
      Message::message(
          Message_fatal( "CLIPPER: Missing Phi for " + Hkl( miller_indices[i] ).format() ) );
    result.push_back( datum.phi() );
  }
  return result;
}

/*! Extracting figure of merit from HKL_data into std::vector<double>.
  \param hkl_data HKL_data containing figure of merit
  \param miller_indices Vector of reflections. */
std::vector<double> GEMMI::extract_figure_of_merit(
    const HKL_data<data64::Phi_fom>& hkl_data, const std::vector<gemmi::Miller>& miller_indices ) {
  std::vector<double> result( miller_indices.size() );
  data64::Phi_fom datum;
  for ( std::size_t i = 0; i < miller_indices.size(); i++ ) {
    if ( !hkl_data.get_data( Hkl( miller_indices[i] ), datum ) )
      Message::message( Message_fatal( "CLIPPER: Failed to retrieve FOM for " +
                                       Hkl( miller_indices[i] ).format() ) );
    if ( datum.missing() )
      Message::message(
          Message_fatal( "CLIPPER: Missing FOM for " + Hkl( miller_indices[i] ).format() ) );
    result.push_back( datum.fom() );
  }
  return result;
}

// compile templates
template void GEMMI::import_xmap<ftype32>( Xmap<ftype32>& xmap,
                                           const gemmi::Ccp4<ftype32>& mapobj );
template void GEMMI::import_xmap<ftype64>( Xmap<ftype64>& xmap,
                                           const gemmi::Ccp4<ftype32>& mapobj );
template void GEMMI::import_nxmap<ftype32>( NXmap<ftype32>& xmap,
                                            const gemmi::Ccp4<ftype32>& mapobj );
template void GEMMI::import_nxmap<ftype64>( NXmap<ftype64>& xmap,
                                            const gemmi::Ccp4<ftype32>& mapobj );

template void GEMMI::export_xmap<ftype32>( const Xmap<ftype32>& xmap,
                                           gemmi::Ccp4<ftype32>& mapobj );
template void GEMMI::export_xmap<ftype64>( const Xmap<ftype64>& xmap,
                                           gemmi::Ccp4<ftype32>& mapobj );
template void GEMMI::export_nxmap<ftype32>( const NXmap<ftype32>& xmap,
                                            gemmi::Ccp4<ftype32>& mapobj, const Cell& unitcell );
template void GEMMI::export_nxmap<ftype64>( const NXmap<ftype64>& xmap,
                                            gemmi::Ccp4<ftype32>& mapobj, const Cell& unitcell );

}  // namespace clipper
