/*! \file gemmi/clipper_gemmi.h
    \author Soon Wen Hoh
    Header file for clipper-gemmi HKL data and map class type conversion
*/

#ifndef CLIPPER_GEMMI
#define CLIPPER_GEMMI

#include "../core/hkl_datatypes.h"
#include "../core/nxmap.h"
#include "../core/xmap.h"
#include <gemmi/asudata.hpp> // AsuData
#include <gemmi/ccp4.hpp>    // Ccp4
#include <gemmi/mtz.hpp>     // Mtz, math, symmetry, unitcell, <algorithm> and <vector>

namespace clipper {
using Miller = gemmi::Miller;
//! Gemmi interface class
/*! This class contains static methods which convert various objects
  between Clipper and Gemmi formats. See https://gemmi.readthedocs.io/en/latest/
*/
class GEMMI {
 public:
  //! Convert Clipper Spacegroup to Gemmi SpaceGroup
  static const ::gemmi::SpaceGroup* spacegroup( const Spacegroup& spgr );
  //! Convert Gemmi SpaceGroup to Clipper Spacegroup
  static Spacegroup spacegroup( const gemmi::SpaceGroup& spgr );
  //! Convert Clipper RTop_orth to Gemmi Transfrom
  static ::gemmi::Transform transform( const RTop_orth& rtop );
  //! Convert Gemmi Transfrom to Clipper RTop_orth
  static RTop_orth transform( const gemmi::Transform& rtop );

  //! Convert Clipper HKL to Gemmi Miller
  static Miller Hkl( const HKL& hkl );
  //! Convert Gemmi Miller to Clipper HKL
  static HKL Hkl( const Miller& hkl );
  //! Convert Clipper HKL to Gemmi Miller
  /*! This version is inlined and uses reinterpret cast, and so
      is very fast, but depends on the binary representation of use
      in the two packages remaining compatible. */
  static inline const Miller& hkl( const HKL& hkl ) {
    return reinterpret_cast<const Miller&>( hkl );
  }
  //! Convert Gemmi Miller to Clipper HKL
  /*! This version is inlined and uses reinterpret cast, and so
      is very fast, but depends on the binary representation of use
      in the two packages remaining compatible. */
  static inline const HKL& hkl( const Miller& hkl ) { return reinterpret_cast<const HKL&>( hkl ); }

  //! Convert Clipper cell to Gemmi cell
  static gemmi::UnitCell cell( const Cell& cell );
  //! Convert Gemmi cell to Clipper cell
  static Cell cell( const gemmi::UnitCell& cell );

  //! Set HKL_info from Gemmi::Mtz
  static HKL_info as_HKL_info( const gemmi::Mtz& mtzobj, double tolerance = 1.e-8,
                               const bool& generate = false );
  //! Conversion of Gemmi unit cell and spacegroup to HKL_info
  static HKL_info as_HKL_info( const gemmi::UnitCell& unit_cell, const gemmi::SpaceGroup& sg,
                               double d_min, double tolerance = 1.e-8,
                               const bool& generate = false );
  //! Conversion of  Gemmi unit cell and spacegroup to HKL_info
  static HKL_info as_HKL_info( const gemmi::UnitCell& unit_cell, const gemmi::SpaceGroup& sg,
                               const std::vector<Miller>& miller_indices,
                               double tolerance = 1.e-8 );

  //! Import data from gemmi's Mtz object into HKL_data object.
  static void import_hkl_data( HKL_data_base& cdata, const gemmi::Mtz& mtzobj,
                               const String mtzpath );

  //! Conversion of F, sigF as HKL_data
  static HKL_data<data64::F_sigF> as_HKL_data_fsigf( HKL_info& hkl_info,
                                                     const std::vector<Miller>& miller_indices,
                                                     const std::vector<double>& data,
                                                     const std::vector<double>& sigmas );

  //! Conversion of complex structure factors.
  static HKL_data<data64::F_phi> as_HKL_data_fphi(
      HKL_info& hkl_info, const gemmi::AsuData<std::complex<double> >& data );
  //! Conversion of vectors of F and Phi values to HKL_data<data64::F_phi>
  static HKL_data<data64::F_phi> as_HKL_data_fphi( HKL_info& hkl_info,
                                                   const std::vector<Miller>& miller_indices,
                                                   const std::vector<double>& data_f,
                                                   const std::vector<double>& data_phi );
  //! Conversion of complex structure factors.
  static std::vector<std::complex<double> > extract_complex(
      const HKL_data<data64::F_phi>& hkl_data, const std::vector<Miller>& miller_indices );

  //! Conversion of Hendrickson-Lattman coefficients
  static HKL_data<data64::ABCD> as_HKL_data( HKL_info& hkl_info,
                                             const std::vector<Miller>& miller_indices,
                                             const std::vector<const gemmi::Mtz::Column*>& data );
  //! Conversion of Hendrickson-Lattman coefficients.
  static std::vector<std::array<double, 4> > extract_hendrickson_lattman(
      const HKL_data<data64::ABCD>& hkl_data, const std::vector<Miller>& miller_indices );

  //! Conversion of centroid phases.
  static std::vector<double> extract_centroid_phases( const HKL_data<data64::Phi_fom>& hkl_data,
                                                      const std::vector<Miller>& miller_indices );
  //! Conversion of figures of merit
  static std::vector<double> extract_figure_of_merit( const HKL_data<data64::Phi_fom>& hkl_data,
                                                      const std::vector<Miller>& miller_indices );

  //! Import map data from gemmi's Ccp4.grid into Xmap
  template <class T>
  static void import_xmap( Xmap<T>& xmap, const gemmi::Ccp4<ftype32>& mapobj );
  //! Import map data from gemmi's Ccp4.grid into NXmap
  template <class T>
  static void import_nxmap( NXmap<T>& nxmap, const gemmi::Ccp4<ftype32>& mapobj );
  //! Export map data to gemmi's Ccp4.grid from Xmap
  template <class T>
  static void export_xmap( const Xmap<T>& xmap, gemmi::Ccp4<ftype32>& mapobj );
  //! Export map data to gemmi's Ccp4.grid from NXmap
  template <class T>
  static void export_nxmap( const NXmap<T>& nxmap, gemmi::Ccp4<ftype32>& mapobj,
                            const Cell& unitcell );

 private:
  /* Internal: Prepare gemmi::Ccp4 map object header.
  mapobj - gemmi Ccp4 map object
  dim - dimension of grid in std::array<int, 3>
  gfms0 - NXSTART, NYSTART, NZSTART; location of first column, row, section
  grid - full grid size
  orderfms - fast medium slow order
  unitcell - unit cell */
  static void prepare_gemmi_header( gemmi::Ccp4<ftype32>& mapobj, std::array<int, 3> dim,
                                    std::array<int, 3> gfms0, std::array<int, 3> grid,
                                    std::array<int, 3> orderfms, Cell unitcell ) {
    mapobj.prepare_ccp4_header_except_mode_and_stats();
    mapobj.set_header_3i32( 1, dim[0], dim[1], dim[2] );
    mapobj.set_header_3i32( 5, gfms0[0], gfms0[1], gfms0[2] );
    mapobj.set_header_3i32( 8, grid[0], grid[1], grid[2] );
    // fast, medium ,slow order
    mapobj.set_header_3i32( 17, orderfms[0], orderfms[1], orderfms[2] );
    mapobj.set_header_float( 11, ( float )unitcell.a() );
    mapobj.set_header_float( 12, ( float )unitcell.b() );
    mapobj.set_header_float( 13, ( float )unitcell.c() );
    mapobj.set_header_float( 14, ( float )unitcell.alpha_deg() );
    mapobj.set_header_float( 15, ( float )unitcell.beta_deg() );
    mapobj.set_header_float( 16, ( float )unitcell.gamma_deg() );
    mapobj.setup( NAN );
    if ( mapobj.grid.spacegroup->number != 1 ) mapobj.grid.symmetrize_max();
    mapobj.update_ccp4_header( 2 );
  }

};  // class GEMMI

// Implementations templates
/*! Import grid data from gemmi::Ccp4::grid into Xmap.
    Run setup(0) on gemmi::Ccp4 map object before importing grid data to Xmap.
  \param xmap The Xmap to be imported.
  \param mapobj The gemmi::Ccp4 map containing the grid data. */
template <class T>
void GEMMI::import_xmap( Xmap<T>& xmap, const gemmi::Ccp4<ftype32>& mapobj ) {
  // check if HKL_info params are already set
  Spacegroup s = xmap.spacegroup();
  Cell c = xmap.cell();
  Grid_sampling r = xmap.grid_sampling();
  std::array<int, 3> orderfms, orderxyz, gfms0, gfms1, dim, grid, g;
  dim = mapobj.header_3i32( 1 );
  orderfms = mapobj.axis_positions();
  gfms0 = mapobj.header_3i32( 5 );
  grid = mapobj.header_3i32( 8 );
  // import any missing params
  if ( s.is_null() ) s = GEMMI::spacegroup( *mapobj.grid.spacegroup );
  if ( c.is_null() ) c = cell( mapobj.grid.unit_cell );
  if ( r.is_null() ) r = Grid_sampling( grid[0], grid[1], grid[2] );
  xmap.init( s, c, r );
  // get grid bound and axis order
  for ( int i = 0; i < 3; i++ ) {
    gfms1[i] = gfms0[i] + dim[i] - 1;
    orderxyz[orderfms[i]] = i;
  }
  // copy map data from gemmi map object
  Xmap_base::Map_reference_coord x( xmap );
  bool zero_orig =
      ( std::all_of( gfms0.cbegin(), gfms0.cend(), []( const int index ) { return index == 0; } ) );
  for ( g[2] = gfms0[2]; g[2] <= gfms1[2]; g[2]++ ) {
    for ( g[1] = gfms0[1]; g[1] <= gfms1[1]; g[1]++ ) {
      for ( g[0] = gfms0[0]; g[0] <= gfms1[0]; g[0]++ ) {
        x.set_coord( Coord_grid( g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]] ) );
        if ( zero_orig )
          xmap[x] = T( mapobj.grid.get_value_q( g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]] ) );
        else
          xmap[x] = T( mapobj.grid.get_value( g[orderxyz[0]] - gfms0[orderxyz[0]],
                                              g[orderxyz[1]] - gfms0[orderxyz[1]],
                                              g[orderxyz[2]] - gfms0[orderxyz[2]] ) );
      }
    }
  }  // finish copy data
}

/*! Import grid data from gemmi::Ccp4::grid into NXmap.
  \param nxmap The NXmap to be imported.
  \param mapobj The gemmi::Ccp4 map containing the grid data. */
template <class T>
void GEMMI::import_nxmap( NXmap<T>& nxmap, const gemmi::Ccp4<ftype32>& mapobj ) {
  std::array<int, 3> orderfms, orderxyz, gfms0, gfms1, dim, g, grid;
  orderxyz = mapobj.axis_positions();  // starts from 0
  dim = mapobj.header_3i32( 1 );
  gfms0 = mapobj.header_3i32( 5 );
  grid = mapobj.header_3i32( 8 );
  Grid_sampling r( grid[0], grid[1], grid[2] );
  for ( int i = 0; i < 3; i++ ) {
    gfms1[i] = gfms0[i] + dim[i] - 1;
  }
  Grid_range grid_map( Coord_grid( gfms0[orderxyz[0]], gfms0[orderxyz[1]], gfms0[orderxyz[2]] ),
                       Coord_grid( gfms1[orderxyz[0]], gfms1[orderxyz[1]], gfms1[orderxyz[2]] ) );
  nxmap.init( cell( mapobj.grid.unit_cell ), r, grid_map );
  // copy map data from gemmi map object
  for ( g[2] = 0; g[2] <= gfms1[2] - gfms0[2]; g[2]++ ) {
    for ( g[1] = 0; g[1] <= gfms1[1] - gfms0[1]; g[1]++ ) {
      for ( g[0] = 0; g[0] <= gfms1[0] - gfms0[0]; g[0]++ ) {
        nxmap.set_data(
            Coord_grid( g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]] ),
            T( mapobj.grid.get_value_q( g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]] ) ) );
      }
    }
  }  // finish copy data
}

/*! Export map data to gemmi's Ccp4::grid from Xmap.
  \param xmap The Xmap to be exported.
  \param mapobj The gemmi::Ccp4 map to hold the grid data.*/
template <class T>
void GEMMI::export_xmap( const Xmap<T>& xmap, gemmi::Ccp4<ftype32>& mapobj ) {
  std::array<int, 3> orderfms, orderxyz, grid, gfms0, gfms1, dim;
  int spg = xmap.spacegroup().descr().spacegroup_number();
  //  use axis order 1,2,3 (fast, medium, slow) for gemmi Ccp4 map
  orderxyz = { 0, 1, 2 };
  orderfms = { 1, 2, 3 };
  // grids
  for ( int i = 0; i < 3; i++ ) {
    grid[i] = xmap.grid_sampling()[i];
    gfms0[orderxyz[i]] = xmap.grid_asu().min()[i];
    gfms1[orderxyz[i]] = xmap.grid_asu().max()[i];
  }
  for ( int i = 0; i < 3; i++ ) dim[i] = gfms1[i] - gfms0[i] + 1;
  // prepare gemmi::Ccp4::grid
  mapobj.grid.spacegroup = GEMMI::spacegroup( xmap.spacegroup() );
  mapobj.grid.set_unit_cell( cell( xmap.cell() ) );
  mapobj.grid.set_size( grid[0], grid[1], grid[2] );
  mapobj.grid.axis_order = gemmi::AxisOrder::Unknown;
  // copy map data
  int g[3];
  Xmap_base::Map_reference_coord x( xmap );
  for ( g[2] = gfms0[2]; g[2] <= gfms1[2]; g[2]++ ) {
    for ( g[1] = gfms0[1]; g[1] <= gfms1[1]; g[1]++ ) {
      for ( g[0] = gfms0[0]; g[0] <= gfms1[0]; g[0]++ ) {
        x.set_coord( Coord_grid( g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]] ) );
        mapobj.grid.set_value( g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]], ftype32( xmap[x] ) );
      }
    }
  }
  prepare_gemmi_header( mapobj, dim, gfms0, grid, orderfms, xmap.cell() );
}

/*! Export map data to gemmi's Ccp4::grid from NXmap.
  \param nxmap The NXmap to be exported.
  \param mapobj The gemmi::Ccp4 map to hold the grid data.
  \param unitcell The unit cell for the map.*/
template <class T>
void GEMMI::export_nxmap( const NXmap<T>& nxmap, gemmi::Ccp4<ftype32>& mapobj,
                          const Cell& unitcell ) {
  std::array<int, 3> grid, gfms0, gfms1, dim;
  std::array<int, 3> orderfms{ 1, 2, 3 };
  // grids
  Coord_frac c0, c1;
  c0 = nxmap.coord_orth( Coord_map( 0, 0, 0 ) ).coord_frac( unitcell );
  c1 = nxmap.coord_orth( Coord_map( nxmap.grid().nu(), nxmap.grid().nv(), nxmap.grid().nw() ) )
           .coord_frac( unitcell );
  Grid_sampling grid_sam =
      Grid_sampling( Util::intr( ftype( nxmap.grid().nu() ) / ( c1.u() - c0.u() ) ),
                     Util::intr( ftype( nxmap.grid().nv() ) / ( c1.v() - c0.v() ) ),
                     Util::intr( ftype( nxmap.grid().nw() ) / ( c1.w() - c0.w() ) ) );
  Coord_grid g0 = c0.coord_grid( grid_sam );
  Coord_grid g1 = g0 + Coord_grid( nxmap.grid() ) - Coord_grid( 1, 1, 1 );
  Grid_range grid_map = Grid_range( g0, g1 );
  for ( int i = 0; i < 3; i++ ) {
    grid[i] = grid_sam[i];
    gfms0[i] = grid_map.min()[i];
    gfms1[i] = grid_map.max()[i];
    dim[i] = gfms1[i] - gfms0[i] + 1;
  }
  // prepare gemmi::Ccp4::grid
  mapobj.grid.spacegroup = GEMMI::spacegroup( Spacegroup::p1() );
  mapobj.grid.set_unit_cell( cell( unitcell ) );
  mapobj.grid.set_size( grid[0], grid[1], grid[2] );
  mapobj.grid.axis_order = gemmi::AxisOrder::Unknown;
  // copy map data
  int g[3];
  NXmap_base::Map_reference_coord i0, iu, iv, iw;
  i0 = NXmap_base::Map_reference_coord( nxmap, g0 );
  for ( iu = i0; iu.coord().u() <= g1.u(); iu.next_u() ) {
    for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() ) {
      for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {
        mapobj.grid.set_value( iw.coord().u(), iw.coord().v(), iw.coord().w(),
                               ftype32( nxmap[iw] ) );
      }
    }
  }
  prepare_gemmi_header( mapobj, dim, gfms0, grid, orderfms, unitcell );
}

}  // namespace clipper

#endif
