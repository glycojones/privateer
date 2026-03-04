/* map_io.cpp: class file for reflection data  map importer */
//C Copyright (C) 2000-2006 Kevin Cowtan and University of York
//L
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA


#include "ccp4_map_io.h"

#include <ccp4/cmaplib.h>


namespace clipper {


/*! Constructing an CCP4MAPfile does nothing except flag the object as not
  attached to any file for either input or output */
CCP4MAPfile::CCP4MAPfile()
{
  mode = NONE;
}


/*! Close any files which were left open. This is particularly
  important since to access the MAP file efficiently, data reads and
  writes are deferred until the file is closed. */
CCP4MAPfile::~CCP4MAPfile()
{
  switch ( mode ) {
  case READ:
    close_read(); break;
  case WRITE:
    close_write(); break;
  case NONE:
    break;
  }
}


/*! The file is opened for reading. This CCP4MAPfile object will remain
  attached to this file until it is closed. Until that occurs, no
  other file may be opened with this object, however another CCP4MAPfile
  object could be used to access another file.
  \param filename_in The input filename or pathname. */
void CCP4MAPfile::open_read( const String filename_in )
{
  if ( mode != NONE )
    Message::message( Message_fatal( "CCP4MAPfile: open_read - File already open" ) );

  filename = filename_in;
  mode = READ;

  // read the file header
  char symop[81];
  int grid[3], orderfms[3], orderxyz[3], dim[3], gfms0[3], gfms1[3];
  float cp[6];
  CMap_io::CMMFile* file = (CMap_io::CMMFile*)CMap_io::ccp4_cmap_open( filename.c_str(), O_RDONLY );
  if ( file == NULL ) Message::message( Message_fatal( "CCP4MAPfile: open_read - File missing or corrupted: "+filename ) );
  CMap_io::ccp4_cmap_get_cell( file, cp );
  CMap_io::ccp4_cmap_get_grid( file, grid );
  CMap_io::ccp4_cmap_get_order( file, orderfms );
  CMap_io::ccp4_cmap_get_dim( file, dim );
  CMap_io::ccp4_cmap_get_origin( file, gfms0 );

  // get grid bound and axis order
  for ( int i = 0; i < 3; i++ ) gfms1[i] = gfms0[i] + dim[i] - 1;
  for ( int i = 0; i < 3; i++ ) orderxyz[orderfms[i]-1] = i;

  // get spacegroup
  String symops;
  for ( int i = 0; i < CMap_io::ccp4_cmap_num_symop( file ); i++ ) {
    CMap_io::ccp4_cmap_seek_symop( file, i, SEEK_SET );
    CMap_io::ccp4_cmap_get_symop( file, symop );
    symops += String( symop ) + ";";
  }
  CMap_io::ccp4_cmap_close( file );

  // store header info
  spacegroup_ = Spacegroup( Spgr_descr( symops, Spgr_descr::Symops ) );
  cell_ = Cell( Cell_descr( cp[0], cp[1], cp[2], cp[3], cp[4], cp[5] ) );
  grid_sam_ = Grid_sampling( grid[0], grid[1], grid[2] );
  grid_map_ = Grid_range(
    Coord_grid( gfms0[orderxyz[0]], gfms0[orderxyz[1]], gfms0[orderxyz[2]] ),
    Coord_grid( gfms1[orderxyz[0]], gfms1[orderxyz[1]], gfms1[orderxyz[2]] ) );
}


/*! Close the file after reading. */
void CCP4MAPfile::close_read()
{
  if ( mode != READ )
    Message::message( Message_fatal( "CCP4MAPfile: no file open for read" ) );

  mode = NONE;
}


/*! The file is opened for writing. This will be a new file, created
  entirely from data from within the program, rather than by extending
  an existing file. Similar restrictions apply as for open_read().
  \param filename_out The output filename or pathname. */
void CCP4MAPfile::open_write( const String filename_out )
{
  if ( mode != NONE )
    Message::message( Message_fatal( "CCP4MAPfile: open_write - File already open" ) );

  filename = filename_out;
  mode = WRITE;
}


/*! Close the file after writing. */
void CCP4MAPfile::close_write()
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "CCP4MAPfile: no file open for write" ) );

  mode = NONE;
}


/*! When writing an NXmap, the cell for the output map must be set
  using set_cell(). Note that the NXmap is rather more general than
  the CCP4 map, since it can take an arbitrary rotation+skew
  matrix. The resulting map will only be sensible if the NXmap grid
  skew matrix reflects the supplied cell. This is not possible in the
  general case. (Alternatively, opening an equivalent map for read and
  then closing it again will also set the cell).
  \param cell The cell description for the output map file. */
void CCP4MAPfile::set_cell( const Cell& cell )
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "CCP4MAPfile: no file open for write" ) );
  cell_ = cell;
}


/*! Get the spacegroup from the MTZ file. \return The spacegroup. */
const Spacegroup& CCP4MAPfile::spacegroup() const
{ return spacegroup_; }


/*! Get the base cell from the MTZ file. \return The cell. */
const Cell& CCP4MAPfile::cell() const
{ return cell_; }


/*! Get the grid sampling from the MTZ file. \return The grid sampling. */
const Grid_sampling& CCP4MAPfile::grid_sampling() const
{ return grid_sam_; }


/*! Import a complete Xmap object. The supplied Xmap object is
  examined, and if any of the parameters (spacegroup, cell, or
  grid_sampling) are unset, then they will be set using values from the
  file. The data is the imported from the file.

  If the spacegroups mismatch, the resulting map will obey its
  spacegroup symmetry, but no expansion will be performed if the file
  has a higher symmetry and only holds an asymmetric unit.

  \param xmap The Xmap to be imported.
*/
template<class T> void CCP4MAPfile::import_xmap( Xmap<T>& xmap ) const
{
  if ( mode != READ )
    Message::message( Message_fatal( "CCP4MAPfile: no file open for read" ) );

  // first check if the HKL_info params are already set
  Spacegroup    s = xmap.spacegroup();
  Cell          c = xmap.cell();
  Grid_sampling r = xmap.grid_sampling();
  // import any missing params
  if ( s.is_null() ) s = spacegroup_;
  if ( c.is_null() ) c = cell_;
  if ( r.is_null() ) r = grid_sam_;
  xmap.init( s, c, r );

  // read the file header
  int orderfms[3], orderxyz[3], dim[3], gfms0[3], gfms1[3];
  CMap_io::CMMFile* file = (CMap_io::CMMFile*)CMap_io::ccp4_cmap_open( filename.c_str(), O_RDONLY );
  if ( file == NULL ) Message::message( Message_fatal( "CCP4MAPfile: import_xmap - File missing or corrupted: "+filename ) );
  CMap_io::ccp4_cmap_get_order( file, orderfms );
  CMap_io::ccp4_cmap_get_dim( file, dim );
  CMap_io::ccp4_cmap_get_origin( file, gfms0 );
  int dmode = CMap_io::ccp4_cmap_get_datamode( file );
  if ( dmode != 0 && dmode != 2 )
    Message::message( Message_fatal( "CCP4CCP4MAPfile: unsupported data mode" ) );

  // get grid bound and axis order
  for ( int i = 0; i < 3; i++ ) gfms1[i] = gfms0[i] + dim[i] - 1;
  for ( int i = 0; i < 3; i++ ) orderxyz[orderfms[i]-1] = i;

  // read the map data
  int n0 = (gfms1[0]-gfms0[0]+1);
  int n1 = n0 * (gfms1[1]-gfms0[1]+1);
  std::vector<float> section( n1 );
  int index, g[3];
  Xmap_base::Map_reference_coord x( xmap );
  for ( g[2] = gfms0[2]; g[2] <= gfms1[2]; g[2]++ ) {
    index = 0;
    CMap_io::ccp4_cmap_read_section( file, &section[0] );
    if ( dmode == 0 ) for ( int i = n1-1; i >= 0; i-- )  // deal with byte maps
      section[i] = float(((signed char*)&section[0])[i]);
    for ( g[1] = gfms0[1]; g[1] <= gfms1[1]; g[1]++ ) {
      for ( g[0] = gfms0[0]; g[0] <= gfms1[0]; g[0]++ ) {
	x.set_coord( Coord_grid( g[orderxyz[0]], g[orderxyz[1]],
				 g[orderxyz[2]] ) );
        xmap[x] = T( section[ index++ ] );
      }
    }
  }

  // done
  CMap_io::ccp4_cmap_close( file );
}


/*! 
  \param xmap The Xmap to be exported.
*/
template<class T> void CCP4MAPfile::export_xmap( const Xmap<T>& xmap )
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "CCP4MAPfile: no file open for write" ) );

  const char* title = "From clipper Xmap                                                               ";
  char symop[80];
  int grid[3], orderfms[3], orderxyz[3], dim[3], gfms0[3], gfms1[3];
  float cp[6];

  int spg = xmap.spacegroup().descr().spacegroup_number();  // spacegroup

  // get axis order
  switch ( spg ) {
  case 1:  case 2:  case 3:  case 4:
  case 10: case 16: case 17: case 18:
  case 20: case 21: case 23:
    orderfms[0] = 2; orderfms[1] = 1; orderfms[2] = 3; break;
  default:
    orderfms[0] = 3; orderfms[1] = 1; orderfms[2] = 2; break;
  }
  for ( int i = 0; i < 3; i++ ) orderxyz[orderfms[i]-1] = i;

  // grids
  for ( int i = 0; i < 3; i++ ) {
    grid[i] = xmap.grid_sampling()[i];
    gfms0[orderxyz[i]] = xmap.grid_asu().min()[i];
    gfms1[orderxyz[i]] = xmap.grid_asu().max()[i];
  }
  Cell_descr cd = xmap.cell().descr();
  cp[0] = cd.a(); cp[3] = cd.alpha_deg();
  cp[1] = cd.b(); cp[4] = cd.beta_deg ();
  cp[2] = cd.c(); cp[5] = cd.gamma_deg();
  for ( int i = 0; i < 3; i++ ) dim[i] = gfms1[i] - gfms0[i] + 1;

  CMap_io::CMMFile* file = (CMap_io::CMMFile*)CMap_io::ccp4_cmap_open( filename.c_str(), O_WRONLY );
  if ( file == NULL ) Message::message( Message_fatal( "CCP4MAPfile: export_xmap - File missing or corrupted: "+filename ) );
  CMap_io::ccp4_cmap_set_cell( file, cp );
  CMap_io::ccp4_cmap_set_grid( file, grid );
  CMap_io::ccp4_cmap_set_order( file, orderfms );
  CMap_io::ccp4_cmap_set_dim( file, dim );
  CMap_io::ccp4_cmap_set_origin( file, gfms0 );

  CMap_io::ccp4_cmap_set_spacegroup( file, spg );
  CMap_io::ccp4_cmap_set_title( file, title );
  CMap_io::ccp4_cmap_set_datamode( file, 2 );

  // write symops
  for ( int i = 0; i < xmap.spacegroup().num_symops(); i++ ) {
    String strop = xmap.spacegroup().symop(i).format();
    for ( int j = 0; j < 80; j++ ) symop[j] = ' ';
    for ( int j = 0; j < strop.length(); j++ ) symop[j] = strop[j];
    CMap_io::ccp4_cmap_set_symop( file, symop );
  }

  // write the map data
  int n0 = (gfms1[0]-gfms0[0]+1);
  int n1 = n0 * (gfms1[1]-gfms0[1]+1);
  std::vector<float> section( n1 );
  int index, g[3];
  Xmap_base::Map_reference_coord x( xmap );
  for ( g[2] = gfms0[2]; g[2] <= gfms1[2]; g[2]++ ) {
    index = 0;
    for ( g[1] = gfms0[1]; g[1] <= gfms1[1]; g[1]++ ) {
      for ( g[0] = gfms0[0]; g[0] <= gfms1[0]; g[0]++ ) {
	x.set_coord( Coord_grid( g[orderxyz[0]], g[orderxyz[1]],
				 g[orderxyz[2]] ) );
        section[ index++ ] = float( xmap[x] );
      }
    }
    CMap_io::ccp4_cmap_write_section( file, &section[0] );
  }

  // done
  CMap_io::ccp4_cmap_close( file );
}

/*! 
  \param nxmap The NXmap to be imported.
*/
template<class T> void CCP4MAPfile::import_nxmap( NXmap<T>& nxmap ) const
{
  if ( mode != READ )
    Message::message( Message_fatal( "CCP4MAPfile: no file open for read" ) );

  nxmap.init( cell_, grid_sam_, grid_map_ );

  // read the file header
  int orderfms[3], orderxyz[3], dim[3], gfms0[3], gfms1[3];
  CMap_io::CMMFile* file = (CMap_io::CMMFile*)CMap_io::ccp4_cmap_open( filename.c_str(), O_RDONLY );
  if ( file == NULL ) Message::message( Message_fatal( "CCP4MAPfile: import_nxmap - File missing or corrupted: "+filename ) );
  CMap_io::ccp4_cmap_get_order( file, orderfms );
  CMap_io::ccp4_cmap_get_dim( file, dim );
  CMap_io::ccp4_cmap_get_origin( file, gfms0 );
  int dmode = CMap_io::ccp4_cmap_get_datamode( file );
  if ( dmode != 0 && dmode != 2 )
    Message::message( Message_fatal( "CCP4CCP4MAPfile: unsupported data mode" ) );

  // get grid bound and axis order
  for ( int i = 0; i < 3; i++ ) gfms1[i] = gfms0[i] + dim[i] - 1;
  for ( int i = 0; i < 3; i++ ) orderxyz[orderfms[i]-1] = i;

  // read the map data
  int n0 = (gfms1[0]-gfms0[0]+1);
  int n1 = n0 * (gfms1[1]-gfms0[1]+1);
  std::vector<float> section( n1 );
  int index, g[3];
  for ( g[2] = 0; g[2] <= gfms1[2]-gfms0[2]; g[2]++ ) {
    index = 0;
    CMap_io::ccp4_cmap_read_section( file, &section[0] );
    if ( dmode == 0 ) for ( int i = n1-1; i >= 0; i-- )  // deal with byte maps
      section[i] = float(((unsigned char*)&section[0])[i]);
    for ( g[1] = 0; g[1] <= gfms1[1]-gfms0[1]; g[1]++ ) {
      for ( g[0] = 0; g[0] <= gfms1[0]-gfms0[0]; g[0]++ ) {
	nxmap.set_data( Coord_grid( g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]] ), T( section[ index++ ] ) );
      }
    }
  }

  // done
  CMap_io::ccp4_cmap_close( file );
}

/*! The cell for the output map must have been set using
  set_cell_descr(). Note that the NXmap is rather more general than
  the CCP4 map, since it can take an arbitrary rotation+skew
  matrix. The resulting map will only be sensible if the NXmap grid
  skew matrix reflects the supplied cell. This is not possible in the
  general case.
  \param nxmap The NXmap to be exported. */
template<class T> void CCP4MAPfile::export_nxmap( const NXmap<T>& nxmap )
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "CCP4MAPfile: no file open for write" ) );

  const char* title = "From clipper NXmap                                                              ";
  const char* symop = "X, Y, Z                                                                         ";
  int grid[3], orderfms[3], orderxyz[3], dim[3], gfms0[3], gfms1[3];
  float cp[6];

  int spg = 1;                                       // spacegroup
  orderfms[0] = 2; orderfms[1] = 1; orderfms[2] = 3; // axis order
  for ( int i = 0; i < 3; i++ ) orderxyz[orderfms[i]-1] = i;

  /* Because CCP4 maps don't allow an arbitrary skew matrix, we have
     to jump through hoops to try and fit an NXmap into one. This code
     will work for those cases where it is possible and produce
     garbage otherwise. Don't even try to understand this unless you
     are pretty smart. */
  // cell (set by user, or from previous map)
  cp[0] = cell_.descr().a(); cp[3] = cell_.descr().alpha_deg();
  cp[1] = cell_.descr().b(); cp[4] = cell_.descr().beta_deg ();
  cp[2] = cell_.descr().c(); cp[5] = cell_.descr().gamma_deg();
  // grid (calculated to fit with cell provided - assume angles match)
  Coord_frac c0, c1;
  c0 = nxmap.coord_orth( Coord_map(0,0,0) ).coord_frac(cell_);
  c1 = nxmap.coord_orth( Coord_map( nxmap.grid().nu(), nxmap.grid().nv(),
				    nxmap.grid().nw() ) ).coord_frac(cell_);
  grid_sam_ =
    Grid_sampling( Util::intr( ftype(nxmap.grid().nu())/(c1.u()-c0.u()) ),
		   Util::intr( ftype(nxmap.grid().nv())/(c1.v()-c0.v()) ),
		   Util::intr( ftype(nxmap.grid().nw())/(c1.w()-c0.w()) ) );
  Coord_grid g0 = c0.coord_grid(grid_sam_);
  Coord_grid g1 = g0 + Coord_grid(nxmap.grid()) - Coord_grid(1,1,1);
  grid_map_ = Grid_range( g0, g1 );

  for ( int i = 0; i < 3; i++ ) {
    grid[i] = grid_sam_[i];
    gfms0[orderxyz[i]] = grid_map_.min()[i];
    gfms1[orderxyz[i]] = grid_map_.max()[i];
  }
  for ( int i = 0; i < 3; i++ ) dim[i] = gfms1[i] - gfms0[i] + 1;

  CMap_io::CMMFile* file = (CMap_io::CMMFile*)CMap_io::ccp4_cmap_open( filename.c_str(), O_WRONLY );
  if ( file == NULL ) Message::message( Message_fatal( "CCP4MAPfile: export_nxmap - File missing or corrupted: "+filename ) );
  CMap_io::ccp4_cmap_set_cell( file, cp );
  CMap_io::ccp4_cmap_set_grid( file, grid );
  CMap_io::ccp4_cmap_set_order( file, orderfms );
  CMap_io::ccp4_cmap_set_dim( file, dim );
  CMap_io::ccp4_cmap_set_origin( file, gfms0 );

  CMap_io::ccp4_cmap_set_spacegroup( file, spg );
  CMap_io::ccp4_cmap_set_title( file, title );
  CMap_io::ccp4_cmap_set_datamode( file, 2 );

  // write symops
  CMap_io::ccp4_cmap_set_symop( file, symop );

  // write the map data
  int n0 = (gfms1[0]-gfms0[0]+1);
  int n1 = n0 * (gfms1[1]-gfms0[1]+1);
  std::vector<float> section( n1 );
  int index, g[3];
  for ( g[2] = 0; g[2] <= gfms1[2]-gfms0[2]; g[2]++ ) {
    index = 0;
    for ( g[1] = 0; g[1] <= gfms1[1]-gfms0[1]; g[1]++ ) {
      for ( g[0] = 0; g[0] <= gfms1[0]-gfms0[0]; g[0]++ ) {
        section[ index++ ] = float( nxmap.get_data( Coord_grid( g[orderxyz[0]], g[orderxyz[1]], g[orderxyz[2]] ) ) );
      }
    }
    CMap_io::ccp4_cmap_write_section( file, &section[0] );
  }

  // done
  CMap_io::ccp4_cmap_close( file );
}


/*! Import a complete Xmap object. The supplied Xmap object is
  examined, and if any of the parameters (spacegroup, cell, or
  grid_sampling) are unset, then they will be set using values from the
  file. The data is the imported from the file.

  If the spacegroups mismatch, the resulting map will obey its
  spacegroup symmetry, but no expansion will be performed if the file
  has a higher symmetry and only holds an asymmetric unit.

  This version performs a check on the symmatry of the inport map
  file. If the import map contains less than a complete ASU, or if it
  contains more than an ASU and symmetry related positions contain
  inconsistent values, then an error is returned.

  This function is only available for float or double maps.

  \param xmap The Xmap to be imported.
  \return The ASU check flag
*/
template<class T> CCP4MAPfile::ASUerror CCP4MAPfile::import_xmap_check_asu( Xmap<T>& xmap, T missing ) const
{
  if ( mode != READ )
    Message::message( Message_fatal( "CCP4MAPfile: no file open for read" ) );

  // first check if the HKL_info params are already set
  Spacegroup    s = xmap.spacegroup();
  Cell          c = xmap.cell();
  Grid_sampling r = xmap.grid_sampling();
  // import any missing params
  if ( s.is_null() ) s = spacegroup_;
  if ( c.is_null() ) c = cell_;
  if ( r.is_null() ) r = grid_sam_;
  xmap.init( s, c, r );

  // set missing flag
  T flag;
  Util::set_null( flag );
  xmap = flag;
  ftype maxerr = 0.0;

  // read the file header
  int orderfms[3], orderxyz[3], dim[3], gfms0[3], gfms1[3];
  CMap_io::CMMFile* file = (CMap_io::CMMFile*)CMap_io::ccp4_cmap_open( filename.c_str(), O_RDONLY );
  if ( file == NULL ) Message::message( Message_fatal( "CCP4MAPfile: import_xmap - File missing or corrupted: "+filename ) );
  CMap_io::ccp4_cmap_get_order( file, orderfms );
  CMap_io::ccp4_cmap_get_dim( file, dim );
  CMap_io::ccp4_cmap_get_origin( file, gfms0 );
  int dmode = CMap_io::ccp4_cmap_get_datamode( file );
  if ( dmode != 0 && dmode != 2 )
    Message::message( Message_fatal( "CCP4CCP4MAPfile: unsupported data mode" ) );

  // get grid bound and axis order
  for ( int i = 0; i < 3; i++ ) gfms1[i] = gfms0[i] + dim[i] - 1;
  for ( int i = 0; i < 3; i++ ) orderxyz[orderfms[i]-1] = i;

  // read the map data
  int n0 = (gfms1[0]-gfms0[0]+1);
  int n1 = n0 * (gfms1[1]-gfms0[1]+1);
  std::vector<float> section( n1 );
  int index, g[3];
  Xmap_base::Map_reference_coord x( xmap );
  for ( g[2] = gfms0[2]; g[2] <= gfms1[2]; g[2]++ ) {
    index = 0;
    CMap_io::ccp4_cmap_read_section( file, &section[0] );
    if ( dmode == 0 ) for ( int i = n1-1; i >= 0; i-- )  // deal with byte maps
      section[i] = float(((unsigned char*)&section[0])[i]);
    for ( g[1] = gfms0[1]; g[1] <= gfms1[1]; g[1]++ ) {
      for ( g[0] = gfms0[0]; g[0] <= gfms1[0]; g[0]++ ) {
	x.set_coord( Coord_grid( g[orderxyz[0]], g[orderxyz[1]],
				 g[orderxyz[2]] ) );
	T oldval = xmap[x];
	T newval = T( section[ index++ ] );
	if ( !Util::is_nan( oldval ) && !Util::is_nan( newval ) ) {
	  maxerr = Util::max( maxerr, fabs( ftype( newval - oldval ) ) );
          newval = Util::max( oldval, newval );
        }
	xmap[x] = newval;
      }
    }
  }

  // done
  CMap_io::ccp4_cmap_close( file );

  // ASU checks
  typedef Xmap_base::Map_reference_index MRI;
  ASUerror asuerr = ASUCORRECT;
  ftype s0(0.0), s1(0.0), s2(0.0);
  for ( MRI ix = xmap.first(); !ix.last(); ix.next() ) {
    T val = xmap[ix];
    if ( Util::is_null( val ) ) {
      asuerr = ASUINCOMPLETE;
      xmap[ix] = missing;
    } else {
      s0 += 1.0;
      s1 += val;
      s2 += val*val;
    }
  }
  if ( s0 > 0.5 ) s2 = sqrt( s2*s0 - s1*s1 ) / s0;
  if ( maxerr > 0.01*s2 ) asuerr = ASUINCONSISTENT;

  return asuerr;
}


// instantiate the template functions for all reasonable types
// xmap import
template void CCP4MAPfile::import_xmap<unsigned char>( Xmap<unsigned char>& xmap ) const;
template void CCP4MAPfile::import_xmap<char>( Xmap<char>& xmap ) const;
template void CCP4MAPfile::import_xmap<unsigned short>( Xmap<unsigned short>& xmap ) const;
template void CCP4MAPfile::import_xmap<short>( Xmap<short>& xmap ) const;
template void CCP4MAPfile::import_xmap<unsigned int>( Xmap<unsigned int>& xmap ) const;
template void CCP4MAPfile::import_xmap<int>( Xmap<int>& xmap ) const;
template void CCP4MAPfile::import_xmap<ftype32>( Xmap<ftype32>& xmap ) const;
template void CCP4MAPfile::import_xmap<ftype64>( Xmap<ftype64>& xmap ) const;
// xmap export
template void CCP4MAPfile::export_xmap<unsigned char>( const Xmap<unsigned char>& xmap );
template void CCP4MAPfile::export_xmap<char>( const Xmap<char>& xmap );
template void CCP4MAPfile::export_xmap<unsigned short>( const Xmap<unsigned short>& xmap );
template void CCP4MAPfile::export_xmap<short>( const Xmap<short>& xmap );
template void CCP4MAPfile::export_xmap<unsigned int>( const Xmap<unsigned int>& xmap );
template void CCP4MAPfile::export_xmap<int>( const Xmap<int>& xmap );
template void CCP4MAPfile::export_xmap<ftype32>( const Xmap<ftype32>& xmap );
template void CCP4MAPfile::export_xmap<ftype64>( const Xmap<ftype64>& xmap );
// nxmap import
template void CCP4MAPfile::import_nxmap<unsigned char>( NXmap<unsigned char>& nxmap ) const;
template void CCP4MAPfile::import_nxmap<char>( NXmap<char>& nxmap ) const;
template void CCP4MAPfile::import_nxmap<unsigned short>( NXmap<unsigned short>& nxmap ) const;
template void CCP4MAPfile::import_nxmap<short>( NXmap<short>& nxmap ) const;
template void CCP4MAPfile::import_nxmap<unsigned int>( NXmap<unsigned int>& nxmap ) const;
template void CCP4MAPfile::import_nxmap<int>( NXmap<int>& nxmap ) const;
template void CCP4MAPfile::import_nxmap<ftype32>( NXmap<ftype32>& nxmap ) const;
template void CCP4MAPfile::import_nxmap<ftype64>( NXmap<ftype64>& nxmap ) const;
// nxmap export
template void CCP4MAPfile::export_nxmap<unsigned char>( const NXmap<unsigned char>& nxmap );
template void CCP4MAPfile::export_nxmap<char>( const NXmap<char>& nxmap );
template void CCP4MAPfile::export_nxmap<unsigned short>( const NXmap<unsigned short>& nxmap );
template void CCP4MAPfile::export_nxmap<short>( const NXmap<short>& nxmap );
template void CCP4MAPfile::export_nxmap<unsigned int>( const NXmap<unsigned int>& nxmap );
template void CCP4MAPfile::export_nxmap<int>( const NXmap<int>& nxmap );
template void CCP4MAPfile::export_nxmap<ftype32>( const NXmap<ftype32>& nxmap );
template void CCP4MAPfile::export_nxmap<ftype64>( const NXmap<ftype64>& nxmap );
// xmap import and check
template CCP4MAPfile::ASUerror CCP4MAPfile::import_xmap_check_asu<ftype32>( Xmap<ftype32>& xmap, ftype32 missing ) const;
template CCP4MAPfile::ASUerror CCP4MAPfile::import_xmap_check_asu<ftype64>( Xmap<ftype64>& xmap, ftype64 missing ) const;

} // namespace clipper
