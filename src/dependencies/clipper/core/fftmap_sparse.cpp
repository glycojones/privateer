/* fftmap_sparse.cpp: implementation file for P1 fft map */
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

/* Modifications by Tristan Croll, 2016-2019:
 *
 * - Native Windows compatibility
 * - Parallelisation/thread-safety for most key steps:
 *   - map_uv() and map_kl() made thread-safe for parallel filling of data from
 *     Xmap
 *   - Parallelised FFTs along each axis
 */

#include "fftmap_sparse.h"

#include "hkl_datatypes.h"

#include <config.h>
#include <pocketfft_hdronly.h>

#include <thread>
#include <future>
#include <chrono>



namespace clipper {


FFTmap_base::FFTtype FFTmap_sparse_p1_base::default_type_ = FFTmap_base::Estimate;


/*! Initialise an FFTmap_sparse_p1_base for a grid.
  \param grid_sam The grid sampling of the unit cell.
  \param type Can be FFTmap_sparse_base::Measure, ::Estimate.
  This parameter is a carry-over from when Clipper relied on FFTW, and is no longer used. */
void FFTmap_sparse_p1_base::init( const Grid_sampling& grid_sam, const int num_threads, const FFTtype type )
{
  num_threads_ = num_threads;
  type_ = type;
  if ( type_ == Default ) type_ = default_type();

  // allocate data
  grid_real_ = grid_sam;
  grid_reci_ = Grid( grid_real_.nu(), grid_real_.nv(), grid_real_.nw()/2+1 );

  // make section maps
  std::complex<ffttype>* fillptr = NULL;
  row_kl.resize( grid_reci_.nv(), grid_reci_.nw(), fillptr );
  ffttype* rfillptr = NULL;
  row_uv.resize( grid_real_.nu(), grid_real_.nv(), rfillptr );
}

FFTmap_sparse_p1_base::~FFTmap_sparse_p1_base()
{
  int u, v, w;
  for ( w = 0; w < grid_reci_.nw(); w++ )
    for ( v = 0; v < grid_reci_.nv(); v++ )
      if ( row_kl( v, w ) != NULL ) delete[] row_kl( v, w );
  for ( v = 0; v < grid_real_.nv(); v++ )
    for ( u = 0; u < grid_real_.nu(); u++ )
      if ( row_uv( u, v ) != NULL ) delete[] row_uv( u, v );
}

ffttype* FFTmap_sparse_p1_base::map_uv( const int& u, const int& v )
{
  row_uv.wait_lock(u,v);
  ffttype* ptr = row_uv( u, v );
  if ( ptr == NULL ) {
    ptr = new ffttype[ grid_real_.nw()+2 ]; // padding for in-place FFT
    const ffttype zero( 0.0 );
    for ( int w = 0; w < grid_real_.nw()+2; w++ ) ptr[w] = zero;
    row_uv( u, v ) = ptr;
  }
  row_uv.unlock(u,v);
  return ptr;
}

std::complex<ffttype>* FFTmap_sparse_p1_base::map_kl( const int& k, const int& l )
{
  row_kl.wait_lock(k,l);
  std::complex<ffttype>* ptr = row_kl( k, l );
  if ( ptr == NULL ) {
    ptr = new std::complex<ffttype>[ grid_reci_.nu() ];
    const std::complex<ffttype> zero( 0.0, 0.0 );
    for ( int u = 0; u < grid_reci_.nu(); u++ ) ptr[u] = zero;
    row_kl( k, l ) = ptr;
  }
  row_kl.unlock(k,l);
  return ptr;
}

/*! For later initialisation: see init() */
FFTmap_sparse_p1_hx::FFTmap_sparse_p1_hx()
{}

/*! Construct an FFTmap_sparse_p1_hx for a given grid.
  \param grid_sam The grid sampling of the unit cell.
  \param type Can be FFTmap_sparse_base::Measure, ::Estimate.
  Measure performs slow precalculation (first time only) to get faster FFT. */
FFTmap_sparse_p1_hx::FFTmap_sparse_p1_hx( const Grid_sampling& grid_sam, const int num_threads, const FFTtype type )
{ init( grid_sam, num_threads, type ); }

/*! \fn void FFTmap_sparse_p1_hx::require_real_data( const Coord_grid& uvw )
  The given Coord_grid will be required in the final map. ( uvw must
  be in grid_sampling() )
  \param uvw The coordinate to require. */

/*! \fn const ffttype& FFTmap_sparse_p1_hx::real_data( const Coord_grid& uvw ) const
  ( uvw must be in grid_sampling(), and have been requested )
  \param uvw The coordinate to get.
  \return The real value at that coordinate. */

/*! Friedel opposites are handled correctly
  \param hkl The HKL to set.
  \param f The complex value to set. */
void FFTmap_sparse_p1_hx::set_hkl( const HKL& hkl, const std::complex<ffttype>& f )
{
  Coord_grid c;
  c = Coord_grid(hkl).unit( grid_real_ );
  if ( c.w() < grid_reci_.nw() )
    cplx_data(c) = f;
  c = Coord_grid(-hkl).unit( grid_real_ );
  if ( c.w() < grid_reci_.nw() )
    cplx_data(c) = std::conj(f);
}

/*! The 'require' functions must have been called first to mark the
  required data in the target space. (Source space requirements are
  inferred automatically). */
void FFTmap_sparse_p1_hx::fft_h_to_x( const ftype& scale )
{
    const int nmax =
      Util::max(Util::max(grid_real_.nu(),grid_real_.nv()),grid_real_.nw());

  int u, v, w;
  ffttype s = ffttype( scale );
  // make ul map
  map_l.clear(); map_l.resize(grid_reci_.nw(), false);
  row_u.clear(); row_u.resize(grid_real_.nu(), false);
  for ( w = 0; w < grid_reci_.nw(); w++ )
    for ( v = 0; v < grid_reci_.nv(); v++ )
      if ( row_kl( v, w ) != NULL ) map_l[w] = true;
  for ( v = 0; v < grid_real_.nv(); v++ )
    for ( u = 0; u < grid_real_.nu(); u++ )
      if ( row_uv( u, v ) != NULL ) row_u[u] = true;

  hu_checkpoints = std::unique_ptr<std::atomic_bool[]>(new std::atomic_bool[num_threads_]);
  kv_checkpoints = std::unique_ptr<std::atomic_bool[]>(new std::atomic_bool[num_threads_]);

  std::vector<std::thread> threads;
  for (size_t i=0; i<num_threads_; ++i)
  {
      hu_checkpoints[i].store(false);
      kv_checkpoints[i].store(false);
  }
  for (size_t i=0; i<num_threads_; ++i)
  {
      threads.push_back(
            std::thread(&FFTmap_sparse_p1_hx::thread_kernel_,
                this, i, nmax, s
            )
      );
  }
  for (auto& t: threads)
    t.join();

}

void FFTmap_sparse_p1_hx::thread_kernel_(size_t thread_num, int nmax, ffttype s)
{

    // prep fft
    std::vector<std::complex<ffttype> > in(nmax+1), out(nmax+1);

    int start, end;
    int layers_per_thread = grid_reci_.nw() / num_threads_ + 1;
    start = layers_per_thread * thread_num;
    end = std::min(start+layers_per_thread, grid_reci_.nw());
    transform_along_hu_(start, end);
    hu_checkpoints[thread_num].exchange(true);
    // Wait for other threads
    for (size_t i=0; i<num_threads_; ++i)
        while (!hu_checkpoints[i].load()) { std::this_thread::sleep_for(std::chrono::microseconds(1)); }

    transform_along_kv_(in, out, start, end, s, nmax);
    kv_checkpoints[thread_num].exchange(true);
    // Wait for other threads
    for (size_t i=0; i<num_threads_; ++i)
        while (!kv_checkpoints[i].load()) { std::this_thread::sleep_for(std::chrono::microseconds(1)); }

    layers_per_thread = grid_real_.nv()/num_threads_ + 1;
    start = layers_per_thread * thread_num;
    end = std::min(start+layers_per_thread, grid_real_.nv());
    transform_along_lw_(start, end);

}




void FFTmap_sparse_p1_hx::transform_along_hu_(const int& start, const int& end)
{
    std::complex<ffttype>* ptr;
    pocketfft::shape_t shape = {(size_t)grid_real_.nu()};
    pocketfft::shape_t axes = {0};
    ffttype fct = 1;
    pocketfft::stride_t stride = {sizeof(std::complex<ffttype>)};
    int u,v,w;
    for ( w = start; w < end; w++ ) {
        for ( v = 0; v < grid_reci_.nv(); v++ ) {
            ptr = row_kl( v, w );
            if ( ptr != NULL )
              pocketfft::c2c( shape, stride, stride, axes, pocketfft::FORWARD, ptr, ptr, fct);
        }
    }
}

void FFTmap_sparse_p1_hx::transform_along_kv_(
    std::vector<std::complex<ffttype> >& in,
    std::vector<std::complex<ffttype> >& out,
    const int& start, const int& end, const ffttype& s, const int& nmax)
{
    ffttype zero_real = 0.0;
    std::complex<ffttype> zero( zero_real, zero_real );
    const std::complex<ffttype>* ptr; ffttype* rptr;
    int u,v,w;
    int hw = grid_real_.nw()/2;
    pocketfft::shape_t shape = {(size_t)grid_real_.nv()};
    pocketfft::shape_t axes = {0};
    ffttype fct = 1;
    pocketfft::stride_t stride = {sizeof(std::complex<ffttype>)};

    for ( w = start; w < end; w++ ) if ( map_l[w] ) {
        for ( u = 0; u < grid_real_.nu(); u++ ) if ( row_u[u] ) {
            for ( v = 0; v < grid_real_.nv(); v++ ) {
                ptr = row_kl( v, w );
                if ( ptr != NULL )
                    in[v] = s * ptr[u];
                else
                    in[v] = zero;
            }
            pocketfft::c2c(shape, stride, stride, axes, pocketfft::FORWARD, in.data(), out.data(), fct);
            for ( v = 0; v < grid_real_.nv(); v++ ) {
                rptr = row_uv( u, v );
                if ( rptr != NULL ) {
                    rptr[w * 2] = out[v].real();
                    rptr[w * 2 + 1] = -out[v].imag();
                }
            }
        }
    }
}

void FFTmap_sparse_p1_hx::transform_along_lw_(const int& start, const int& end)
{
    ffttype* rptr;
    std::vector<ffttype> out(grid_real_.nw());
    int u,v,w;
    pocketfft::shape_t shape = {(size_t)grid_real_.nw()};
    pocketfft::shape_t axes = {0};
    ffttype fct = 1;
    pocketfft::stride_t stride_in = {sizeof(std::complex<ffttype>)};
    pocketfft::stride_t stride_out = {sizeof(ffttype)};
    if (true)
    {
        for ( v = start; v < end; v++ ) {
            for ( u = 0; u < grid_real_.nu(); u++ ) {
                rptr = row_uv( u, v );
                if (rptr != NULL) {
                  pocketfft::c2r(shape, stride_in, stride_out, axes, pocketfft::BACKWARD, (std::complex<ffttype>*)rptr, out.data(), fct);
                  for (size_t i=0; i<grid_real_.nw(); ++i)
                    *rptr++ = out[i];
                }
            }
        }
    }
}


/*! For later initialisation: see init() */
FFTmap_sparse_p1_xh::FFTmap_sparse_p1_xh()
{}

/*! Construct an FFTmap_sparse_p1_xh for a given grid.
  \param grid_sam The grid sampling of the unit cell.
  \param type Can be FFTmap_sparse_base::Measure, ::Estimate.
  Measure performs slow precalculation (first time only) to get faster FFT. */
FFTmap_sparse_p1_xh::FFTmap_sparse_p1_xh( const Grid_sampling& grid_sam, const int num_threads, const FFTtype type )
{ init( grid_sam, num_threads, type ); }

/*! Friedel opposites are handled correctly
  \param hkl The HKL required. */
void FFTmap_sparse_p1_xh::require_hkl( const HKL& hkl )
{
  Coord_grid c = Coord_grid(hkl).unit( grid_real_ );
  if ( c.w() < grid_reci_.nw() )
    map_kl( c.v(), c.w() );
  else
    map_kl( ( grid_real_.nv() - c.v() ) % grid_real_.nv(),
        ( grid_real_.nw() - c.w() ) % grid_real_.nw() );
}

/*! Friedel opposites are handled correctly
  \param hkl The required. */
const std::complex<ffttype> FFTmap_sparse_p1_xh::get_hkl( const HKL& hkl ) const
{
  Coord_grid c = Coord_grid(hkl).unit( grid_real_ );
  if ( c.w() < grid_reci_.nw() )
    return row_kl( c.v(), c.w() )[ c.u() ];
  else
    return std::conj( row_kl( ( grid_real_.nv() - c.v() ) % grid_real_.nv(),
                  ( grid_real_.nw() - c.w() ) % grid_real_.nw() )
                    [ ( grid_real_.nu() - c.u() ) % grid_real_.nu() ] );
}

/*! \fn void FFTmap_sparse_p1_xh::require_cplx_data( const Coord_grid& hkl )
  The given Coord_grid will be required in the final reflections.
  ( uvw must be in grid_reci() )
  \param uvw The coordinate to require. */

/*! \fn const std::complex<ffttype>& FFTmap_sparse_p1_xh::cplx_data( const Coord_grid& hkl ) const
  ( hkl must be in grid_reci(), and have been requested )
  \param uvw The coordinate to get.
  \return The complex value at that coordinate. */

/*! \fn ffttype& FFTmap_sparse_p1_xh::real_data( const Coord_grid& uvw )
  ( uvw must be in grid_real() ) */

/*! The 'require' functions must have been called first to mark the
  required data in the target space. (Source space requirements are
  inferred automatically). */
void FFTmap_sparse_p1_xh::fft_x_to_h( const ftype& scale )
{

    ffttype s = ffttype( scale ) / grid_real_.size();

  const int nmax =
    Util::max(Util::max(grid_real_.nu(),grid_real_.nv()),grid_real_.nw());

    // make ul map
    map_l.clear(); map_l.resize(grid_reci_.nw(), false);
    row_u.clear(); row_u.resize(grid_real_.nu(), false);

    int u, v, w;

    for ( w = 0; w < grid_reci_.nw(); w++ )
      for ( v = 0; v < grid_reci_.nv(); v++ )
        if ( row_kl( v, w ) != NULL ) map_l[w] = true;
    for ( v = 0; v < grid_real_.nv(); v++ )
      for ( u = 0; u < grid_real_.nu(); u++ )
        if ( row_uv( u, v ) != NULL ) row_u[u] = true;

  lw_checkpoints = std::unique_ptr<std::atomic_bool[]>(new std::atomic_bool[num_threads_]);
  kv_checkpoints = std::unique_ptr<std::atomic_bool[]>(new std::atomic_bool[num_threads_]);

  std::vector<std::thread> threads;
  for (size_t i=0; i<num_threads_; ++i)
  {
      lw_checkpoints[i].exchange(false);
      kv_checkpoints[i].exchange(false);
  }
  for (size_t i=0; i<num_threads_; ++i)
  {
      threads.push_back(
          std::thread(&FFTmap_sparse_p1_xh::thread_kernel_,
              this, i, nmax, s
          )
      );
  }
  for (auto& t: threads)
    t.join();

}

void FFTmap_sparse_p1_xh::thread_kernel_(size_t thread_num,
    int nmax, ffttype s)
{

    std::vector<std::complex<ffttype> > in(nmax+1), out(nmax+1);

    int start, end;
    int layers_per_thread = grid_real_.nv() / num_threads_ + 1;
    start = layers_per_thread * thread_num;
    end = std::min(start+layers_per_thread, grid_real_.nv());
    transform_along_lw_(start, end, nmax);
    lw_checkpoints[thread_num].exchange(true);
    // Wait for other threads

    for (size_t i=0; i<num_threads_; ++i)
        while (!lw_checkpoints[i].load()) { std::this_thread::sleep_for(std::chrono::microseconds(1)); }

    layers_per_thread = grid_reci_.nw() / num_threads_ + 1;
    start = layers_per_thread * thread_num;
    end = std::min(start+layers_per_thread, grid_reci_.nw());
    transform_along_kv_(in, out, start, end, s, nmax);
    kv_checkpoints[thread_num].exchange(true);
    // Wait for other threads

    for (size_t i=0; i<num_threads_; ++i)
        while (!kv_checkpoints[i].load()) { std::this_thread::sleep_for(std::chrono::microseconds(1)); }

    transform_along_hu_(start, end, nmax);
}

void FFTmap_sparse_p1_xh::transform_along_lw_(
    const int& start, const int& end, const int& nmax)
{
    ffttype* rptr;
    int u,v,w;
    std::vector<std::complex<ffttype> > out(grid_reci_.nw());

    pocketfft::shape_t shape = {(size_t)grid_real_.nw()};
    pocketfft::shape_t axes = {0};
    ffttype fct = 1;
    pocketfft::stride_t stride_in = {sizeof(ffttype)};
    pocketfft::stride_t stride_out = {sizeof(std::complex<ffttype>)};
    for ( v = start; v < end; v++ ) {
        for ( u = 0; u < grid_real_.nu(); u++ ) {
            rptr = row_uv( u, v );
            if ( rptr != NULL )
            {
              pocketfft::r2c(shape, stride_in, stride_out, axes, pocketfft::FORWARD, rptr, out.data(), fct);
              for (const auto& o: out)
              {
                *rptr++ = o.real();
                *rptr++ = o.imag();
              }
            }
                
        }
    }
}

void FFTmap_sparse_p1_xh::transform_along_kv_(
    std::vector<std::complex<ffttype> >& in,
    std::vector<std::complex<ffttype> >& out,
    const int& start, const int& end,
    const ffttype& s, const int& nmax)
{
    std::complex<ffttype>* ptr; const ffttype* rptr;
    int u,v,w;
    int hw = grid_real_.nw()/2;

    ffttype zero_real = 0.0;
    std::complex<ffttype> zero(zero_real, zero_real);

    pocketfft::shape_t shape = {(size_t)grid_real_.nv()};
    pocketfft::shape_t axes = {0};
    ffttype fct = 1;
    pocketfft::stride_t stride = {sizeof(std::complex<ffttype>)};


    for ( w = start; w < end; w++ ) if ( map_l[w] ) {
        for ( u = 0; u < grid_real_.nu(); u++ ) if ( row_u[u] ) {
            for ( v = 0; v < grid_real_.nv(); v++ ) {
                rptr = row_uv( u, v );
                if ( rptr != NULL ) {
                    if (w != hw)
                        in[v] = std::complex<ffttype>(rptr[w * 2], -rptr[w * 2 + 1]);
                    else
                        in[v] = std::complex<ffttype>(rptr[w * 2], zero_real);
                } else {
                    in[v] = zero;
                }
            }
            pocketfft::c2c(shape, stride, stride, axes, pocketfft::BACKWARD, in.data(), out.data(), fct );
            for ( v = 0; v < grid_real_.nv(); v++ ) {
                ptr = row_kl( v, w );
                if ( ptr != NULL ) ptr[u] = s * out[v];
            }
        }
    }

}

void FFTmap_sparse_p1_xh::transform_along_hu_(
    const int& start, const int& end, const int& nmax)
{
    std::complex<ffttype>* ptr;
    std::vector<std::complex<ffttype> > out(nmax+1);
    int u,v,w;

    pocketfft::shape_t shape = {(size_t)grid_real_.nu()};
    pocketfft::shape_t axes = {0};
    ffttype fct = 1;
    pocketfft::stride_t stride = {sizeof(std::complex<ffttype>)};

    for ( w = start; w < end; w++ ) {
        for ( v = 0; v < grid_reci_.nv(); v++ ) {
            ptr = row_kl( v, w );
            if ( ptr != NULL )
              pocketfft::c2c(shape, stride, stride, axes, pocketfft::BACKWARD, ptr, ptr, fct);
        }
    }

}

} // namespace clipper
