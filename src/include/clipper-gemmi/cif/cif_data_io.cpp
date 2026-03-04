/* cif_io.cpp: class file for reflection data  cif importer               */
//c Copyright (C) 2000-2006 Paul Emsley, Kevin Cowtan and University of York
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


#include "cif_data_io.h"

extern "C" {
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
}
#include <vector>
#include <fstream>
#include <gemmi/cif.hpp>


namespace clipper {

/*! Get values into array from specified columns.
  \param value array of values to set
  \param datasize number of elements in array
  \param p_item pointer to gemmi::cif::Item
  \param irow row position of data
  \param icol array of column positions of data
  \return The value >0 if str is null, ? or .;0 otherwise. */
int get_real(xtype value[], const int datasize, const gemmi::cif::Item *p_item, const int irow, const int icol[]) {
  int ret_val = 0;
  for (int i = 0; i < datasize; i++) {
    value[i] = 0.0;
    std::string strval = p_item->loop.val(irow, icol[i]);
    if (!gemmi::cif::is_null(strval)) {
      value[i] = clipper::String(strval).f64();
    } else {
      ret_val += 1;
    }
  }
  return ret_val;
}

/*! Get hkl index from cif table.
	\param hkl HKL index to set
	\param p_item pointer to gemmi::cif::Item
	\param irow row position of data
	\param icol array of column positions for H,K,L
	\return The value >0 if str is null, ? or .;0 otherwise. */
int get_hkl(HKL &hkl, const gemmi::cif::Item *p_item, const int irow, const int icol[]){
  int ret_val =0;
  for (int i=0;i<3;i++){
    hkl[i] = 0;
    std::string strval = p_item->loop.val(irow, icol[i]);
    if(!gemmi::cif::is_null(strval)) {
      hkl[i] = clipper::String(strval).i();
    } else {
      ret_val += 1;
    }
  }
  return ret_val;
}

/*! Constructing an CIFfile does nothing except flag the object as not
  attached to any file for either input or output */
CIFfile::CIFfile()
{
   mode = NONE;

   // The
   // problem was that we were trying to import into f_phi_i in
   // close_read() when f_phi_i had not been initialised by calling
   // import_hkl_data for it (because it comes from a function that
   // will calculate the phases from the model itself).  So what do we
   // do?  We will set flags in the CIFfile object that mean that the
   // HKL_data is not to be read into.  That value is null.  When we
   // call import_hkl_data, those values get modified to be pointers
   // to HKL_data.  We will add tests for non-nullness to close_read
   // before we import into those HKL_data.
   f_sigf_i = NULL;
   f_phi_i  = NULL;
   rfree_i  = NULL;
   d_sigd_i     = NULL;
   ABCD_i       = NULL;
   I_sigI_i     = NULL;
   I_sigI_ano_i = NULL;
   f_sigf_ano_i = NULL;
   phi_fom_i    = NULL;
   
   // set these when we find them....
   clipper_cell_set_flag = 0; 
   clipper_reso_set_flag = 0; 
   clipper_symm_set_flag = 0; 
}

/*! Close any files which were left open. This is particularly
  important since to access the CIF file efficiently, data reads and
  writes are deferred until the file is closed. */
CIFfile::~CIFfile()
{
  switch ( mode ) {
  case READ:
     
     // No way do we want to close_read() here! If try to get
     // import_hkl_info() and that fails with a Message fatal, then we
     // pass the cif reading code, import_hkl_data() and close_read()
     // code, because if import_hkl_data() fails those functions are
     // useless.  So, on distruction of the CIFfile, we don't want to
     // close_read() - which may well crash in data_import() or some
     // such.
     // 
     // close_read();
     break;
//   case WRITE:
//     close_write(); break;
  case NONE:
    break;
  }
}


/*! The file is opened for reading. This CIFfile object will remain
  attached to this file until it is closed. Until that occurs, no
  other file may be opened with this object, however another CIFfile
  object could be used to access another file.
  \param filename_in The input filename or pathname.
*/
void CIFfile::open_read( const String filename_in )
{
  if ( mode != NONE )
    Message::message( Message_fatal( "CIFfile: open_read - File already open" ) );

  // open the cif
  f_sigf_i = NULL;
  f_phi_i  = NULL;
  rfree_i  = NULL;
  
  filename = filename_in;

  try{
    doc_.clear();
    doc_ = gemmi::cif::read_file(filename.c_str());
  }
  catch (const std::exception &err) {
    std::string mess = "CIFfile: open_read  - : \n ";
    mess += err.what();
    mess += " \n";
    Message::message( Message_fatal(mess) );
  }

  mode = READ;
  filename = filename_in;

  set_cell_symm_reso(filename);
  if ( clipper_cell_set_flag && !clipper_reso_set_flag ) {
    resolution_ = resolution( cell_ );
    clipper_reso_set_flag = 1;
  }
  if ( clipper_cell_set_flag && clipper_reso_set_flag ) {
    hkl_sampling_ = clipper::HKL_sampling( cell_, resolution_ );
  }
}

/*! Return true if cif file contains phases. */
bool CIFfile::contains_phases_p() const {
  int ret_val = 0;
  if (!doc_.blocks.empty()){
    for(size_t i = 0; i<doc_.blocks.size(); i++) {
      if (doc_.blocks[i].has_tag("_refln.phase_calc")) { 
        ret_val = 1;
        break; // should we check all? or just one block?
      }
    }
  }
  return ret_val;
}


/*! Close the file after reading. This command also actually fills in
  the data in any HKL_data structures which have been marked for
  import.

  Note that we attempt to read in calculated structure factors too.
  These rely on the tags "F_calc" or "F_calc_au" and "phase_calc".  

  http://pdb.rutgers.edu/mmcif/dictionaries/cif_mm.dic/Categories/refln.html

  It is quite usual then for the HKL_data vector of type F_sigF to
  have a different size to the HKL_data for the F_phi's (often this
  will be zero).  I hope that this will not be a problem.

  We test for f_phi_i being non-null before we import data into it -
  CIFfiles can be used with or without reading calculated phases.

  Note to self: how about we make a function of a CIFfile that says
  whether or not it contains phases...

  Note to self: th-is text need to be properly marked up in doxygen format.
*/
void CIFfile::close_read()
{
  if ( mode != READ )
     Message::message( Message_fatal( "CIFfile: no file open for read" ) );

  // make sure the data list is sized
  if ( f_sigf_i     != NULL ) f_sigf_i->update();
  if (  f_phi_i     != NULL )  f_phi_i->update();
  if (  rfree_i     != NULL )  rfree_i->update();
  if ( d_sigd_i     != NULL ) d_sigd_i->update();
  if (   ABCD_i     != NULL )   ABCD_i->update();
  if ( I_sigI_i     != NULL ) I_sigI_i->update();
  if ( I_sigI_ano_i != NULL ) I_sigI_ano_i->update();
  if ( f_sigf_ano_i != NULL ) f_sigf_ano_i->update();

  int ret_val = 0;
  int n_calc_data = 0; 

  // read the data from the CIF file
  // stat cif_file_name.c_str() here, make sure it exists, is readable.

  int ierr;
  int ierr_f;
  int ierr_calc;
  int ierr_rfree_flag;
  int ierr_anom_flag;
  int ierr_intensity_flag;
  int ierr_ABCD_flag;
  
  if (!doc_.blocks.empty()) { 
    for (gemmi::cif::Block &data : doc_.blocks) { //loop through data/blocks
      if (!data.has_mmcif_category("_refln.")) {
        String warn = "CIFfile: close_read  - :\n ";
        warn += data.name + " : cannot find _refln. block.";
        Message::message(Message_warn(warn));
        continue;
      }
      ret_val = 1; // success:
      HKL hkl;
      int hkl_pos[3] = {-1,-1,-1};
      int rfree_pos;
      int FsigF_pos[2], FPhi_pos[2], anomF_pos[4];
      int anomD_pos[2], IsigI_pos[2], anomI_pos[4], ABCD_pos[4];
      xtype x1[2], x_anom[5], x_hl[4];
      const gemmi::cif::Item *p_item = data.find_loop_item("_refln.index_h");

      if (p_item->loop.length() != 0) { 
        // get column positions for hkl and data needed
        if (data.has_tag("_refln.index_h")) {
          hkl_pos[0] = p_item->loop.find_tag("_refln.index_h");
          hkl_pos[1] = p_item->loop.find_tag("_refln.index_k");
          hkl_pos[2] = p_item->loop.find_tag("_refln.index_l");
        }
        if (f_sigf_i) {
          FsigF_pos[0] = -1;
          if (data.has_tag("_refln.F_meas") &&
              data.has_tag("_refln.F_meas_sigma")) {
            FsigF_pos[0] = p_item->loop.find_tag("_refln.F_meas");
            FsigF_pos[1] = p_item->loop.find_tag("_refln.F_meas_sigma");
          } else if (f_sigf_i && data.has_tag("_refln.F_meas_au") &&
                     data.has_tag("_refln.F_meas_sigma_au")) {
            FsigF_pos[0] = p_item->loop.find_tag("_refln.F_meas_au");
            FsigF_pos[1] = p_item->loop.find_tag("_refln.F_meas_sigma_au");
          }
        }
        if (f_phi_i) {
          FPhi_pos[0] = -1;
          if (data.has_tag("_refln.F_calc") &&
              data.has_tag("_refln.phase_calc")) {
            FPhi_pos[0] = p_item->loop.find_tag("_refln.F_calc");
            FPhi_pos[1] = p_item->loop.find_tag("_refln.phase_calc");
          } else if (f_phi_i && data.has_tag("_refln.F_calc_au") &&
                     data.has_tag("_refln.phase_calc")) {
            FPhi_pos[0] = p_item->loop.find_tag("_refln.F_calc_au");
            FPhi_pos[1] = p_item->loop.find_tag("_refln.phase_calc");
          }
        }
        if (rfree_i) {
          rfree_pos = -1;
          if (data.has_tag("_refln.status"))
          rfree_pos = p_item->loop.find_tag("_refln.status");
        }

        if (f_sigf_ano_i) {
          anomF_pos[0] = -1;
          if (data.has_tag("_refln.pdbx_F_plus") &&
              data.has_tag("_refln.pdbx_F_plus_sigma") &&
              data.has_tag("_refln.pdbx_F_minus") &&
              data.has_tag("_refln.pdbx_F_minus_sigma")) {
            anomF_pos[0] = p_item->loop.find_tag("_refln.pdbx_F_plus");
            anomF_pos[1] = p_item->loop.find_tag("_refln.pdbx_F_plus_sigma");
            anomF_pos[2] = p_item->loop.find_tag("_refln.pdbx_F_minus");
            anomF_pos[3] = p_item->loop.find_tag("_refln.pdbx_F_minus_sigma");
          }
        }
        if (d_sigd_i) {
          anomD_pos[0] = -1;
          if (data.has_tag("_refln.pdbx_anom_difference") &&
              data.has_tag("_refln.pdbx_anom_difference_sigma")) {
            anomD_pos[0] = p_item->loop.find_tag("_refln.pdbx_anom_difference");
            anomD_pos[1] =
                p_item->loop.find_tag("_refln.pdbx_anom_difference_sigma");
          }
        }
        if (I_sigI_i) {
          IsigI_pos[0] = -1;
          if (data.has_tag("_refln.intensity_meas") &&
              data.has_tag("_refln.intensity_sigma")) {
            IsigI_pos[0] = p_item->loop.find_tag("_refln.intensity_meas");
            IsigI_pos[1] = p_item->loop.find_tag("_refln.intensity_sigma");
          }
        }
        if (I_sigI_ano_i) {
          anomI_pos[0] = -1;
          if (data.has_tag("_refln.pdbx_I_plus") &&
              data.has_tag("_refln.pdbx_I_plus_sigma") &&
              data.has_tag("_refln.pdbx_I_minus") &&
              data.has_tag("_refln.pdbx_I_minus_sigma")) {
            anomI_pos[0] = p_item->loop.find_tag("_refln.pdbx_I_plus");
            anomI_pos[1] = p_item->loop.find_tag("_refln.pdbx_I_plus_sigma");
            anomI_pos[2] = p_item->loop.find_tag("_refln.pdbx_I_minus");
            anomI_pos[3] = p_item->loop.find_tag("_refln.pdbx_I_minus_sigma");
          }
        }
        if (ABCD_i) {
          ABCD_pos[0] = -1;
          if (data.has_tag("_refln.pdbx_HLA") &&
              data.has_tag("_refln.pdbx_HLB") &&
              data.has_tag("_refln.pdbx_HLC") &&
              data.has_tag("_refln.pdbx_HLD")) {
            ABCD_pos[0] = p_item->loop.find_tag("_refln.pdbx_HLA");
            ABCD_pos[1] = p_item->loop.find_tag("_refln.pdbx_HLB");
            ABCD_pos[2] = p_item->loop.find_tag("_refln.pdbx_HLC");
            ABCD_pos[3] = p_item->loop.find_tag("_refln.pdbx_HLD");
          }
        }
        for (size_t j = 0; j < p_item->loop.length(); j++) {
          // Get HKL indices
          ierr = get_hkl(hkl, p_item, j, hkl_pos);
          if (!ierr) {
            // get F_meas or F_meas_au (arbitrary units)
            if (f_sigf_i != NULL && FsigF_pos[0] != -1) {
              ierr_f =
                  get_real(x1, f_sigf_i->data_size(), p_item, j, FsigF_pos);
              if (x1[0] < -0.9e10 || x1[1] < -0.9e10) ierr_f++;
              if (!ierr_f) {
                if (!f_sigf_i->is_null()) {
                  f_sigf_i->data_import(hkl, x1);
                }
              }
            }
            // get F_calc and Phi
            if (f_phi_i != NULL && FPhi_pos[0] != -1) {
              ierr_calc =
                  get_real(x1, f_phi_i->data_size(), p_item, j, FPhi_pos);
              if (x1[0] < -0.9e10 || x1[1] < -0.9e10) ierr_calc++;
              if (!ierr_calc) {
                if (!f_phi_i->is_null()) {
                  f_phi_i->data_import(hkl, x1);
                  n_calc_data++;
                }
              }
            }
            // get RFree flag
            if (rfree_i != NULL && rfree_pos != -1) {
              char s = p_item->loop.val(j, rfree_pos)[0];
              x1[0] = -1;
              // could be also "x"; not observed, just listed
              if (s == 'o')
                x1[0] = 1;
              else if (s == 'f')
                x1[0] = 0;
              //
              // Every output reflection gets one of these, even if
              // it is set to -1.  Is that the correct thing to do?
              // I have looked at over 200 recent sf mmCIFs and if
              // the contain conventional data then they all have
              // the status flag.
              //
              // However, SFS from EM data do not
              //
              rfree_i->data_import(hkl, x1);
            }
            // Anomalous Fs
            //  (only add if all 4 F+ and sigF+ and F- and sigF- are
            //  present), that is, don't add an F if a sigma is
            //  missing.
            if (f_sigf_ano_i != NULL && anomF_pos[0] != -1) {
              ierr_anom_flag = get_real(x_anom, f_sigf_ano_i->data_size() - 1,
                                        p_item, j, anomF_pos);
              if (!ierr_anom_flag) {
                x_anom[4] = 1.0; // no covarience in cif files, hack a value.
                if (!f_sigf_ano_i->is_null())
                  f_sigf_ano_i->data_import(hkl, x_anom);
              }
            }
            // Anomalous Differences (on F, presumably)
            if (d_sigd_i != NULL && anomD_pos[0] != -1) {
              ierr_anom_flag =
                  get_real(x1, d_sigd_i->data_size(), p_item, j, anomD_pos);
              if (!ierr_anom_flag) {
                if (!d_sigd_i->is_null())
                  d_sigd_i->data_import(hkl, x1);
              }
            }
            // Intensities
            if (I_sigI_i != NULL && IsigI_pos[0] != -1) {
              ierr_intensity_flag =
                  get_real(x1, I_sigI_i->data_size(), p_item, j, IsigI_pos);
              if (!ierr_intensity_flag) {
                if (!I_sigI_i->is_null())
                  I_sigI_i->data_import(hkl, x1);
              }
            }
            // Anomalous Intensities
            if (I_sigI_ano_i != NULL && anomI_pos[0] != -1) {
              ierr_intensity_flag = get_real(
                  x_anom, I_sigI_ano_i->data_size() - 1, p_item, j, anomI_pos);
              if (!ierr_intensity_flag) {
                x_anom[4] = 1.0; // no covariance in cif files, hack a value.
                if (!I_sigI_ano_i->is_null())
                  I_sigI_ano_i->data_import(hkl, x_anom);
              }
            }
            // Hendrickson Lattman coefficients (not many files have these)
            if (ABCD_i != NULL && ABCD_pos[0] != -1) {
              ierr_ABCD_flag =
                  get_real(x_hl, ABCD_i->data_size(), p_item, j, ABCD_pos);
              if (!ierr_ABCD_flag) {
                if (!ABCD_i->is_null())
                  ABCD_i->data_import(hkl, x_hl);
              }
            }
          }
        }
      }
    }
  }
  mode = NONE;
}

/*! Get the spacegroup from the MTZ file. \return The spacegroup. */
const Spacegroup& CIFfile::spacegroup() const
{ return space_group; }

/*! Get the base cell from the MTZ file. \return The cell. */
const Cell& CIFfile::cell() const
{ return cell_; }

/*! Get the resolution limit from the MTZ file. \return The resolution. */
const Resolution& CIFfile::resolution() const
{ return resolution_; }

/*! Get the HKL sampling from the MTZ file. \return The hkl_sampling. */
const HKL_sampling& CIFfile::hkl_sampling() const
{ return hkl_sampling_; }

/*! Get the resolution limit from the CIF file.
  Since a CIF file does not contain cell information, a Cell object
  must be supplied, which will be used to determine the resultion.
  The result is the resolution determined by the most extreme
  reflection in the file.
  \param cell Cell object
  \return The resolution. */
Resolution CIFfile::resolution( const Cell& cell ) const
{
  if ( mode != READ )
    Message::message( Message_fatal( "CIFfile: resolution - no file open for read" ) );

  HKL hkl;
  int h, k, l;
  ftype slim = 0.0;
  if (!doc_.blocks.empty()) { 
    int ierr;
    for(size_t i = 0; i<doc_.blocks.size(); i++) {
      if (doc_.blocks[i].has_tag("_refln.index_h")) { 
        const gemmi::cif::Item *p_item = doc_.blocks[i].find_loop_item("_refln.index_h");
        int hklcol[3] = {0,1,2};
        hklcol[0] = p_item->loop.find_tag("_refln.index_h"); 
	      hklcol[1] = p_item->loop.find_tag("_refln.index_k"); 
	      hklcol[2] = p_item->loop.find_tag("_refln.index_l");
	      for (size_t r = 0; r < p_item->loop.length(); r++ ) {
          ierr = get_hkl(hkl, p_item, r, hklcol);
          if (!ierr) {
		        slim = Util::max( slim, hkl.invresolsq(cell) );
	        } 
        }
	    }
	  }
	}
  return Resolution( 1.0/sqrt(slim) );
}

/*! Import the list of reflection HKLs from an CIF file into an
  HKL_info object. At the start of the routine try to determine the
  space group, cell and resolution. If the resolution limit was found
  and if the resolution limit of the HKL_info object is lower than the
  limit of the file, any excess reflections will be rejected, as will
  any systematic absences or duplicates.  If the resolution is not
  found then the resolution will be determined from the input hkl data
  and \parm target will be init()ed before returning from this function.
  \param target The HKL_info object to be initialised. */
void CIFfile::import_hkl_info( HKL_info& target )
{
  std::vector<HKL> hkls;
  if ( mode != READ )
    Message::message( Message_fatal( "CIFfile: import_hkl_info - no file open for read"+filename ) );

  std::string new_str = filename; 
  int icell = 0;
  if ( clipper_cell_set_flag && clipper_symm_set_flag ){
    icell = 1;
  } else { int icell = set_cell_symm_reso(new_str); }
  
  if (! icell) {
    if (! clipper_cell_set_flag)
	    Message::message( Message_fatal( "CIFfile: import_hkl_info - error getting cell "+filename ) );
     
    if (! clipper_symm_set_flag)
	    Message::message( Message_fatal( "CIFfile: import_hkl_info - error getting symm "+filename ) );
  } else {
    // we have the cell and symmetry, proceed.
    if (! clipper_reso_set_flag)
	    resolution_.init(2.0); // just a dummy value, we will set it to
			// the right value once we know it after
			// we have read in all the reflection
			// hkls.  If clipper_reso_set_flag *was*
			// set, then we have initialised
			// resolution already.
       
    // import any missing params
    target.init( space_group, cell_, resolution_ );

    // read the reflections from the cif
    ftype slim = target.resolution().invresolsq_limit();
    ftype tmp_lim = 0.0; 
  
    if (!doc_.blocks.empty()) { 
      int ierr;
      HKL hkl;
      int hklcol[3] = {-1, -1, -1};
      for (gemmi::cif::Block &data : doc_.blocks) {
        if (data.has_mmcif_category("_refln.")){
          const gemmi::cif::Item *p_item = data.find_loop_item("_refln.index_h");
          hklcol[0] = p_item->loop.find_tag("_refln.index_h");
          hklcol[1] = p_item->loop.find_tag("_refln.index_k");
          hklcol[2] = p_item->loop.find_tag("_refln.index_l");
          for (size_t j=0; j<p_item->loop.length(); j++) {
            ierr = get_hkl(hkl, p_item, j, hklcol);
            if (!ierr) {
              if (clipper_reso_set_flag) {
                if (hkl.invresolsq(target.cell()) < slim ) {
                  hkls.push_back(hkl);
                }
              } else { // resolution had not been set
                if (hkl.invresolsq(target.cell()) > tmp_lim) {
                  tmp_lim = hkl.invresolsq(target.cell());
                }
              }
            }
          }
        }  
      }
    }
    // Now we can initialise target properly, if we OK reading the hkls.
    if (!clipper_reso_set_flag) {
	    if (tmp_lim > 0.0) { // the starting value
	      resolution_.init( 1/sqrt(tmp_lim)); // 210170 // just above
	      target.init( space_group, cell_, resolution_ );
	      std::cout << "Resolution limit set to " << resolution_.limit() << std::endl; 
	    } else {
	      std::cout << "Disaster couldn't set resolution" << std::endl;
	    }
    }
  }
  // Quiet! Silence is the clipper way :)
  // std::cout << "import_hkl_info read " << hkls.size() << " hkls" << std::endl; 
  target.add_hkl_list( hkls );
}

/*! Import data from an CIF file into an HKL_data object.

  This routine does not actually read any data, but rather marks the
  data to be read when the file is closed.

  The data to be read (F_sigF or Phi_fom) will be selected based on
  the type of the HKL_data object.

  \param cdata The HKL_data object into which data is to be imported. */
void CIFfile::import_hkl_data( HKL_data_base& cdata )
{
  if ( mode != READ )
    Message::message( Message_fatal( "CIFfile: import_hkl_data - no file open for read" ) );

  if ( cdata.is_null() ) cdata.init( space_group, cell_, hkl_sampling_ );

  if      ( cdata.type() == data32::F_sigF::type() )        f_sigf_i = &cdata;
  else if ( cdata.type() == data32::F_phi::type()  )         f_phi_i = &cdata;
  else if ( cdata.type() == data32::I_sigI::type() )        I_sigI_i = &cdata;
  else if ( cdata.type() == data32::Flag::type()   )         rfree_i = &cdata;
  else if ( cdata.type() == data32::Phi_fom::type())       phi_fom_i = &cdata;
  else if ( cdata.type() == data32::F_sigF_ano::type()) f_sigf_ano_i = &cdata;
  else if ( cdata.type() == data32::I_sigI_ano::type()) I_sigI_ano_i = &cdata;
  else if ( cdata.type() == data32::ABCD::type())             ABCD_i = &cdata;
  else if ( cdata.type() == data32::D_sigD::type())         d_sigd_i = &cdata;
  else {
     clipper::String m = "CIFfile: import_hkl_data error";
     m += " - data must be F_sigF, F_phi or Flag";
     Message::message( Message_fatal( m ) );
  }
}

/*! Set cell, symmetry and resolution if found in CIF file
    Return non-zero if we find a cell and symmetry and can set them.
    \param cif_file_name The input filename or path 
    \return non-zero if cell and symmetry are found and set. */ 
int CIFfile::set_cell_symm_reso_by_cif(std::string cif_file_name)
{
  if (!doc_.blocks.empty()) {   
    for(gemmi::cif::Block &data : doc_.blocks) {
      // get cell paramaters
      if (data.has_mmcif_category("_cell.")) {
         const std::string *a = data.find_value("_cell.length_a");
         const std::string *b = data.find_value("_cell.length_b");
         const std::string *c = data.find_value("_cell.length_c");
         const std::string *alpha = data.find_value("_cell.angle_alpha");
         const std::string *beta = data.find_value("_cell.angle_beta");
         const std::string *gamma = data.find_value("_cell.angle_gamma");
         if (!gemmi::cif::is_null(*a) && !gemmi::cif::is_null(*b) && !gemmi::cif::is_null(*c) &&
             !gemmi::cif::is_null(*alpha) && !gemmi::cif::is_null(*beta) && !gemmi::cif::is_null(*gamma)) {
            // set clipper cell
            clipper_cell_set_flag = 1;
            cell_.init(clipper::Cell_descr(String(*a).f64(), String(*b).f64(), String(*c).f64(), String(*alpha).f64(),
                                           String(*beta).f64(), String(*gamma).f64()));
            // std::cout << "got cell from cif: "
            // << cell_.format() << std::endl;
         }
      }

	    // Try reading symmetry construction
      if (data.has_tag("_symmetry.space_group_name_H-M")) {
         const std::string *spgHM = data.find_value("_symmetry.space_group_name_H-M");
         if (!gemmi::cif::is_null(*spgHM)) {
            space_group.init(clipper::Spgr_descr(*spgHM));
            clipper_symm_set_flag = 1;
            // 			std::cout << "INFO space_group from symmetry in cif: "
            // 				  << space_group.descr().symbol_hm() << std::endl;
         }
      }

	    // Have another try for symmetry: e.g. from shelx cif files:
	    // (such files tell us the symmetry operators, so we can get
	    // the space group from those rather than the name using a
	    // clipper function).
      if (data.has_tag("_symmetry_equiv.pos_as_xyz")) {
        std::string symmetry_ops("");
        for (std::string &ops : data.find_loop("_symmetry_equiv.pos_as_xyz")) {
          symmetry_ops += ops + " ; ";
        }
        if (symmetry_ops != "") {
          clipper_symm_set_flag = 1;
          space_group.init(clipper::Spgr_descr(symmetry_ops));
        }
      }
	 
	    // Reflection meta data:
	    if ( data.has_tag("_reflns.d_resolution_high") ) {
	      std::string dres_high = *data.find_value("_refln.d_resolution_high");
	      if (!gemmi::cif::is_null(dres_high)) { 
	        xtype reso = clipper::String(dres_high).f64();
		      clipper_reso_set_flag = 1;
		      resolution_.init(reso);
          // std::cout << "got resolution from cif: "
          //    << resolution_.limit() << std::endl;
		    }
	    }
	  }
	}
  return clipper_symm_set_flag && clipper_cell_set_flag; 
}

/*! Return non-zero if we find a cell and symmetry and can set them.
  \param cif_file_name Cif file name
  \return non-zero if cell and symmetry are found and set. */
int CIFfile::set_cell_symm_reso(std::string cif_file_name) {

   if (set_cell_symm_reso_by_cif(cif_file_name)) {
      return 1;
   } else {
      return set_cell_symm_reso_by_kludge(cif_file_name);
   }
}

/*! Find cell and symmetry through the lines in the Cif file.
  \param cif_file_name The input filename or path
  \return non-zero if the cell and symmetry are found and set
*/
int CIFfile::set_cell_symm_reso_by_kludge(std::string cif_file_name) { 
   
   // Cif files from the RCSB do not contain cell and symmetry
   // information in cif format (now how totally crap is that? - but I
   // won't go into that now).
   // 
   // However, they do provide us with comment lines that look like a
   // conventional pdb file.  e.g:
   // 
   // #CRYST1  114.300  114.300  155.200  90.00  90.00 120.00 P 61 2 2     12
   // #REMARK 290 SYMMETRY OPERATORS FOR SPACE GROUP: P 61                   
   // 
   // So the plan is the open the file conventionally and grab the
   // cell and symmetry then after we have got to the "#CRYST1" field
   // or find "loop_", we close up.  Notice how we don't grab the
   // space group symbol, since this is fraught with parsing
   // difficulties, instead we grab the symmetry operators and use
   // clipper cleverness to determine the space group.
   // 
   // Return the success status. True is success, 0 failure.
   //
   //
      
   // std::string line; 
   char word[800];
   short int read_symm_flag = 0;   // 2 is good for read
   short int cell_coming_flag = 0; // 1 is good for read
   short int local_found_symm = 0; 
   short int local_found_cell = 0; 
   //short int local_found_reso = 0; 
   std::string clipper_symm_string(""); 
   std::vector<double> cell_bits; 

   std::ifstream from;
   float reso; 
   short int read_reso_flag = 0; 

   from.open(cif_file_name.c_str()); 

   while (from >> &word[0]) { 
 
      if (read_reso_flag == 1) { 

	 // We can do reso immediately.  The others (cell and
	 // symmetry) we set flags and do them at the end.
	 // 
	 char **endptr = new (char *); 
	 reso = strtod(word, endptr); 
	 if (endptr != NULL &&  *endptr != word) { 
	    resolution_.init(reso); 
	    std::cout << " Found reso: " << resolution_.limit() << std::endl; 
	    clipper_reso_set_flag = 1; 
	 }
	 delete endptr;
	 endptr = 0;
	 read_reso_flag = 0; 
      } 

      if (! strncmp(word, "RESOLUTION.", 11)) { 
	 read_reso_flag = 1; 
      }

      if (! strncmp(word, "NNNMMM",6)) { 
	 read_symm_flag = 1; 
      }

      if (! strncmp(word, "OPERATOR",6)) { 
	 read_symm_flag++;
      }

      if (! strncmp(word, "WHERE",5)) { 
	 read_symm_flag =0;
      }

      if (read_symm_flag == 2) { 
	 
	 if (strchr(word,'X') && 
	     strchr(word,'Y') && 
	     strchr(word,'Z') ) { 

	    clipper_symm_string += word;
	    clipper_symm_string += ";";
	    local_found_symm = 1; 
	 }
      }

      // position dependent code (should be before CRYST1 comparison)
      if (cell_coming_flag == 1) { 
	 cell_bits.push_back(atof(word));
	 local_found_cell = 1; 
      } 

      if (! strncmp(word, "#CRYST1", 7)) {
	 cell_coming_flag = 1; 
      }

      if (cell_bits.size() == 6)
	 break; 

      if (! strncmp(word, "_refln.index_h",14)) break; 
   }
      
   from.close(); 

   if (clipper_symm_set_flag == 0) {
      if (local_found_symm) {
	 space_group.init(clipper::Spgr_descr(clipper_symm_string));
	 clipper_symm_set_flag = 1; 
	 std::cout << " Symm: " << space_group.descr().symbol_hm() << std::endl;
      }
   }
   
   if (clipper_cell_set_flag == 0) {
      if (local_found_cell) {
	 if (cell_bits.size() == 6 && clipper_symm_string != "") { 
	    if ( !clipper::Util::is_nan(cell_bits[0]) && !clipper::Util::is_nan(cell_bits[1]) && 
	         !clipper::Util::is_nan(cell_bits[2]) && !clipper::Util::is_nan(cell_bits[3]) && 
		 !clipper::Util::is_nan(cell_bits[4]) && !clipper::Util::is_nan(cell_bits[5]) ) { 
								      
	       cell_.init(clipper::Cell_descr(cell_bits[0], cell_bits[1],
					      cell_bits[2], cell_bits[3],
					      cell_bits[4], cell_bits[5]));
	       std::cout << cell_.format() << std::endl;
	       clipper_cell_set_flag = 1;
	    }
	 }
      }
   }

   return clipper_cell_set_flag && clipper_symm_set_flag; 
}
 

} // namespace clipper
