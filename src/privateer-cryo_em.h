/*! \file privateer-cryo_em.h
  Header file for handling cryo-EM maps and models */

// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2019 Jon Agirre
// York Structural Biology Laboratory
// Department of Chemistry
// University of York
// mailto: jon.agirre@york.ac.uk
//
// Work funded by The Royal Society
// University Research Fellowship
// award UF160039

#ifndef PRIVATEER_CRYO_EM_H_INCLUDED
#define PRIVATEER_CRYO_EM_H_INCLUDED

#include <clipper/clipper.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>
#include <clipper/clipper-contrib.h>

typedef clipper::HKL_data_base::HKL_reference_index HRI;

namespace privateer
{
  namespace cryo_em
  {
    void read_cryoem_map  ( clipper::String const pathname, clipper::HKL_info& hklinfo, clipper::Xmap<double>& output_map, clipper::CCP4MAPfile& mrcin, float const resolution_value );

    void initialize_dummy_fobs(clipper::HKL_data<clipper::data32::F_sigF>& fobs, clipper::HKL_data<clipper::data32::F_phi>& fc_cryoem_obs);

    void calculate_sfcs_of_fc_maps ( clipper::HKL_data<clipper::data32::F_phi>& fc_all_cryoem_data, clipper::HKL_data<clipper::data32::F_phi>& fc_ligands_only_cryoem_data, clipper::Atom_list& allAtoms, clipper::Atom_list& ligandAtoms);

    bool generate_output_map_coefficients (clipper::HKL_data<clipper::data32::F_phi>& difference_coefficients, clipper::HKL_data<clipper::data32::F_phi>& fc_cryoem_obs, clipper::HKL_data<clipper::data32::F_phi>& fc_all_cryoem_data, clipper::HKL_info& hklinfo);
    

    std::pair<double, double> calculate_rscc  ( clipper::Xmap<double> &experimental_map,
                                                clipper::Xmap<double> &fc_map, 
                                                clipper::Xmap<double> &mask,
                                                clipper::HKL_info& hklinfo,
                                                clipper::Grid_sampling& mygrid,
                                                clipper::Coord_orth& origin,
                                                clipper::Coord_orth& destination );

    void write_cryoem_map ( clipper::String const pathname, 
                            clipper::Xmap<float> const &input_map );
  }
}

#endif