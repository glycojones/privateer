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

namespace privateer
{
  namespace cryo_em
  {
    void read_cryoem_map  ( clipper::String const pathname, clipper::MiniMol& mmol, clipper::HKL_info& hklinfo, clipper::Xmap<double>& output_map, clipper::CCP4MAPfile& mrcin, float const resolution_value );

    void calculate_fc_map ( clipper::NXmap<float> const &reference_map,
                            clipper::NXmap<float> &output_map );

    void mask_from_model  ( std::vector <clipper::MMonomer> &input_model,
                                            clipper::NXmap<float> &output_mask,
                                            float radius=2.5 );

    float calculate_rscc  ( clipper::NXmap<float> const &experimental_map,
                                            clipper::NXmap<float> const &fc_map,
                                            clipper::NXmap<float> const &mask,
                                            clipper::Coord_orth origin,
                                            clipper::Coord_orth destination );

    void write_cryoem_map ( clipper::String const pathname, 
                            clipper::NXmap<float> const &input_map );
  }
}

#endif