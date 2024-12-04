// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York


#ifndef PRIVATEER_XRAY_H_INCLUDED
#define PRIVATEER_XRAY_H_INCLUDED

#include <clipper/clipper.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>
#include <clipper/clipper-contrib.h>


namespace privateer
{
  namespace xray
  {
        void read_xray_map ( clipper::String const pathname, clipper::String const input_model_path, clipper::MiniMol& mmol, clipper::HKL_info& hklinfo, clipper::CCP4MTZfile& mtzin );
        void initialize_experimental_dataset(clipper::CCP4MTZfile& mtzin, clipper::CCP4MTZfile& ampmtzin, clipper::String const input_column_fobs, clipper::HKL_data<clipper::data32::F_sigF>& fobs, clipper::HKL_info& hklinfo, clipper::MTZcrystal& opxtal, clipper::MTZdataset& opdset, clipper::String const input_reflections_mtz);
        std::pair<double, double> calculate_rscc  ( clipper::Xmap<float>& sigmaa_all_map,
                                                clipper::Xmap<float>& sigmaa_omit_fd,
                                                clipper::Xmap<float> &ligandmap, 
                                                clipper::Xmap<float> &mask,
                                                clipper::HKL_info& hklinfo,
                                                clipper::Grid_sampling& mygrid,
                                                clipper::Coord_orth& origin,
                                                clipper::Coord_orth& destination,
                                                bool useSigmaa);
  }
}

#endif