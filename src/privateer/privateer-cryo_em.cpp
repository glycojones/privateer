
/*! \file privateer-cryo_em.cpp
  Contains code for handling cryo-EM maps and models */

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

#include "privateer-cryo_em.h"

void privateer::cryo_em::read_cryoem_map  ( clipper::String const pathname, clipper::MiniMol& mmol, clipper::HKL_info& hklinfo, clipper::Xmap<double>& output_map, clipper::CCP4MAPfile& mrcin, float const resolution_value )
{
    std::cout << "Reading " << pathname.trim().c_str() << "... ";
    fflush(0);
    try
    {
        mrcin.open_read( pathname.trim() );
        mrcin.import_xmap( output_map );
        mrcin.close_read();

        std::cout << "done." << std::endl;

        clipper::Resolution resolution(resolution_value);
        hklinfo = clipper::HKL_info(output_map.spacegroup(), output_map.cell(), resolution, true);
    }
    catch (...)
    {
        throw; // hand control to the caller so we can abort
    }
}

void privateer::cryo_em::initialize_dummy_fobs(clipper::HKL_data<clipper::data32::F_sigF>& fobs, clipper::HKL_data<clipper::data32::F_phi>& fc_cryoem_obs)
{
    int iteration = 0;
    for (HRI ih = fc_cryoem_obs.first(); !ih.last(); ih.next() ) // we want to use all available reflections
        {
            if(!fc_cryoem_obs[ih].missing())
            {
                fobs[ih].f() = fc_cryoem_obs[ih].f();
                // fobs[ih].sigf() = fobs[ih].f();
                fobs[ih].sigf() = 1;
                iteration++;
                if(iteration == 420 || iteration == 6743) std::cout << "fobs[" << ih.index() << "].f() = " << fobs[ih].f() << "\tfobs[" << ih.index() << "].sigf() = " << fobs[ih].sigf() << std::endl;
                // std::cout << "fobs[" << ih.index() << "].f() = " << fobs[ih].f() << "\tfobs[" << ih.index() << "].sigf() = " << fobs[ih].sigf() << std::endl;
            }
        }
}

void privateer::cryo_em::mask_from_model  ( std::vector <clipper::MMonomer> &input_model,
                                            clipper::NXmap<float> &output_mask,
                                            float radius )
{

}

void privateer::cryo_em::calculate_fc_map ( clipper::NXmap<float> const &reference_map, clipper::NXmap<float> &output_map )
{

}

float privateer::cryo_em::calculate_rscc  ( clipper::NXmap<float> const &experimental_map,
                                            clipper::NXmap<float> const &fc_map,
                                            clipper::NXmap<float> const &mask,
                                            clipper::Coord_orth origin,
                                            clipper::Coord_orth destination )
{
  clipper::HKL_info hklinfo = clipper::HKL_info(); //placeholder
  clipper::Grid_sampling mygrid = clipper::Grid_sampling(); //placeholder

  double mean_rho_exp, mean_rho_calc, num, den1, den2, corr_coeff;
  mean_rho_calc = mean_rho_exp = num = den1 = den2 = corr_coeff = 0.0;

  int n_points = 0;
  clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;

  double accum = 0.0;

  // calculation of the mean densities of the two input maps

  // i0 = clipper::Xmap_base::Map_reference_coord( experimental_map, origin.coord_frac(hklinfo.cell()).coord_grid(mygrid) );

  // for ( iu = i0; iu.coord().u() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).u(); iu.next_u() )
  //     for ( iv = iu; iv.coord().v() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).v(); iv.next_v() )
  //         for ( iw = iv; iw.coord().w() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).w(); iw.next_w() )
  //         {
  //             if ( mask[iw] == 1.0)
  //             {
  //                 mean_rho_calc = mean_rho_calc + fc_map[iw];
  //                 mean_rho_exp = mean_rho_exp + experimental_map[iw];
  //                 n_points++;
  //             }
  //         }

  // accum = mean_rho_exp / ms.std_dev();
  // accum /= n_points;

  // mean_rho_calc = mean_rho_calc / n_points;
  // mean_rho_exp = mean_rho_exp / n_points;

  // i0 = clipper::Xmap_base::Map_reference_coord( experimental_map, origin.coord_frac(hklinfo.cell()).coord_grid(mygrid) );

  // for ( iu = i0; iu.coord().u() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).u(); iu.next_u() )
  //   for ( iv = iu; iv.coord().v() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).v(); iv.next_v() )
  //       for ( iw = iv; iw.coord().w() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).w(); iw.next_w() )
  //       {
  //           if ( mask[iw] == 1.0)
  //           {
  //             num = num + (experimental_map[iw] - mean_rho_exp) * (fc_map[iw] - mean_rho_calc);
  //             den1 = den1 + pow((experimental_map[iw] - mean_rho_exp),2);
  //             den2 = den2 + pow((fc_map[iw] - mean_rho_calc),2);
  //           }
  //       }

  // corr_coeff = num / (sqrt(den1) * sqrt(den2));

  return corr_coeff;
}

void privateer::cryo_em::write_cryoem_map ( clipper::String const pathname, clipper::NXmap<float> const &input_map )
{
  clipper::CCP4MAPfile file;
  try
  {
    file.open_write( pathname );
    file.export_nxmap( input_map );
    file.close_write();
  }
  catch (...)
  {
    throw; // hand control to the caller so we can abort
  }
}