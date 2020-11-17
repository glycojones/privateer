
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

void privateer::cryo_em::read_cryoem_map  ( clipper::String const pathname, clipper::HKL_info& hklinfo, clipper::Xmap<double>& output_map, clipper::CCP4MAPfile& mrcin, float const resolution_value )
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
                // fobs[ih].sigf() = fc_cryoem_obs[ih].f();
                // fobs[ih].sigf() = 0;
                fobs[ih].sigf() = 1;
                iteration++;
                if(iteration == 420 || iteration == 6743) std::cout << "fobs[" << ih.index() << "].f() = " << fobs[ih].f() << "\tfobs[" << ih.index() << "].sigf() = " << fobs[ih].sigf() << std::endl;
                // std::cout << "fobs[" << ih.index() << "].f() = " << fobs[ih].f() << "\tfobs[" << ih.index() << "].sigf() = " << fobs[ih].sigf() << std::endl;
            }
        }
}


void privateer::cryo_em::calculate_sfcs_of_fc_maps ( clipper::HKL_data<clipper::data32::F_phi>& fc_all_cryoem_data, clipper::HKL_data<clipper::data32::F_phi>& fc_ligands_only_cryoem_data, clipper::Atom_list& allAtoms, clipper::Atom_list& ligandAtoms) //calculated ligandmap here, lacks atom list of ligands. Replace reference_map with direct object of fc_ligands_bsc
{
  clipper::SFcalc_iso_fft<float> sfcligands;
  clipper::SFcalc_iso_fft<float> sfcall;

  try
    {
    #pragma omp parallel sections
      {
      #pragma omp section
                  sfcligands( fc_ligands_only_cryoem_data,  ligandAtoms );
      #pragma omp section
                  sfcall( fc_all_cryoem_data, allAtoms );
      }
    }
    catch ( ... ) 
    {
      std::cout << "\nThe input file has unrecognised atoms. Might cause unexpected results...\n";  // this causes clipper to freak out, so better remove those unknowns
    }
}

bool privateer::cryo_em::generate_output_map_coefficients (clipper::HKL_data<clipper::data32::F_phi>& difference_coefficients, clipper::HKL_data<clipper::data32::F_phi>& fc_cryoem_obs, clipper::HKL_data<clipper::data32::F_phi>& fc_all_cryoem_data, clipper::HKL_info& hklinfo)
{
  clipper::data32::F_phi fo, twofo, fc, fzero(0.0, 0.0);

  const int n_scl_param = 20;

  // scale and difference data
  std::vector<clipper::ftype> params( n_scl_param, 1.0 );
  clipper::BasisFn_spline basis( hklinfo, n_scl_param, 2.0 );
  clipper::TargetFn_scaleF1F2<clipper::data32::F_phi, clipper::data32::F_phi> scaling_target( fc_all_cryoem_data, fc_cryoem_obs );
  clipper::ResolutionFn scalefn( hklinfo, basis, scaling_target, params );
  
  for (HRI ih = hklinfo.first(); !ih.last(); ih.next() )
  {
    clipper::data32::F_phi temp(fc_all_cryoem_data[ih].f(), fc_all_cryoem_data[ih].phi());
    temp.scale( sqrt( scalefn.f(ih) ) ); // scale

    fo = clipper::data32::F_phi( fc_cryoem_obs[ih].f(), fc_cryoem_obs[ih].phi() );
    fc = clipper::data32::F_phi( temp.f(), temp.phi() );

    if ( !fc_all_cryoem_data[ih].missing() ) 
    {
      if ( !fc_cryoem_obs[ih].missing() )
      {
        difference_coefficients[ih] = fo - fc;
      }
      else
      {
        difference_coefficients[ih] = fzero;
      }
    }
    else
    {
      if ( !fc_cryoem_obs[ih].missing() )
      {
        difference_coefficients[ih] = fzero;
      }
      else 
      {
        difference_coefficients[ih] = fzero;
      }
    }
  }

  return true;
}

std::pair<double, double> privateer::cryo_em::calculate_rscc  ( clipper::Xmap<double>& experimental_map,
                                            clipper::Xmap<double>& fc_map, // equivalent to lignadmap in xray implementation
                                            clipper::Xmap<double>& mask,
                                            clipper::HKL_info& hklinfo,
                                            clipper::Grid_sampling& mygrid,
                                            clipper::Coord_orth& origin,
                                            clipper::Coord_orth& destination )
{
  clipper::Map_stats ms;


  ms = clipper::Map_stats(experimental_map);

  double mean_rho_exp, mean_rho_calc, num, den1, den2, corr_coeff;
  mean_rho_calc = mean_rho_exp = num = den1 = den2 = corr_coeff = 0.0;

  int n_points = 0;
  clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;

  double accum = 0.0;

  // calculation of the mean densities of the two input maps

  i0 = clipper::Xmap_base::Map_reference_coord( experimental_map, origin.coord_frac(hklinfo.cell()).coord_grid(mygrid) );

  for ( iu = i0; iu.coord().u() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).u(); iu.next_u() )
      for ( iv = iu; iv.coord().v() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).v(); iv.next_v() )
          for ( iw = iv; iw.coord().w() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).w(); iw.next_w() )
          {
              if ( mask[iw] == 1.0)
              {
                  mean_rho_calc = mean_rho_calc + fc_map[iw];
                  mean_rho_exp = mean_rho_exp + experimental_map[iw];
                  n_points++;
              }
          }


  accum = mean_rho_exp / ms.std_dev();
  accum /= n_points;

  mean_rho_calc = mean_rho_calc / n_points;
  mean_rho_exp = mean_rho_exp / n_points;


  i0 = clipper::Xmap_base::Map_reference_coord( experimental_map, origin.coord_frac(hklinfo.cell()).coord_grid(mygrid) );

  for ( iu = i0; iu.coord().u() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).u(); iu.next_u() )
    for ( iv = iu; iv.coord().v() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).v(); iv.next_v() )
        for ( iw = iv; iw.coord().w() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).w(); iw.next_w() )
        {
            if ( mask[iw] == 1.0)
            {

              num = num + (experimental_map[iw] - mean_rho_exp) * (fc_map[iw] - mean_rho_calc);
              den1 = den1 + pow((experimental_map[iw] - mean_rho_exp),2);
              den2 = den2 + pow((fc_map[iw] - mean_rho_calc),2);
            }
        }

  corr_coeff = num / (sqrt(den1) * sqrt(den2));


  return std::make_pair(corr_coeff, accum);
}

  void privateer::cryo_em::write_cryoem_map ( clipper::String const pathname, clipper::Xmap<float> const &input_map )
  {
    clipper::CCP4MAPfile file;
    try
    {
      file.open_write( pathname );
      file.export_xmap( input_map );
      file.close_write();
    }
    catch (...)
    {
      throw; // hand control to the caller so we can abort
    }
}