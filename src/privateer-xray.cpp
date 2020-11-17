
/*! \file privateer-xray.cpp
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

#include "privateer-xray.h"


void privateer::xray::read_xray_map ( clipper::String const pathname, clipper::String const input_model_path, clipper::MiniMol& mmol, clipper::HKL_info& hklinfo, clipper::CCP4MTZfile& mtzin)
{
    std::cout << "Reading " << pathname.trim().c_str() << "... ";
    fflush(0);
    
    mtzin.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
    mtzin.open_read( pathname.trim() );


    std::cout << "done." << std::endl;

    try // we could be in trouble should the MTZ file have no cell parameters
    {
        mtzin.import_hkl_info( hklinfo );     // read spacegroup, cell, resolution, HKL's
    }
    catch (...)
    {
        std::cout << "\nReading cell and spacegroup parameters from the CRYST1 card in ";
        std::cout << input_model_path << ":\n Spacegroup (" << mmol.spacegroup().spacegroup_number() << ")\n" << mmol.cell().format() << "\n\n" ;

        clipper::Resolution myRes(0.96);
        hklinfo = clipper::HKL_info( mmol.spacegroup(), mmol.cell(), myRes, true);
    }
}

void privateer::xray::initialize_experimental_dataset(clipper::CCP4MTZfile& mtzin, clipper::CCP4MTZfile& ampmtzin, clipper::String const input_column_fobs, clipper::HKL_data<clipper::data32::F_sigF>& fobs, clipper::HKL_info& hklinfo, clipper::MTZcrystal& opxtal, clipper::MTZdataset& opdset, clipper::String const input_reflections_mtz)
{
    bool notFound = true;
    // initialize_experimental_dataset 
    std::vector<clipper::String> mtzColumns;
    mtzColumns = mtzin.column_labels();

    if (input_column_fobs != "NONE")
    {
        std::cout << "MTZ file supplied. Using " << input_column_fobs << "...\n";
        mtzin.import_hkl_data( fobs, "*/*/["+ input_column_fobs+"]" );
        mtzin.import_crystal(opxtal, input_column_fobs);
        mtzin.import_dataset(opdset, input_column_fobs);
        mtzin.close_read();
        notFound = false;
    }
    else
    {
        for ( int i = 0 ; ((i < mtzColumns.size()) && notFound) ; i++ )
        {
            if (mtzColumns[i].find("/FOBS ") != -1)
            {
                std::cout << "\nMTZ file supplied. Using FOBS & SIGFOBS...\n";
                mtzin.import_hkl_data( fobs, "*/*/[FOBS,SIGFOBS]" );
                mtzin.import_crystal(opxtal, "*/*/[FOBS,SIGFOBS]" );
                mtzin.import_dataset(opdset, "*/*/[FOBS,SIGFOBS]" );
                mtzin.close_read();
                notFound = false;
            }
            else if (mtzColumns[i].find("/FP ") != -1)
            {
                std::cout << "\nMTZ file supplied. Using FP & SIGFP...\n";
                mtzin.import_hkl_data( fobs, "*/*/[FP,SIGFP]" );
                mtzin.import_crystal(opxtal, "*/*/[FP,SIGFP]" );
                mtzin.import_dataset(opdset, "*/*/[FP,SIGFP]" );
                mtzin.close_read();
                notFound = false;
            }
            else if (mtzColumns[i].find("/FOSC ") != -1)
            {
                std::cout << "\nMTZ file supplied. Using FOSC & SIGFOSC...\n";
                mtzin.import_hkl_data( fobs, "*/*/[FOSC,SIGFOSC]" );
                mtzin.import_crystal(opxtal, "*/*/[FOSC,SIGFOSC]" );
                mtzin.import_dataset(opdset, "*/*/[FOSC,SIGFOSC]" );
                mtzin.close_read();
                notFound = false;
            }
            else if (mtzColumns[i].find("/F-obs ") != -1)
            {
                std::cout << "\nMTZ file supplied. Using F-obs & SIGF-obs...\n";
                mtzin.import_hkl_data( fobs, "*/*/[F-obs,SIGF-obs]" );
                mtzin.import_crystal(opxtal, "*/*/[F-obs,SIGF-obs]" );
                mtzin.import_dataset(opdset, "*/*/[F-obs,SIGF-obs]" );
                mtzin.close_read();
                notFound = false;
            }
            else if (mtzColumns[i].find("/F ") != -1)
            {
                std::cout << "\nMTZ file supplied. Using F & SIGF...\n";
                mtzin.import_hkl_data( fobs, "*/*/[F,SIGF]" );
                mtzin.import_crystal(opxtal, "*/*/[F,SIGF]" );
                mtzin.import_dataset(opdset, "*/*/[F,SIGF]" );
                mtzin.close_read();
                notFound = false;
            }
        }
    }

    if (notFound)
    {
        std::cout << "\nNo suitable amplitudes have been found in the MTZ file!\n\nSummoning ctruncate in case what we have are intensities...\n\n";

        char cmd[100];
        int exitCodeCTruncate;


        sprintf(cmd, "$CBIN/ctruncate -hklin %s -mtzout amplitudes.mtz -colin '/*/*/[I,SIGI]'", input_reflections_mtz.c_str());

        exitCodeCTruncate = system(cmd);

        // For future developer: Because I relocated this code from privateer.cpp to this file, I then couldn't return EXIT_FAILURE. If this thing is even called to begin with, then pass int exitCodeCTruncate by reference as a variable to this function as a fix to whatever bug may appear.
        if (exitCodeCTruncate != EXIT_SUCCESS)
            std::cout << "CTruncate exited with a failure, Privateer is likely to fail subsequently." << std::endl;
        else
        {
            std::cout << "\nReading output from ctruncate...\n" << "Previous hklinfo: " << hklinfo.cell().format() << std::endl;
            std::cout << " " << hklinfo.spacegroup().spacegroup_number() << " " << hklinfo.num_reflections() << "\n";
    

            ampmtzin.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
            ampmtzin.open_read( "amplitudes.mtz" );   // open file, no security checks
            ampmtzin.import_hkl_data( fobs, "*/*/[F,SIGF]" );
            ampmtzin.import_crystal(opxtal, "*/*/[F,SIGF]" );
            ampmtzin.import_dataset(opdset, "*/*/[F,SIGF]" );
            ampmtzin.close_read();

            std::cout << "\nPresent hklinfo: " << hklinfo.cell().format() << " " << hklinfo.spacegroup().spacegroup_number() << " " << hklinfo.num_reflections() << "\n";
        }
    }
}