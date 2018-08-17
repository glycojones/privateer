
// Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
//    -> The name combines Mark Knopfler's 'Privateering' cracking album & Dr Cowtan's tradition of using maritime acronyms
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2017 Jon Agirre & Kevin Cowtan
// York Structural Biology Laboratory
// The University of York
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk
//
// Targets supported: structures containing known or pre-defined furanose and pyranose rings
//
// The program, initially conceived for simple conformational validation, has evolved into a sort of
// Swiss-army knife, supporting re-refinement of sugars and graphical analysis of glycans through SVG graphics
//
// Some of it's features are:
//
//   - Cremer-pople puckering parameters for pyranose and furanose rings
//   - Comparison between detected and minimal-energy conformation
//   - Glycochemical validation: Privateer will detect wrong links and nomenclature
//   - Linkage torsion calculation and annotation
//   - Annotated SVG graphics files for N-, O- and S-glycans
//   - Generation of dictionaries containing aperiodic torsions (for re-refinement)
//   - Python interface for connecting to other programs (CCP4i2, CCP4mg and Coot)
//   - Coot-friendly output through Python and Scheme scripts
//   - SFcalc with bulk solvent correction from model (internal)
//   - Fo-Fc sigmaa weighted omit maps (output as map file or map coeffs)
//   - Real Space Correlation coefficient + score (console and XML)
//   - Mean B-factors, <mFo> and other metrics

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iostream>
#include "privateer-lib.h"
#include "clipper-glyco.h"
#include "blobs.h"
#include <clipper/clipper.h>
#include <clipper/clipper-cif.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <clipper/contrib/sfcalc_obs.h>
#include <clipper/minimol/minimol_utils.h>



clipper::String program_version = "MKIV_a";
using clipper::data32::F_sigF;
using clipper::data32::F_phi;
using clipper::data32::Phi_fom;
using clipper::data32::Flag;
typedef clipper::HKL_data_base::HKL_reference_index HRI;
using namespace std;


int main(int argc, char** argv)
{

    CCP4Program prog( "Privateer", program_version.c_str(), "$Date: 2015/07/10" );

    prog.set_termination_message( "Failed" );

    std::cout << "\nCopyright 2013-2016 Jon Agirre, Kevin Cowtan and The University of York." << std::endl  ;
    std::cout << "\nWith contributions from Haroldas Bagdonas, Marcin Wojdyr and Jacob Sorensen." << std::endl  ;
    std::cout << "\n\nPlease reference these articles: "<< std::endl ;
    std::cout << "\n  'Privateer: software for the conformational validation of carbohydrate structures'";
    std::cout << "\n   Agirre J, Fernandez-Iglesias J, Rovira C, Davies GJ, Wilson KS and Cowtan KD. (2015) Nature Structural & Molecular Biology 22 (11), 833-834." << std::endl;
    std::cout << "\n  'Carbohydrate anomalies in the PDB'";
    std::cout << "\n   Agirre J, Davies GJ, Wilson KS and Cowtan KD. (2015) Nature Chemical Biology 11 (5), 303." << std::endl << std::endl;

    std::size_t found = program_version.find("_"); // if found, this indicates that we are running a test version

    if ( found != std::string::npos )
        std::cout << "   \nWARNING: TEST VERSION - DO NOT USE IN A PRODUCTION ENVIRONMENT.\n" << std::endl;

    clipper::HKL_info hklinfo; // allocate space for the hkl metadata
    clipper::CIFfile cifin;
    clipper::CCP4MTZfile mtzin, ampmtzin;
    clipper::String ippdb    = "NONE";
    clipper::String ipcol_fo = "NONE";
    clipper::String ipsfcif  = "NONE";
    clipper::String ipmmcif  = "NONE";
    clipper::String ipcode   = "XXX";
    clipper::String opfile   = "privateer-hklout.mtz";
    clipper::String title    = "generic title";
    clipper::String ipmtz    = "NONE";
    clipper::String expsys   = "undefined";
    clipper::String validation_string = "";
    std::vector<clipper::String> validation_options;
    clipper::data::sugar_database_entry external_validation;
    bool useSigmaa = false;
    bool vertical(false), original(false), invert(false);
    int n_refln = 1000;
    int n_param = 20;
    bool useMTZ = false;
    bool batch = false;
    bool noMaps = false;
    bool allSugars = true;
    bool showGeom = false;
    bool check_unmodelled = false;
    float ipradius = 2.5;    // default value, punishing enough!
    FILE *output;
    bool output_mtz = false;
    std::vector < clipper::MGlycan > list_of_glycans;
    clipper::CCP4MTZfile opmtz_best, opmtz_omit;
    clipper::MTZcrystal opxtal;
    clipper::MTZdataset opdset;
    clipper::MGlycology mgl;

    // command input
    CCP4CommandInput args( argc, argv, true );
    int arg = 0;
    while ( ++arg < args.size() )
    {
        if ( args[arg] == "-title" )
        {
            if ( ++arg < args.size() )
                title = args[arg];
        }
        else if ( args[arg] == "-colin-fo" )
        {
            if  (++arg < args.size() )
                ipcol_fo = args[arg];
        }
        else if ( args[arg] == "-mtzout" )
        {
            if ( ++arg < args.size() )
            {
                opfile = args[arg];
                output_mtz = true;
            }
        }
        else if ( args[arg] == "-pdbin" )
        {
            if ( ++arg < args.size() )
                ippdb = args[arg];
        }
        else if ( args[arg] == "-mtzin" )
        {
            if ( ++arg < args.size() )
            {
                useMTZ = true;
                ipmtz = args[arg];
            }
        }
        else if ( args[arg] == "-cifin" )
        {
            if ( ++arg < args.size() )
                ipsfcif = args[arg];
        }
        else if ( args[arg] == "-list" )
        {
            privateer::util::print_supported_code_list();
            prog.set_termination_message( "Success" );
            return 0;
        }
        else if ( args[arg] == "-expression" )
        {
            if ( ++arg < args.size() )
                expsys = args[arg];
            if ( expsys.trim() != "undefined" &&
                 expsys.trim() != "fungal" &&
                 expsys.trim() != "yeast" &&
                 expsys.trim() != "plant" &&
                 expsys.trim() != "insect" &&
                 expsys.trim() != "mammalian" &&
                 expsys.trim() != "human" )
            {
                std::cout << std::endl << std::endl << "Error: " << expsys << " is not a supported expression system" << std::endl << std::endl;
                privateer::util::print_usage();
                prog.set_termination_message( "Failed" );
                return 1;
            }
        }
        else if ( args[arg] == "-codein" )
        {
            if ( ++arg < args.size() )
            {
                ipcode = clipper::String(args[arg]);
                allSugars = false;

                if (ipcode.trim().length() != 3)
                {
                    std::cout << std::endl << std::endl << "Error: the sugar code must be three characters long, e.g. GLC."
                    << "\nPlease refer to the Chemical Component Dictionary (http://www.wwpdb.org/ccd.html).\n"
                    << "Alternatively, use the -list option to get a full list of supported codes\nExiting..." << std::endl << std::endl;
                    prog.set_termination_message( "Failed" );
                    return 1;
                }
            }
        }
        else if ( args[arg] == "-mode" )
        {
            if ( ++arg < args.size() )
                if (clipper::String(args[arg]) == "ccp4i2")
                    batch = true;
        }
        else if ( args[arg] == "-vertical" )
        {
            vertical = true;
        }
        else if ( args[arg] == "-essentials" )
        {
            original = true;
        }
        else if ( args[arg] == "-invert" )
        {
            invert = true;
        }
        else if ( args[arg] == "-radiusin" )
        {
            if ( ++arg < args.size() )
            {
                ipradius = clipper::String(args[arg]).f();
                if (ipradius < 1.0)
                {
                    std::cout << "\n\nMask radius is too small!" << std::endl << std::endl;
                    prog.set_termination_message( "Failed" );
                    return 1;
                }
            }
        }
        else if ( args[arg] == "-valstring" )
        {
            if ( ++arg < args.size() )
            {
                validation_string = args[arg];

                validation_options = validation_string.split(",");
                if ( privateer::util::compute_and_print_external_validation ( validation_options, external_validation ) )
                {
                    prog.set_termination_message ( "Failed" );
                    return -1;
                }
            }
        }
        else if ( args[arg] == "-showgeom" )
            showGeom = true;

        else if ( args[arg] == "-nomaps" )
            noMaps = true;

        else if ( args[arg] == "-check-unmodelled" )
            check_unmodelled = true;

        else
        {
            std::cout << "\nUnrecognised:\t" << args[arg] << std::endl;
            args.clear();
        }
    }

    if ( ipcode != "XXX" )
    {
        if ((( !clipper::MSugar::search_database(ipcode.c_str())) && (clipper::MDisaccharide::search_disaccharides(ipcode.c_str())==-1)) && (validation_options.size() == 0) )
        {
            std::cout << "\n\nError: no internal validation data found for " << ipcode << std::endl;
            std::cout << "\nYou can provide external validation data with -valstring <data>" << std::endl;
            std::cout << "\n\tAccepted format: SUG,O5/C1/C2/C3/C4/C5,A,D,4c1\n";
            std::cout << "\tThree-letter code, ring atoms separated by /, anomer, handedness, expected conformation.\n" << std::endl;
            prog.set_termination_message("Failed");
            return 1;
        }
    }


    if (batch)
    {
        output = fopen("validation_data-privateer","w");
        if (NULL == output)
        {
            std::cout << std::endl << "Error: unable to create output file. Please check directory permissions." << std::endl;
            prog.set_termination_message( "Failed" );
            return 1;
        }
    }

    if ( (ippdb != "NONE") && ((ipsfcif == "NONE") && (ipmtz == "NONE")) )
        noMaps = true;

    if ( (ippdb == "NONE") || ((ipsfcif == "NONE") && (ipmtz == "NONE") && (noMaps == false)) )
    {
        privateer::util::print_usage();
        prog.set_termination_message( "Failed" );
        return 1;
    }

    clipper::MMDBfile mfile;
    clipper::MiniMol mmol;

    privateer::util::read_coordinate_file ( mfile, mmol, ippdb, batch);
    int pos_slash = ippdb.rfind("/");



    // Fast mode, no maps nor correlation calculations

    if ( noMaps )
    {
        clipper::Atom_list mainAtoms;
        clipper::Atom_list ligandAtoms;
        clipper::Atom_list allAtoms;
        std:vector< std::string > enable_torsions_for;

        if (!batch)
        {
            std::cout << std::endl << "Analysing carbohydrates... \n";
            fflush(0);
        }

        std::vector < std::pair <clipper::String , clipper::MSugar> > ligandList; // we store the Chain ID and create an MSugar to be scored
        std::vector < clipper::MMonomer > sugarList; // store the original MMonomer

        const clipper::MAtomNonBond& manb = clipper::MAtomNonBond( mmol, 1.0 ); // was 1.0

        mgl = clipper::MGlycology(mmol, manb, expsys);

        list_of_glycans = mgl.get_list_of_glycans();

        if ( !batch ) std::cout << std::endl << "Number of detected glycosylations: " << list_of_glycans.size();

        if ( list_of_glycans.size() > 0 )
        {
            clipper::String current_chain = "" ;

            for (int i = 0; i < list_of_glycans.size() ; i++ )
            {
                if ( current_chain != list_of_glycans[i].get_chain() )
                {
                    current_chain = list_of_glycans[i].get_chain();
                    std::cout << std::endl << std::endl << "Chain " << current_chain[0] << std::endl << "-------" << std::endl ;
                }
                std::cout << std::endl << list_of_glycans[i].print_linear ( true, false, true ) << std::endl;
                privateer::glycoplot::Plot plot(vertical, original, list_of_glycans[i].get_root_by_name(), invert, true);
                plot.plot_glycan ( list_of_glycans[i] );
                std::ostringstream os;
                os << list_of_glycans[i].get_root_for_filename() << ".svg";
                plot.write_to_file ( os.str() );
                plot.delete_shapes();
            }
        }

        if ( !batch ) std::cout << "\n\nDetailed validation data" << std::endl;
        if ( !batch ) std::cout << "------------------------" << std::endl;

        // erase ligand atoms from the model and then calculate phases using
        // the omitted model, effectively computing an omit map

        for ( int p = 0; p < mmol.size(); p++ )
        {
            for ( int m = 0; m < mmol[p].size(); m++ )
            {
                if (allSugars)
                {
                    if ( clipper::MDisaccharide::search_disaccharides(mmol[p][m].type().c_str()) != -1 ) // treat disaccharide
                    {
                        clipper::MDisaccharide md(mmol, manb, mmol[p][m] );
                        sugarList.push_back ( mmol[p][m] );
                        sugarList.push_back ( mmol[p][m] );
                        clipper::String id = mmol[p].id();
                        id.resize(1);

                        ligandList.push_back ( std::pair < clipper::String, clipper::MSugar> (id, md.get_first_sugar()));
                        ligandList.push_back ( std::pair < clipper::String, clipper::MSugar> (id, md.get_second_sugar()));
                    }
                    else if ( !clipper::MSugar::search_database(mmol[p][m].type().c_str()) ) // true if strings are different
                    {
                        for (int id = 0; id < mmol[p][m].size(); id++ )
                        {
                            mainAtoms.push_back(mmol[p][m][id]); // cycle through atoms and copy them
                            allAtoms.push_back(mmol[p][m][id]);
                        }
                    }
                    else // it's one of the sugars contained in the database
                    {
                        clipper::MSugar msug, msug_b;

                        std::vector <char> conformers = privateer::util::number_of_conformers(mmol[p][m]);

                        #ifdef DUMP
                            std::cout << "number of alternate conformations: " << conformers.size() << std::endl;
                        #endif

                        int n_conf = conformers.size();

                        if ( n_conf > 0 )
                        {
                            if ( n_conf == 1 )
                                msug   = clipper::MSugar(mmol, mmol[p][m], manb, conformers[0]);
                            else
                            {
                                msug   = clipper::MSugar(mmol, mmol[p][m], manb, conformers[0]);
                                msug_b = clipper::MSugar(mmol, mmol[p][m], manb, conformers[1]);
                            }

                        }
                        else
                        {
                            msug = clipper::MSugar(mmol, mmol[p][m], manb);
                        }

                        sugarList.push_back(mmol[p][m]);
                        clipper::String id = mmol[p].id();
                        id.resize(1);

                        #ifdef DUMP
                            std::cout << "Looking at chain " << id << std::endl;
                        #endif

                        ligandList.push_back(std::pair<clipper::String, clipper::MSugar> (id, msug));
                        // add both conformers if the current monomer contains more than one
                        if ( n_conf == 2 )
                            ligandList.push_back(std::pair<clipper::String, clipper::MSugar> (id, msug_b));

                        if ( msug.type_of_sugar() == "unsupported" )
                        {
                            std::cout << std::endl << "Error: at least one of the sugars in the supplied PDB file is missing expected atoms"
                                      << std::endl << std::endl;
                            prog.set_termination_message( "Unrecoverable error" );
                            return 0;
                        }

                        for (int id = 0; id < mmol[p][m].size(); id++ )
                        {
                            ligandAtoms.push_back(mmol[p][m][id]);  // add the ligand atoms to a second array
                            allAtoms.push_back(mmol[p][m][id]);
                        }
                    }

                }

                else
                {
                    if ( strncmp( mmol[p][m].type().c_str(), ipcode.trim().c_str(), 3 )) // true if strings are different
                    {
                        for (int id = 0; id < mmol[p][m].size(); id++ )
                        {
                            mainAtoms.push_back(mmol[p][m][id]); // cycle through atoms and copy them
                            allAtoms.push_back(mmol[p][m][id]);
                        }
                    }
                    else // it's the one sugar we're looking to omit
                    {
                        if ( validation_options.size() > 0 )
                        {
                            const clipper::MSugar msug ( mmol, mmol[p][m], manb, external_validation );

                            sugarList.push_back(mmol[p][m]);
                            ligandList.push_back(std::pair<clipper::String, clipper::MSugar> (mmol[p].id().trim(), msug));
                        }
                        else
                        {
                            const clipper::MSugar msug(mmol, mmol[p][m], manb);

                            sugarList.push_back(mmol[p][m]);
                            ligandList.push_back(std::pair<clipper::String, clipper::MSugar> (mmol[p].id().trim(), msug));
                        }

                        for (int id = 0; id < mmol[p][m].size(); id++ )
                        {
                            ligandAtoms.push_back(mmol[p][m][id]);  // add the ligand atoms to a second array
                            allAtoms.push_back(mmol[p][m][id]);
                        }
                    }
                }
            }
        }

        if (!batch) printf("\nPDB \t    Sugar   \t  Q  \t Phi  \tTheta \t   Detected type   \tCnf\t<Bfac>\tCtx\t Ok?");
        if (!batch && showGeom) printf("\tBond lengths, angles and torsions, reported clockwise with in-ring oxygen as first vertex");
        if (!batch) printf("\n----\t------------\t-----\t------\t------\t-------------------\t---\t------\t---\t-----");
        if (!batch && showGeom) printf("\t------------------------------------------------------------------------------------------------------------");
        if (!batch) printf("\n");

        for (int index = 0; index < ligandList.size(); index++)
        {
            float x, y, z, maxX, maxY, maxZ, minX, minY, minZ;
            x = y = z = 0.0;
            maxX = maxY = maxZ = -999999.0;
            minX = minY = minZ =  999999.0;

            for (int natom = 0; natom < sugarList[index].size(); natom++)
            {
                if(sugarList[index][natom].coord_orth().x() > maxX) maxX=sugarList[index][natom].coord_orth().x(); // calculation of the sugar centre
                if(sugarList[index][natom].coord_orth().y() > maxY) maxY=sugarList[index][natom].coord_orth().y();
                if(sugarList[index][natom].coord_orth().z() > maxZ) maxZ=sugarList[index][natom].coord_orth().z();
                if(sugarList[index][natom].coord_orth().x() < minX) minX=sugarList[index][natom].coord_orth().x();
                if(sugarList[index][natom].coord_orth().y() < minY) minY=sugarList[index][natom].coord_orth().y();
                if(sugarList[index][natom].coord_orth().z() < minZ) minZ=sugarList[index][natom].coord_orth().z();
            }

            x = minX + ((maxX - minX)/2);
            y = minY + ((maxY - minY)/2);
            z = minZ + ((maxZ - minZ)/2);

            if (batch)
            {
                fprintf(output, "%c%c%c%c\t%s-",ippdb[1+pos_slash],ippdb[2+pos_slash],ippdb[3+pos_slash],ippdb[4+pos_slash], ligandList[index].second.type().c_str());
                fprintf(output, "%s-%s   ", ligandList[index].first.c_str(), ligandList[index].second.id().trim().c_str());
            }
            else
            {
                printf("%c%c%c%c\t%s-",ippdb[1+pos_slash],ippdb[2+pos_slash],ippdb[3+pos_slash],ippdb[4+pos_slash], ligandList[index].second.type().c_str());
                std::cout << ligandList[index].first << "-" << ligandList[index].second.id().trim() << "  ";
            }

            if (batch)
            {
                std::vector<clipper::ftype> cpParams(10, 0);
                cpParams = ligandList[index].second.cremer_pople_params();
                fprintf(output,"\t%1.3f\t%3.2f\t",cpParams[0],cpParams[1]);  // output cremer-pople parameters
                if ( cpParams[2] == -1 ) fprintf ( output, " --  \t" ); else fprintf ( output, "%3.2f\t", cpParams[2] );
                fprintf(output,"%s\t", ligandList[index].second.type_of_sugar().c_str()); // output the type of sugar, e.g. alpha-D-aldopyranose
                fprintf(output,"%s\t", ligandList[index].second.conformation_name().c_str()); // output a 3 letter code for the conformation

                float bfac = 0.0;
                for (int i=0; i < ligandList[index].second.size(); i++)
                    bfac+=ligandList[index].second[i].u_iso();
                bfac /= ligandList[index].second.size();
                bfac  = clipper::Util::u2b(bfac);

                fprintf ( output, "%3.2f", bfac );     // output <bfactor>


                std::vector < clipper::MGlycan > list_of_glycans = mgl.get_list_of_glycans();
                bool found_in_tree = false;

                for ( int i = 0 ; i < list_of_glycans.size() ; i++ )
                {
                    std::vector < clipper::MSugar > list_of_sugars = list_of_glycans[i].get_sugars();

                    for ( int j = 0 ; j < list_of_sugars.size() ; j++ )
                    {
                        if ( list_of_sugars[j].id().trim() == ligandList[index].second.id().trim() )
                        {
                            if ( list_of_glycans[i].get_type() == "n-glycan" )
                            {
                                ligandList[index].second.set_context ( "n-glycan" );
                                fprintf ( output, "\t(n) " );
                            }
                            else
                            {
                                ligandList[index].second.set_context ( "o-glycan" );
                                fprintf ( output, "\t(o) " );
                            }
                            found_in_tree = true;
                            break;
                        }
                    }

                    if ( found_in_tree ) break;
                }

                if ( !found_in_tree )
                {
                    ligandList[index].second.set_context ( "ligand" );
                    fprintf ( output, "\t(l) ");
                }


                if (ligandList[index].second.in_database(ligandList[index].second.type().trim()))
                {
                    if ((ligandList[index].second.ring_members().size() == 6 ))
                    {
                        if (ligandList[index].second.is_sane())
                        {
                            if ( ! ligandList[index].second.ok_with_conformation() )
                            {
                                fprintf(output, "\tcheck");
                            }
                            else
                                fprintf(output, "\tyes");
                        }
                        else fprintf (output, "\tno");
                    }
                    else
                        if (ligandList[index].second.is_sane())
                            fprintf(output, "\tyes");
                        else
                        {
                            fprintf(output, "\tno");
                        }
                }
                else fprintf(output, "\tunk");

                if ( ! ligandList[index].second.ok_with_conformation () )
                    enable_torsions_for.push_back (ligandList[index].second.type().trim());

                bool occupancy_check = false;
                std::vector<clipper::MAtom> ringcomponents = ligandList[index].second.ring_members();

                for ( int i = 0 ; i < ringcomponents.size() ; i++ )
                {
                    if (privateer::util::get_altconformation(ringcomponents[i]) != ' ')
                        occupancy_check = true;
                }

                if (showGeom)
                {
                    std::vector<clipper::ftype> rangles = ligandList[index].second.ring_angles();
                    std::vector<clipper::ftype> rbonds  = ligandList[index].second.ring_bonds();
                    std::vector<clipper::ftype> rtorsions = ligandList[index].second.ring_torsions();

                    for (int i = 0 ; i < ligandList[index].second.ring_members().size(); i++ )
                        fprintf(output, "\t%1.2f", rbonds[i]);

                    for (int i = 0 ; i < ligandList[index].second.ring_members().size(); i++ )
                        fprintf(output, "\t%3.1f", rangles[i]);

                    for (int i = 0 ; i < ligandList[index].second.ring_members().size(); i++ )
                        fprintf(output, "\t%3.1f", rtorsions[i]);
                }

                if (occupancy_check)
                    fprintf(output, " (*)");

                fprintf(output, "\n");

            }
            else
            {
                std::vector<clipper::ftype> cpParams(10, 0);
                cpParams = ligandList[index].second.cremer_pople_params();
                printf("\t%1.3f\t%3.2f\t",cpParams[0],cpParams[1]);  // output cremer-pople parameters
                if ( cpParams[2] == -1 ) printf ( " --  \t" ); else printf ( "%3.2f\t", cpParams[2] );
                printf("%s\t", ligandList[index].second.type_of_sugar().c_str()); // output the type of sugar, e.g. alpha-D-aldopyranose
                printf("%s\t", ligandList[index].second.conformation_name().c_str()); // output a 3 letter code for the conformation

                float bfac = 0.0;

                for (int i=0; i < ligandList[index].second.size(); i++)
                    bfac+=ligandList[index].second[i].u_iso();

                bfac /= ligandList[index].second.size();
                bfac  = clipper::Util::u2b(bfac);
                printf ( "%3.2f", bfac ); // output <Bfactor>


                std::vector < clipper::MGlycan > list_of_glycans = mgl.get_list_of_glycans();
                bool found_in_tree = false;

                for ( int i = 0 ; i < list_of_glycans.size() ; i++ )
                {
                    std::vector < clipper::MSugar > list_of_sugars = list_of_glycans[i].get_sugars();

                    for ( int j = 0 ; j < list_of_sugars.size() ; j++ )
                    {
                        if ( list_of_sugars[j].id().trim() == ligandList[index].second.id().trim() )
                        {
                            if ( list_of_glycans[i].get_type() == "n-glycan" )
                            {
                                ligandList[index].second.set_context ( "n-glycan" );
                                std::cout << "\t(n) ";
                            }
                            else
                            {
                                ligandList[index].second.set_context ( "o-glycan" );
                                std::cout << "\t(o) ";
                            }
                            found_in_tree = true;
                            break;
                        }
                    }

                    if ( found_in_tree ) break;
                }

                if ( !found_in_tree )
                {
                    ligandList[index].second.set_context ( "ligand" );
                    std::cout << "\t(l) ";
                }


                if (ligandList[index].second.in_database(ligandList[index].second.type().trim()))
                {
                    if ((ligandList[index].second.ring_members().size() == 6 ))
                    {
                        if (ligandList[index].second.is_sane())
                        {
                            if ( ! ligandList[index].second.ok_with_conformation () )
                            {
                                printf("\tcheck");
                            }
                            else printf("\tyes");
                        }
                        else
                            printf ("\tno");
                    }
                    else
                        if (ligandList[index].second.is_sane())
                            printf("\tyes");
                        else
                        {
                            printf("\tno");
                        }
                }
                else printf("\tunk");

                if ( ! ligandList[index].second.ok_with_conformation () )
                    enable_torsions_for.push_back (ligandList[index].second.type().trim());

                bool occupancy_check = false;
                std::vector<clipper::MAtom> ringcomponents = ligandList[index].second.ring_members();

                for ( int i = 0 ; i < ringcomponents.size() ; i++ )
                {
                    if (privateer::util::get_altconformation(ringcomponents[i]) != ' ')
                        occupancy_check = true;
                }

                if (showGeom)
                {
                    std::vector<clipper::ftype> rangles = ligandList[index].second.ring_angles();
                    std::vector<clipper::ftype> rbonds  = ligandList[index].second.ring_bonds();
                    std::vector<clipper::ftype> rtorsions = ligandList[index].second.ring_torsions();


                    for (int i = 0 ; i < ligandList[index].second.ring_members().size(); i++ )
                        printf("\t%1.2f", rbonds[i]);
                    for (int i = 0 ; i < ligandList[index].second.ring_members().size(); i++ )
                        printf("\t%3.1f", rangles[i]);
                    for (int i = 0 ; i < ligandList[index].second.ring_members().size(); i++ )
                        printf("\t%3.1f", rtorsions[i]);
                }

                if (occupancy_check)
                    std::cout << " (*)";

                std::cout << std::endl;
            }
        }

        if (!batch)
        {
            std::cout << "\nPartially occupied monosaccharides, if any, are marked with an asterisk (*)";
            std::cout << std::endl << std::endl;
        }
        else
            fclose(output);

        std::fstream of_scm; std::fstream of_py;
        of_scm.open("privateer-results.scm", std::fstream::out);
        of_py.open("privateer-results.py", std::fstream::out);
        privateer::coot::insert_coot_prologue_scheme ( of_scm );
        privateer::coot::insert_coot_prologue_python ( of_py );

        clipper::String all_MapName, dif_MapName, omit_dif_MapName;
        all_MapName = ""; dif_MapName = ""; omit_dif_MapName = "";

        privateer::coot::insert_coot_files_loadup_scheme (of_scm, ippdb, all_MapName, dif_MapName, omit_dif_MapName, batch );
        privateer::coot::insert_coot_files_loadup_python (of_py , ippdb, all_MapName, dif_MapName, omit_dif_MapName, batch );

        int n_geom = 0, n_anomer = 0, n_config = 0, n_pucker = 0, n_conf = 0;

        int sugar_count = 0;

        for (int k = 0 ; k < ligandList.size() ; k++)
        {
            int n_errors = 0;

            clipper::String diagnostic = ligandList[k].second.type().trim() + "/" + ligandList[k].first + "/" + ligandList[k].second.id().trim() + ": " ;

            if (ligandList[k].second.is_supported() )
            {
                if ( ! ligandList[k].second.is_sane() )
                {
                    sugar_count++;

                    if (!ligandList[k].second.ok_with_anomer())
                    {
                        if ( n_errors > 0 )
                            diagnostic.append(", wrong anomer");
                        else
                            diagnostic.append("wrong anomer");

                        n_errors++; n_anomer++;
                    }

                    if (!ligandList[k].second.ok_with_chirality())
                    {
                        if ( n_errors > 0 )
                            diagnostic.append(", wrong configuration at " + ligandList[k].second.configurational_carbon().name().trim());
                        else
                            diagnostic.append("wrong configuration at " + ligandList[k].second.configurational_carbon().name().trim());

                        n_errors++; n_config++;
                    }

                    if (!ligandList[k].second.ok_with_puckering())
                    {
                        if ( n_errors > 0 )
                            diagnostic.append(", Q=" + clipper::String(ligandList[k].second.puckering_amplitude()) );
                        else
                            diagnostic.append("Q=" + clipper::String( ligandList[k].second.puckering_amplitude() ));

                        n_errors++; n_pucker++;
                    }

                    if (!ligandList[k].second.ok_with_conformation() )
                    {
                        if ( n_errors > 0 )
                            diagnostic.append(", conformation (" + ligandList[k].second.conformation_name() + clipper::String(") might be mistaken") );
                        else
                            diagnostic.append("conformation (" + ligandList[k].second.conformation_name() + clipper::String(") might be mistaken"));

                        n_errors++; n_conf++;
                    }

                    privateer::coot::insert_coot_go_to_sugar_scheme ( of_scm, ligandList[k].second.ring_centre(), diagnostic );
                    privateer::coot::insert_coot_go_to_sugar_python ( of_py, ligandList[k].second.ring_centre(), diagnostic );
                }
                else // sugar is sane, but still need to check higher-energy conformations
                {
                    if ( !ligandList[k].second.ok_with_conformation() )
                    {
                        if ( n_errors > 0 )
                            diagnostic.append(", conformation (" + ligandList[k].second.conformation_name() + clipper::String(") might be mistaken") );
                        else
                            diagnostic.append("conformation (" + ligandList[k].second.conformation_name() + clipper::String(") might be mistaken"));

                        n_errors++;
                        n_conf++;
                        sugar_count++;

                        privateer::coot::insert_coot_go_to_sugar_scheme ( of_scm, ligandList[k].second.ring_centre(), diagnostic );
                        privateer::coot::insert_coot_go_to_sugar_python ( of_py, ligandList[k].second.ring_centre(), diagnostic );
                    }
                }
            }
        }

        clipper::String status_msg = "Blue map: 2mFo-DFc. Pink map: omit mFo-DFc. Torsion restraints have been enabled.";

        privateer::coot::insert_coot_epilogue_scheme ( of_scm );
        privateer::coot::insert_coot_epilogue_python ( of_py );

        privateer::coot::insert_coot_statusbar_text_scheme ( of_scm, status_msg );
        privateer::coot::insert_coot_statusbar_text_python ( of_py, status_msg );

        of_scm.close();
        of_py.close();

        std::cout << "SUMMARY: " << std::endl << std::endl ;
        std::cout << "   Wrong anomer: " << n_anomer << std::endl;
        std::cout << "   Wrong configuration: " << n_config << std::endl;
        std::cout << "   Unphysical puckering amplitude: " << n_pucker << std::endl;
        std::cout << "   In higher-energy conformations: " << n_conf << std::endl;
        std::cout << std::endl;
        std::cout << "   Privateer has identified " << n_anomer + n_config + n_pucker + n_conf;
        std::cout << " issues, with " << sugar_count << " of " << ligandList.size() << " sugars affected." << std::endl;

        if ( enable_torsions_for.size() > 0 )
        {
            privateer::util::write_refmac_keywords ( enable_torsions_for );
            privateer::util::write_libraries( enable_torsions_for );
        }

        prog.set_termination_message( "Normal termination" );
        system("touch scored");
        return 0;

    }


    // Full analysis, slower but much more complete //

    std::vector< std::string > enable_torsions_for;


    if (useMTZ)
    {
        // grab the structure factors from an MTZ file

        if (!batch)
        {
            std::cout << "Reading " << ipmtz.trim().c_str() << "... ";
            fflush(0);
        }
        mtzin.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
        mtzin.open_read( ipmtz.trim() );

        if (!batch)
            std::cout << "done." << std::endl;

        try // we could be in trouble should the MTZ file have no cell parameters
        {
            mtzin.import_hkl_info( hklinfo );     // read spacegroup, cell, resolution, HKL's
        }
        catch (...)
        {
            if (!batch)
            {
                std::cout << "\nReading cell and spacegroup parameters from the CRYST1 card in ";
                std::cout << ippdb << ":\n Spacegroup (" << mmol.spacegroup().spacegroup_number() << ")\n" << mmol.cell().format() << "\n\n" ;
            }

            clipper::Resolution myRes(0.96);
            hklinfo = clipper::HKL_info( mmol.spacegroup(), mmol.cell(), myRes, true);
        }
    }
    else
    {   // assume CIF file format instead
        if (!batch)
            std::cout << "Opening CIF file: " << ipsfcif.trim() << std::endl;

        cifin.open_read( ipsfcif.trim().c_str() );

        try // we could be in trouble should the CIF file not have cell parameters
        {
            cifin.import_hkl_info( hklinfo );     // read spacegroup, cell, resolution, HKL's
        }
        catch (...)
        {
            if (!batch) std::cout << "\nReading cell and spacegroup parameters from the CRYST1 card in " << ippdb;
            std::cout << ":\n Spacegroup (" << mmol.spacegroup().spacegroup_number() << ")\n" << mmol.cell().format() << "\n\n";

            clipper::Resolution myRes(0.96);
            hklinfo = clipper::HKL_info( mmol.spacegroup(), mmol.cell(), myRes, true);
        }
    }

    clipper::HKL_data<clipper::data32::F_sigF> fobs ( hklinfo );            // allocate space for F and sigF
    clipper::HKL_data<clipper::data32::F_sigF> fobs_scaled ( hklinfo );     // allocate space for scaled F and sigF
    clipper::HKL_data<clipper::data32::F_phi> fc_omit_bsc ( hklinfo );      // allocate space for the omit calc data with bsc
    clipper::HKL_data<clipper::data32::F_phi> fc_all_bsc ( hklinfo );       // allocate space for the whole calculated data with bsc
    clipper::HKL_data<clipper::data32::F_phi> fc_ligands_bsc( hklinfo );    // allocate space for the ligand calculated data

    bool notFound = true;

    // now scan for some of the most used column ID's

    if (useMTZ)
    {
        std::vector<clipper::String> mtzColumns;
        mtzColumns = mtzin.column_labels();

        if (ipcol_fo != "NONE")
        {
            if (!batch) std::cout << "MTZ file supplied. Using " << ipcol_fo << "...\n";
            mtzin.import_hkl_data( fobs, "*/*/["+ ipcol_fo+"]" );
            mtzin.import_crystal(opxtal, ipcol_fo);
            mtzin.import_dataset(opdset, ipcol_fo);
            mtzin.close_read();
            notFound = false;
        }
        else
        {
            for ( int i = 0 ; ((i < mtzColumns.size()) && notFound) ; i++ )
            {
                if (mtzColumns[i].find("/FOBS ") != -1)
                {
                    if (!batch) std::cout << "\nMTZ file supplied. Using FOBS & SIGFOBS...\n";
                    mtzin.import_hkl_data( fobs, "*/*/[FOBS,SIGFOBS]" );
                    mtzin.import_crystal(opxtal, "*/*/[FOBS,SIGFOBS]" );
                    mtzin.import_dataset(opdset, "*/*/[FOBS,SIGFOBS]" );
                    mtzin.close_read();
                    notFound = false;
                }
                else if (mtzColumns[i].find("/FP ") != -1)
                {
                    if (!batch) std::cout << "\nMTZ file supplied. Using FP & SIGFP...\n";
                    mtzin.import_hkl_data( fobs, "*/*/[FP,SIGFP]" );
                    mtzin.import_crystal(opxtal, "*/*/[FP,SIGFP]" );
                    mtzin.import_dataset(opdset, "*/*/[FP,SIGFP]" );
                    mtzin.close_read();
                    notFound = false;
                }
                else if (mtzColumns[i].find("/FOSC ") != -1)
                {
                    if (!batch) std::cout << "\nMTZ file supplied. Using FOSC & SIGFOSC...\n";
                    mtzin.import_hkl_data( fobs, "*/*/[FOSC,SIGFOSC]" );
                    mtzin.import_crystal(opxtal, "*/*/[FOSC,SIGFOSC]" );
                    mtzin.import_dataset(opdset, "*/*/[FOSC,SIGFOSC]" );
                    mtzin.close_read();
                    notFound = false;
                }
                else if (mtzColumns[i].find("/F-obs ") != -1)
                {
                    if (!batch) std::cout << "\nMTZ file supplied. Using F-obs & SIGF-obs...\n";
                    mtzin.import_hkl_data( fobs, "*/*/[F-obs,SIGF-obs]" );
                    mtzin.import_crystal(opxtal, "*/*/[F-obs,SIGF-obs]" );
                    mtzin.import_dataset(opdset, "*/*/[F-obs,SIGF-obs]" );
                    mtzin.close_read();
                    notFound = false;
                }
                else if (mtzColumns[i].find("/F ") != -1)
                {
                    if (!batch) std::cout << "\nMTZ file supplied. Using F & SIGF...\n";
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
            if (!batch)
                std::cout << "\nNo suitable amplitudes have been found in the MTZ file!\n\nSummoning ctruncate in case what we have are intensities...\n\n";

            char cmd[80];
            int exitCodeCTruncate;

            if (!batch)
                sprintf(cmd, "$CBIN/ctruncate -hklin %s -mtzout amplitudes.mtz -colin '/*/*/[I,SIGI]'", ipmtz.c_str());
            else
                sprintf(cmd, "$CBIN/ctruncate -hklin %s -mtzout amplitudes.mtz -colin '/*/*/[I,SIGI]' >/dev/null", ipmtz.c_str());

            exitCodeCTruncate = system(cmd);

            if (exitCodeCTruncate != EXIT_SUCCESS)
                return EXIT_FAILURE;

            else
            {
                if (!batch)
                {
                    std::cout << "\nReading output from ctruncate...\n" << "Previous hklinfo: " << hklinfo.cell().format() << std::endl;
                    std::cout << " " << hklinfo.spacegroup().spacegroup_number() << " " << hklinfo.num_reflections() << "\n";
                }

                ampmtzin.set_column_label_mode( clipper::CCP4MTZfile::Legacy );
                ampmtzin.open_read( "amplitudes.mtz" );   // open file, no security checks
                ampmtzin.import_hkl_data( fobs, "*/*/[F,SIGF]" );
                ampmtzin.import_crystal(opxtal, "*/*/[F,SIGF]" );
                ampmtzin.import_dataset(opdset, "*/*/[F,SIGF]" );
                ampmtzin.close_read();

                if (!batch) std::cout << "\nPresent hklinfo: " << hklinfo.cell().format() << " " << hklinfo.spacegroup().spacegroup_number() << " " << hklinfo.num_reflections() << "\n";
            }
        }
    }

    else
    {
        try
        {
            cifin.import_hkl_data( fobs );
            cifin.close_read();
        }
        catch (...)
        {
            if (!batch) std::cout << "\nNo suitable amplitudes found in the CIF file!\n\nExiting...\n\n";
            return 1;
        }
    }

    fobs_scaled = fobs;


        if ( !batch )
        {

            std::vector<std::vector<GlycosylationMonomerMatch> > glycosylationFromSequence = get_matching_monomer_positions(ippdb);
            clipper::MiniMol modelRemovedWaters = get_model_without_waters(ippdb);

            clipper::Atom_list withoutWaterModelAtomList = modelRemovedWaters.atom_list();

            clipper::Grid_sampling fakegrid( hklinfo.spacegroup(), hklinfo.cell(), hklinfo.resolution() );

            clipper::Xmap<float> two_times_observed_minus_one_model(mmol.spacegroup(), mmol.cell(), fakegrid);
            clipper::Xmap<float> observed_minus_model(mmol.spacegroup(), mmol.cell(), fakegrid);

            bool test = privateer::util::calculate_sigmaa_maps (withoutWaterModelAtomList,
                                                                fobs,
                                                                two_times_observed_minus_one_model,
                                                                observed_minus_model,
                                                                n_refln,
                                                                n_param);

           std::cout << "Test says: ";
           // << std::boolalpha << test << std::endl;

        }

    if (!batch)
    {
        std::cout << std::endl << " " << fobs.num_obs() << " reflections have been loaded";
        std::cout << std::endl << std::endl << " Resolution " << hklinfo.resolution().limit() << "Å" << std::endl << hklinfo.cell().format() << std::endl;
    }

<<<<<<< HEAD
    if ( check_unmodelled )
    {

      std::cout << "Scanning a waterless difference map for unmodelled glycosylation..." << std::endl;

      clipper::Grid_sampling mygrid( hklinfo.spacegroup(), hklinfo.cell(), hklinfo.resolution() );
      clipper::Xmap<float> sigmaa_all_map( hklinfo.spacegroup(), hklinfo.cell(), mygrid );
      clipper::Xmap<float> sigmaa_dif_map( hklinfo.spacegroup(), hklinfo.cell(), mygrid );
      bool no_errors = privateer::util::calculate_sigmaa_maps ( mmol.atom_list(),
                                                                fobs,
                                                                sigmaa_all_map,
                                                                sigmaa_dif_map );
    }
=======
>>>>>>> e3c93a73c6ca6f378c8f765ca8273574092045e2

    clipper::Atom_list mainAtoms;
    clipper::Atom_list ligandAtoms;
    clipper::Atom_list allAtoms;

    std::vector<std::pair< clipper::String , clipper::MSugar> > ligandList; // we store the Chain ID and create an MSugar to be scored
    std::vector<clipper::MMonomer> sugarList; // store the original MMonomer

    const clipper::MAtomNonBond& manb = clipper::MAtomNonBond( mmol, 1.0 ); // was 1.0

    mgl = clipper::MGlycology(mmol, manb, expsys);

    list_of_glycans = mgl.get_list_of_glycans();



    if ( !batch )
    {
        std::cout << std::endl << "Number of detected glycosylations: " << list_of_glycans.size();

        if ( list_of_glycans.size() > 0 )
        {
            clipper::String current_chain = "" ;

            for (int i = 0; i < list_of_glycans.size() ; i++ )
            {
                if ( current_chain != list_of_glycans[i].get_chain() )
                {
                    current_chain = list_of_glycans[i].get_chain();
                    std::cout << std::endl << std::endl << "Chain " << current_chain[0] << std::endl << "-------" << std::endl ;
                }
                std::cout << std::endl << list_of_glycans[i].print_linear ( true, false, true ) << std::endl;
                privateer::glycoplot::Plot plot(vertical, original, list_of_glycans[i].get_root_by_name(), invert, true);
                plot.plot_glycan ( list_of_glycans[i] );
                std::ostringstream os;
                os << list_of_glycans[i].get_root_for_filename() << ".svg";
                plot.write_to_file ( os.str() );
            }
        }

        std::cout << std::endl << "Analysing carbohydrates... ";
        fflush(0);
    }

    else
    {
        for (int i = 0; i < list_of_glycans.size() ; i++ )
        {
            privateer::glycoplot::Plot plot(vertical, original, list_of_glycans[i].get_root_by_name(), invert, true);
            plot.plot_glycan ( list_of_glycans[i] );
            std::ostringstream os;
            os << list_of_glycans[i].get_root_for_filename() << ".svg";
            plot.write_to_file ( os.str() );
        }

    }

    // erase ligand atoms from the model and then calculate phases

    for ( int p = 0; p < mmol.size(); p++ )
    {
        for ( int m = 0; m < mmol[p].size(); m++ )
        {
            if (allSugars)
            {
                if ( clipper::MDisaccharide::search_disaccharides(mmol[p][m].type().c_str()) != -1 ) // treat disaccharide
                {
                    clipper::MDisaccharide md(mmol, manb, mmol[p][m] );
                    sugarList.push_back ( mmol[p][m] );
                    sugarList.push_back ( mmol[p][m] );
                    clipper::String id = mmol[p].id();
                    id.resize(1);

                    ligandList.push_back ( std::pair < clipper::String, clipper::MSugar> (id, md.get_first_sugar()));
                    ligandList.push_back ( std::pair < clipper::String, clipper::MSugar> (id, md.get_second_sugar()));

                    if (( md.get_first_sugar().type_of_sugar() == "unsupported" ) || ( md.get_second_sugar().type_of_sugar() == "unsupported" ) )
                    {
                        std::cout << std::endl;
                        std::cout << "Error: strangely, at least one of the sugars in the supplied PDB file is missing required atoms. Stopping..." << std::endl << std::endl;
                        prog.set_termination_message( "Unrecoverable error" );
                        return 0;
                    }

                    for (int id = 0; id < mmol[p][m].size(); id++ )
                    {
                        ligandAtoms.push_back(mmol[p][m][id]);  // add the ligand atoms to a second array
                        allAtoms.push_back(mmol[p][m][id]);
                    }

                }
                else if ( !clipper::MSugar::search_database(mmol[p][m].type().c_str()) )
                {
                    for (int id = 0; id < mmol[p][m].size(); id++ )
                    {
                        mainAtoms.push_back(mmol[p][m][id]); // cycle through atoms and copy them
                        allAtoms.push_back(mmol[p][m][id]);
                    }
                }
                else // it's one of the sugars contained in the database
                {
                    clipper::MSugar msug, msug_b;

                    std::vector <char> conformers = privateer::util::number_of_conformers(mmol[p][m]);

                    #ifdef DUMP
                        std::cout << "number of alternate conformations: " << conformers.size() << std::endl;
                    #endif

                    int n_conf = conformers.size();

                    if ( n_conf > 0 )
                    {
                        if ( n_conf == 1 )
                            msug   = clipper::MSugar(mmol, mmol[p][m], manb, conformers[0]);
                        else
                        {
                            msug   = clipper::MSugar(mmol, mmol[p][m], manb, conformers[0]);
                            msug_b = clipper::MSugar(mmol, mmol[p][m], manb, conformers[1]);
                        }

                    }
                    else
                    {
                        msug = clipper::MSugar(mmol, mmol[p][m], manb);
                    }

                    sugarList.push_back(mmol[p][m]);
                    clipper::String id = mmol[p].id();
                    id.resize(1);

                    ligandList.push_back(std::pair<clipper::String, clipper::MSugar> (id, msug));
                    // add both conformers if the current monomer contains more than one
                    if ( n_conf == 2 )
                        ligandList.push_back(std::pair<clipper::String, clipper::MSugar> (id, msug_b));

                    if ( msug.type_of_sugar() == "unsupported" )
                    {
                        std::cout << std::endl << "Error: at least one of the sugars in the supplied PDB file is missing expected atoms"
                                  << std::endl << std::endl;
                        prog.set_termination_message( "Unrecoverable error" );
                        return 0;
                    }

                    for (int id = 0; id < mmol[p][m].size(); id++ )
                    {
                        ligandAtoms.push_back(mmol[p][m][id]);  // add the ligand atoms to a second array
                        allAtoms.push_back(mmol[p][m][id]);
                    }
                }

            }

            else
            {
                if ( strncmp( mmol[p][m].type().c_str(), ipcode.trim().c_str(), 3 )) // true if strings are different
                {
                    for (int id = 0; id < mmol[p][m].size(); id++ )
                    {
                        mainAtoms.push_back(mmol[p][m][id]); // cycle through atoms and copy them
                        allAtoms.push_back(mmol[p][m][id]);
                    }
                }
                else // it's the one sugar we're looking to omit
                {
                    if ( validation_options.size() > 0 )
                    {
                        const clipper::MSugar msug ( mmol, mmol[p][m], manb, external_validation );

                        sugarList.push_back(mmol[p][m]);
                        ligandList.push_back(std::pair<clipper::String, clipper::MSugar> (mmol[p].id().trim(), msug));
                    }
                    else
                    {
                        const clipper::MSugar msug(mmol, mmol[p][m], manb);

                        sugarList.push_back(mmol[p][m]);
                        ligandList.push_back(std::pair<clipper::String, clipper::MSugar> (mmol[p].id().trim(), msug));

                    }
                    for (int id = 0; id < mmol[p][m].size(); id++ )
                    {
                        ligandAtoms.push_back(mmol[p][m][id]);  // add the ligand atoms to a second array
                        allAtoms.push_back(mmol[p][m][id]);
                    }
                }
            }
        }
    }

    if (!batch) std::cout << "done.\nCalculating structure factors with bulk solvent correction... "; fflush(0);

    // calculate structure factors

    clipper::SFcalc_obs_bulk<float> sfcbligands;
    clipper::SFcalc_obs_bulk<float> sfcb;
    clipper::SFcalc_obs_bulk<float> sfcball;

    try
    {   // calculate structure factors with bulk solvent correction
#pragma omp parallel sections
        {
#pragma omp section
            sfcbligands( fc_ligands_bsc, fobs, ligandAtoms ); // was fobs_scaled
#pragma omp section
            sfcb( fc_omit_bsc, fobs, mainAtoms );  // calculation of omit SF with bulk solvent correction
#pragma omp section
            sfcball( fc_all_bsc, fobs, allAtoms ); // calculation of SF with bulk solvent correction
        }
    }
    catch ( ... )
    {
        if (!batch) std::cout << "\nThe input file has unrecognised atoms. Might cause unexpected results...\n";  // this causes clipper to freak out, so better remove those unknowns
    }

    fc_ligands_bsc[0].set_null();
    fc_omit_bsc[0].set_null();
    fc_all_bsc[0].set_null();

    if (!batch)
    {
        std::cout << "done." << std::endl << "Computing 2mFo-DFc, mFo-DFc and mFo-omit_DFc maps... ";
        fflush(0);
    }

    clipper::Grid_sampling mygrid( hklinfo.spacegroup(), hklinfo.cell(), hklinfo.resolution() );  // define grid
    clipper::Xmap<float> sigmaa_all_map( hklinfo.spacegroup(), hklinfo.cell(), mygrid );          // define sigmaa best map
    clipper::Xmap<float> sigmaa_dif_map( hklinfo.spacegroup(), hklinfo.cell(), mygrid );          // define sigmaa diff  map
    clipper::Xmap<float> sigmaa_omit_fd( hklinfo.spacegroup(), hklinfo.cell(), mygrid );          // define sigmaa omit diff map
    clipper::Xmap<float> ligandmap( hklinfo.spacegroup(), hklinfo.cell(), mygrid );

    // scale data and flag R-free

    HRI ih;
    clipper::HKL_data<Flag> flag( hklinfo );     // same flag for both calculations, omit absent reflections
    clipper::SFscale_aniso<float> sfscale;

#pragma omp parallel sections
    {
#pragma omp section
        {
            sfscale( fobs_scaled, fc_all_bsc );  // anisotropic scaling of Fobs. We scale Fobs to Fcalc instead of scaling our 3 Fcalcs to Fobs
        }
#pragma omp section
        {
            for ( ih = flag.first(); !ih.last(); ih.next() ) // we want to use all available reflections
            {
                if ( !fobs_scaled[ih].missing() ) flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
                else flag[ih].flag() = clipper::SFweight_spline<float>::NONE;
            }
        }
    }

    double FobsFcalcSum = 0.0;
    double FobsFcalcAllSum = 0.0;
    double FobsSum = 0.0;

    clipper::HKL_data<F_phi> fb_omit( hklinfo ); // new variables for omit sigmaa weighting calculation
    clipper::HKL_data<F_phi> fd_omit( hklinfo );
    clipper::HKL_data<Phi_fom> phiw_omit( hklinfo );
    clipper::HKL_data<F_phi> fb_all( hklinfo ); // variables for all atom sigmaa weighting calculation
    clipper::HKL_data<F_phi> fd_all( hklinfo );
    clipper::HKL_data<Phi_fom> phiw_all( hklinfo );

    // now do sigmaa calc
#pragma omp parallel sections
    {
#pragma omp section
        {
            clipper::SFweight_spline<float> sfw_omit (n_refln, n_param );
            sfw_omit( fb_omit, fd_omit, phiw_omit, fobs_scaled, fc_omit_bsc, flag ); // sigmaa omit
        }

#pragma omp section
        {
            clipper::SFweight_spline<float> sfw_all( n_refln, n_param );
            sfw_all( fb_all,  fd_all,  phiw_all,  fobs_scaled, fc_all_bsc,  flag ); // sigmaa all atoms
        }
    }

    // fb:          output best map coefficients
    // fd:          output difference map coefficients
    // phiw:        output phase and fom
    // fobs_scaled: input observed structure factors, previously scaled
    // fc_omit_bsc: input calculated omit data, with bulk solvent correction
    // fc_all_bsc:  input calculated data, with bsc

    std::vector<double> params( n_param, 2.0 );
    clipper::BasisFn_spline wrk_basis( hklinfo, n_param, 2.0 );

    clipper::TargetFn_scaleF1F2<F_phi,F_sigF> wrk_target_omit( fc_omit_bsc, fobs_scaled ); // was just fobs
    clipper::TargetFn_scaleF1F2<F_phi,F_sigF> wrk_target_all ( fc_all_bsc, fobs_scaled );
    clipper::ResolutionFn wrk_scale_omit( hklinfo, wrk_basis, wrk_target_omit, params );
    clipper::ResolutionFn wrk_scale_all ( hklinfo, wrk_basis, wrk_target_all,  params );

    double Fo, Fc_all, Fc_omit;

#pragma omp parallel sections
    {
#pragma omp section
        sigmaa_all_map.fft_from( fb_all );  // calculate the maps
#pragma omp section
        sigmaa_dif_map.fft_from( fd_all );
#pragma omp section
        sigmaa_omit_fd.fft_from( fd_omit );
#pragma omp section
        ligandmap.fft_from( fc_ligands_bsc );       // this is the map that will serve as Fc map for the RSCC calculation
#pragma omp section
        for ( HRI ih = fobs_scaled.first(); !ih.last(); ih.next() )
        {
            if ( !fobs_scaled[ih].missing() )
            {
                Fo = fobs_scaled[ih].f();
                Fc_all = sqrt ( wrk_scale_all.f(ih) ) * fc_all_bsc[ih].f() ;
                Fc_omit = sqrt ( wrk_scale_omit.f(ih) ) * fc_omit_bsc[ih].f() ;
                FobsFcalcSum += fabs( Fo - Fc_omit); // R factor calculation
                FobsFcalcAllSum += fabs( Fo- Fc_all);
                FobsSum += Fo;
            }
        }
    }

    if (!batch)
        std::cout << "done." << std::endl;

    if ( output_mtz )
    {
        if (!batch)
        {
            std::cout << "Writing map coefficients to " << opfile << "... ";
            fflush(0);
        }
        if (useMTZ)
        {
            clipper::CCP4MTZfile mtzout;
            mtzout.open_append(ipmtz, opfile );
            mtzout.export_hkl_data( fb_all, "*/*/BEST" );
            mtzout.export_hkl_data( fd_all, "*/*/DIFF" );
            mtzout.export_hkl_data( fd_omit,"*/*/OMIT");
            mtzout.close_append();

            if (!batch)
                std::cout << "done" << std::endl;
        }
        else
            if (!batch)
                std::cout << "skipped. You must supply an input MTZ from which columns can be read and transferred to the output MTZ." << std::endl << std::endl;
    }

    if (batch) // create miniMTZ files for ccp4i2
    {
        std::cout << opxtal.crystal_name() << " " << opxtal.project_name() << " " << opdset.dataset_name() << " " << opdset.wavelength() << std::endl;
        clipper::String path = "/" + opxtal.crystal_name() + "/" + opdset.dataset_name() + "/[F,PHI]";

        opmtz_best.open_write ( "FPHIOUT.mtz" );
        opmtz_best.export_hkl_info ( hklinfo );
        opmtz_best.export_crystal ( opxtal, path );
        opmtz_best.export_dataset ( opdset, path );
        opmtz_best.export_hkl_data ( fb_all, path );
        opmtz_best.close_write ();

        opmtz_omit.open_write ( "OMITFPHIOUT.mtz" );
        opmtz_omit.export_hkl_info ( hklinfo );
        opmtz_omit.export_crystal ( opxtal, path );
        opmtz_omit.export_dataset ( opdset, path );
        opmtz_omit.export_hkl_data ( fd_omit, path );
        opmtz_omit.close_write ();
    }



    if (!batch)
        printf("\n R-all = %1.3f  R-omit = %1.3f\n", (FobsFcalcAllSum / FobsSum), (FobsFcalcSum / FobsSum));

    if (!batch)
        if (((FobsFcalcAllSum / FobsSum)*10) > hklinfo.resolution().limit() + 0.6)
            std::cout << " Warning: R-work is unusually high. Please ensure that your PDB file contains full B-factors instead of residuals after TLS refinement!" << std::endl;

    float difference = (FobsFcalcSum / FobsSum) - (FobsFcalcAllSum / FobsSum);

    if (( difference > 0.15 ) || (clipper::Util::is_nan((FobsFcalcSum / FobsSum))))
    {
        useSigmaa = true;

        if (!batch)
            std::cout << std::endl << " The studied portions of the model account for a very significant part of the data. Calculating RSCC against a regular 2mFo-DFc map" << std::endl;
    }

    if (!batch)
    {
        std::cout << "\nWriting maps to disk... ";
        fflush(0);
    }

    clipper::CCP4MAPfile sigmaa_all_MapOut;
    clipper::CCP4MAPfile sigmaa_dif_MapOut;
    clipper::CCP4MAPfile sigmaa_omit_fd_MapOut;

    if (allSugars)
        ipcode = "all";

    clipper::Map_stats ms;

    if (useSigmaa)
        ms = clipper::Map_stats(sigmaa_all_map);
    else
        ms = clipper::Map_stats(sigmaa_omit_fd);

    if (!batch)
    {
#pragma omp parallel sections
        {
#pragma omp section
            {
                sigmaa_all_MapOut.open_write( "sigmaa_best.map" );      // write maps
                sigmaa_all_MapOut.export_xmap( sigmaa_all_map );
                sigmaa_all_MapOut.close_write();
            }
#pragma omp section
            {
                sigmaa_dif_MapOut.open_write( "sigmaa_diff.map" );
                sigmaa_dif_MapOut.export_xmap( sigmaa_dif_map );
                sigmaa_dif_MapOut.close_write();
            }
#pragma omp section
            {
                sigmaa_omit_fd_MapOut.open_write( "sigmaa_omit.map" );
                sigmaa_omit_fd_MapOut.export_xmap( sigmaa_omit_fd );
                sigmaa_omit_fd_MapOut.close_write();
            }
        }

        std::cout << "done." << std::endl;
        std::cout << "\n\nDetailed validation data" << std::endl;
        std::cout << "------------------------" << std::endl;
    }

    if (!batch)
        printf("\nPDB \t    Sugar   \tRsln\t  Q  \t Phi  \tTheta \tRSCC\t   Detected type   \tCnf\t<mFo>\t<Bfac>\tCtx\t Ok?");
    if (!batch && showGeom)
        printf("\tBond lengths, angles and torsions, reported clockwise with in-ring oxygen as first vertex");
    if (!batch)
        printf("\n----\t------------\t----\t-----\t------\t------\t----\t-------------------\t---\t-----\t------\t---\t-----");
    if (!batch && showGeom)
        printf("\t------------------------------------------------------------------------------------------------------------");
    if (!batch)
        printf("\n");

    for (int index = 0; index < ligandList.size(); index++)
    {
        float x,y,z,maxX,maxY,maxZ,minX,minY,minZ;
        x=y=z=0.0;
        maxX=maxY=maxZ=-999999.0;
        minX=minY=minZ=999999.0;

        for (int natom = 0; natom < sugarList[index].size(); natom++)
        {
            if(sugarList[index][natom].coord_orth().x() > maxX) maxX=sugarList[index][natom].coord_orth().x(); // calculation of the sugar centre
            if(sugarList[index][natom].coord_orth().y() > maxY) maxY=sugarList[index][natom].coord_orth().y();
            if(sugarList[index][natom].coord_orth().z() > maxZ) maxZ=sugarList[index][natom].coord_orth().z();
            if(sugarList[index][natom].coord_orth().x() < minX) minX=sugarList[index][natom].coord_orth().x();
            if(sugarList[index][natom].coord_orth().y() < minY) minY=sugarList[index][natom].coord_orth().y();
            if(sugarList[index][natom].coord_orth().z() < minZ) minZ=sugarList[index][natom].coord_orth().z();
        }

        x = minX + ((maxX - minX)/2);
        y = minY + ((maxY - minY)/2);
        z = minZ + ((maxZ - minZ)/2);

        if (batch)
        {
            fprintf(output, "%c%c%c%c\t%s-",ippdb[1+pos_slash],ippdb[2+pos_slash],ippdb[3+pos_slash],ippdb[4+pos_slash], ligandList[index].second.type().trim().c_str());
            fprintf(output, "%s-%s   ", ligandList[index].first.c_str(), ligandList[index].second.id().trim().c_str());
        }
        else
        {
            printf("%c%c%c%c\t%s-",ippdb[1+pos_slash],ippdb[2+pos_slash],ippdb[3+pos_slash],ippdb[4+pos_slash], ligandList[index].second.type().c_str());
            std::cout << ligandList[index].first << "-" << ligandList[index].second.id().trim() << "  ";
        }

        // now calculate the correlation between the weighted experimental & calculated maps
        // maps are scanned only inside a sphere containing the sugar for performance reasons,
        // although RSCC and <RMS> are restricted to a mask surrounding the model

        double meanDensityExp, meanDensityCalc, num, den1, den2, corr_coeff;
        meanDensityCalc = meanDensityExp = num = den1 = den2 = corr_coeff = 0.0;

        int n_points = 0;

        //////// mask calculation //////////

        clipper::Xmap<float> mask( hklinfo.spacegroup(), hklinfo.cell(), mygrid );

        clipper::EDcalc_mask<float> masker( ipradius );
        masker(mask, sugarList[index].atom_list());

        ////////////////////////////////////

        clipper::Coord_orth origin(minX-2,minY-2,minZ-2);
        clipper::Coord_orth destination(maxX+2,maxY+2,maxZ+2);

        clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;

        double accum = 0.0;

        // calculation of the mean densities of the calc (ligandmap) and weighted obs (sigmaamap) maps

        std::vector<clipper::Xmap_base::Map_reference_coord> buffer_coord;

        if (useSigmaa)
            i0 = clipper::Xmap_base::Map_reference_coord( sigmaa_all_map, origin.coord_frac(hklinfo.cell()).coord_grid(mygrid) );
        else
            i0 = clipper::Xmap_base::Map_reference_coord( sigmaa_omit_fd, origin.coord_frac(hklinfo.cell()).coord_grid(mygrid) );

        for ( iu = i0; iu.coord().u() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).u(); iu.next_u() )
            for ( iv = iu; iv.coord().v() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).v(); iv.next_v() )
                for ( iw = iv; iw.coord().w() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).w(); iw.next_w() )
                {
                    if ( mask[iw] == 1.0)
                    {
                        meanDensityCalc = meanDensityCalc + ligandmap[iw];

                        if (useSigmaa)
                            meanDensityExp = meanDensityExp + sigmaa_all_map[iw];
                        else
                            meanDensityExp = meanDensityExp + sigmaa_omit_fd[iw];

                        n_points++;
                    }
                }

        accum = meanDensityExp / ms.std_dev();
        accum /= n_points;

        meanDensityCalc = meanDensityCalc / n_points;
        meanDensityExp = meanDensityExp / n_points;

        // calculation of the correlation coefficient between calc (ligandmap) and weighted obs (sigmaamap) maps

        if (useSigmaa)
            i0 = clipper::Xmap_base::Map_reference_coord( sigmaa_all_map, origin.coord_frac(hklinfo.cell()).coord_grid(mygrid) );
        else
            i0 = clipper::Xmap_base::Map_reference_coord( sigmaa_omit_fd, origin.coord_frac(hklinfo.cell()).coord_grid(mygrid) );

        for ( iu = i0; iu.coord().u() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).u(); iu.next_u() )
            for ( iv = iu; iv.coord().v() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).v(); iv.next_v() )
                for ( iw = iv; iw.coord().w() <= destination.coord_frac(hklinfo.cell()).coord_grid(mygrid).w(); iw.next_w() )
                {
                    if ( mask[iw] == 1.0)
                    {
                        if (useSigmaa)
                        {
                            num = num + (sigmaa_all_map[iw] - meanDensityExp) * (ligandmap[iw] - meanDensityCalc);
                            den1 = den1 + pow((sigmaa_all_map[iw] - meanDensityExp),2);
                            den2 = den2 + pow((ligandmap[iw] - meanDensityCalc),2);
                        }
                        else
                        {
                            num = num + (sigmaa_omit_fd[iw] - meanDensityExp) * (ligandmap[iw] - meanDensityCalc);
                            den1 = den1 + pow((sigmaa_omit_fd[iw] - meanDensityExp),2);
                            den2 = den2 + pow((ligandmap[iw] - meanDensityCalc),2);
                        }
                    }
                }

        corr_coeff = num / (sqrt(den1) * sqrt(den2));

        ///////////// here we deal with the sugar /////////////

        if (batch)
        {
            std::vector<clipper::ftype> cpParams(10, 0);
            cpParams = ligandList[index].second.cremer_pople_params();
            fprintf(output,"\t%1.2f\t%1.3f\t%3.2f\t",hklinfo.resolution().limit(),cpParams[0],cpParams[1] );    // output cremer-pople parameters
            if ( cpParams[2] == -1 ) fprintf ( output, " --  \t" ); else fprintf ( output, "%3.2f\t", cpParams[2] );
            fprintf(output,"%1.2f\t", corr_coeff);                                              // output RSCC and data resolution
            fprintf(output,"%s\t", ligandList[index].second.type_of_sugar().c_str());           // output the type of sugar, e.g. alpha-D-aldopyranose
            fprintf(output,"%s\t", ligandList[index].second.conformation_name().c_str());       // output a 3 letter code for the conformation
            fprintf(output,"%1.3f \t", accum);
            ligandList[index].second.set_rscc ( corr_coeff );

            float bfac = ligandList[index].second.get_bfactor ();

            fprintf ( output, "%3.2f", bfac ); // output <bfactor>

            std::vector < clipper::MGlycan > list_of_glycans = mgl.get_list_of_glycans();
            bool found_in_tree = false;

            for ( int i = 0 ; i < list_of_glycans.size() ; i++ )
            {
                std::vector < clipper::MSugar > list_of_sugars = list_of_glycans[i].get_sugars();

                for ( int j = 0 ; j < list_of_sugars.size() ; j++ )
                {
                    if ( list_of_sugars[j].id().trim() == ligandList[index].second.id().trim() )
                    {
                        if ( list_of_glycans[i].get_type() == "n-glycan" )
                        {
                            ligandList[index].second.set_context ( "n-glycan" );
                            fprintf ( output, "\t(n)");
                        }
                        else
                        {
                            ligandList[index].second.set_context ( "o-glycan" );
                            fprintf ( output, "\t(o)");
                        }
                        found_in_tree = true;
                        break;
                    }
                }
                if ( found_in_tree ) break;
            }

            if ( !found_in_tree )
            {
                ligandList[index].second.set_context ( "ligand" );
                fprintf ( output, "\t(l)");
            }

            if (ligandList[index].second.in_database(ligandList[index].second.type().trim()))
            {
                if ((ligandList[index].second.ring_members().size() == 6 ))
                {
                    if (ligandList[index].second.is_sane())
                    {
                        if ( ! ligandList[index].second.ok_with_conformation () )
                        {
                            fprintf(output, "\tcheck");
                        }
                        else fprintf(output, "\tyes");
                    }
                    else
                        fprintf (output, "\tno");
                }
                else
                    if (ligandList[index].second.is_sane())
                        fprintf(output, "\tyes");
                    else
                    {
                        fprintf(output, "\tno");
                    }
            }
            else
                fprintf(output, "\tunk");

            if ( ! ligandList[index].second.ok_with_conformation () )
                enable_torsions_for.push_back (ligandList[index].second.type().trim());

            bool occupancy_check = false;
            std::vector<clipper::MAtom> ringcomponents = ligandList[index].second.ring_members();

            for ( int i = 0 ; i < ringcomponents.size() ; i++ )
                if (privateer::util::get_altconformation(ringcomponents[i]) != ' ')
                    occupancy_check = true;


            if (showGeom)
            {
                std::vector<clipper::ftype> rangles = ligandList[index].second.ring_angles();
                std::vector<clipper::ftype> rbonds  = ligandList[index].second.ring_bonds();
                std::vector<clipper::ftype> rtorsions = ligandList[index].second.ring_torsions();

                for (int i = 0 ; i < ligandList[index].second.ring_members().size(); i++ )
                    fprintf(output, "\t%1.2f", rbonds[i]);

                for (int i = 0 ; i < ligandList[index].second.ring_members().size(); i++ )
                    fprintf(output, "\t%3.1f", rangles[i]);

                for (int i = 0 ; i < ligandList[index].second.ring_members().size(); i++ )
                    fprintf(output, "\t%3.1f", rtorsions[i]);
            }

            if (occupancy_check)
                fprintf(output, " (*)");

            fprintf(output, "\n");

        }
        else
        {
            std::vector<clipper::ftype> cpParams(10, 0);
            cpParams = ligandList[index].second.cremer_pople_params();
            printf("\t%1.2f\t%1.3f\t%3.2f\t",hklinfo.resolution().limit(),cpParams[0],cpParams[1]);             // output cremer-pople parameters
            if ( cpParams[2] == -1 ) printf ( " --  \t" ); else printf ( "%3.2f\t", cpParams[2] );
            printf("%1.2f\t", corr_coeff);                                                                                              // output RSCC and data resolution
            printf("%s\t", ligandList[index].second.type_of_sugar().c_str());                   // output the type of sugar, e.g. alpha-D-aldopyranose
            printf("%s\t", ligandList[index].second.conformation_name().c_str());               // output a 3 letter code for the conformation
            printf("%1.3f \t", accum);                                                                                                  // output <mFo>
            ligandList[index].second.set_rscc ( corr_coeff );

            float bfac = 0.0;

            for (int i=0; i < ligandList[index].second.size(); i++)
                bfac+=ligandList[index].second[i].u_iso();

            bfac /= ligandList[index].second.size();
            bfac  = clipper::Util::u2b(bfac);

            printf ( "%3.2f", bfac );                 // output <Bfactor>

            std::vector < clipper::MGlycan > list_of_glycans = mgl.get_list_of_glycans();
            bool found_in_tree = false;

            for ( int i = 0 ; i < list_of_glycans.size() ; i++ )
            {
                std::vector < clipper::MSugar > list_of_sugars = list_of_glycans[i].get_sugars();

                for ( int j = 0 ; j < list_of_sugars.size() ; j++ )
                {
                    if ( list_of_sugars[j].id().trim() == ligandList[index].second.id().trim() )
                    {
                        if ( list_of_glycans[i].get_type() == "n-glycan" )
                        {
                            ligandList[index].second.set_context ( "n-glycan" );
                            std::cout << "\t(n) ";
                        }
                        else
                        {
                            std::cout << "\t(o) ";
                            ligandList[index].second.set_context ( "o-glycan" );
                        }
                        found_in_tree = true;
                        break;
                    }
                }
                if ( found_in_tree ) break;
            }

            if ( !found_in_tree )
            {
                ligandList[index].second.set_context ( "ligand" );
                std::cout << "\t(l) ";
            }

            if (ligandList[index].second.in_database(ligandList[index].second.type().trim()))
            {
                if ((ligandList[index].second.ring_members().size() == 6 ))
                {
                    if (ligandList[index].second.is_sane())
                    {
                        if ( ! ligandList[index].second.ok_with_conformation () )
                            printf("\tcheck");
                        else
                            printf("\tyes");
                    }
                    else
                        printf ("\tno");
                }
                else
                    if (ligandList[index].second.is_sane())
                        printf("\tyes");
                    else printf("\tno");
            }
            else
                printf("\tunk");

            if ( ! ligandList[index].second.ok_with_conformation () )
                enable_torsions_for.push_back (ligandList[index].second.type().trim());

            bool occupancy_check = false;
            std::vector<clipper::MAtom> ringcomponents = ligandList[index].second.ring_members();

            for ( int i = 0 ; i < ringcomponents.size() ; i++ )
                if (privateer::util::get_altconformation(ringcomponents[i]) != ' ')
                    occupancy_check = true;


            if (showGeom)
            {
                std::vector<clipper::ftype> rangles = ligandList[index].second.ring_angles();
                std::vector<clipper::ftype> rbonds  = ligandList[index].second.ring_bonds();
                std::vector<clipper::ftype> rtorsions = ligandList[index].second.ring_torsions();

                for (int i = 0 ; i < ligandList[index].second.ring_members().size(); i++ )
                    printf("\t%1.2f", rbonds[i]);

                for (int i = 0 ; i < ligandList[index].second.ring_members().size(); i++ )
                    printf("\t%3.1f", rangles[i]);

                for (int i = 0 ; i < ligandList[index].second.ring_members().size(); i++ )
                    printf("\t%3.1f", rtorsions[i]);
            }

            if (occupancy_check)
                std::cout << " (*)";

            std::cout << std::endl;
        }
    }

    if (!batch)
    {
        printf("\n");
        std::cout << "\nPartially occupied monosaccharides, if any, are marked with an asterisk (*)";
        std::cout << std::endl << std::endl;
    }
    else
        fclose(output);

    std::fstream of_scm; std::fstream of_py;
    of_scm.open("privateer-results.scm", std::fstream::out);
    of_py.open("privateer-results.py", std::fstream::out);

    privateer::coot::insert_coot_prologue_scheme ( of_scm );
    privateer::coot::insert_coot_prologue_python ( of_py );
    privateer::coot::insert_coot_files_loadup_scheme (of_scm, ippdb, "sigmaa_best.map", "sigmaa_diff.map", "sigmaa_omit.map", batch );
    privateer::coot::insert_coot_files_loadup_python (of_py,  ippdb, "sigmaa_best.map", "sigmaa_diff.map", "sigmaa_omit.map", batch );

    int n_geom = 0, n_anomer = 0, n_config = 0, n_pucker = 0, n_conf = 0;

    int sugar_count = 0;

    for (int k = 0 ; k < ligandList.size() ; k++)
    {
        int n_errors = 0;

        clipper::String diagnostic = ligandList[k].second.type().trim() + "/" + ligandList[k].first + "/" + ligandList[k].second.id().trim();
        clipper::String report = "";
        clipper::String sugarRSCC = clipper::String( ligandList[k].second.get_rscc() );
        sugarRSCC.resize(4);
        diagnostic.append ( ": RSCC=" + sugarRSCC + ", " );

        if (ligandList[k].second.is_supported() )
        {
            if ( ! ligandList[k].second.is_sane() )
            {
                sugar_count++;

                if (!ligandList[k].second.ok_with_anomer())
                {
                    if ( n_errors > 0 )
                    {
                        diagnostic.append(", wrong anomer");
                        report.append(", wrong anomer");
                    }
                    else
                    {
                        diagnostic.append("wrong anomer");
                        report.append("Wrong anomer");
                    }

                    n_errors++;
                    n_anomer++;
                }

                if (!ligandList[k].second.ok_with_chirality())
                {
                    if ( n_errors > 0 )
                    {
                        diagnostic.append(", wrong configuration at " + ligandList[k].second.configurational_carbon().name().trim());
                        report.append(", wrong configuration at " + ligandList[k].second.configurational_carbon().name().trim());
                    }
                    else
                    {
                        diagnostic.append("wrong configuration at " + ligandList[k].second.configurational_carbon().name().trim());
                        report.append("Wrong configuration at " + ligandList[k].second.configurational_carbon().name().trim());
                    }
                    n_errors++;
                    n_config++;
                }

                if ( !ligandList[k].second.ok_with_puckering() )
                {
                    if ( n_errors > 0 )
                    {
                        diagnostic.append(", Q=" + clipper::String(ligandList[k].second.puckering_amplitude()) );
                        report.append(", Q=" + clipper::String(ligandList[k].second.puckering_amplitude()) );
                    }
                    else
                    {
                        diagnostic.append("Q=" + clipper::String( ligandList[k].second.puckering_amplitude() ));
                        report.append("Q=" + clipper::String(ligandList[k].second.puckering_amplitude()) );
                    }

                    n_errors++;
                    n_pucker++;
                }

                if ( !ligandList[k].second.ok_with_conformation() )
                {
                    if ( n_errors > 0 )
                    {
                        diagnostic.append(", conformation (" + ligandList[k].second.conformation_name() + clipper::String(") might be mistaken") );
                        report.append(", conformation might be mistaken");
                    }
                    else
                    {
                        diagnostic.append("conformation (" + ligandList[k].second.conformation_name() + clipper::String(") might be mistaken"));
                        report.append("Conformation might be mistaken");
                    }

                    n_errors++;
                    n_conf++;
                }

                privateer::coot::insert_coot_go_to_sugar_scheme ( of_scm, ligandList[k].second.ring_centre(), diagnostic );
                privateer::coot::insert_coot_go_to_sugar_python ( of_py, ligandList[k].second.ring_centre(), diagnostic );
            }
            else // sugar is sane, but still need to check higher-energy conformations
            {
                if ( !ligandList[k].second.ok_with_conformation() )
                {
                    if ( n_errors > 0 )
                    {
                        diagnostic.append(", conformation (" + ligandList[k].second.conformation_name() + clipper::String(") might be mistaken") );
                        report.append(", conformation might be mistaken");
                    }
                    else
                    {
                        diagnostic.append("conformation (" + ligandList[k].second.conformation_name() + clipper::String(") might be mistaken"));
                        report.append("Conformation might be mistaken");
                    }

                    n_errors++;
                    n_conf++;

                    sugar_count++;
                    privateer::coot::insert_coot_go_to_sugar_scheme ( of_scm, ligandList[k].second.ring_centre(), diagnostic );
                    privateer::coot::insert_coot_go_to_sugar_python ( of_py, ligandList[k].second.ring_centre(), diagnostic );
                }

            }
        }

        if ( report == "" )
            report = "Ok";

        ligandList[k].second.set_diagnostic ( report );

    }

    privateer::coot::insert_coot_epilogue_scheme ( of_scm );
    privateer::coot::insert_coot_epilogue_python ( of_py );
    of_scm.close();
    of_py.close();

    std::cout << "SUMMARY: " << std::endl << std::endl ;
    std::cout << "   Wrong anomer: " << n_anomer << std::endl;
    std::cout << "   Wrong configuration: " << n_config << std::endl;
    std::cout << "   Unphysical puckering amplitude: " << n_pucker << std::endl;
    std::cout << "   In higher-energy conformations: " << n_conf << std::endl;
    std::cout << std::endl;
    std::cout << "   Privateer has identified " << n_anomer + n_config + n_pucker + n_conf;
    std::cout << " issues, with " << sugar_count << " of " << ligandList.size() << " sugars affected." << std::endl;

    privateer::util::print_XML(ligandList, list_of_glycans, ippdb);

    if ( enable_torsions_for.size() > 0 )
    {
        privateer::util::write_refmac_keywords ( enable_torsions_for );
        privateer::util::write_libraries( enable_torsions_for );
    }

    prog.set_termination_message( "Normal termination" );
    system("touch scored");
    return 0;
}
