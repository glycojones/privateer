
// Clipper app to score a target cyclic carbohydrate
// See below for a complete description of the calculations
// 2013 Jon Agirre & Kevin Cowtan @ University of York
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk
//
// Targets supported: structures containing furanose and pyranose rings
//
//  Currently, the application calculates:
//   - SFcalc from model (internal)
//   - Fo-Fc sigmaa weighted omit maps (output)
//   - Fo-Fc omit maps for comparison (output)
//   - Fc maps for RSCC calculation (internal)
//   - Correlation coefficient + score
//   - Cremer-pople puckering parameters (q2,phi2) + q3 for pyranose rings, (q2,phi2) for furanose rings

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
#include <clipper/clipper.h>
#include <clipper/clipper-cif.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <clipper/contrib/sfcalc_obs.h>
#include <clipper/minimol/minimol_utils.h>


using clipper::data32::F_sigF;
using clipper::data32::F_phi;
using clipper::data32::Phi_fom;
using clipper::data32::Flag;
typedef clipper::HKL_data_base::HKL_reference_index HRI;
using namespace std;


// begin helper function declarations //

void print_usage ();

void print_supported_code_list ();

void print_XML ( std::vector < std::pair < clipper::String, clipper::MSugar > >, std::vector < clipper::MGlycan > list_of_glycans, clipper::String pdbname );

bool compute_and_print_external_validation ( const std::vector<clipper::String> validation_options, clipper::data::sugar_database_entry& external_validation );

bool read_coordinate_file (clipper::MMDBfile& mfile, clipper::MiniMol& mmol, clipper::String& ippdb, bool batch);

char get_altconformation(clipper::MAtom ma);

void write_refmac_keywords ( std::vector < clipper::String > code_list );

// end helper function declarations //



int main(int argc, char** argv)
{
    
#ifndef SVN_REV
    clipper::String program_version = "MKII";
    CCP4Program prog( "privateer-validate", program_version.c_str(), "$Date: 2015/07/10" );
#else
    clipper::String program_version = "MKII-";
    program_version.append(SVN_REV);
    CCP4Program prog( "prval", program_version.c_str(), "$Date: 2015/07/10" );
#endif
    
    prog.set_termination_message( "Failed" );
    
    std::cout << std::endl << "Copyright 2013-2015 Jon Agirre, Kevin Cowtan and The University of York." << std::endl;
    
    clipper::HKL_info hklinfo; // allocate space for the hkl metadata
    clipper::CIFfile cifin;
    clipper::CCP4MTZfile mtzin, ampmtzin;
    clipper::String ippdb    = "NONE";
    clipper::String ipcol_fo = "NONE";
    clipper::String ipcif    = "NONE";
    clipper::String ipcode   = "XXX";
    clipper::String opfile   = "privateer-validate.mtz";
    clipper::String title    = "generic title";
    clipper::String ipmtz    = "NONE";
    clipper::String validation_string = "";
    std::vector<clipper::String> validation_options;
    clipper::data::sugar_database_entry external_validation;
    bool useSigmaa = false;
    int n_refln = 1000;
    int n_param = 20;
    bool useMTZ = false;
    bool batch = false;
    bool noMaps = false;
    bool allSugars = true;
    bool showGeom = false;
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
                ipcif = args[arg];
        }
        else if ( args[arg] == "-list" )
        {
            print_supported_code_list();
            prog.set_termination_message( "Success" );
            return 0;
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
                if ( compute_and_print_external_validation ( validation_options, external_validation ) )
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
            std::cout << std::endl << std::endl << "Error: no internal validation data found for " << ipcode << std::endl;
            std::cout << "\nYou can provide external validation data with -valstring <data>" << std::endl
            << "\n\tAccepted format: SUG,O5/C1/C2/C3/C4/C5,A,D,4c1\n"
            << "\tThree-letter code, ring atoms separated by /, anomer, handedness, expected conformation.\n" << std::endl;
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
    
    if ( (ippdb != "NONE") && ((ipcif == "NONE") && (ipmtz == "NONE")) )
        noMaps = true;
    
    if ( (ippdb == "NONE") || ((ipcif == "NONE") && (ipmtz == "NONE") && (noMaps == false)) )
    {
        print_usage();
        prog.set_termination_message( "Failed" );
        return 1;
    }
    
    clipper::MMDBfile mfile;
    clipper::MiniMol mmol;
    
    read_coordinate_file ( mfile, mmol, ippdb, batch);
    int pos_slash = ippdb.rfind("/");
    
    
    
    
    
    
    // Fast mode, no maps nor correlation calculations
    
    if ( noMaps )
    {
        clipper::Atom_list mainAtoms;
        clipper::Atom_list ligandAtoms;
        clipper::Atom_list allAtoms;
        std:vector<clipper::String> enable_torsions_for;
        
        if (!batch)
        {
            std::cout << std::endl << "Analysing carbohydrates... \n";
            fflush(0);
        }
        
        std::vector < std::pair <clipper::String , clipper::MSugar> > ligandList; // we store the Chain ID and create an MSugar to be scored
        std::vector < clipper::MMonomer > sugarList; // store the original MMonomer
        
        const clipper::MAtomNonBond& manb = clipper::MAtomNonBond( mmol, 1.0 ); // was 1.0
        
        mgl = clipper::MGlycology(mmol, manb);
        
        list_of_glycans = mgl.get_list_of_glycans();
        
        if ( !batch ) std::cout << std::endl << "Number of detected glycans: " << list_of_glycans.size();
        
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
            }
        }
        
        std::cout << "\n\nDetailed validation data" << std::endl;
        std::cout << "------------------------" << std::endl;
        
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
                        clipper::MSugar msug(mmol, mmol[p][m], manb);
                        
                        sugarList.push_back(mmol[p][m]);
                        clipper::String id = mmol[p].id();
                        id.resize(1);
                        ligandList.push_back(std::pair<clipper::String, clipper::MSugar> (id, msug));
                        
                        if ( msug.type_of_sugar() == "unsupported" )
                        {
                            std::cout << std::endl << "Error: at least one of the sugars in the supplied PDB file is missing required atoms" << std::endl << std::endl;
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
                    else // it's the sugar we're looking to omit
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
        
        if (!batch) printf("\nPDB \t    Sugar   \t  Q  \t Phi  \tTheta \t   Detected type   \tCnf\t<Bfac>\tBonds\tAngles\tCtx\t Ok?");
        if (!batch && showGeom) printf("\tBond lengths, angles and torsions, reported clockwise with in-ring oxygen as first vertex");
        if (!batch) printf("\n----\t------------\t-----\t------\t------\t-------------------\t---\t------\t-----\t------\t---\t-----");
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
                fprintf(output, "%s-%s  ", ligandList[index].first.c_str(), ligandList[index].second.id().trim().c_str());
            }
            else
            {
                printf("%c%c%c%c\t%s-",ippdb[1+pos_slash],ippdb[2+pos_slash],ippdb[3+pos_slash],ippdb[4+pos_slash], ligandList[index].second.type().c_str());
                std::cout << ligandList[index].first << "-" << ligandList[index].second.id().trim();
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
                
                fprintf(output,"%3.2f\t%1.3f\t%1.3f", bfac, ligandList[index].second.ring_bond_rmsd(), ligandList[index].second.ring_angle_rmsd());     // output <bfactor>
                
                
                std::vector < clipper::MGlycan > list_of_glycans = mgl.get_list_of_glycans();
                bool found_in_tree = false;
                
                for ( int i = 0 ; i < list_of_glycans.size() ; i++ )
                {
                    std::vector < clipper::MSugar > list_of_sugars = list_of_glycans[i].get_sugars();
                    
                    for ( int j = 0 ; j < list_of_sugars.size() ; j++ )
                    {
                        if ( list_of_sugars[j].id().trim() == ligandList[index].second.id().trim() )
                        {
                            if ( list_of_glycans[i].get_kind_of_glycan() == "n-glycan" )
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
                    if (get_altconformation(ringcomponents[i]) != ' ')
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
                printf("%3.2f\t%1.3f\t%1.3f", bfac, ligandList[index].second.ring_bond_rmsd(), ligandList[index].second.ring_angle_rmsd()); // output <Bfactor>
                
                
                std::vector < clipper::MGlycan > list_of_glycans = mgl.get_list_of_glycans();
                bool found_in_tree = false;
                
                for ( int i = 0 ; i < list_of_glycans.size() ; i++ )
                {
                    std::vector < clipper::MSugar > list_of_sugars = list_of_glycans[i].get_sugars();
                    
                    for ( int j = 0 ; j < list_of_sugars.size() ; j++ )
                    {
                        if ( list_of_sugars[j].id().trim() == ligandList[index].second.id().trim() )
                        {
                            if ( list_of_glycans[i].get_kind_of_glycan() == "n-glycan" )
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
                    if (get_altconformation(ringcomponents[i]) != ' ')
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
            std::cout << "\nThe results for those monosaccharides marked with (*) correspond to the first of at least two possible conformations, each with occupancy < 1.0";
            std::cout << std::endl << std::endl;
        }
        else
            fclose(output);
        
        std::fstream of_scm; std::fstream of_py;
        of_scm.open("privateer-results.scm", std::fstream::out);
        of_py.open("privateer-results.py", std::fstream::out);
        privateer::insert_coot_prologue_scheme ( of_scm );
        privateer::insert_coot_prologue_python ( of_py );
        
        clipper::String all_MapName, dif_MapName, omit_dif_MapName;
        all_MapName = ""; dif_MapName = ""; omit_dif_MapName = "";
        
        privateer::insert_coot_files_loadup_scheme (of_scm, ippdb, all_MapName, dif_MapName, omit_dif_MapName, batch );
        privateer::insert_coot_files_loadup_python (of_py , ippdb, all_MapName, dif_MapName, omit_dif_MapName, batch );
        
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
                    
                    if ( (!ligandList[k].second.ok_with_ring()) || (!ligandList[k].second.ok_with_bonds_rmsd()) || (!ligandList[k].second.ok_with_angles_rmsd()) )
                    {
                        diagnostic.append("check geometry");
                        n_errors++; n_geom++;
                    }
                    
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
                    
                    privateer::insert_coot_go_to_sugar_scheme ( of_scm, ligandList[k].second.ring_centre(), diagnostic );
                    privateer::insert_coot_go_to_sugar_python ( of_py, ligandList[k].second.ring_centre(), diagnostic );
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
                        
                        privateer::insert_coot_go_to_sugar_scheme ( of_scm, ligandList[k].second.ring_centre(), diagnostic );
                        privateer::insert_coot_go_to_sugar_python ( of_py, ligandList[k].second.ring_centre(), diagnostic );
                    }
                }
            }
        }
        
        clipper::String status_msg = "Blue map: 2mFo-DFc. Pink map: omit mFo-DFc. Torsion restraints have been enabled.";
        
        privateer::insert_coot_epilogue_scheme ( of_scm );
        privateer::insert_coot_epilogue_python ( of_py );
        
        privateer::insert_coot_statusbar_text_scheme ( of_scm, status_msg );
        privateer::insert_coot_statusbar_text_python ( of_py, status_msg );
        
        of_scm.close();
        of_py.close();
        
        std::cout << "SUMMARY: " << std::endl << std::endl ;
        std::cout << "   Check ring geometry: " << n_geom << std::endl;
        std::cout << "   Wrong anomer: " << n_anomer << std::endl;
        std::cout << "   Wrong configuration: " << n_config << std::endl;
        std::cout << "   Unphysical puckering amplitude: " << n_pucker << std::endl;
        std::cout << "   In higher-energy conformations: " << n_conf << std::endl;
        std::cout << std::endl;
        std::cout << "   Privateer-validate has identified " << n_geom + n_anomer + n_config + n_pucker + n_conf;
        std::cout << " issues, with " << sugar_count << " of " << ligandList.size() << " sugars affected." << std::endl;
        
        if ( enable_torsions_for.size() > 0 )
            write_refmac_keywords ( enable_torsions_for );
        
        prog.set_termination_message( "Normal termination" );
        system("touch scored");
        return 0;
        
    }
    
    
    
    
    
    
    
    std::vector<clipper::String> enable_torsions_for;
    
    // Full analysis, slower but much more complete //
    
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
            std::cout << "Opening CIF file: " << ipcif.trim() << std::endl;
        
        cifin.open_read( ipcif.trim().c_str() );
        
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
    
    if (!batch)
    {
        std::cout << std::endl << " " << fobs.num_obs() << " reflections have been loaded";
        std::cout << std::endl << std::endl << " Resolution " << hklinfo.resolution().limit() << "Å" << std::endl << hklinfo.cell().format() << std::endl;
    }
    
    clipper::Atom_list mainAtoms;
    clipper::Atom_list ligandAtoms;
    clipper::Atom_list allAtoms;
    
    std::vector<std::pair< clipper::String , clipper::MSugar> > ligandList; // we store the Chain ID and create an MSugar to be scored
    std::vector<clipper::MMonomer> sugarList; // store the original MMonomer
    
    const clipper::MAtomNonBond& manb = clipper::MAtomNonBond( mmol, 1.0 ); // was 1.0
    
    mgl = clipper::MGlycology(mmol, manb);
    
    list_of_glycans = mgl.get_list_of_glycans();
    
    if ( !batch )
    {
        std::cout << std::endl << "Number of detected glycans: " << list_of_glycans.size();
        
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
            }
        }
        
        std::cout << std::endl << "Analysing carbohydrates... ";
        fflush(0);
    }
    
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
                    clipper::MSugar msug(mmol, mmol[p][m], manb);
                    
                    sugarList.push_back(mmol[p][m]);
                    
                    clipper::String id = mmol[p].id();
                    id.resize(1);
                    ligandList.push_back(std::pair<clipper::String, clipper::MSugar> (id, msug));
                    
                    if ( msug.type_of_sugar() == "unsupported" )
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
                else // it's the ligand we're looking to omit
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
    
    char all_MapName[40];
    
    if ( batch )
        sprintf( all_MapName, "sigmaa_best.map\0" );
    else
        sprintf ( all_MapName, "%c%c%c%c_sigmaa_best_%s.map\0",ippdb[1+pos_slash],ippdb[2+pos_slash],ippdb[3+pos_slash],ippdb[4+pos_slash], ipcode.trim().c_str());
    
    char dif_MapName[40];
    
    if ( batch )
        sprintf ( dif_MapName, "sigmaa_diff.map\0" );
    else
        sprintf(dif_MapName, "%c%c%c%c_sigmaa_diff_%s.map\0",ippdb[1+pos_slash],ippdb[2+pos_slash],ippdb[3+pos_slash],ippdb[4+pos_slash], ipcode.trim().c_str());
    
    char omit_dif_MapName[40];
    
    if ( batch )
        sprintf ( omit_dif_MapName, "sigmaa_diff_omit.map\0" );
    else
        sprintf(omit_dif_MapName, "%c%c%c%c_sigmaa_omit_diff_%s.map",ippdb[1+pos_slash],ippdb[2+pos_slash],ippdb[3+pos_slash],ippdb[4+pos_slash], ipcode.trim().c_str());
    
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
                sigmaa_all_MapOut.open_write( all_MapName );      // write maps
                sigmaa_all_MapOut.export_xmap( sigmaa_all_map );
                sigmaa_all_MapOut.close_write();
            }
#pragma omp section
            {
                sigmaa_dif_MapOut.open_write( dif_MapName );
                sigmaa_dif_MapOut.export_xmap( sigmaa_dif_map );
                sigmaa_dif_MapOut.close_write();
            }
#pragma omp section
            {
                sigmaa_omit_fd_MapOut.open_write( omit_dif_MapName );
                sigmaa_omit_fd_MapOut.export_xmap( sigmaa_omit_fd );
                sigmaa_omit_fd_MapOut.close_write();
            }
        }
        
        std::cout << "done." << std::endl;
        std::cout << "\n\nDetailed validation data" << std::endl;
        std::cout << "------------------------" << std::endl;
    }
    
    if (!batch)
        printf("\nPDB \t    Sugar   \tRsln\t  Q  \t Phi  \tTheta \tRSCC\t   Detected type   \tCnf\t<mFo>\t<Bfac>\tBonds\tAngles\tCtx\t Ok?");
    if (!batch && showGeom)
        printf("\tBond lengths, angles and torsions, reported clockwise with in-ring oxygen as first vertex");
    if (!batch)
        printf("\n----\t------------\t----\t-----\t------\t------\t----\t-------------------\t---\t-----\t------\t-----\t------\t---\t-----");
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
            fprintf(output, "%s-%s  ", ligandList[index].first.c_str(), ligandList[index].second.id().trim().c_str());
        }
        else
        {
            printf("%c%c%c%c\t%s-",ippdb[1+pos_slash],ippdb[2+pos_slash],ippdb[3+pos_slash],ippdb[4+pos_slash], ligandList[index].second.type().c_str());
            std::cout << ligandList[index].first << "-" << ligandList[index].second.id().trim();
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
            
            fprintf(output,"%3.2f\t%1.3f\t%1.3f", bfac, ligandList[index].second.ring_bond_rmsd(), ligandList[index].second.ring_angle_rmsd()); // output <bfactor>
            
            std::vector < clipper::MGlycan > list_of_glycans = mgl.get_list_of_glycans();
            bool found_in_tree = false;
            
            for ( int i = 0 ; i < list_of_glycans.size() ; i++ )
            {
                std::vector < clipper::MSugar > list_of_sugars = list_of_glycans[i].get_sugars();
                
                for ( int j = 0 ; j < list_of_sugars.size() ; j++ )
                {
                    if ( list_of_sugars[j].id().trim() == ligandList[index].second.id().trim() )
                    {
                        if ( list_of_glycans[i].get_kind_of_glycan() == "n-glycan" )
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
                if (get_altconformation(ringcomponents[i]) != ' ')
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
            
            printf("%3.2f\t%1.3f\t%1.3f", bfac, ligandList[index].second.ring_bond_rmsd(), ligandList[index].second.ring_angle_rmsd());                 // output <Bfactor>
            
            std::vector < clipper::MGlycan > list_of_glycans = mgl.get_list_of_glycans();
            bool found_in_tree = false;
            
            for ( int i = 0 ; i < list_of_glycans.size() ; i++ )
            {
                std::vector < clipper::MSugar > list_of_sugars = list_of_glycans[i].get_sugars();
                
                for ( int j = 0 ; j < list_of_sugars.size() ; j++ )
                {
                    if ( list_of_sugars[j].id().trim() == ligandList[index].second.id().trim() )
                    {
                        if ( list_of_glycans[i].get_kind_of_glycan() == "n-glycan" )
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
                if (get_altconformation(ringcomponents[i]) != ' ')
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
        std::cout << "The results for those monosaccharides marked with (*) correspond to the first of at least two possible conformations, each with occupancy < 1.0";
        std::cout << std::endl << std::endl;
    }
    else
        fclose(output);
    
    std::fstream of_scm; std::fstream of_py;
    of_scm.open("privateer-results.scm", std::fstream::out);
    of_py.open("privateer-results.py", std::fstream::out);
    
    privateer::insert_coot_prologue_scheme ( of_scm );
    privateer::insert_coot_prologue_python ( of_py );
    privateer::insert_coot_files_loadup_scheme (of_scm, ippdb, all_MapName, dif_MapName, omit_dif_MapName, batch );
    privateer::insert_coot_files_loadup_python (of_py, ippdb, all_MapName, dif_MapName, omit_dif_MapName, batch );
    
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
                
                if ( (!ligandList[k].second.ok_with_ring()) || (!ligandList[k].second.ok_with_bonds_rmsd()) || (!ligandList[k].second.ok_with_angles_rmsd()) )
                {
                    diagnostic.append("check geometry");
                    report.append("Check geometry");
                    n_errors++;
                    n_geom++;
                }
                
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
                
                privateer::insert_coot_go_to_sugar_scheme ( of_scm, ligandList[k].second.ring_centre(), diagnostic );
                privateer::insert_coot_go_to_sugar_python ( of_py, ligandList[k].second.ring_centre(), diagnostic );
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
                    privateer::insert_coot_go_to_sugar_scheme ( of_scm, ligandList[k].second.ring_centre(), diagnostic );
                    privateer::insert_coot_go_to_sugar_python ( of_py, ligandList[k].second.ring_centre(), diagnostic );
                }
                
            }
        }
        
        if ( report == "" )
            report = "Ok";
        
        ligandList[k].second.set_diagnostic ( report );
        
    }
    
    privateer::insert_coot_epilogue_scheme ( of_scm );
    privateer::insert_coot_epilogue_python ( of_py );
    of_scm.close();
    of_py.close();
    
    std::cout << "SUMMARY: " << std::endl << std::endl ;
    std::cout << "   Check ring geometry: " << n_geom << std::endl;
    std::cout << "   Wrong anomer: " << n_anomer << std::endl;
    std::cout << "   Wrong configuration: " << n_config << std::endl;
    std::cout << "   Unphysical puckering amplitude: " << n_pucker << std::endl;
    std::cout << "   In higher-energy conformations: " << n_conf << std::endl;
    std::cout << std::endl;
    std::cout << "   Privateer-validate has identified " << n_geom + n_anomer + n_config + n_pucker + n_conf;
    std::cout << " issues, with " << sugar_count << " of " << ligandList.size() << " sugars affected." << std::endl;
    
    print_XML(ligandList, list_of_glycans, ippdb);
    
    if ( enable_torsions_for.size() > 0 )
        write_refmac_keywords ( enable_torsions_for );
    
    prog.set_termination_message( "Normal termination" );
    system("touch scored");
    return 0;
}










// ############# helper functions ############## //





char get_altconformation(clipper::MAtom ma)
{
    clipper::String identifier = ma.id();
    if (identifier.size() > 5)
        return identifier[5];
    else
        return ' ';     // The alternate conformation code is the fifth character in the complete identificator.
}                       // We will return a blank character if there is no code present or if it is, but is blank






void print_XML ( std::vector < std::pair < clipper::String, clipper::MSugar > > sugarList, std::vector < clipper::MGlycan > list_of_glycans, clipper::String pdbname )
{
    
    std::fstream of_xml;
    
    of_xml.open("program.xml", std::fstream::out);
    
    of_xml << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
    of_xml << "<?xml-stylesheet type=\"text/xsl\" href=\"program.xml\"?>\n";
    of_xml << "<xsl:stylesheet id=\"stylesheet\" version=\"1.0\" xmlns:xsl=\"http://www.w3.org/1999/XSL/Transform\">\n";
    of_xml << "<xsl:value-of select=\".\" disable-output-escaping = \"yes\"/>\n" ;
    
    of_xml << "<PrivateerResult>\n";
    of_xml << "  <ValidationData>\n";
    
    for (int i = 0; i < sugarList.size() ; i++ )
    {
        if ( ( sugarList[i].second.ring_cardinality() == 6 ) || ( sugarList[i].second.ring_cardinality() == 5 ) )
        {
            std::vector<clipper::ftype> cpParams = sugarList[i].second.cremer_pople_params();
            
            clipper::String sugarRSCC = clipper::String ( sugarList[i].second.get_rscc() );
            sugarRSCC.resize(4);
            
            clipper::String puckering_amplitude = clipper::String ( sugarList[i].second.puckering_amplitude() );
            puckering_amplitude.resize(5);
            
            if (sugarList[i].second.ring_cardinality() == 6 )
                of_xml << "    <Pyranose>\n";
            else
                of_xml << "    <Furanose>\n";
            
            of_xml << "      <SugarPDB>"            << pdbname                                          << "</SugarPDB>\n"          ;
            of_xml << "      <SugarName>"           << sugarList[i].second.type()                       << "</SugarName>\n"         ;
            of_xml << "      <SugarChain>"          << sugarList[i].first                               << "</SugarChain>\n"        ;
            of_xml << "      <SugarQ>"              << puckering_amplitude                              << "</SugarQ>\n"            ;
            
            of_xml << "      <SugarPhi>"            << cpParams[1]                                      << "</SugarPhi>\n"          ;
            
            if (sugarList[i].second.ring_cardinality() == 6 )
                of_xml << "      <SugarTheta>"          << cpParams[2]                                  << "</SugarTheta>\n"        ;
            
            of_xml << "      <SugarAnomer>"         << sugarList[i].second.anomer()                     << "</SugarAnomer>\n"       ;
            of_xml << "      <SugarHand>"           << sugarList[i].second.handedness()                 << "</SugarHand>\n"         ;
            of_xml << "      <SugarConformation><![CDATA["   << sugarList[i].second.conformation_name_iupac()          << "]]></SugarConformation>\n" ;
            of_xml << "      <SugarBondRMSD>"       << sugarList[i].second.ring_bond_rmsd()             << "</SugarBondRMSD>\n"     ;
            of_xml << "      <SugarAngleRMSD>"      << sugarList[i].second.ring_angle_rmsd()            << "</SugarAngleRMSD>\n"    ;
            of_xml << "      <SugarBFactor>"        << sugarList[i].second.get_bfactor()                << "</SugarBFactor>\n"      ;
            of_xml << "      <SugarRSCC>"           << sugarRSCC                                        << "</SugarRSCC>\n"         ;     
            of_xml << "      <SugarDiagnostic>"     << sugarList[i].second.get_diagnostic()             << "</SugarDiagnostic>\n"   ;
            of_xml << "      <SugarContext>"        << sugarList[i].second.get_context()                << "</SugarContext>\n"      ;
            
            if (sugarList[i].second.ring_cardinality() == 6 ) 
                of_xml << "    </Pyranose>\n\n";
            else
                of_xml << "    </Furanose>\n\n";
        }
    }
    
    for ( int i = 0 ; i < list_of_glycans.size() ; i++ )
    {
        of_xml << "    <Glycan>\n" ;
        of_xml << "     <GlycanPDB>"        << pdbname                                           << "</GlycanPDB>\n"                                         ;
        of_xml << "     <GlycanType>"       << list_of_glycans[i].get_type()                     << "</GlycanType>\n"                                       ;
        of_xml << "     <GlycanRoot>"       << list_of_glycans[i].get_root().type().trim() + list_of_glycans[i].get_root().id().trim() << "</GlycanRoot>\n" ;
        of_xml << "     <GlycanChain>"      << list_of_glycans[i].get_chain()                    << "</GlycanChain>\n"                                      ;
        of_xml << "     <GlycanText><![CDATA["<< list_of_glycans[i].print_linear ( false, true, true )    << "]]></GlycanText>\n"                           ;
        of_xml << "     <GlycanSVG>"        << "placeholder.svg"                                 << "</GlycanSVG>\n"                                        ;
        of_xml << "    </Glycan>\n" ;    
    }
    
    of_xml << "  </ValidationData>\n";
    of_xml << "</PrivateerResult>\n";
    of_xml << "</xsl:stylesheet>\n\n";
    
    of_xml.close(); 
}






void print_usage ( )
{
    std::cout << "Usage: privateer-validate" << std::endl << std::endl;
    std::cout << "\t-pdbin <.pdb>\t\t\tCOMPULSORY: input PDB file to examine" << std::endl;
    std::cout << "\t-cifin <.cif> OR -mtzin <.mtz>\tCIF file containing F's or MTZ file containing I's or F's\n";
    std::cout << "\t\t\t\t\tIf intensities are supplied, Privateer will run 'ctruncate'" << std::endl;
    std::cout << "\t\t\t\t\tin order to convert intensities to amplitudes" << std::endl;
    std::cout << "\t\t\t\t\tObservations are required for RSCC calculation and map output" << std::endl;
    std::cout << "\t-mtzout <.mtz>\t\t\tOutput best and difference map coefficients to MTZ files" << std::endl;
    std::cout << "\t-colin-fo\t\t\tColumns containing F & SIGF, e.g. FOBS,SIGFOBS" << std::endl;
    std::cout << "\t\t\t\t\tIf not supplied, Privateer will try to guess the path" << std::endl;
    std::cout << "\t-codein <3-letter code>\t\tA 3-letter code for the target sugar" << std::endl;
    std::cout << "\t\t\t\t\tIf not supplied, Privateer will check all known sugar codes" << std::endl;
    std::cout << "\t-valstring <options>\t\tUse external validation options" << std::endl;
    std::cout << "\t\t\t\t\tExample: SUG,O5/C1/C2/C3/C4/C5,A,D,4c1" << std::endl;
    std::cout << "\t\t\t\t\tThree-letter code, ring atoms, anomer, handedness, expected conformation" << std::endl;
    std::cout << "\t-showgeom\t\t\tRing bond lengths, angles and torsions are reported clockwise" << std::endl;
    std::cout << "\t\t\t\t\tExample: first bond, angle and torsion for an aldopyranose:" << std::endl;
    std::cout << "\t\t\t\t\t\t (O5-C1), (C5-O5-C1), (C5-O5-C1-C2)" << std::endl;
    std::cout << "\t-radiusin <value>\t\tA radius for the calculation of a mask around the target sugar" << std::endl;
    std::cout << "\t\t\t\t\tDefault value is 1.5 angstroems" << std::endl;
    std::cout << "\t-list\t\t\t\tProduces a list of space-separated supported 3-letter codes and stops" << std::endl;
    std::cout << "\t-mode <normal|ccp4i2>\t\tOptional: mode of operation. Defaults to normal" << std::endl;
    std::cout << "\t\t\t\t\tccp4i2 mode produces XML and tabbed-value files" << std::endl << std::endl;
    std::cout << "\tThe program will also produce a visual checklist with the conflicting sugar models in the form\n";
    std::cout << "\tof Scheme and Python scripts for use with Coot" << std::endl;
    std::cout << "\n\tTo use them: 'coot --script privateer-results.scm' or 'coot --script privateer-results.py'\n";
    std::cout << "\tBlue map: 2mFo-DFc map. Pink map: mFo-DFc omit map. RSCC is computed against this map" << std::endl;

}






bool compute_and_print_external_validation ( const std::vector<clipper::String> validation_options, clipper::data::sugar_database_entry& external_validation )
{
    std::vector<clipper::String> atoms_buffer = validation_options[1].split("/");
    
    // required format: SUG,O5/C1/C2/C3/C4/C5,A,D,4c1
    
    if (( validation_options.size() != 5 ) || ( atoms_buffer.size() > 6 ) || ( atoms_buffer.size() < 5 ) )
    {
        std::cout << "\n\nCannot interpret the supplied validation options correctly\n\n"
        << "\tAccepted format: SUG,O5/C1/C2/C3/C4/C5,A,D,4c1\n"
        << "\tThree-letter code, ring atoms separated by /, anomer, handedness, expected conformation."
        << std::endl << std::endl;
        return true;
    }
    else
    {
        external_validation.name_short       = validation_options[0];
        external_validation.anomer           = validation_options[2];
        external_validation.handedness       = validation_options[3];
        external_validation.ref_conformation = validation_options[4];
        external_validation.name_long        = "Unidentified sugar" ;
        
        for ( int index = 0; index < atoms_buffer.size() ; index++ )
            external_validation.ring_atoms += atoms_buffer[index] + " ";
        
        if (( external_validation.ref_conformation == "4c1" ) || ( external_validation.ref_conformation == "1c4" ))
        {
            external_validation.ref_puckering   = 0.540;  // these are knowledge-based values
            external_validation.ref_bonds_rmsd  = 0.010;
            external_validation.ref_angles_rmsd = 1.500;
        }
        else
        {
            external_validation.ref_puckering   = 0.640;
            external_validation.ref_bonds_rmsd  = 0.020;
            external_validation.ref_angles_rmsd = 4.000;
        }
        std::cout << std::endl       << "External validation data supplied" << std::endl
        << "---------------------------------" << std::endl
        << "Name: "        << external_validation.name_short << std::endl
        << "Anomer: "      << external_validation.anomer << std::endl
        << "Handedness: "  << external_validation.handedness << std::endl
        << "Conformation: "<< external_validation.ref_conformation << std::endl
        << "Ring atoms: "  << external_validation.ring_atoms << std::endl << std::endl;
    }
    return false;
}







bool read_coordinate_file (clipper::MMDBfile& mfile, clipper::MiniMol& mmol, clipper::String& ippdb, bool batch)
{
    if (!batch)
    {
        std::cout << std::endl << "Reading " << ippdb.trim().c_str() << "... ";
        fflush(0);
    }
    
    const int mmdbflags = mmdb::MMDBF_IgnoreBlankLines | mmdb::MMDBF_IgnoreDuplSeqNum | mmdb::MMDBF_IgnoreNonCoorPDBErrors | mmdb::MMDBF_IgnoreRemarks | mmdb::MMDBF_EnforceUniqueChainID;
    mfile.SetFlag( mmdbflags );
    
    mfile.read_file( ippdb.trim() );
    mfile.import_minimol( mmol );
    
    
    if (!batch)
        std::cout << "done." << std::endl;
    
    if ( mmol.cell().is_null() )
    {
        std::cout << std::endl << " Spacegroup/cell information is missing from the PDB file." << std::endl;
        std::cout << " Privateer-validate will still run, but may miss any important contacts described by crystallographic symmetry." << std::endl << std::endl;
        mmol.init ( clipper::Spacegroup::p1(), clipper::Cell(clipper::Cell_descr ( 300, 300, 300, 90, 90, 90 )) );
        return false;
    }
    
    return true;

}







void print_supported_code_list ()
{
    std::cout << std::endl << std::endl << "Printing list of supported three-letter codes" << std::endl ;
    std::cout << std::endl << "Pyranoses: " << std::endl;;
    
    int n_printed = 0;
    
    for ( int i = 0 ; i < clipper::data::sugar_database_size ; i++ )
    {
        if ( clipper::data::sugar_database[i].ring_atoms.split(" ").size() == 6 )
        {
            if ( (n_printed % 22) == 0 )
                std::cout << std::endl;
            std::cout << clipper::data::sugar_database[i].name_short << " " ;
            n_printed++;
        }
    }
    std::cout << std::endl << std::endl ;
    std::cout << std::endl << "Furanoses: " << std::endl;;
    
    n_printed = 0;
    
    for ( int i = 0 ; i < clipper::data::sugar_database_size ; i++ )
    {
        if ( clipper::data::sugar_database[i].ring_atoms.split(" ").size() == 5 )
        {
            if ( (n_printed % 22) == 0 )
                std::cout << std::endl;
            std::cout << clipper::data::sugar_database[i].name_short << " " ;
            n_printed++;
        }
    }
    std::cout << std::endl << std::endl ;
    
    n_printed = 0;
    
    std::cout << std::endl << "Disaccharides: " << std::endl;
    
    for ( int i = 0 ; i < clipper::data::disaccharide_database_size ; i++ )
    {
        if ( (n_printed % 22) == 0 )
            std::cout << std::endl;
        std::cout << clipper::data::disaccharide_database[i].name_short << " " ;
        n_printed++;
    }
    std::cout << std::endl << std::endl ;

}






void write_refmac_keywords ( std::vector < clipper::String > code_list )
{
    std::fstream of_keywords;
    of_keywords.open("keywords_refmac5.txt", std::fstream::out);

    // remove replicas

    std::sort ( code_list.begin(), code_list.end() );
    std::vector< clipper::String>::iterator last = std::unique(code_list.begin(), code_list.end());
    code_list.erase ( last, code_list.end() );
    
    // output unique keywords
    
    for ( int index = 0; index < code_list.size() ; index++ )
    {
        of_keywords << "restr tors include resi " << code_list[index] << std::endl;
        if ( index == 0 )
            std::cout << std::endl << std::endl << "The keywords_refmac5.txt file has been generated. This activates torsion restraints in refmac5 for " << code_list[index];
        else if ( index == code_list.size() -1 )
            std::cout << " and " << code_list[index] << std::endl;
        else
            std::cout << ", " << code_list[index];
    }
    std::cout << std::endl << std::endl;
    of_keywords.close();
}





