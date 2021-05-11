// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#include "privateer-pyanalysis.h"

using namespace pybind11::literals;

// #define DUMP 1
// #define DBG std::cout << "[" << __FUNCTION__ << "] - "

///////////////////////////////////////////////// Class GlycosylationComposition ////////////////////////////////////////////////////////////////////

void privateer::pyanalysis::GlycosylationComposition::read_from_file( std::string path_to_model_file, std::string expression_system ) {

    if(path_to_model_file == "undefined")
    {
        throw std::invalid_argument( "No path was provided for model file input! Aborting." );
    }

    this->path_to_model_file = path_to_model_file;
    this->expression_system = expression_system;
    
    clipper::MiniMol mmol;
    clipper::MMDBfile mfile;
    clipper::String path_to_model_file_clipper = path_to_model_file;

    privateer::util::read_coordinate_file_mtz(mfile, mmol, path_to_model_file_clipper, true);

    this->mgl = clipper::MGlycology(mmol, expression_system);

    initialize_summary_of_detected_glycans(mgl);
    
    std::vector<clipper::MGlycan> list_of_glycans = mgl.get_list_of_glycans();
    auto list = pybind11::list();
    for(int i = 0; i < list_of_glycans.size(); i++)
    {
        auto glycanObject = privateer::pyanalysis::GlycanStructure(mgl, i, *this);
        list.append(glycanObject);
    }
    this->glycans = list;
}

void privateer::pyanalysis::GlycosylationComposition::initialize_summary_of_detected_glycans( clipper::MGlycology& mglObject )
{
    std::vector<clipper::MGlycan> list_of_glycans = mglObject.get_list_of_glycans();

    this->numberOfGlycanChains = list_of_glycans.size();

    auto list = pybind11::list();
    for(int i = 0; i < list_of_glycans.size(); i++)
    {
        std::string wurcsNotation = list_of_glycans[i].generate_wurcs();
        std::string kindOfGlycan = list_of_glycans[i].get_type();
        
        list_of_glycans[i].get_root_by_name();
        std::string proteinResidue = list_of_glycans[i].get_root().first.type().trim();
        std::string proteinResidueID = list_of_glycans[i].get_root().first.id().trim();
        std::string proteinChainID = list_of_glycans[i].get_chain().substr(0,1);

        auto rootSummary = pybind11::dict ("ProteinResidueType"_a=proteinResidue, "ProteinResidueID"_a=std::stoi(proteinResidueID), "ProteinChainID"_a=proteinChainID);
        
        std::vector<float> torsions = list_of_glycans[i].get_glycosylation_torsions();
        auto protein_glycan_linkage_torsion = pybind11::dict ("Phi"_a=torsions[0], "Psi"_a=torsions[1]);
        
        auto dict = pybind11::dict ("GlycanID"_a=i, "WURCS"_a=wurcsNotation, "GlycosylationType"_a=kindOfGlycan, "RootInfo"_a=rootSummary, "ProteinGlycanLinkageTorsion"_a=protein_glycan_linkage_torsion);
        list.append(dict);
    }
    this->glycosylationSummary = list;
}


privateer::pyanalysis::GlycanStructure privateer::pyanalysis::GlycosylationComposition::get_glycan(const int glycanID)
{
    if(glycanID >= numberOfGlycanChains || glycanID < 0)
    {
        throw std::invalid_argument( "Provided ID is out of bounds and exceeds/inceeds number of glycans detected in the model. \nInput: " + std::to_string(glycanID) + "\tPermitted Range: [0-" + std::to_string(numberOfGlycanChains - 1) + "]");
    }

    auto glycanObject = glycans[glycanID].cast<privateer::pyanalysis::GlycanStructure>();
    return glycanObject;
}
///////////////////////////////////////////////// Class GlycosylationComposition END ////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////// Class GlycanStructure ////////////////////////////////////////////////////////////////////
void privateer::pyanalysis::GlycanStructure::pyinit( const clipper::MGlycology& mgl, const int glycanID, privateer::pyanalysis::GlycosylationComposition& parentGlycosylationComposition)
{
    this->parentGlycosylation = parentGlycosylationComposition;
    std::vector<clipper::MGlycan> list_of_glycans = mgl.get_list_of_glycans();

    clipper::MGlycan inputGlycan = list_of_glycans[glycanID];

    this->glycan = inputGlycan;
    this->sugars_in_glycan = inputGlycan.get_sugars();
    
    this->glycanID = glycanID;
    this->numberOfSugars = inputGlycan.number_of_nodes();
    this->glycanWURCS = inputGlycan.generate_wurcs();

    auto uniqueMonosaccharidesPyList = pybind11::list();
    std::vector<std::string> uniqueMonosaccharidesVector = inputGlycan.obtain_unique_residue_codes();
    for(int i = 0; i < uniqueMonosaccharidesVector.size(); i++)
    {
        uniqueMonosaccharidesPyList.append(uniqueMonosaccharidesVector[i]);
    }
    this->uniqueMonosaccharides = uniqueMonosaccharidesPyList;
    this->numberOfGlycosidicBonds = inputGlycan.obtain_total_number_of_glycosidic_bonds();
    this->glycosylationType = inputGlycan.get_type();
    
    std::string proteinResidue = inputGlycan.get_root().first.type().trim();
    std::string proteinResidueID = inputGlycan.get_root().first.id().trim();
    std::string proteinChainID = inputGlycan.get_chain().substr(0,1);

    auto rootSummary = pybind11::dict("ProteinResidueType"_a=proteinResidue, "ProteinResidueID"_a=std::stoi(proteinResidueID), "ProteinChainID"_a=proteinChainID);
    this->rootSummary = rootSummary;

    std::vector<float> torsions = inputGlycan.get_glycosylation_torsions();
    auto protein_glycan_linkage_torsion = pybind11::dict("Phi"_a=torsions[0], "Psi"_a=torsions[1]);
    this->protein_glycan_linkage_torsion = protein_glycan_linkage_torsion;

    std::vector<clipper::MSugar> list_of_sugars = inputGlycan.get_sugars();

    auto list = pybind11::list();
    for(int i = 0; i < list_of_sugars.size(); i++)
    {
        auto sugarObject = privateer::pyanalysis::CarbohydrateStructure(inputGlycan, i, glycanID, parentGlycosylation, *this);
        list.append(sugarObject);
    }
    this->sugars = list;

    initialize_summary_of_glycan();
}

void privateer::pyanalysis::GlycanStructure::initialize_summary_of_glycan( )
{
    auto dict = pybind11::dict ("GlycanID"_a=glycanID, "WURCS"_a=glycanWURCS, "UniqueMonosaccharides"_a=uniqueMonosaccharides, "TotalSugars"_a=numberOfSugars, "NumberOfGlycosidicBonds"_a=numberOfGlycosidicBonds, "GlycosylationType"_a=glycosylationType, "RootInfo"_a=rootSummary, "ProteinGlycanLinkageTorsion"_a=protein_glycan_linkage_torsion);

    this->glycanSummary = dict;
}

privateer::pyanalysis::CarbohydrateStructure privateer::pyanalysis::GlycanStructure::get_monosaccharide(const int sugarID)
{
    if(sugarID >= numberOfSugars || sugarID < 0)
    {
        throw std::invalid_argument( "Provided ID is out of bounds and exceeds/inceeds number of sugars detected in the glycan. \nInput: " + std::to_string(sugarID) + "\tPermitted Range: [0-" + std::to_string(numberOfSugars - 1) + "]");
    }

    auto sugarObject = sugars[sugarID].cast<privateer::pyanalysis::CarbohydrateStructure>();
    return sugarObject;
}
///////////////////////////////////////////////// Class GlycanStructure END ////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////// Class CarbohydrateStructure  ////////////////////////////////////////////////////////////////////
void privateer::pyanalysis::CarbohydrateStructure::pyinit( clipper::MGlycan& mglycan, const int sugarID, const int glycanID, privateer::pyanalysis::GlycosylationComposition& parentGlycosylationComposition, privateer::pyanalysis::GlycanStructure& parentGlycanStructure)
{
    this->parentGlycosylation = parentGlycosylationComposition;
    this->parentGlycanStructure = parentGlycanStructure;
    // std::vector<clipper::MGlycan> list_of_glycans = mgl.get_list_of_glycans();
    std::vector<clipper::MSugar> list_of_sugars = mglycan.get_sugars();

    clipper::MSugar inputSugar = list_of_sugars[sugarID];
    
    this->parentGlycan = mglycan;
    this->sugar = inputSugar;
    this->sugarID = sugarID;
    this->glycanID = glycanID;
    this->sugar_conformation_code = inputSugar.conformation_code();
    this->sugar_conformation_name = clipper::data::conformational_landscape[sugar_conformation_code];
    this->sugar_conformation_name_iupac = clipper::data::iupac_conformational_landscape[sugar_conformation_code];
    this->sugar_puckering_amplitude = inputSugar.puckering_amplitude();
    this->sugar_anomer = inputSugar.anomer();
    this->sugar_handedness = inputSugar.handedness();
    this->sugar_denomination = inputSugar.type_of_sugar();
    this->sugar_ring_cardinality = inputSugar.ring_cardinality();

    std::vector<clipper::ftype> cremer_pople_params_vector = inputSugar.cremer_pople_params();
    auto sugar_cremer_pople_params = pybind11::list();
    for(int i = 0; i < cremer_pople_params_vector.size(); i++)
    {
        float currentParam = cremer_pople_params_vector[i];
        sugar_cremer_pople_params.append(currentParam);
    }
    this->sugar_cremer_pople_params=sugar_cremer_pople_params;

    this->sugar_sane = inputSugar.is_sane();

    std::string sugardiagnostic;
    if(inputSugar.is_sane())
    {
        if(!inputSugar.ok_with_conformation())
        {
            sugardiagnostic = "check";
        }
        else
            sugardiagnostic = "yes";
    }
    else
        sugardiagnostic = "no";

    this->privateer_diagnostic = sugardiagnostic;

    this->sugar_name_full = inputSugar.full_name();
    this->sugar_name_short = inputSugar.short_name();
    this->sugar_pdb_id = std::stoi(inputSugar.id());
    this->sugar_pdb_chain = mglycan.get_chain().substr(0,1);
    this->sugar_type = inputSugar.full_type();

    std::vector<clipper::ftype> sugar_ring_angles_vector = inputSugar.ring_angles();
    auto sugar_ring_angles = pybind11::list();
    for(int i = 0; i < sugar_ring_angles_vector.size(); i++)
    {
        float currentParam = sugar_ring_angles_vector[i];
        sugar_ring_angles.append(currentParam);
    }
    this->sugar_ring_angles=sugar_ring_angles;

    std::vector<clipper::ftype> sugar_ring_bonds_vector = inputSugar.ring_bonds();
    auto sugar_ring_bonds = pybind11::list();
    for(int i = 0; i < sugar_ring_bonds_vector.size(); i++)
    {
        float currentParam = sugar_ring_bonds_vector[i];
        sugar_ring_bonds.append(currentParam);
    }
    this->sugar_ring_bonds=sugar_ring_bonds;

    std::vector<clipper::ftype> sugar_ring_torsions_vector = inputSugar.ring_torsions();
    auto sugar_ring_torsion = pybind11::list();
    for(int i = 0; i < sugar_ring_torsions_vector.size(); i++)
    {
        float currentParam = sugar_ring_torsions_vector[i];
        sugar_ring_torsion.append(currentParam);
    }
    this->sugar_ring_torsion=sugar_ring_torsion;

    this->sugar_ring_bond_rmsd=inputSugar.ring_bond_rmsd();
    this->sugar_ring_angle_rmsd=inputSugar.ring_angle_rmsd();
    this->sugar_bfactor=inputSugar.get_bfactor();
    this->sugar_supported=inputSugar.is_supported();
    this->sugar_diag_ring=inputSugar.ok_with_ring();
    this->sugar_diag_bonds_rmsd=inputSugar.ok_with_bonds_rmsd();
    this->sugar_diag_angles_rmsd=inputSugar.ok_with_angles_rmsd();
    this->sugar_diag_anomer=inputSugar.ok_with_anomer();
    this->sugar_diag_chirality=inputSugar.ok_with_chirality();
    this->sugar_diag_conformation=inputSugar.ok_with_conformation();
    this->sugar_diag_puckering=inputSugar.ok_with_puckering();
    this->sugar_context=mglycan.get_type();



    initialize_summary_of_sugar();
}


void privateer::pyanalysis::CarbohydrateStructure::initialize_summary_of_sugar( )
{
    auto dict = pybind11::dict ("sugarID"_a=sugarID, "glycanID"_a=glycanID, "sugar_pdb_id"_a=sugar_pdb_id, "sugar_pdb_chain"_a=sugar_pdb_chain, "sugar_name_full"_a=sugar_name_full, "sugar_name_short"_a=sugar_name_short, "sugar_anomer"_a=sugar_anomer, "is_sane"_a=sugar_sane);

    this->sugarSummary = dict;
}

///////////////////////////////////////////////// Class CarbohydrateStructure END ////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////// Class XRayData ////////////////////////////////////////////////////////////////////

void privateer::pyanalysis::XRayData::read_from_file( std::string& path_to_mtz_file, std::string& path_to_model_file, std::string& input_column_fobs_user, float ipradius, int nThreads) {

    privateer::thread_pool pool(0);

    bool useSigmaa = false;
    int detectedThreads = std::thread::hardware_concurrency();
    bool useParallelism = true;
    bool showGeom = false;

    
    if(nThreads < 0)
        nThreads = detectedThreads;
    else if(nThreads < 2 && nThreads > -1)
    {
        useParallelism = false;
    }
    else if(nThreads > detectedThreads)
    {
        std::cout << "Error: More cores/threads were inputted as an argument, than detected on the system." 
        << "\n\tNumber of Available Cores/Threads detected on the system: " << detectedThreads 
        << "\n\tNumber of Cores/Threads requested via input argument: " << nThreads << "." << std::endl;

        throw std::invalid_argument( "Number of inputted threads exceed the number of detected threads." );
    }

    if(useParallelism)
    {
        std::cout << std::endl << "THREADING: Resizing and initiating a thread pool object with " << nThreads << " threads..." << std::endl;
        pool.resize(nThreads);
        std::cout << "THREADING: Successfully initialized a thread pool with " << pool.size() << " threads!" << std::endl << std::endl;
    }

    if(path_to_model_file == "undefined" || path_to_model_file == "")
    {
        throw std::invalid_argument( "No path was provided for model file input! Aborting." );
    }

    if(path_to_mtz_file == "undefined" || path_to_mtz_file == "")
    {
        throw std::invalid_argument( "No path was provided for .mtz file input! Aborting." );
    }

    clipper::String input_column_fobs_clipper = input_column_fobs_user;
    this->input_column_fobs = input_column_fobs_clipper;
    clipper::String path_to_model_file_clipper = path_to_model_file;
    this->path_to_model_file=path_to_model_file;
    clipper::String path_to_mtz_file_clipper = path_to_mtz_file;
    clipper::MMDBfile mfile;
    clipper::CCP4MTZfile mtzin;
    clipper::CCP4MTZfile ampmtzin;
    clipper::MTZcrystal opxtal;
    clipper::MTZdataset opdset;
    clipper::HKL_info hklinfo;

    privateer::util::read_coordinate_file_mtz(mfile, mmol, path_to_model_file_clipper, true);
    int pos_slash = path_to_model_file.rfind("/");
    privateer::xray::read_xray_map( path_to_mtz_file_clipper, path_to_model_file_clipper, mmol, hklinfo, mtzin );

    clipper::HKL_data<clipper::data32::F_sigF> fobs;            // allocate space for F and sigF
    clipper::HKL_data<clipper::data32::F_sigF> fobs_scaled;     // allocate space for scaled F and sigF
    clipper::HKL_data<clipper::data32::F_phi> fc_omit_bsc;      // allocate space for the omit calc data with bsc
    clipper::HKL_data<clipper::data32::F_phi> fc_all_bsc;       // allocate space for the whole calculated data with bsc
    clipper::HKL_data<clipper::data32::F_phi> fc_ligands_bsc;    // allocate space for the ligand calculated data

    fobs = clipper::HKL_data<clipper::data32::F_sigF> ( hklinfo );
    fobs_scaled = clipper::HKL_data<clipper::data32::F_sigF> ( hklinfo );
    fc_omit_bsc = clipper::HKL_data<clipper::data32::F_phi> ( hklinfo );
    fc_all_bsc = clipper::HKL_data<clipper::data32::F_phi> ( hklinfo );
    fc_ligands_bsc = clipper::HKL_data<clipper::data32::F_phi> ( hklinfo );
    privateer::xray::initialize_experimental_dataset( mtzin, ampmtzin, input_column_fobs, fobs, hklinfo, opxtal, opdset, path_to_mtz_file_clipper);

    std::cout << std::endl << " " << fobs.num_obs() << " reflections have been loaded";
        std::cout << std::endl << std::endl << " Resolution " << hklinfo.resolution().limit() << "Å" << std::endl << hklinfo.cell().format() << std::endl;

    this->hklinfo = hklinfo;
    fobs_scaled = fobs;

    clipper::Atom_list mainAtoms;
    clipper::Atom_list ligandAtoms;
    clipper::Atom_list allAtoms;

    std::vector<std::pair< clipper::String , clipper::MSugar> > ligandList; // we store the Chain ID and create an MSugar to be scored
    std::vector<clipper::MMonomer> sugarList; // store the original MMonomer

    const clipper::MAtomNonBond& manb = clipper::MAtomNonBond( mmol, 1.0 ); // was 1.0

    clipper::MGlycology mgl = clipper::MGlycology(mmol, manb, "undefined");

    for ( int p = 0; p < mmol.size(); p++ )
    {
        for ( int m = 0; m < mmol[p].size(); m++ )
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
                    throw std::invalid_argument( "Error: strangely, at least one of the sugars in the supplied PDB file is missing required atoms. Stopping..." );
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

                // #if DUMP
                //     std::cout << "number of alternate conformations: " << conformers.size() << std::endl;
                // #endif

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
                {
                    ligandList.push_back(std::pair<clipper::String, clipper::MSugar> (id, msug_b));
                    sugarList.push_back(mmol[p][m]);
                }

                if ( msug.type_of_sugar() == "unsupported" )
                {
                    throw std::invalid_argument( "Error: strangely, at least one of the sugars in the supplied PDB file is missing required atoms. Stopping..." );
                }

                for (int id = 0; id < mmol[p][m].size(); id++ )
                {
                    ligandAtoms.push_back(mmol[p][m][id]);  // add the ligand atoms to a second array
                    allAtoms.push_back(mmol[p][m][id]);
                }
            }
        }
    }

    clipper::SFcalc_obs_bulk<float> sfcbligands;
    clipper::SFcalc_obs_bulk<float> sfcb;
    clipper::SFcalc_obs_bulk<float> sfcball;

    try
    {   // calculate structure factors with bulk solvent correction
        if(useParallelism)
        {
            pool.push([&sfcbligands, &fc_ligands_bsc, &fobs, &ligandAtoms](int id)
            { 
                #if DUMP
                    std::cout << std::endl;
                    DBG << "Calculating sfcbligands from Thread ID: " << id << '.' << std::endl;
                #endif
                
                sfcbligands( fc_ligands_bsc, fobs, ligandAtoms );
            });
    
            pool.push([&sfcb, &fc_omit_bsc, &fobs, &mainAtoms](int id)
            { 
                #if DUMP
                    DBG << "Calculating sfcb from Thread ID: " << id << '.' << std::endl;
                #endif
                
                sfcb( fc_omit_bsc, fobs, mainAtoms );
            });

            pool.push([&sfcball, &fc_all_bsc, &fobs, &allAtoms](int id)
            { 
                #if DUMP
                    DBG << "Calculating sfcball from Thread ID: " << id << '.' << std::endl;
                #endif
                
                sfcball( fc_all_bsc, fobs, allAtoms );
            });
        }
        else
        {
            sfcbligands( fc_ligands_bsc, fobs, ligandAtoms ); // was fobs_scaled

            sfcb( fc_omit_bsc, fobs, mainAtoms );  // calculation of omit SF with bulk solvent correction

            sfcball( fc_all_bsc, fobs, allAtoms ); // calculation of SF with bulk solvent correction
        }
    }
    catch ( ... )
    {
        std::cout << "\nThe input file has unrecognised atoms. Might cause unexpected results...\n";  // this causes clipper to freak out, so better remove those unknowns
    }

    if (useParallelism)
        {
            #if DUMP
                std::cout << std::endl;
                DBG << "Number of jobs in the queue: " << pool.n_remaining_jobs() << std::endl;
            #endif
            
            std::this_thread::sleep_for(std::chrono::seconds(1));

            while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                pool.sync();
            
            #if DUMP
                DBG << "Number of jobs in the queue: " << pool.n_remaining_jobs() << " after sync operation!" << std::endl;
            #endif
        }


    fc_ligands_bsc[0].set_null();
    fc_omit_bsc[0].set_null();
    fc_all_bsc[0].set_null();


    clipper::Grid_sampling mygrid( hklinfo.spacegroup(), hklinfo.cell(), hklinfo.resolution() );  // define grid
    clipper::Xmap<float> sigmaa_all_map( hklinfo.spacegroup(), hklinfo.cell(), mygrid );          // define sigmaa best map
    clipper::Xmap<float> sigmaa_dif_map( hklinfo.spacegroup(), hklinfo.cell(), mygrid );          // define sigmaa diff  map
    clipper::Xmap<float> sigmaa_omit_fd( hklinfo.spacegroup(), hklinfo.cell(), mygrid );          // define sigmaa omit diff map
    clipper::Xmap<float> ligandmap( hklinfo.spacegroup(), hklinfo.cell(), mygrid );

    // scale data and flag R-free

    HRI ih;
    clipper::HKL_data<clipper::data32::Flag> flag( hklinfo );     // same flag for both calculations, omit absent reflections
    clipper::SFscale_aniso<float> sfscale;


    if(useParallelism)
        {
            pool.push([&sfscale, &fobs_scaled, &fc_all_bsc](int id)
            { 
                #if DUMP
                    std::cout << std::endl;
                    DBG << "Calculating sfscale from Thread ID: " << id << '.' << std::endl;
                #endif
                
                sfscale( fobs_scaled, fc_all_bsc );
            });

            pool.push([&fobs_scaled, &flag, &ih](int id)
            { 
                #if DUMP
                    DBG << "Calculating flag[ih].flag() from Thread ID: " << id << '.' << std::endl;
                #endif

                for ( ih = flag.first(); !ih.last(); ih.next() ) // we want to use all available reflections
                {
                    if ( !fobs_scaled[ih].missing() ) flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
                    else flag[ih].flag() = clipper::SFweight_spline<float>::NONE;
                }
            });
        }
        else
        {
            sfscale( fobs_scaled, fc_all_bsc );  // anisotropic scaling of Fobs. We scale Fobs to Fcalc instead of scaling our 3 Fcalcs to Fobs

            for ( ih = flag.first(); !ih.last(); ih.next() ) // we want to use all available reflections
            {
                if ( !fobs_scaled[ih].missing() ) flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
                else flag[ih].flag() = clipper::SFweight_spline<float>::NONE;
            }
        }

    int n_refln = 1000;
    int n_param = 20;
    double FobsFcalcSum = 0.0;
    double FobsFcalcAllSum = 0.0;
    double FobsSum = 0.0;

    clipper::HKL_data<clipper::data32::F_phi> fb_omit( hklinfo ); // new variables for omit sigmaa weighting calculation
    clipper::HKL_data<clipper::data32::F_phi> fd_omit( hklinfo );
    clipper::HKL_data<clipper::data32::Phi_fom> phiw_omit( hklinfo );
    clipper::HKL_data<clipper::data32::F_phi> fb_all( hklinfo ); // variables for all atom sigmaa weighting calculation
    clipper::HKL_data<clipper::data32::F_phi> fd_all( hklinfo );
    clipper::HKL_data<clipper::data32::Phi_fom> phiw_all( hklinfo );

    if (useParallelism)
    {
        #if DUMP
            std::cout << std::endl;
            DBG << "Number of jobs in the queue: " << pool.n_remaining_jobs() << std::endl;
        #endif
        
        std::this_thread::sleep_for(std::chrono::seconds(1));

        while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
            pool.sync();
        
        #if DUMP
            DBG << "Number of jobs in the queue: " << pool.n_remaining_jobs() << " after sync operation!" << std::endl;
        #endif
    }


    if(useParallelism)
    {
        pool.push([&fb_omit, &fd_omit, &phiw_omit, &fobs_scaled, &fc_omit_bsc, &flag, n_refln, n_param](int id)
        { 
            #if DUMP
                std::cout << std::endl;
                DBG << "Calculating sfw_omit from Thread ID: " << id << '.' << std::endl;
            #endif
            
            clipper::SFweight_spline<float> sfw_omit (n_refln, n_param );
            sfw_omit( fb_omit, fd_omit, phiw_omit, fobs_scaled, fc_omit_bsc, flag );
        });

        pool.push([&fb_all, &fd_all, &phiw_all, &fobs_scaled, &fc_all_bsc, &flag, n_refln, n_param](int id)
        { 
            #if DUMP
                DBG << "Calculating sfw_all from Thread ID: " << id << '.' << std::endl;
            #endif
            
            clipper::SFweight_spline<float> sfw_all( n_refln, n_param );
            sfw_all( fb_all,  fd_all,  phiw_all,  fobs_scaled, fc_all_bsc,  flag );
        });
    }
    else
    {
        std::cout << "Before calculating sfw omit and sfw all" << std::endl;
        clipper::SFweight_spline<float> sfw_omit (n_refln, n_param );
        sfw_omit( fb_omit, fd_omit, phiw_omit, fobs_scaled, fc_omit_bsc, flag );

        clipper::SFweight_spline<float> sfw_all( n_refln, n_param );
        sfw_all( fb_all,  fd_all,  phiw_all,  fobs_scaled, fc_all_bsc,  flag );
        std::cout << "After calculating sfw omit and sfw all" << std::endl;
    }
    
    // fb:          output best map coefficients
    // fd:          output difference map coefficients
    // phiw:        output phase and fom
    // fobs_scaled: input observed structure factors, previously scaled
    // fc_omit_bsc: input calculated omit data, with bulk solvent correction
    // fc_all_bsc:  input calculated data, with bsc


    if (useParallelism)
    {
        #if DUMP
            std::cout << std::endl;
            DBG << "Number of jobs in the queue: " << pool.n_remaining_jobs() << std::endl;
        #endif
        
        std::this_thread::sleep_for(std::chrono::seconds(1));

        while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
        {
            pool.sync();
        }
            
        
        #if DUMP
            DBG << "Number of jobs in the queue: " << pool.n_remaining_jobs() << " after sync operation!" << std::endl;
        #endif
    }

    std::vector<double> params( n_param, 2.0 );
    clipper::BasisFn_spline wrk_basis( hklinfo, n_param, 2.0 );

    clipper::TargetFn_scaleF1F2<clipper::data32::F_phi,clipper::data32::F_sigF> wrk_target_omit( fc_omit_bsc, fobs_scaled ); // was just fobs
    clipper::TargetFn_scaleF1F2<clipper::data32::F_phi,clipper::data32::F_sigF> wrk_target_all ( fc_all_bsc, fobs_scaled );
    clipper::ResolutionFn wrk_scale_omit( hklinfo, wrk_basis, wrk_target_omit, params );
    clipper::ResolutionFn wrk_scale_all ( hklinfo, wrk_basis, wrk_target_all,  params );

    double Fo, Fc_all, Fc_omit;

    if(useParallelism)
    {
        pool.push([&sigmaa_all_map, &fb_all](int id)
        { 
            #if DUMP
                std::cout << std::endl;
                DBG << "Calculating sigmaa_all_map from Thread ID: " << id << '.' << std::endl;
            #endif
            
            sigmaa_all_map.fft_from(fb_all);
        });

        pool.push([&sigmaa_dif_map, &fd_all](int id)
        { 
            #if DUMP
                DBG << "Calculating sigmaa_dif_map from Thread ID: " << id << '.' << std::endl;
            #endif
            
            sigmaa_dif_map.fft_from(fd_all);
        });     

        pool.push([&sigmaa_omit_fd, &fd_omit](int id)
        { 
            #if DUMP
                DBG << "Calculating sigmaa_omit_fd from Thread ID: " << id << '.' << std::endl;
            #endif
            
            sigmaa_omit_fd.fft_from(fd_omit);
        });  

        pool.push([&ligandmap, &fc_ligands_bsc](int id)
        { 
            ligandmap.fft_from(fc_ligands_bsc);
            
            #if DUMP
                DBG << "Calculating ligandmap from Thread ID: " << id << '.' << std::endl;
            #endif
        });

        pool.push([&fobs_scaled, &Fo, &Fc_all, &Fc_omit, &wrk_scale_all, &fc_all_bsc, &wrk_scale_omit, &fc_omit_bsc, &FobsFcalcSum, &FobsFcalcAllSum, &FobsSum](int id)
        { 
            #if DUMP
                DBG << "Calculating FobsSum, FobsFCalcSum and FobsFcalcAllSum from Thread ID: " << id << '.' << std::endl;
            #endif

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
        });
    }
    else
    {
        sigmaa_all_map.fft_from( fb_all );  // calculate the maps

        sigmaa_dif_map.fft_from( fd_all );

        sigmaa_omit_fd.fft_from( fd_omit );

        ligandmap.fft_from( fc_ligands_bsc );       // this is the map that will serve as Fc map for the RSCC calculation

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

    if (useParallelism)
    {
        #if DUMP
            std::cout << std::endl;
            DBG << "Number of jobs in the queue: " << pool.n_remaining_jobs() << std::endl;
        #endif
        
        std::this_thread::sleep_for(std::chrono::seconds(1));

        while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                pool.sync();
        
        #if DUMP
            DBG << "Number of jobs in the queue: " << pool.n_remaining_jobs() << " after sync operation!" << std::endl;
        #endif
    }

    if (((FobsFcalcAllSum / FobsSum)*10) > hklinfo.resolution().limit() + 0.6)
        std::cout << " Warning: R-work is unusually high. Please ensure that your PDB file contains full B-factors instead of residuals after TLS refinement!" << std::endl;

    float difference = (FobsFcalcSum / FobsSum) - (FobsFcalcAllSum / FobsSum);

    if (( difference > 0.15 ) || (clipper::Util::is_nan((FobsFcalcSum / FobsSum))))
    {
        useSigmaa = true;

        std::cout << std::endl << " The studied portions of the model account for a very significant part of the data. Calculating RSCC against a regular 2mFo-DFc map" << std::endl;
    }

    if(useParallelism)
    {
        int processedMonomers = 0;
        for (int index = 0; index < ligandList.size(); index++)
        {
            // if(pool.n_remaining_jobs() >= (pool.n_idle() - 1))
            //     pool.sync();

            if(pool.n_idle() == 0)
                pool.greedy_sync();

            pool.push([&sugarList, &path_to_model_file_clipper, &ligandList, &hklinfo, &mygrid, &sigmaa_all_map, &sigmaa_omit_fd, &ligandmap, &mgl, showGeom, ipradius, pos_slash, index, useSigmaa](int id)
            { 
                #if DUMP
                    std::cout << std::endl;
                    DBG << "Calculating RSCC score from Thread ID: " << id << " for nth " << index << " index out of " << ligandList.size() << " total indices." << std::endl;
                #endif
                
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

                // now calculate the correlation between the weighted experimental & calculated maps
                // maps are scanned only inside a sphere containing the sugar for performance reasons,
                // although RSCC and <RMS> are restricted to a mask surrounding the model

                //////// mask calculation //////////
                clipper::Xmap<float> mask( hklinfo.spacegroup(), hklinfo.cell(), mygrid );
                clipper::EDcalc_mask<float> masker( ipradius );
                masker(mask, sugarList[index].atom_list());

                ////////////////////////////////////

                clipper::Coord_orth origin(minX-2,minY-2,minZ-2);
                clipper::Coord_orth destination(maxX+2,maxY+2,maxZ+2);


                double accum = 0.0;
                double corr_coeff = 0.0;
                std::pair<double, double> rscc_and_accum;

                rscc_and_accum = privateer::xray::calculate_rscc(sigmaa_all_map, sigmaa_omit_fd, ligandmap, mask, hklinfo, mygrid, origin, destination, useSigmaa);
                corr_coeff = rscc_and_accum.first;
                accum = rscc_and_accum.second;

                ligandList[index].second.set_rscc ( corr_coeff );
                ligandList[index].second.set_accum_score ( accum );
                // calculation of the mean densities of the calc (ligandmap) and weighted obs (sigmaamap) maps

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
                            }
                            else if ( list_of_glycans[i].get_type() == "c-glycan" )
                            {
                                ligandList[index].second.set_context ( "c-glycan" );
                            }
                            else if ( list_of_glycans[i].get_type() == "o-glycan" )
                            {
                                ligandList[index].second.set_context ( "o-glycan" );
                            }
                            else if ( list_of_glycans[i].get_type() == "s-glycan" )
                            {
                                ligandList[index].second.set_context ( "s-glycan" );
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
                }


                bool occupancy_check = false;
                std::vector<clipper::MAtom> ringcomponents = ligandList[index].second.ring_members();

                for ( int i = 0 ; i < ringcomponents.size() ; i++ )
                    if (privateer::util::get_altconformation(ringcomponents[i]) != ' ')
                        occupancy_check = true;

                ligandList[index].second.set_occupancy_check ( occupancy_check );
                
            });
            processedMonomers++;

            #if DUMP
                std::cout << std::endl;
                DBG << "Processed " << processedMonomers << "/" << ligandList.size() << " monomers..." << std::endl;
            #endif

        }
        #if DUMP
            std::cout << std::endl;
            DBG << "Number of jobs in the queue: " << pool.n_remaining_jobs() << std::endl;
        #endif

        std::this_thread::sleep_for(std::chrono::seconds(1));
        
        while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
            pool.sync();
        
        #if DUMP
            DBG << "Number of jobs in the queue: " << pool.n_remaining_jobs() << " after sync operation!" << std::endl;
        #endif
    }
    else
    {
        int processedMonomers = 0;
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

            // now calculate the correlation between the weighted experimental & calculated maps
            // maps are scanned only inside a sphere containing the sugar for performance reasons,
            // although RSCC and <RMS> are restricted to a mask surrounding the model

            //////// mask calculation //////////

            clipper::Xmap<float> mask( hklinfo.spacegroup(), hklinfo.cell(), mygrid );

            clipper::EDcalc_mask<float> masker( ipradius );
            masker(mask, sugarList[index].atom_list());

            ////////////////////////////////////

            clipper::Coord_orth origin(minX-2,minY-2,minZ-2);
            clipper::Coord_orth destination(maxX+2,maxY+2,maxZ+2);

            double accum = 0.0;
            double corr_coeff = 0.0;
            std::pair<double, double> rscc_and_accum;

            rscc_and_accum = privateer::xray::calculate_rscc(sigmaa_all_map, sigmaa_omit_fd, ligandmap, mask, hklinfo, mygrid, origin, destination, useSigmaa);

            corr_coeff = rscc_and_accum.first;
            accum = rscc_and_accum.second;

            ligandList[index].second.set_rscc ( corr_coeff );
            ligandList[index].second.set_accum_score ( accum );
            // calculation of the mean densities of the calc (ligandmap) and weighted obs (sigmaamap) maps

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
                        }
                        else if ( list_of_glycans[i].get_type() == "c-glycan" )
                        {
                            ligandList[index].second.set_context ( "c-glycan" );
                        }
                        else if ( list_of_glycans[i].get_type() == "o-glycan" )
                        {
                            ligandList[index].second.set_context ( "o-glycan" );
                        }
                        else if ( list_of_glycans[i].get_type() == "s-glycan" )
                        {
                            ligandList[index].second.set_context ( "s-glycan" );
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
            }

        
            bool occupancy_check = false;
            std::vector<clipper::MAtom> ringcomponents = ligandList[index].second.ring_members();

            for ( int i = 0 ; i < ringcomponents.size() ; i++ )
                if (privateer::util::get_altconformation(ringcomponents[i]) != ' ')
                    occupancy_check = true;

            ligandList[index].second.set_occupancy_check ( occupancy_check );

            processedMonomers++;

            #if DUMP
                std::cout << std::endl;
                DBG << "Processed " << processedMonomers << "/" << ligandList.size() << " monomers..." << std::endl;
            #endif

        }
    }

        

    this->finalLigandList = ligandList;
    this->sugar_summary_of_experimental_data = generate_sugar_experimental_data_summary(finalLigandList);
}

// private methods //
pybind11::list privateer::pyanalysis::XRayData::generate_sugar_experimental_data_summary(std::vector<std::pair< clipper::String , clipper::MSugar>>& finalLigandList)
{
    pybind11::list output; 

    for(int index = 0; index < finalLigandList.size(); index++)
    {
        std::string ccdCode = finalLigandList[index].second.type().c_str();
        std::string chainID = finalLigandList[index].first;
        int sugarID = index;
        int pdbID = std::stoi(finalLigandList[index].second.id().trim());
        float RSCC = finalLigandList[index].second.get_rscc();
        float accum = finalLigandList[index].second.get_accum();

        std::string sugardiagnostic;
        if(finalLigandList[index].second.is_sane())
        {
            if(!finalLigandList[index].second.ok_with_conformation())
            {
                sugardiagnostic = "check";
            }
            else
                sugardiagnostic = "yes";
        }
        else
            sugardiagnostic = "no";

        auto currentSugar = pybind11::dict ("three_letter_code"_a=ccdCode, "Chain"_a=chainID, "sugar_index_internal"_a=sugarID, "PDB_ID"_a=pdbID, "RSCC"_a=RSCC, "mFo"_a=accum, "privateer_diagnostic"_a=sugardiagnostic);
        output.append(currentSugar);
    }
    return output;
}
// private methods end //

//Needs to return a list of RSCC and accumm values. Maybe a pybind11::list of dicts - ccdCode, chainID, sugarID, pdbID, RSCC, accum. 
//Need a private method that acts like privateer::util::print_monosaccharide_summary_python, but creates this python object instead. 


///////////////////////////////////////////////// Class XrayData END ////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////// PYBIND11 BINDING DEFINITIONS ////////////////////////////////////////////////////////////////////
namespace py = pybind11;
namespace pa = privateer::pyanalysis;

void init_pyanalysis(py::module& m)
{
    py::class_<pa::GlycosylationComposition>(m, "GlycosylationComposition")
        .def(py::init<>())
        .def(py::init<std::string&, std::string>(), py::arg("path_to_model_file")="undefined", py::arg("expression_system")="undefined")
        .def("get_path_of_model_file_used",  &pa::GlycosylationComposition::get_path_of_model_file_used)
        .def("get_expression_system_used",  &pa::GlycosylationComposition::get_expression_system_used)
        .def("get_number_of_glycan_chains_detected",  &pa::GlycosylationComposition::get_number_of_glycan_chains_detected)
        .def("get_summary_of_detected_glycans",  &pa::GlycosylationComposition::get_summary_of_detected_glycans)
        .def("get_glycan",  &pa::GlycosylationComposition::get_glycan);

    py::class_<pa::GlycanStructure>(m, "GlycanStructure")
        .def(py::init<>())
        .def(py::init<const clipper::MGlycology&, const int, privateer::pyanalysis::GlycosylationComposition&>())
        .def("get_glycan_id", &pa::GlycanStructure::get_glycan_id)
        .def("get_total_number_of_sugars", &pa::GlycanStructure::get_total_number_of_sugars)
        .def("get_wurcs_notation", &pa::GlycanStructure::get_wurcs_notation)
        .def("get_unique_monosaccharide_codes", &pa::GlycanStructure::get_unique_monosaccharide_codes)
        .def("get_total_of_glycosidic_bonds", &pa::GlycanStructure::get_total_of_glycosidic_bonds)
        .def("get_glycosylation_type", &pa::GlycanStructure::get_glycosylation_type)
        .def("get_root_info", &pa::GlycanStructure::get_root_info)
        .def("get_protein_glycan_linkage_torsions", &pa::GlycanStructure::get_protein_glycan_linkage_torsions)
        .def("get_glycan_summary", &pa::GlycanStructure::get_glycan_summary)
        .def("get_monosaccharide", &pa::GlycanStructure::get_monosaccharide)
        .def("get_all_monosaccharides", &pa::GlycanStructure::get_all_monosaccharides);

    py::class_<pa::CarbohydrateStructure>(m, "CarbohydrateStructure")
        .def(py::init<>())
        .def(py::init<clipper::MGlycan&, const int, const int, privateer::pyanalysis::GlycosylationComposition&, privateer::pyanalysis::GlycanStructure&>())
        .def(py::self == py::self)
        .def("get_sugar_summary", &pa::CarbohydrateStructure::get_sugar_summary)
        .def("get_sugar_id", &pa::CarbohydrateStructure::get_sugar_id)
        .def("get_glycan_id", &pa::CarbohydrateStructure::get_glycan_id)
        .def("get_sugar_pdb_id", &pa::CarbohydrateStructure::get_sugar_pdb_id)
        .def("get_sugar_pdb_chain", &pa::CarbohydrateStructure::get_sugar_pdb_chain)
        .def("get_conformation_name", &pa::CarbohydrateStructure::get_conformation_name)
        .def("get_conformation_name_iupac", &pa::CarbohydrateStructure::get_conformation_name_iupac)
        .def("get_puckering_amplitude", &pa::CarbohydrateStructure::get_puckering_amplitude)
        .def("get_anomer", &pa::CarbohydrateStructure::get_anomer)
        .def("get_handedness", &pa::CarbohydrateStructure::get_handedness)
        .def("get_denomination", &pa::CarbohydrateStructure::get_denomination)
        .def("get_ring_cardinality", &pa::CarbohydrateStructure::get_ring_cardinality)
        .def("get_cremer_pople_params", &pa::CarbohydrateStructure::get_cremer_pople_params)
        .def("is_sane", &pa::CarbohydrateStructure::is_sane)
        .def("get_privateer_diagnostic", &pa::CarbohydrateStructure::get_privateer_diagnostic)
        .def("get_name_full", &pa::CarbohydrateStructure::get_name_full)
        .def("get_name_short", &pa::CarbohydrateStructure::get_name_short)
        .def("get_type", &pa::CarbohydrateStructure::get_type)
        .def("get_ring_angles", &pa::CarbohydrateStructure::get_ring_angles)
        .def("get_ring_bonds", &pa::CarbohydrateStructure::get_ring_bonds)
        .def("get_ring_torsion", &pa::CarbohydrateStructure::get_ring_torsion)
        .def("get_ring_bond_rmsd", &pa::CarbohydrateStructure::get_ring_bond_rmsd)
        .def("get_ring_angle_rmsd", &pa::CarbohydrateStructure::get_ring_angle_rmsd)
        .def("get_bfactor", &pa::CarbohydrateStructure::get_bfactor)
        .def("is_supported", &pa::CarbohydrateStructure::is_supported)
        .def("ok_with_ring", &pa::CarbohydrateStructure::ok_with_ring)
        .def("ok_with_bonds_rmsd", &pa::CarbohydrateStructure::ok_with_bonds_rmsd)
        .def("ok_with_angles_rmsd", &pa::CarbohydrateStructure::ok_with_angles_rmsd)
        .def("ok_with_anomer", &pa::CarbohydrateStructure::ok_with_anomer)
        .def("ok_with_chirality", &pa::CarbohydrateStructure::ok_with_chirality)
        .def("ok_with_conformation", &pa::CarbohydrateStructure::ok_with_conformation)
        .def("ok_with_puckering", &pa::CarbohydrateStructure::ok_with_puckering)
        .def("get_glycosylation_context", &pa::CarbohydrateStructure::get_glycosylation_context);

    py::class_<pa::XRayData>(m, "XRayData")
        .def(py::init<>())
        .def(py::init<std::string&, std::string&, std::string&, float, int>(), py::arg("path_to_mtz_file")="undefined", py::arg("path_to_model_file")="undefined", py::arg("input_column_fobs_user")="NONE", py::arg("ipradius")=2.5, py::arg("nThreads")=-1)
        .def("get_sugar_summary_with_experimental_data", &pa::XRayData::get_sugar_summary_with_experimental_data)
        .def("print_cpp_console_output_summary", &pa::XRayData::print_cpp_console_output_summary);
}

///////////////////////////////////////////////// PYBIND11 BINDING DEFINITIONS END////////////////////////////////////////////////////////////////////