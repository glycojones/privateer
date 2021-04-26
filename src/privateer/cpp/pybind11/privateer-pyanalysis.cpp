// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#include "privateer-pyanalysis.h"

using namespace pybind11::literals;


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

    auto glycanObject = privateer::pyanalysis::GlycanStructure(mgl, glycanID);
    
    return glycanObject;
}
///////////////////////////////////////////////// Class GlycosylationComposition END ////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////// Class GlycanStructure ////////////////////////////////////////////////////////////////////
void privateer::pyanalysis::GlycanStructure::pyinit( const clipper::MGlycology& mgl, const int glycanID)
{
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
        auto sugarObject = privateer::pyanalysis::CarbohydrateStructure(inputGlycan, i, glycanID);
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

    auto sugarObject = privateer::pyanalysis::CarbohydrateStructure(glycan, sugarID, glycanID);
    
    return sugarObject;
}
///////////////////////////////////////////////// Class GlycanStructure END ////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////// Class CarbohydrateStructure  ////////////////////////////////////////////////////////////////////
void privateer::pyanalysis::CarbohydrateStructure::pyinit( clipper::MGlycan& mglycan, const int sugarID, const int glycanID )
{
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

void privateer::pyanalysis::XRayData::read_from_file( std::string& path_to_mtz_file, std::string& path_to_model_file, std::string& input_column_fobs_user) {
    
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
    clipper::String path_to_mtz_file_clipper = path_to_mtz_file;
    clipper::MMDBfile mfile;
    clipper::CCP4MTZfile ampmtzin;
    clipper::MTZcrystal opxtal;
    clipper::MTZdataset opdset;

    clipper::HKL_data<clipper::data32::F_sigF> fobs;            // allocate space for F and sigF
    clipper::HKL_data<clipper::data32::F_sigF> fobs_scaled;     // allocate space for scaled F and sigF
    clipper::HKL_data<clipper::data32::F_phi> fc_omit_bsc;      // allocate space for the omit calc data with bsc
    clipper::HKL_data<clipper::data32::F_phi> fc_all_bsc;       // allocate space for the whole calculated data with bsc
    clipper::HKL_data<clipper::data32::F_phi> fc_ligands_bsc;    // allocate space for the ligand calculated data

    privateer::util::read_coordinate_file_mtz(mfile, mmol, path_to_model_file_clipper, true);
    privateer::xray::read_xray_map( path_to_mtz_file_clipper, path_to_model_file_clipper, mmol, hklinfo, mtzin );

    fobs = clipper::HKL_data<clipper::data32::F_sigF> ( hklinfo );
    fobs_scaled = clipper::HKL_data<clipper::data32::F_sigF> ( hklinfo );
    fc_omit_bsc = clipper::HKL_data<clipper::data32::F_phi> ( hklinfo );
    fc_all_bsc = clipper::HKL_data<clipper::data32::F_phi> ( hklinfo );
    fc_ligands_bsc = clipper::HKL_data<clipper::data32::F_phi> ( hklinfo );

    std::cout << "Value of input_column_fobs_user = " << input_column_fobs_user << std::endl;
    std::cout << "Value of input_column_fobs = " << input_column_fobs << std::endl;
    privateer::xray::initialize_experimental_dataset( mtzin, ampmtzin, input_column_fobs, fobs, hklinfo, opxtal, opdset, path_to_mtz_file_clipper);

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
            // if (allSugars)
            // {
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

            // }

            // else
            // {
            //     if ( strncmp( mmol[p][m].type().c_str(), input_ccd_code.trim().c_str(), 3 )) // true if strings are different
            //     {
            //         for (int id = 0; id < mmol[p][m].size(); id++ )
            //         {
            //             mainAtoms.push_back(mmol[p][m][id]); // cycle through atoms and copy them
            //             allAtoms.push_back(mmol[p][m][id]);
            //         }
            //     }
            //     else // it's the one sugar we're looking to omit
            //     {
            //         if ( input_validation_options.size() > 0 )
            //         {
            //             const clipper::MSugar msug ( mmol, mmol[p][m], manb, external_validation );

            //             sugarList.push_back(mmol[p][m]);
            //             ligandList.push_back(std::pair<clipper::String, clipper::MSugar> (mmol[p].id().trim(), msug));
            //         }
            //         else
            //         {
            //             const clipper::MSugar msug(mmol, mmol[p][m], manb);

            //             sugarList.push_back(mmol[p][m]);
            //             ligandList.push_back(std::pair<clipper::String, clipper::MSugar> (mmol[p].id().trim(), msug));

            //         }
            //         for (int id = 0; id < mmol[p][m].size(); id++ )
            //         {
            //             ligandAtoms.push_back(mmol[p][m][id]);  // add the ligand atoms to a second array
            //             allAtoms.push_back(mmol[p][m][id]);
            //         }
            //     }
            // }
        }
    }

    std::cout << std::endl << " " << fobs.num_obs() << " reflections have been loaded";
            std::cout << std::endl << std::endl << " Resolution " << hklinfo.resolution().limit() << "Å" << std::endl << hklinfo.cell().format() << std::endl;
}

///////////////////////////////////////////////// Class GlycosylationComposition END ////////////////////////////////////////////////////////////////////


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
        .def(py::init<const clipper::MGlycology&, const int>())
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
        .def(py::init<clipper::MGlycan&, const int, const int>())
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
        .def(py::init<std::string&, std::string&, std::string&>(), py::arg("path_to_mtz_file")="undefined", py::arg("path_to_model_file")="undefined", py::arg("input_column_fobs_user")="NONE");
}

///////////////////////////////////////////////// PYBIND11 BINDING DEFINITIONS END////////////////////////////////////////////////////////////////////