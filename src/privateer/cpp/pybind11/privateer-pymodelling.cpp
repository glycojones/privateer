// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#include "privateer-pymodelling.h"

#define DBG std::cout << "[" << __FUNCTION__ << "] - "


privateer::pymodelling::Builder::Builder(std::string& path_to_receiving_model_file, bool enable_user_messages) 
{
    this->trim_donor_when_clashes_detected = false;
    this->remove_donor_when_clashes_detected = false;
    this->ANY_search_policy = false;
    this->enable_user_messages = enable_user_messages;
    this->debug_output = false;
    this->nThreads = -1;
    this->import_receiving_model_only(path_to_receiving_model_file);
};

privateer::pymodelling::Builder::Builder(std::string& path_to_receiving_model_file, std::string& path_to_donor_model_file, int nThreads, bool trim_donor_when_clashes_detected, bool remove_donor_when_clashes_detected, bool ANY_search_policy, bool enable_user_messages, bool debug_output) 
{
    this->trim_donor_when_clashes_detected = trim_donor_when_clashes_detected;
    this->remove_donor_when_clashes_detected = remove_donor_when_clashes_detected;
    this->ANY_search_policy = ANY_search_policy;
    this->enable_user_messages = enable_user_messages;
    this->debug_output = debug_output;
    this->nThreads = nThreads;
    this->read_from_file ( path_to_receiving_model_file, path_to_donor_model_file, nThreads, trim_donor_when_clashes_detected, remove_donor_when_clashes_detected, ANY_search_policy, enable_user_messages, debug_output );
};

std::string privateer::pymodelling::Builder::convert_three_letter_code_to_single_letter(std::string three_letter_code)
{
    std::unordered_map<std::string, std::string> code_conversion = 
    {
        {"ALA", "A"}, {"ARG", "R"}, {"ASN", "N"}, {"ASP", "D"},
        {"CYS", "C"}, {"GLN", "Q"}, {"GLU", "E"}, {"GLY", "G"},
        {"HIS", "H"}, {"ILE", "I"}, {"LEU", "L"}, {"LYS", "K"},
        {"MET", "M"}, {"PHE", "F"}, {"PRO", "P"}, {"SER", "S"},
        {"THR", "T"}, {"TRP", "W"}, {"TYR", "Y"}, {"VAL", "V"},
        {"HYP", "P"}, {"LYZ", "K"}, {"SEP", "S"}
    };

    std::unordered_map<std::string, std::string>::const_iterator result = code_conversion.find(three_letter_code);

    if ( result == code_conversion.end() )
        return three_letter_code + "_";
    else
        return result->second;
}

std::string privateer::pymodelling::Builder::convert_three_letter_code_to_single_letter_for_grafter(std::string three_letter_code)
{
    std::unordered_map<std::string, std::string> code_conversion = 
    {
        {"ALA", "A"}, {"ARG", "R"}, {"ASN", "N"}, {"ASP", "D"},
        {"CYS", "C"}, {"GLN", "Q"}, {"GLU", "E"}, {"GLY", "G"},
        {"HIS", "H"}, {"ILE", "I"}, {"LEU", "L"}, {"LYS", "K"},
        {"MET", "M"}, {"PHE", "F"}, {"PRO", "P"}, {"SER", "S"},
        {"THR", "T"}, {"TRP", "W"}, {"TYR", "Y"}, {"VAL", "V"},
        {"HYP", "P"}, {"LYZ", "K"}, {"SEP", "S"}
    };

    std::unordered_map<std::string, std::string>::const_iterator result = code_conversion.find(three_letter_code);

    if ( result == code_conversion.end() )
        return "_";
    else
        return result->second;
}

void privateer::pymodelling::Builder::import_receiving_model_only( std::string& path_to_receiving_model_file) {

    if(path_to_receiving_model_file == "undefined")
    {
        throw std::invalid_argument( "No path was provided for model file input! Aborting." );
    }

    this->path_to_receiving_model = path_to_receiving_model_file;
    this->path_to_donor_model = "undefined";
    
    clipper::MiniMol receiver_mmol;
    clipper::MMDBfile receiver_mfile;
    clipper::String path_to_model_file_clipper_receiver = path_to_receiving_model_file;

    privateer::util::read_coordinate_file_mtz(receiver_mfile, receiver_mmol, path_to_model_file_clipper_receiver, true);

    if(enable_user_messages && !debug_output)
    {
        std::cout << "Importing receiver model ONLY with the path of: " << path_to_model_file_clipper_receiver << std::endl;
        std::cout << "Number of chains detected in receiver model: " << receiver_mmol.size() << std::endl;
        std::cout << "Donor model was not imported!" << std::endl;
    }


    this->imported_receiving_model = receiver_mmol;
    this->imported_receiving_model_seq_info = get_protein_sequence_information(imported_receiving_model);
}

void privateer::pymodelling::Builder::read_from_file( std::string& path_to_receiving_model_file, std::string& path_to_donor_model_file, int nThreads, bool trim_donor_when_clashes_detected, bool remove_donor_when_clashes_detected, bool ANY_search_policy, bool enable_user_messages, bool debug_output ) {

    if(path_to_receiving_model_file == "undefined")
    {
        throw std::invalid_argument( "No path was provided for model file input! Aborting." );
    }

    if(path_to_donor_model_file == "undefined")
    {
        throw std::invalid_argument( "No path was provided for model file input! Aborting." );
    }


    int detectedThreads = std::thread::hardware_concurrency();
    bool useParallelism = true;
    
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
        std::cout << std::endl << "THREADING: Privateer will use " << nThreads << " of the " << detectedThreads << " available threads..." << std::endl;
    }

    this->path_to_receiving_model = path_to_receiving_model_file;
    this->path_to_donor_model = path_to_donor_model_file;
    
    clipper::MiniMol receiver_mmol;
    clipper::MiniMol donor_mmol;
    clipper::MMDBfile receiver_mfile;
    clipper::MMDBfile donor_mfile;
    clipper::String path_to_model_file_clipper_receiver = path_to_receiving_model_file;
    clipper::String path_to_model_file_clipper_donor = path_to_donor_model_file;

    privateer::util::read_coordinate_file_mtz(receiver_mfile, receiver_mmol, path_to_model_file_clipper_receiver, true);
    if(debug_output)
    {
        DBG << "Imported receiver model with the path of: " << path_to_model_file_clipper_receiver << std::endl;
        DBG << "Number of chains detected in receiver model: " << receiver_mmol.size() << std::endl;
    }
    
    if(enable_user_messages && !debug_output)
    {
        std::cout << "Imported receiver model with the path of: " << path_to_model_file_clipper_receiver << std::endl;
        std::cout << "Number of chains detected in receiver model: " << receiver_mmol.size() << std::endl;
    }

    privateer::util::read_coordinate_file_mtz(donor_mfile, donor_mmol, path_to_model_file_clipper_donor, true);
    if(debug_output)
    {
        DBG << "Imported donor model with the path of: " << path_to_model_file_clipper_donor << std::endl;
        DBG << "Number of chains detected in donor model: " << donor_mmol.size() << std::endl;
    }
    
    if(enable_user_messages && !debug_output)
    {
        std::cout << "Imported donor model with the path of: " << path_to_model_file_clipper_donor << std::endl;
        std::cout << "Number of chains detected in donor model: " << donor_mmol.size() << std::endl;
    }


    this->imported_receiving_model = receiver_mmol;
    this->imported_donor_model = donor_mmol;
    this->imported_receiving_model_seq_info = get_protein_sequence_information(imported_receiving_model);
    
    if(debug_output)
        DBG << "Generated a pybind11 list with receiver sequence info. Length of output list: " << imported_receiving_model_seq_info.size() << std::endl;

    privateer::pyanalysis::GlycosylationComposition glycosylation(path_to_donor_model, "undefined", true, false);
    this->glycan_summary_donor = glycosylation.get_summary_of_detected_glycans();
    
    if(debug_output)
        DBG << "Generated a pybind11 list with summary of glycans detected in donor. Length of output list: " << glycan_summary_donor.size() << std::endl;

    this->grafter = privateer::modelling::Grafter(imported_receiving_model, imported_donor_model, nThreads, trim_donor_when_clashes_detected, remove_donor_when_clashes_detected, enable_user_messages, debug_output);
    
}

// return single letter code also
pybind11::list privateer::pymodelling::Builder::get_protein_sequence_information( clipper::MiniMol& input ) 
{
    pybind11::list output;
    for(int i = 0; i < input.size(); i++)
    {
        pybind11::list residues;
        std::string complete_polymer_seq;
        for(int j = 0; j < input[i].size(); j++)
        {
            std::string residueType = input[i][j].type().trim();
            std::string residueCode = convert_three_letter_code_to_single_letter_for_grafter(residueType);
            if (j + 1 == input[i].size() && residueCode.back() == '_')
                residueCode.pop_back();
            int residueSeqnum = input[i][j].seqnum();
            auto residue_dict = pybind11::dict("index"_a=j, "residueType"_a=residueType, "residueCode"_a=residueCode, "residueSeqnum"_a=residueSeqnum);
            complete_polymer_seq += residueCode;
            residues.append(residue_dict);
        }
        std::string chainID = input[i].id().trim();
        auto mpolymer_dict = pybind11::dict("index"_a = i, "ChainID"_a=chainID, "Sequence"_a=complete_polymer_seq, "Residues"_a=residues);
        output.append(mpolymer_dict);
    }
    return output;
}

void privateer::pymodelling::Builder::graft_glycan_to_receiver(int mglycanindex, int receiver_chain_index, int received_residue_index)
{
    std::vector<clipper::MGlycan> donor_glycans = grafter.get_donor_glycans();
    clipper::MGlycan glycan_to_graft = donor_glycans[mglycanindex];
    
    if(debug_output)
    {
        DBG << "Grafting glycan with the WURCS of: " << glycan_to_graft.generate_wurcs() << std::endl;
    }
    
    if(enable_user_messages && !debug_output)
    {
        std::cout << "Grafting glycan with the WURCS of:\n" << glycan_to_graft.generate_wurcs() << std::endl;
    }

    clipper::MMonomer receiver_monomer = imported_receiving_model[receiver_chain_index][received_residue_index];
    
    if(debug_output)
    {
        DBG << "Grafting input glycan to: " << imported_receiving_model[receiver_chain_index].id().trim() << "/" << imported_receiving_model[receiver_chain_index][received_residue_index].type().trim() << "-" << imported_receiving_model[receiver_chain_index][received_residue_index].id().trim() << "\tseqnum = " << imported_receiving_model[receiver_chain_index][received_residue_index].seqnum() << std::endl;
    }
    
    if(enable_user_messages && !debug_output)
    {
        std::cout << "Grafting input glycan to: " << imported_receiving_model[receiver_chain_index].id().trim() << "/" << imported_receiving_model[receiver_chain_index][received_residue_index].type().trim() << "-" << imported_receiving_model[receiver_chain_index][received_residue_index].id().trim() << std::endl;
    }

    if(userValuesChanged)
    {
        grafter.setPhi(userPhi);
        grafter.setPsi(userPsi);
        grafter.setPhi_error(userPhi_error);
        grafter.setPhi_error(userPsi_error);
        grafter.setIteration_step(userIteration_step);
    }

    grafter.graft_mpolymer_to_receiving_model(glycan_to_graft, receiver_monomer, imported_receiving_model[receiver_chain_index].id().trim(), ANY_search_policy);
    
    std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MMonomer, clipper::String> > > grafted_clashes = grafter.get_grafted_clashes();

    double currentAverageDistanceBetweenResidues = 0.0;
    double totalAveragesOfPerResidueDistances = 0.0;

    auto glycanClashes = pybind11::list();
    for(int i = 0; i < grafted_clashes.size(); i++)
    {
        int n_protein_side_chain_atoms = grafted_clashes[i].first.first.size();
        int n_sugar_atoms = grafted_clashes[i].second.first.size();
        double averageDistance = 0.0;
        double totalDistance = 0.0;
        int totalIterations = 0;

        for(int sugarAtom = 0; sugarAtom < n_sugar_atoms; sugarAtom++)
        {
            for(int proteinResidueAtom = 0; proteinResidueAtom < n_protein_side_chain_atoms; proteinResidueAtom++)
            {
                double currentDistance = clipper::Coord_orth::length(grafted_clashes[i].first.first[proteinResidueAtom].coord_orth(), grafted_clashes[i].second.first[sugarAtom].coord_orth());
                totalDistance = totalDistance + currentDistance;
                totalIterations++;
            }
        }

        averageDistance = totalDistance / totalIterations;
        totalAveragesOfPerResidueDistances = totalAveragesOfPerResidueDistances + averageDistance;

        std::string protein_chain = grafted_clashes[i].first.second;
        std::string protein_residue_id = grafted_clashes[i].first.first.id().trim();
        std::string protein_residue_type = grafted_clashes[i].first.first.type().trim();

        std::string sugar_chain = grafted_clashes[i].second.second;
        std::string sugar_residue_id = grafted_clashes[i].second.first.id().trim();
        std::string sugar_residue_type = grafted_clashes[i].second.first.type().trim();

        auto protein_residue_dict = pybind11::dict("Chain"_a=protein_chain, "ResidueID"_a=protein_residue_id, "ResidueType"_a=protein_residue_type);
        auto sugar_residue_dict = pybind11::dict("Chain"_a=sugar_chain, "ResidueID"_a=sugar_residue_id, "ResidueType"_a=sugar_residue_type);
        auto clash_dict = pybind11::dict("index"_a=i, "AvgAtomicDistance"_a=averageDistance, "ProteinSideChainResidue"_a=protein_residue_dict, "GraftedSugar"_a=sugar_residue_dict);

        glycanClashes.append(clash_dict);
    }

    if(!grafted_clashes.empty())
    {
        currentAverageDistanceBetweenResidues = totalAveragesOfPerResidueDistances / grafted_clashes.size();
        std::string graftedGlycanWURCS = grafter.get_grafted_glycan().generate_wurcs();
        std::string proteinChainID = imported_receiving_model[receiver_chain_index].id().trim();
        std::string receivermonomertype = receiver_monomer.type().trim();
        std::string receivermonomerID = receiver_monomer.id().trim();
        
        std::string grafted_glycan_chainID = grafter.get_grafted_glycan_chainID();
        std::string rootSugarType = grafter.get_grafted_glycan().get_node(0).get_sugar().type().trim();
        std::string rootSugarPDBID = grafter.get_grafted_glycan().get_node(0).get_sugar().id().trim();

        float graftedPhi = grafter.getGraftedPhi();
        float graftedPsi = grafter.getGraftedPsi();

        auto glycosidic_linkage_torsions = pybind11::dict("Phi"_a = graftedPhi, "Psi"_a=graftedPsi);

        bool graftStatus = grafter.get_graft_status();
        
        auto glycan_clash_dict = pybind11::dict("index"_a=summary_of_grafted_glycans.size(), "graftedGlycanWURCS"_a=graftedGlycanWURCS, "receiver_chain_index"_a=receiver_chain_index, "receiving_protein_residue_chain_PDBID"_a=proteinChainID, "receiver_residue_index"_a=received_residue_index, "receiving_protein_residue_monomer_PDBID"_a=receivermonomerID, "receiving_protein_residue_monomer_type"_a=receivermonomertype, "glycan_grafted_as_chainID"_a=grafted_glycan_chainID, "donor_glycan_root_PDBID"_a=rootSugarPDBID, "donor_glycan_root_type"_a=rootSugarType, "glycosidic_linkage_torsions"_a=glycosidic_linkage_torsions, "AvgTotalAtomicDistance"_a=currentAverageDistanceBetweenResidues, "ClashingResidues"_a=glycanClashes, "GraftStatus"_a=graftStatus);
        summary_of_grafted_glycans.append(glycan_clash_dict);
        

    }
    else
    {
        std::string graftedGlycanWURCS = grafter.get_grafted_glycan().generate_wurcs();
        std::string proteinChainID = imported_receiving_model[receiver_chain_index].id().trim();
        std::string receivermonomertype = receiver_monomer.type().trim();
        std::string receivermonomerID = receiver_monomer.id().trim();
        
        std::string grafted_glycan_chainID = grafter.get_grafted_glycan_chainID();
        std::string rootSugarType = grafter.get_grafted_glycan().get_node(0).get_sugar().type().trim();
        std::string rootSugarPDBID = grafter.get_grafted_glycan().get_node(0).get_sugar().id().trim();

        float graftedPhi = grafter.getGraftedPhi();
        float graftedPsi = grafter.getGraftedPsi();

        auto glycosidic_linkage_torsions = pybind11::dict("Phi"_a = graftedPhi, "Psi"_a=graftedPsi);

        bool graftStatus = grafter.get_graft_status();
        
        auto glycan_clash_dict = pybind11::dict("index"_a=summary_of_grafted_glycans.size(), "graftedGlycanWURCS"_a=graftedGlycanWURCS, "receiver_chain_index"_a=receiver_chain_index, "receiving_protein_residue_chain_PDBID"_a=proteinChainID, "receiver_residue_index"_a=received_residue_index, "receiving_protein_residue_monomer_PDBID"_a=receivermonomerID, "receiving_protein_residue_monomer_type"_a=receivermonomertype, "glycan_grafted_as_chainID"_a=grafted_glycan_chainID, "donor_glycan_root_PDBID"_a=rootSugarPDBID, "donor_glycan_root_type"_a=rootSugarType, "glycosidic_linkage_torsions"_a=glycosidic_linkage_torsions, "AvgTotalAtomicDistance"_a=pybind11::cast<pybind11::none>(Py_None), "ClashingResidues"_a=glycanClashes, "GraftStatus"_a=graftStatus);
        summary_of_grafted_glycans.append(glycan_clash_dict);
    }
    
    this->export_model = grafter.get_final_receiving_model();
}

void privateer::pymodelling::Builder::export_grafted_model(std::string& output_path)
{
    clipper::MMDBfile pdbfile;
    pdbfile.export_minimol( export_model );
    std::cout << "Attempting to write PDB file with grafted glycans in the following path: " << output_path << std::endl;
    pdbfile.write_file( output_path );
}




///////////////////////////////////////////////// PYBIND11 BINDING DEFINITIONS ////////////////////////////////////////////////////////////////////
namespace py = pybind11;
namespace pm = privateer::pymodelling;

void init_pymodelling(py::module& m)
{
    py::class_<pm::Builder>(m, "Builder")
        .def(py::init<>())
        .def(py::init<std::string&, bool>(), py::arg("path_to_receiving_model_file")="undefined", py::arg("enable_user_messages")=true)
        .def(py::init<std::string&, std::string&, int, bool, bool, bool, bool, bool>(), py::arg("path_to_receiving_model_file")="undefined", py::arg("path_to_donor_model")="undefined", py::arg("nThreads")=-1, py::arg("trim_donor_when_clashes_detected")=true, py::arg("remove_donor_when_clashes_detected")=true, py::arg("ANY_search_policy")=true, py::arg("enable_user_messages")=true, py::arg("debug_output")=false)
        .def("get_path_of_receiving_model_file_used", &pm::Builder::get_path_of_receiving_model_file_used)
        .def("get_path_of_donor_model_file_used", &pm::Builder::get_path_of_donor_model_file_used)
        .def("get_receiving_model_sequence_info", &pm::Builder::get_receiving_model_sequence_info)
        .def("get_glycan_summary_from_donor", &pm::Builder::get_glycan_summary_from_donor)
        .def("get_summary_of_grafted_glycans", &pm::Builder::get_summary_of_grafted_glycans)
        .def("graft_glycan_to_receiver", &pm::Builder::graft_glycan_to_receiver)
        .def("export_grafted_model", &pm::Builder::export_grafted_model)
        .def("set_phi", &pm::Builder::setPhi)
        .def("set_psi", &pm::Builder::setPsi)
        .def("set_phi_error", &pm::Builder::setPhi_error)
        .def("set_psi_error", &pm::Builder::setPsi_error)
        .def("set_iterator_step", &pm::Builder::setIteration_step)
        ;
}

///////////////////////////////////////////////// PYBIND11 BINDING DEFINITIONS END////////////////////////////////////////////////////////////////////

