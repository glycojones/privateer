// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#include "privateer-pymodelling.h"

#define DBG std::cout << "[" << __FUNCTION__ << "] - "




privateer::pymodelling::Builder::Builder(std::string& path_to_receiving_model_file, std::string& path_to_donor_model_file, bool ANY_search_policy, bool enable_user_messages, bool debug_output) 
{
    this->ANY_search_policy = ANY_search_policy;
    this->enable_user_messages = enable_user_messages;
    this->debug_output = debug_output;
    this->read_from_file ( path_to_receiving_model_file, path_to_donor_model_file, ANY_search_policy, enable_user_messages, debug_output );
};

void privateer::pymodelling::Builder::read_from_file( std::string& path_to_receiving_model_file, std::string& path_to_donor_model_file, bool ANY_search_policy, bool enable_user_messages, bool debug_output ) {

    if(path_to_receiving_model == "undefined")
    {
        throw std::invalid_argument( "No path was provided for model file input! Aborting." );
    }

    if(path_to_donor_model == "undefined")
    {
        throw std::invalid_argument( "No path was provided for model file input! Aborting." );
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

    privateer::pyanalysis::GlycosylationComposition glycosylation(path_to_donor_model, "undefined", false);
    this->glycan_summary_donor = glycosylation.get_summary_of_detected_glycans();
    
    if(debug_output)
        DBG << "Generated a pybind11 list with summary of glycans detected in donor. Length of output list: " << glycan_summary_donor.size() << std::endl;
}

pybind11::list privateer::pymodelling::Builder::get_protein_sequence_information( clipper::MiniMol& input ) 
{
    pybind11::list output;
    for(int i = 0; i < input.size(); i++)
    {
        pybind11::list residues;
        for(int j = 0; j < input[i].size(); j++)
        {
            std::string residueType = input[i][j].type().trim();
            int residueSeqnum = input[i][j].seqnum();
            auto residue_dict = pybind11::dict("index"_a=j, "residueType"_a=residueType, "residueSeqnum"_a=residueSeqnum);
            residues.append(residue_dict);
        }
        std::string chainID = input[i].id().trim();
        auto mpolymer_dict = pybind11::dict("index"_a = i, "ChainID"_a=chainID, "Residues"_a=residues);
        output.append(mpolymer_dict);
    }
    return output;
}

void privateer::pymodelling::Builder::graft_glycan_to_receiver(int mglycanindex, int receiver_chain_index, int received_residue_index)
{
    privateer::modelling::Grafter Grafter(imported_receiving_model, imported_donor_model, enable_user_messages, debug_output);

    std::vector<clipper::MGlycan> donor_glycans = Grafter.get_donor_glycans();
    clipper::MGlycan glycan_to_graft = donor_glycans[mglycanindex];
    
    if(debug_output)
    {
        DBG << "Grafting glycan with the WURCS of: " << glycan_to_graft.generate_wurcs() << std::endl;
    }
    
    if(enable_user_messages && !debug_output)
    {
        std::cout << "Grafting glycan with the WURCS of:\n" << glycan_to_graft.generate_wurcs() << std::endl;
    }

    clipper::MPolymer converted_mglycan = Grafter.convert_mglycan_to_mpolymer(glycan_to_graft);

    if(debug_output)
        DBG << "Converted MGlycan to MPolymer. MPolymer.size(): " << converted_mglycan.size() << std::endl;

    clipper::MMonomer receiver_monomer = imported_receiving_model[receiver_chain_index][received_residue_index];
    
    if(debug_output)
    {
        DBG << "Grafting input glycan to: " << imported_receiving_model[receiver_chain_index].id().trim() << "/" << imported_receiving_model[receiver_chain_index][received_residue_index].type().trim() << "-" << imported_receiving_model[receiver_chain_index][received_residue_index].id().trim() << "\tseqnum = " << imported_receiving_model[receiver_chain_index][received_residue_index].seqnum() << std::endl;
    }
    
    if(enable_user_messages && !debug_output)
    {
        std::cout << "Grafting input glycan to: " << imported_receiving_model[receiver_chain_index].id().trim() << "/" << imported_receiving_model[receiver_chain_index][received_residue_index].type().trim() << "-" << imported_receiving_model[receiver_chain_index][received_residue_index].id().trim() << std::endl;
    }

    int receiver_atom_index = Grafter.lookup_protein_backbone_glycosylation_database(receiver_monomer.type().trim());

    clipper::String residue_name;       
    clipper::String connected_atom;          
    clipper::String vector_origin;
    clipper::String vector_target;
    clipper::ftype Phi;        
    clipper::ftype Psi; 

    if(receiver_atom_index != -1)
    {
        residue_name = privateer::modelling::backbone_instructions[receiver_atom_index].residue_name;
        connected_atom = privateer::modelling::backbone_instructions[receiver_atom_index].connected_atom;
        vector_origin = privateer::modelling::backbone_instructions[receiver_atom_index].vector_origin;
        vector_target = privateer::modelling::backbone_instructions[receiver_atom_index].vector_target;
        Phi = privateer::modelling::backbone_instructions[receiver_atom_index].Phi;
        Psi = privateer::modelling::backbone_instructions[receiver_atom_index].Psi;

        if(enable_user_messages && !debug_output)
            std::cout << "Successfully located " << residue_name << " instructions. Will connect glycan to " << connected_atom << " with " << vector_origin << " and " << vector_target << " used to generate rotation-translation matrix. Phi = " << Phi << " Psi = " << Psi << std::endl;

        if(debug_output)
            DBG << "Successfully located " << residue_name << " instructions. Will connect glycan to " << connected_atom << " with " << vector_origin << " and " << vector_target << " used to generate rotation-translation matrix. Phi = " << Phi << " Psi = " << Psi << std::endl;

    }
    else
    {
        DBG << "Unable to locate " << receiver_monomer.type().trim() << " in protein_backbone_glycosylation_instruction_set from receiving model! ...Aborting..." << std::endl;
        throw std::invalid_argument( "Unable to generate instructions for input monomer from input receiver." );
    }

    clipper::String glycan_type = "non-ideal";
    for(int atom = 0; atom < converted_mglycan[0].size(); atom++)
    {
        if(converted_mglycan[0][atom].id().trim() == "O1")
        {
            glycan_type = "ideal";
        }
    }

    clipper::String sugar_connection_atom;
    clipper::String sugar_coplanar_atom;        
    clipper::String sugar_supporting_atom;

    int glycan_grafting_type = Grafter.lookup_glycan_type_glycosylation_database(glycan_type);
    if(glycan_grafting_type != -1)
    {
        sugar_connection_atom = privateer::modelling::sugar_instructions[glycan_grafting_type].connection_atom;
        sugar_coplanar_atom = privateer::modelling::sugar_instructions[glycan_grafting_type].coplanar_atom;
        sugar_supporting_atom = privateer::modelling::sugar_instructions[glycan_grafting_type].supporting_atom;

        if(enable_user_messages && !debug_output)
            std::cout << "Successfully located " << glycan_type << " glycan grafting instructions. Will connect glycan via " << sugar_connection_atom << " with " << sugar_coplanar_atom << " and " << sugar_supporting_atom << " used to generate rotation-translation matrix." << std::endl;

        if(debug_output)
            DBG << "Successfully located " << glycan_type << " glycan grafting instructions. Will connect glycan via " << sugar_connection_atom << " with " << sugar_coplanar_atom << " and " << sugar_supporting_atom << " used to generate rotation-translation matrix." << std::endl;
    }
    else
    {
        DBG << "Unable to locate " << glycan_type << " in sugar_attachment_instruction_set from donor model! ...Aborting..." << std::endl;
        throw std::invalid_argument( "Unable to generate instructions from donor model from input receiver." );
    }

    clipper::Coord_orth protein_connecting_target;
    clipper::Coord_orth protein_vector_origin;
    clipper::Coord_orth protein_vector_target;
    clipper::Coord_orth sugar_connecting_target;
    clipper::Coord_orth sugar_coplanar_target;
    clipper::Coord_orth sugar_supporting_target;
    if(ANY_search_policy)
    {
        protein_connecting_target = receiver_monomer.find(connected_atom, clipper::MM::ANY).coord_orth(); // ND2
        protein_vector_origin = receiver_monomer.find(vector_origin, clipper::MM::ANY).coord_orth(); // CB
        protein_vector_target = receiver_monomer.find(vector_target, clipper::MM::ANY).coord_orth(); // CB
        sugar_connecting_target = converted_mglycan[0].find(sugar_connection_atom, clipper::MM::ANY).coord_orth(); // C1
        sugar_coplanar_target = converted_mglycan[0].find(sugar_coplanar_atom, clipper::MM::ANY).coord_orth(); // C1 - this used for vector
        sugar_supporting_target = converted_mglycan[0].find(sugar_supporting_atom, clipper::MM::ANY).coord_orth(); // O1 - used for vector
    }
    else
    {
        protein_connecting_target = receiver_monomer.find(connected_atom, clipper::MM::UNIQUE).coord_orth(); // ND2
        protein_vector_origin = receiver_monomer.find(vector_origin, clipper::MM::UNIQUE).coord_orth(); // CB
        protein_vector_target = receiver_monomer.find(vector_target, clipper::MM::UNIQUE).coord_orth(); // CB
        sugar_connecting_target = converted_mglycan[0].find(sugar_connection_atom, clipper::MM::UNIQUE).coord_orth(); // C1
        sugar_coplanar_target = converted_mglycan[0].find(sugar_coplanar_atom, clipper::MM::UNIQUE).coord_orth(); // C1 - this used for vector
        sugar_supporting_target = converted_mglycan[0].find(sugar_supporting_atom, clipper::MM::UNIQUE).coord_orth(); // O1 - used for vector
    }


    clipper::Coord_orth potential_C1_position = get_glycan_target_point(protein_connecting_target, protein_vector_origin, protein_vector_target, 1.5);
    clipper::Vec3<clipper::ftype> sugarVector((sugar_supporting_target.x()-sugar_coplanar_target.x()),(sugar_supporting_target.y()-sugar_coplanar_target.y()), (sugar_supporting_target.z()-sugar_coplanar_target.z()));
    clipper::Vec3<clipper::ftype> proteinTargetVector((potential_C1_position.x()-protein_connecting_target .x()),(potential_C1_position.y()-protein_connecting_target .y()), (potential_C1_position.z()-protein_connecting_target .z()));
    
    // clipper::Coord_orth source(sugarVector);
    // std::vector<clipper::Coord_orth> sourceVector;
    // sourceVector.push_back(source);

    // clipper::Coord_orth target(proteinTargetVector);
    // std::vector<clipper::Coord_orth> targetVector;
    // targetVector.push_back(target);

    clipper::Coord_orth source(sugar_connecting_target);
    std::vector<clipper::Coord_orth> sourceVector;
    sourceVector.push_back(source);

    clipper::Coord_orth target(potential_C1_position);
    std::vector<clipper::Coord_orth> targetVector;
    targetVector.push_back(target);
    
    
    clipper::RTop_orth relocator(sourceVector, targetVector);

    if(enable_user_messages && !debug_output)
    {
        std::cout << "Coordinates of " << connected_atom << ":\t" << protein_connecting_target.format() << std::endl;
        std::cout << "Coordinates of " << vector_origin << ":\t" << protein_vector_origin.format() << std::endl;
        std::cout << "Coordinates of " << vector_target << ":\t" << protein_vector_target.format() << std::endl;
        std::cout << "Coordinates of " << sugar_connection_atom << ":\t" << sugar_connecting_target.format() << std::endl;
        std::cout << "Coordinates of " << sugar_coplanar_atom << ":\t" << sugar_coplanar_target.format() << std::endl;
        std::cout << "Coordinates of " << sugar_supporting_atom << ":\t" << sugar_supporting_target.format() << std::endl;
        std::cout << "Coordinates of potential C1 position:\t" << potential_C1_position.format() << std::endl;
        std::cout << "Vector values of sugarVector:\t" << sugarVector.format() << std::endl;
        std::cout << "Vector values of proteinTargetVector:\t" << proteinTargetVector.format() << std::endl;
        std::cout << "Coordinates of source:\t" << source.format() << std::endl;
        std::cout << "Coordinates of target:\t" << target.format() << std::endl;
    }

    if(debug_output)
    {
        DBG << "Coordinates of " << connected_atom << ":\t" << protein_connecting_target.format() << std::endl;
        DBG << "Coordinates of " << vector_origin << ":\t" << protein_vector_origin.format() << std::endl;
        DBG << "Coordinates of " << vector_target << ":\t" << protein_vector_target.format() << std::endl;
        DBG << "Coordinates of " << sugar_connection_atom << ":\t" << sugar_connecting_target.format() << std::endl;
        DBG << "Coordinates of " << sugar_coplanar_atom << ":\t" << sugar_coplanar_target.format() << std::endl;
        DBG << "Coordinates of " << sugar_supporting_atom << ":\t" << sugar_supporting_target.format() << std::endl;
        DBG << "Coordinates of potential C1 position:\t" << potential_C1_position.format() << std::endl;
        DBG << "Vector values of sugarVector:\t" << sugarVector.format() << std::endl;
        DBG << "Vector values of proteinTargetVector:\t" << proteinTargetVector.format() << std::endl;
        DBG << "Coordinates of source:\t" << source.format() << std::endl;
        DBG << "Coordinates of target:\t" << target.format() << std::endl;
    }

    Grafter.graft_mpolymer_to_receiving_model(relocator, converted_mglycan);
    this->imported_receiving_model = Grafter.get_receiving_model();
}

void privateer::pymodelling::Builder::export_grafted_model(std::string& output_path)
{
    clipper::MMDBfile pdbfile;
    pdbfile.export_minimol( imported_receiving_model );
    pdbfile.write_file( output_path );
}



clipper::Coord_orth privateer::pymodelling::Builder::get_glycan_target_point(clipper::Coord_orth& connecting_atom, clipper::Coord_orth& vector_origin, clipper::Coord_orth& vector_target, float vectorShiftDistance)
{
	clipper::Coord_orth coord; 

	clipper::Vec3<clipper::ftype> baseVector((vector_target.x()-vector_origin.x()),(vector_target.y()-vector_origin.y()), (vector_target.z()-vector_origin.z()));
	// Create a 1A unit vector out of baseVector, to be used later in vector shifting
	clipper::Vec3<clipper::ftype> unitVector = baseVector.unit();

	// Obtain coordinates in the middle of suspected glycan density via 5A vector shift. This is the nearest glycan bonded via ND2 atom to ASN residue.
	coord = clipper::Coord_orth( (connecting_atom.x()+(unitVector[0]*vectorShiftDistance)), (connecting_atom.y()+(unitVector[1]*vectorShiftDistance)), (connecting_atom.z()+(unitVector[2]*vectorShiftDistance)) );

	return coord;
}




///////////////////////////////////////////////// PYBIND11 BINDING DEFINITIONS ////////////////////////////////////////////////////////////////////
namespace py = pybind11;
namespace pm = privateer::pymodelling;

void init_pymodelling(py::module& m)
{
    py::class_<pm::Builder>(m, "Builder")
        .def(py::init<>())
        .def(py::init<std::string&, std::string&, bool, bool, bool>(), py::arg("path_to_receiving_model_file")="undefined", py::arg("path_to_donor_model")="undefined", py::arg("ANY_search_policy")=true, py::arg("enable_user_messages")=true, py::arg("debug_output")=false)
        .def("get_path_of_receiving_model_file_used", &pm::Builder::get_path_of_receiving_model_file_used)
        .def("get_path_of_donor_model_file_used", &pm::Builder::get_path_of_donor_model_file_used)
        .def("get_receiving_model_sequence_info", &pm::Builder::get_receiving_model_sequence_info)
        .def("get_glycan_summary_from_donor", &pm::Builder::get_glycan_summary_from_donor)
        .def("graft_glycan_to_receiver", &pm::Builder::graft_glycan_to_receiver)
        .def("export_grafted_model", &pm::Builder::export_grafted_model);
}

///////////////////////////////////////////////// PYBIND11 BINDING DEFINITIONS END////////////////////////////////////////////////////////////////////

