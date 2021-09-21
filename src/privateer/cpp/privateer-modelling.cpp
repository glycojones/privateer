// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#include "privateer-modelling.h"

#define DBG std::cout << "[" << __FUNCTION__ << "] - "

bool compareMGlycanNode(clipper::MGlycan::Node& node_left, clipper::MGlycan::Node& node_right) {
    return node_left.get_sugar().id().trim() < node_right.get_sugar().id().trim();
}

namespace privateer
{
    namespace modelling 
    {
        const protein_sidechain_glycosylation backbone_instructions[] =
        {
            { "ASN", "ND2", "CG", "CB", -97.5, 178.0, 25.0, 25.0 },
            { "ASP", "OD2", "CG", "CB", -97.5, 178.0, 25.0, 25.0 }
        };
        const int backbone_instructions_size = sizeof( backbone_instructions ) / sizeof( backbone_instructions[0] );

        const sugar_attachment sugar_instructions[] =
        {
            { "ideal", "C1", "O1", "" },   
            { "non-ideal", "C1", "C4", "O4" }
        };
        const int sugar_instructions_size = sizeof( sugar_instructions ) / sizeof( sugar_instructions[0] );
    
        Grafter::Grafter(clipper::MiniMol receiving_model, clipper::MiniMol donor_model, int nThreads, bool trim_donor_when_clashes_detected, bool enable_user_messages, bool debug_output)
        {
            this->nThreads = nThreads;
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
            this->useParallelism = useParallelism;
            this->enable_user_messages = enable_user_messages;
            this->debug_output = debug_output;
            this->receiving_model = receiving_model;
            this->export_model = receiving_model;
            this->trim_donor_when_clashes_detected = trim_donor_when_clashes_detected;
            clipper::MGlycology donor_model_mgl = clipper::MGlycology(donor_model, debug_output, "undefined");
            std::vector<clipper::MGlycan> list_of_glycans_from_donor = donor_model_mgl.get_list_of_glycans();

            if(debug_output)
                DBG << "Number of glycans detected in donor model: " << list_of_glycans_from_donor.size() << std::endl;
            
            if(enable_user_messages && !debug_output)
                std::cout << "Number of glycans detected in donor model: " << list_of_glycans_from_donor.size() << std::endl;
                
            this->donor_glycans=list_of_glycans_from_donor;

            this->numDonorGlycansDetected=list_of_glycans_from_donor.size();
        }

        clipper::MPolymer Grafter::convert_mglycan_to_mpolymer(clipper::MGlycan input)
        {
            clipper::MPolymer output;
            int numNodes = input.number_of_nodes();
            // std::vector<clipper::MGlycan::Node> input.
            
            for(int i = 0; i < numNodes; i++)
            {
                clipper::MGlycan::Node currentNode = input.get_node(i);
                clipper::MSugar currentSugar = currentNode.get_sugar();
                clipper::MMonomer convertedMSugar = currentSugar;
                output.insert(convertedMSugar);
            }

            return output;
        }

        clipper::Coord_orth Grafter::get_glycan_target_point(clipper::Coord_orth connecting_atom, clipper::Coord_orth vector_origin, clipper::Coord_orth vector_target, float vectorShiftDistance)
        {
            clipper::Coord_orth coord; 

            clipper::Vec3<clipper::ftype> baseVector((vector_target.x()-vector_origin.x()),(vector_target.y()-vector_origin.y()), (vector_target.z()-vector_origin.z()));
            // Create a 1A unit vector out of baseVector, to be used later in vector shifting
            clipper::Vec3<clipper::ftype> unitVector = baseVector.unit();

            // Obtain coordinates in the middle of suspected glycan density via 5A vector shift. This is the nearest glycan bonded via ND2 atom to ASN residue.
            coord = clipper::Coord_orth( (connecting_atom.x()+(unitVector[0]*vectorShiftDistance)), (connecting_atom.y()+(unitVector[1]*vectorShiftDistance)), (connecting_atom.z()+(unitVector[2]*vectorShiftDistance)) );

            return coord;
        }

        // Need to come up with a method that would put the dummy O1 atom in correct position for both types of anomers. Currently this does not work.
        clipper::Coord_orth Grafter::get_dummy_O1_position(clipper::MAtom& sugar_connection_target, clipper::MAtom& sugar_vector_point_alpha_target, clipper::MAtom& sugar_vector_point_bravo_target, clipper::String residue_name)
        {
            clipper::Coord_orth coord; 
            
            clipper::Vec3<clipper::ftype> baseVector((sugar_vector_point_bravo_target.coord_orth().x()-sugar_vector_point_alpha_target.coord_orth().x()),(sugar_vector_point_bravo_target.coord_orth().y()-sugar_vector_point_alpha_target.coord_orth().y()), (sugar_vector_point_bravo_target.coord_orth().z()-sugar_vector_point_alpha_target.coord_orth().z()));

            coord = clipper::Coord_orth(sugar_connection_target.coord_orth().x() + baseVector[0], sugar_connection_target.coord_orth().y() + baseVector[1], sugar_connection_target.coord_orth().z() + baseVector[2]);

            return coord;
        }

        void Grafter::overlay_mglycan_via_atom(clipper::Coord_orth target, clipper::Coord_orth origin, clipper::MPolymer& converted_mglycan)
        {
            clipper::Mat33<clipper::ftype> identity_matrix;
            identity_matrix = identity_matrix.identity();

            clipper::Vec3<clipper::ftype> origin_to_target_vector(  target.x() - origin.x(), 
                                                                    target.y() - origin.y(), 
                                                                    target.z() - origin.z());
            

            clipper::RTop_orth target_origin_overlayer(identity_matrix, origin_to_target_vector);
            converted_mglycan.transform(target_origin_overlayer);
        }

        void Grafter::graft_mpolymer_to_receiving_model(clipper::MGlycan& glycan_to_graft, clipper::MMonomer& input_protein_side_chain_residue, clipper::String root_chain_id, bool ANY_search_policy)
        { 

            int receiver_atom_index = lookup_protein_backbone_glycosylation_database(input_protein_side_chain_residue.type().trim());

            clipper::String residue_name;       
            clipper::String connected_atom;          
            clipper::String vector_point_alpha;
            clipper::String vector_point_bravo;
            clipper::ftype targetPhi;        
            clipper::ftype targetPsi;
            clipper::ftype Phi_error;        
            clipper::ftype Psi_error; 

            if(receiver_atom_index != -1)
            {
                residue_name = privateer::modelling::backbone_instructions[receiver_atom_index].residue_name;
                connected_atom = privateer::modelling::backbone_instructions[receiver_atom_index].connected_atom;
                vector_point_alpha = privateer::modelling::backbone_instructions[receiver_atom_index].vector_point_alpha;
                vector_point_bravo = privateer::modelling::backbone_instructions[receiver_atom_index].vector_point_bravo;
                targetPhi = privateer::modelling::backbone_instructions[receiver_atom_index].Phi;
                targetPsi = privateer::modelling::backbone_instructions[receiver_atom_index].Psi;
                Phi_error = privateer::modelling::backbone_instructions[receiver_atom_index].Phi_error;
                Psi_error = privateer::modelling::backbone_instructions[receiver_atom_index].Psi_error;

                if(enable_user_messages && !debug_output)
                    std::cout << "Successfully located " << residue_name << " instructions. Will connect glycan to " << connected_atom << " with " << vector_point_alpha << " and " << vector_point_bravo << " used to generate rotation-translation matrix. Target Phi = " << targetPhi << ", target Psi = " << targetPsi << std::endl;

                if(debug_output)
                    DBG << "Successfully located " << residue_name << " instructions. Will connect glycan to " << connected_atom << " with " << vector_point_alpha << " and " << vector_point_bravo << " used to generate rotation-translation matrix. Target Phi = " << targetPhi << " target Psi = " << targetPsi << std::endl;

            }
            else
            {
                DBG << "Unable to locate " << input_protein_side_chain_residue.type().trim() << " in protein_backbone_glycosylation_instruction_set from receiving model! ...Aborting..." << std::endl;
                throw std::invalid_argument( "Unable to generate instructions for input monomer from input receiver." );
            }

            const clipper::String chainids = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
            std::vector<clipper::String> used_chain_ids;
            for(int i = 0; i < export_model.size(); i++)
            {
                clipper::String currentChainID = export_model[i].id().trim();
                used_chain_ids.push_back(currentChainID);
            }

            char new_chain_id = chainids[used_chain_ids.size()];
            std::string tempchainid(1, new_chain_id);
            clipper::String chainID(tempchainid);
            this->grafted_glycan_chainID = chainID;

            clipper::MPolymer converted_mglycan = convert_mglycan_to_mpolymer(glycan_to_graft);

            if(debug_output)
                DBG << "Converted MGlycan to MPolymer. MPolymer.size(): " << converted_mglycan.size() << std::endl;
 
            clipper::String glycan_type = "non-ideal";
            for(int atom = 0; atom < converted_mglycan[0].size(); atom++)
            {
                if(converted_mglycan[0][atom].id().trim() == "O1")
                {
                    glycan_type = "ideal";
                }
            }

            if(glycan_type == "non-ideal")
            {
                std::cout << "!!! THIS METHOD OF GRAFTING IS CURRENTLY UNOPTIMIZED. THIS WILL BE CORRECTED IN DUE FUTURE. PLEASE USE THE SUPPLIED GLYCAN BLOCKS AS THEY CONTAIN O1 ATOMS !!!" << std::endl;
                std::cout << "!!! THIS METHOD OF GRAFTING IS CURRENTLY UNOPTIMIZED. THIS WILL BE CORRECTED IN DUE FUTURE. PLEASE USE THE SUPPLIED GLYCAN BLOCKS AS THEY CONTAIN O1 ATOMS !!!" << std::endl;
                std::cout << "!!! THIS METHOD OF GRAFTING IS CURRENTLY UNOPTIMIZED. THIS WILL BE CORRECTED IN DUE FUTURE. PLEASE USE THE SUPPLIED GLYCAN BLOCKS AS THEY CONTAIN O1 ATOMS !!!" << std::endl;
                std::cout << "!!! THIS METHOD OF GRAFTING IS CURRENTLY UNOPTIMIZED. THIS WILL BE CORRECTED IN DUE FUTURE. PLEASE USE THE SUPPLIED GLYCAN BLOCKS AS THEY CONTAIN O1 ATOMS !!!" << std::endl;
                std::cout << "!!! THIS METHOD OF GRAFTING IS CURRENTLY UNOPTIMIZED. THIS WILL BE CORRECTED IN DUE FUTURE. PLEASE USE THE SUPPLIED GLYCAN BLOCKS AS THEY CONTAIN O1 ATOMS !!!" << std::endl;
                std::cout << "!!! THIS METHOD OF GRAFTING IS CURRENTLY UNOPTIMIZED. THIS WILL BE CORRECTED IN DUE FUTURE. PLEASE USE THE SUPPLIED GLYCAN BLOCKS AS THEY CONTAIN O1 ATOMS !!!" << std::endl;
                converted_mglycan.set_id(chainID);
                clipper::MM::MODE search_policy;
                if(ANY_search_policy)
                    search_policy = clipper::MM::MODE::ANY;
                else
                    search_policy = clipper::MM::MODE::UNIQUE;

                clipper::String sugar_connection_atom;

                int glycan_grafting_type = lookup_glycan_type_glycosylation_database(glycan_type);
                if(glycan_grafting_type != -1)
                {
                    sugar_connection_atom = privateer::modelling::sugar_instructions[glycan_grafting_type].connection_atom;

                    if(enable_user_messages && !debug_output)
                        std::cout << "Successfully located " << glycan_type << " glycan grafting instructions. Will connect glycan via " << sugar_connection_atom << " used to generate rotation-translation matrix." << std::endl;

                    if(debug_output)
                        DBG << "Successfully located " << glycan_type << " glycan grafting instructions. Will connect glycan using " << sugar_connection_atom << " used to translate in vicinity of " << connected_atom << " atom." << std::endl;
                }
                else
                {
                    DBG << "Unable to locate " << glycan_type << " in sugar_attachment_instruction_set from donor model! ...Aborting..." << std::endl;
                    throw std::invalid_argument( "Unable to generate instructions from donor model from input receiver." );
                }

                clipper::MAtom protein_connecting_target;
                clipper::MAtom protein_vector_point_alpha;
                clipper::MAtom protein_vector_point_bravo;
                clipper::MAtom sugar_connection_target;

                protein_connecting_target = input_protein_side_chain_residue.find(connected_atom, search_policy); // ND2
                protein_vector_point_alpha = input_protein_side_chain_residue.find(vector_point_alpha, search_policy); // CB
                protein_vector_point_bravo = input_protein_side_chain_residue.find(vector_point_bravo, search_policy); // CG
                sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, search_policy); // C1

                clipper::Coord_orth potential_C1_position = get_glycan_target_point(protein_connecting_target.coord_orth(), protein_vector_point_bravo.coord_orth(), protein_vector_point_alpha.coord_orth(), 1.50);
                
                overlay_mglycan_via_atom(potential_C1_position, sugar_connection_target.coord_orth(), converted_mglycan);

                // Update the coordinates of C1 and O1, get O5 coords from the glycan.
                clipper::String ring_oxygen_name = glycan_to_graft.get_sugars()[0].ring_members()[0].id().trim();
                clipper::MAtom ring_oxygen;

                sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, search_policy); // C1
                ring_oxygen = converted_mglycan[0].find(ring_oxygen_name, search_policy); // O5

                double currentPsiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(sugar_connection_target.coord_orth(), protein_connecting_target.coord_orth(), protein_vector_point_alpha.coord_orth(), protein_vector_point_bravo.coord_orth()));
                std::vector<std::pair<clipper::MAtom, std::string>> psiTorsionAtoms = { std::make_pair(sugar_connection_target, "sugar"), std::make_pair(protein_connecting_target, "protein"), std::make_pair(protein_vector_point_alpha, "protein"), std::make_pair(protein_vector_point_bravo, "protein") };
                
                clipper::Coord_orth psiDirection = protein_vector_point_alpha.coord_orth() - protein_connecting_target.coord_orth(); // CG(origin_shift) - ND2(base)
                rotate_mglycan_until_torsion_angle_fulfilled(converted_mglycan, input_protein_side_chain_residue, psiDirection, protein_vector_point_alpha.coord_orth(), psiTorsionAtoms, targetPsi, debug_output);

                sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, search_policy); // C1
                ring_oxygen = converted_mglycan[0].find(ring_oxygen_name, search_policy); // O5

                currentPsiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(sugar_connection_target.coord_orth(), protein_connecting_target.coord_orth(), protein_vector_point_alpha.coord_orth(), protein_vector_point_bravo.coord_orth()));
                
                if(debug_output)
                    DBG << "Psi value after rotation: " << currentPsiTorsionAngle << std::endl;

                double currentPhiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(ring_oxygen.coord_orth(), sugar_connection_target.coord_orth(), protein_connecting_target.coord_orth(), protein_vector_point_alpha.coord_orth()));
                std::vector<std::pair<clipper::MAtom, std::string>> phiTorsionAtoms = { std::make_pair(ring_oxygen, "sugar"), std::make_pair(sugar_connection_target, "sugar"), std::make_pair(protein_connecting_target, "protein"), std::make_pair(protein_vector_point_alpha, "protein") };

                clipper::Coord_orth phiDirection = protein_connecting_target.coord_orth() - sugar_connection_target.coord_orth(); // ND2(origin_shift) - C1(base)
                rotate_mglycan_until_torsion_angle_fulfilled(converted_mglycan, input_protein_side_chain_residue, phiDirection, protein_connecting_target.coord_orth(), phiTorsionAtoms, targetPhi, debug_output);

                sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, search_policy); // C1
                ring_oxygen = converted_mglycan[0].find(ring_oxygen_name, search_policy); // O5

                currentPhiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(ring_oxygen.coord_orth(), sugar_connection_target.coord_orth(), protein_connecting_target.coord_orth(), protein_vector_point_alpha.coord_orth()));
                
                if(debug_output)
                    DBG << "Phi value after rotation: " << currentPhiTorsionAngle << std::endl;

                clipper::MiniMol tmp_clash_model = export_model;
                tmp_clash_model.insert(converted_mglycan);

                std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > clashes_after_initial_grafting = check_for_clashes(tmp_clash_model, converted_mglycan[0], input_protein_side_chain_residue, root_chain_id, chainID);
                this->clashes = clashes_after_initial_grafting;

                if(debug_output)
                    DBG << "clashes_after_initial_grafting.size(): " << clashes_after_initial_grafting.size() << std::endl;

                if(clashes_after_initial_grafting.size() > 0)
                {
                    tmp_clash_model = export_model;
                    sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, search_policy); // C1
                    ring_oxygen = converted_mglycan[0].find(ring_oxygen_name, search_policy); // O5
                    
                    if(enable_user_messages)
                    {
                        std::cout << "Unfortunately, the grafted glycan resulted in " << clashes_after_initial_grafting.size() << " clashes with protein backbone. Attempting to eliminate or minimize clashes through further manipulation of Phi/Psi torsion angles." << std::endl;
                        for(int i = 0; i < clashes_after_initial_grafting.size(); i++)
                        {
                            std::cout << i+1 << "/" << clashes_after_initial_grafting.size() << " initial clash: " << clashes_after_initial_grafting[i].first.second << "/" << clashes_after_initial_grafting[i].first.first.id().trim() << "-" << clashes_after_initial_grafting[i].first.first.type().trim() << "\t\t\tSugar: " << clashes_after_initial_grafting[i].second.second << "/" << clashes_after_initial_grafting[i].second.first.id().trim() << "-" << clashes_after_initial_grafting[i].second.first.type().trim() << std::endl;
                        }
                    }
                    
                    clipper::MPolymer minimized_clashes_glycan;
                    if(useParallelism)
                        minimized_clashes_glycan = rotate_mglycan_until_clashes_are_minimized_parallelized(tmp_clash_model, converted_mglycan, input_protein_side_chain_residue, phiTorsionAtoms, psiTorsionAtoms, Phi_error, Psi_error, root_chain_id, chainID, debug_output);
                    else
                        minimized_clashes_glycan = rotate_mglycan_until_clashes_are_minimized_singlethreaded(tmp_clash_model, converted_mglycan, input_protein_side_chain_residue, phiTorsionAtoms, psiTorsionAtoms, Phi_error, Psi_error, root_chain_id, chainID, debug_output);
                    
                    sugar_connection_target = minimized_clashes_glycan[0].find(sugar_connection_atom, search_policy); // C1
                    ring_oxygen = minimized_clashes_glycan[0].find(ring_oxygen_name, search_policy); // O5
                    
                    tmp_clash_model = export_model;
                    tmp_clash_model.insert(minimized_clashes_glycan);

                    
                    std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > clashes_after_manipulation = check_for_clashes(tmp_clash_model, converted_mglycan[0], input_protein_side_chain_residue, root_chain_id, chainID);
                
                    if(enable_user_messages)
                    {
                        currentPhiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(ring_oxygen.coord_orth(), sugar_connection_target.coord_orth(), protein_connecting_target.coord_orth(), protein_vector_point_alpha.coord_orth()));
                        currentPsiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(sugar_connection_target.coord_orth(), protein_connecting_target.coord_orth(), protein_vector_point_alpha.coord_orth(), protein_vector_point_bravo.coord_orth()));
                        if(enable_user_messages)
                        {
                            std::cout << "After manipulating Phi/Psi torsion angles further, the following were found to minimize clashes the most, Psi: " << currentPsiTorsionAngle << "\t\tPhi: " << currentPhiTorsionAngle << std::endl; 
                            std::cout << "Managed to eliminate " << clashes_after_initial_grafting.size() - clashes_after_manipulation.size() << " clashes. Remaining clashes: " << std::endl;
                        }
                        for(int i = 0; i < clashes_after_manipulation.size(); i++)
                        {
                            std::cout << i+1 << "/" << clashes_after_manipulation.size() << " remaining clash: " << clashes_after_manipulation[i].first.second << "/" << clashes_after_manipulation[i].first.first.id().trim() << "-" << clashes_after_manipulation[i].first.first.type().trim() << "\t\t\tSugar: " << clashes_after_manipulation[i].second.second << "/" << clashes_after_manipulation[i].second.first.id().trim() << "-" << clashes_after_manipulation[i].second.first.type().trim() << std::endl;
                        }
                    }
                    this->clashes = clashes_after_manipulation;
                    converted_mglycan = minimized_clashes_glycan;
                }

                if(enable_user_messages && !debug_output)
                    std::cout << "Grafting input glycan with Chain ID of: " << chainID << " to " << residue_name << "-" << input_protein_side_chain_residue.id().trim() << std::endl;
                
                if(debug_output)
                    DBG << "Grafting input glycan with Chain ID of: " << chainID << " to " << residue_name << "-" << input_protein_side_chain_residue.id().trim() << std::endl;
                
                this->grafted_glycan = glycan_to_graft;
                export_model.insert(converted_mglycan);

                if(trim_donor_when_clashes_detected)
                {
                    std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > clashes_after_grafting = check_for_clashes(export_model, converted_mglycan[0], input_protein_side_chain_residue, root_chain_id, chainID);
                    clipper::MiniMol clash_free = trim_graft_until_no_clashes_left(export_model, glycan_to_graft, std::make_pair(input_protein_side_chain_residue, root_chain_id), chainID, clashes_after_grafting);
                    export_model = clash_free;
                }
                
                if(enable_user_messages && !debug_output)
                    std::cout << "Glycan has been grafted!" << std::endl;
                
                if(debug_output)
                    DBG << "Glycan has been grafted!" << std::endl;
            }
            else if(glycan_type == "ideal")
            {
                converted_mglycan.set_id(chainID);
                clipper::MM::MODE search_policy;
                if(ANY_search_policy)
                    search_policy = clipper::MM::MODE::ANY;
                else
                    search_policy = clipper::MM::MODE::UNIQUE;

                clipper::String sugar_connection_atom;
                clipper::String sugar_vector_point;

                int glycan_grafting_type = lookup_glycan_type_glycosylation_database(glycan_type);
                if(glycan_grafting_type != -1)
                {
                    sugar_connection_atom = privateer::modelling::sugar_instructions[glycan_grafting_type].connection_atom;
                    sugar_vector_point = privateer::modelling::sugar_instructions[glycan_grafting_type].vector_point_alpha;

                    if(enable_user_messages && !debug_output)
                        std::cout << "Successfully located " << glycan_type << " glycan grafting instructions. Will connect glycan via " << sugar_connection_atom << " used to generate rotation-translation matrix." << std::endl;

                    if(debug_output)
                        DBG << "Successfully located " << glycan_type << " glycan grafting instructions. Will connect glycan using " << sugar_vector_point << " used to translate onto " << connected_atom << " atom." << std::endl;
                }
                else
                {
                    DBG << "Unable to locate " << glycan_type << " in sugar_attachment_instruction_set from donor model! ...Aborting..." << std::endl;
                    throw std::invalid_argument( "Unable to generate instructions from donor model from input receiver." );
                }

                clipper::MAtom protein_connecting_target;
                clipper::MAtom protein_vector_point_alpha;
                clipper::MAtom protein_vector_point_bravo;
                clipper::MAtom sugar_connection_target;
                clipper::MAtom sugar_vector_point_target;

                protein_connecting_target = input_protein_side_chain_residue.find(connected_atom, search_policy); // ND2
                protein_vector_point_alpha = input_protein_side_chain_residue.find(vector_point_alpha, search_policy); // CB
                protein_vector_point_bravo = input_protein_side_chain_residue.find(vector_point_bravo, search_policy); // CG
                sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, search_policy); // C1
                sugar_vector_point_target = converted_mglycan[0].find(sugar_vector_point, search_policy); // O1


                overlay_mglycan_via_atom(protein_connecting_target.coord_orth(), sugar_vector_point_target.coord_orth(), converted_mglycan);

                // Update the coordinates of C1 and O1, get O5 coords from the glycan.
                clipper::String ring_oxygen_name = glycan_to_graft.get_sugars()[0].ring_members()[0].id().trim();
                clipper::MAtom ring_oxygen;

                sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, search_policy); // C1
                sugar_vector_point_target = converted_mglycan[0].find(sugar_vector_point, search_policy); // O1
                ring_oxygen = converted_mglycan[0].find(ring_oxygen_name, search_policy); // O5

                double currentPsiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(sugar_connection_target.coord_orth(), protein_connecting_target.coord_orth(), protein_vector_point_alpha.coord_orth(), protein_vector_point_bravo.coord_orth()));
                std::vector<std::pair<clipper::MAtom, std::string>> psiTorsionAtoms = { std::make_pair(sugar_connection_target, "sugar"), std::make_pair(protein_connecting_target, "protein"), std::make_pair(protein_vector_point_alpha, "protein"), std::make_pair(protein_vector_point_bravo, "protein") };
                
                clipper::Coord_orth psiDirection = protein_vector_point_alpha.coord_orth() - protein_connecting_target.coord_orth(); // CG(origin_shift) - ND2(base)
                rotate_mglycan_until_torsion_angle_fulfilled(converted_mglycan, input_protein_side_chain_residue, psiDirection, protein_vector_point_alpha.coord_orth(), psiTorsionAtoms, targetPsi, debug_output);

                sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, search_policy); // C1
                sugar_vector_point_target = converted_mglycan[0].find(sugar_vector_point, search_policy); // O1
                ring_oxygen = converted_mglycan[0].find(ring_oxygen_name, search_policy); // O5

                currentPsiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(sugar_connection_target.coord_orth(), protein_connecting_target.coord_orth(), protein_vector_point_alpha.coord_orth(), protein_vector_point_bravo.coord_orth()));
                
                if(debug_output)
                    DBG << "Psi value after rotation: " << currentPsiTorsionAngle << std::endl;

                double currentPhiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(ring_oxygen.coord_orth(), sugar_connection_target.coord_orth(), protein_connecting_target.coord_orth(), protein_vector_point_alpha.coord_orth()));
                std::vector<std::pair<clipper::MAtom, std::string>> phiTorsionAtoms = { std::make_pair(ring_oxygen, "sugar"), std::make_pair(sugar_connection_target, "sugar"), std::make_pair(protein_connecting_target, "protein"), std::make_pair(protein_vector_point_alpha, "protein") };

                clipper::Coord_orth phiDirection = protein_connecting_target.coord_orth() - sugar_connection_target.coord_orth(); // ND2(origin_shift) - C1(base)
                rotate_mglycan_until_torsion_angle_fulfilled(converted_mglycan, input_protein_side_chain_residue, phiDirection, protein_connecting_target.coord_orth(), phiTorsionAtoms, targetPhi, debug_output);

                sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, search_policy); // C1
                sugar_vector_point_target = converted_mglycan[0].find(sugar_vector_point, search_policy); // O1
                ring_oxygen = converted_mglycan[0].find(ring_oxygen_name, search_policy); // O5

                currentPhiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(ring_oxygen.coord_orth(), sugar_connection_target.coord_orth(), protein_connecting_target.coord_orth(), protein_vector_point_alpha.coord_orth()));
                
                if(debug_output)
                    DBG << "Phi value after rotation: " << currentPhiTorsionAngle << std::endl;

                clipper::MiniMol tmp_clash_model = export_model;
                tmp_clash_model.insert(converted_mglycan);

                std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > clashes_after_initial_grafting = check_for_clashes(tmp_clash_model, converted_mglycan[0], input_protein_side_chain_residue, root_chain_id, chainID);
                this->clashes = clashes_after_initial_grafting;

                if(debug_output)
                    DBG << "clashes_after_initial_grafting.size(): " << clashes_after_initial_grafting.size() << std::endl;

                if(clashes_after_initial_grafting.size() > 0)
                {
                    tmp_clash_model = export_model;
                    sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, search_policy); // C1
                    sugar_vector_point_target = converted_mglycan[0].find(sugar_vector_point, search_policy); // O1
                    ring_oxygen = converted_mglycan[0].find(ring_oxygen_name, search_policy); // O5
                    
                    if(enable_user_messages)
                    {
                        std::cout << "Unfortunately, the grafted glycan resulted in " << clashes_after_initial_grafting.size() << " clashes with protein backbone. Attempting to eliminate or minimize clashes through further manipulation of Phi/Psi torsion angles." << std::endl;
                        for(int i = 0; i < clashes_after_initial_grafting.size(); i++)
                        {
                            std::cout << i+1 << "/" << clashes_after_initial_grafting.size() << " initial clash: " << clashes_after_initial_grafting[i].first.second << "/" << clashes_after_initial_grafting[i].first.first.id().trim() << "-" << clashes_after_initial_grafting[i].first.first.type().trim() << "\t\t\tSugar: " << clashes_after_initial_grafting[i].second.second << "/" << clashes_after_initial_grafting[i].second.first.id().trim() << "-" << clashes_after_initial_grafting[i].second.first.type().trim() << std::endl;
                        }
                    }

                    clipper::MPolymer minimized_clashes_glycan;
                    if(useParallelism)
                        minimized_clashes_glycan = rotate_mglycan_until_clashes_are_minimized_parallelized(tmp_clash_model, converted_mglycan, input_protein_side_chain_residue, phiTorsionAtoms, psiTorsionAtoms, Phi_error, Psi_error, root_chain_id, chainID, debug_output);
                    else
                        minimized_clashes_glycan = rotate_mglycan_until_clashes_are_minimized_singlethreaded(tmp_clash_model, converted_mglycan, input_protein_side_chain_residue, phiTorsionAtoms, psiTorsionAtoms, Phi_error, Psi_error, root_chain_id, chainID, debug_output);
                    
                    sugar_connection_target = minimized_clashes_glycan[0].find(sugar_connection_atom, search_policy); // C1
                    sugar_vector_point_target = minimized_clashes_glycan[0].find(sugar_vector_point, search_policy); // O1
                    ring_oxygen = minimized_clashes_glycan[0].find(ring_oxygen_name, search_policy); // O5
                    
                    tmp_clash_model = export_model;
                    tmp_clash_model.insert(minimized_clashes_glycan);

                    
                    std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > clashes_after_manipulation = check_for_clashes(tmp_clash_model, converted_mglycan[0], input_protein_side_chain_residue, root_chain_id, chainID);
                
                    if(enable_user_messages)
                    {
                        currentPhiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(ring_oxygen.coord_orth(), sugar_connection_target.coord_orth(), protein_connecting_target.coord_orth(), protein_vector_point_alpha.coord_orth()));
                        currentPsiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(sugar_connection_target.coord_orth(), protein_connecting_target.coord_orth(), protein_vector_point_alpha.coord_orth(), protein_vector_point_bravo.coord_orth()));
                        if(enable_user_messages)
                        {
                            std::cout << "After manipulating Phi/Psi torsion angles further, the following were found to minimize clashes the most, Psi: " << currentPsiTorsionAngle << "\t\tPhi: " << currentPhiTorsionAngle << std::endl; 
                            std::cout << "Managed to eliminate " << clashes_after_initial_grafting.size() - clashes_after_manipulation.size() << " clashes. Remaining clashes: " << std::endl;
                        }
                        for(int i = 0; i < clashes_after_manipulation.size(); i++)
                        {
                            std::cout << i+1 << "/" << clashes_after_manipulation.size() << " remaining clash: " << clashes_after_manipulation[i].first.second << "/" << clashes_after_manipulation[i].first.first.id().trim() << "-" << clashes_after_manipulation[i].first.first.type().trim() << "\t\t\tSugar: " << clashes_after_manipulation[i].second.second << "/" << clashes_after_manipulation[i].second.first.id().trim() << "-" << clashes_after_manipulation[i].second.first.type().trim() << std::endl;
                        }
                    }
                    this->clashes = clashes_after_manipulation;
                    converted_mglycan = minimized_clashes_glycan;
                }

                bool first_sugar_has_hydrogens = check_if_residue_has_hydrogens(converted_mglycan[0]);

                if(first_sugar_has_hydrogens)
                {
                    converted_mglycan = delete_atom_from_mglycan(converted_mglycan, sugar_vector_point_target);

                    clipper::String vector_point_target_hydrogen_name = "H" + sugar_vector_point.substr(1,1);

                    if(debug_output)
                        DBG << "Attempting to delete the following H atom: '" << vector_point_target_hydrogen_name << "'" << std::endl;

                    clipper::MAtom vector_point_target_hydrogen = converted_mglycan[0].find(vector_point_target_hydrogen_name, search_policy);

                    converted_mglycan = delete_atom_from_mglycan(converted_mglycan, vector_point_target_hydrogen);

                    if(enable_user_messages && !debug_output)
                        std::cout << "Grafting input glycan with Chain ID of: " << chainID << " to " << residue_name << "-" << input_protein_side_chain_residue.id().trim() << std::endl;
                    
                    if(debug_output)
                        DBG << "Grafting input glycan with Chain ID of: " << chainID << " to " << residue_name << "-" << input_protein_side_chain_residue.id().trim() << std::endl;
                    
                    this->grafted_glycan = glycan_to_graft;
                    export_model.insert(converted_mglycan);

                    if(trim_donor_when_clashes_detected)
                    {
                        std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > clashes_after_grafting = check_for_clashes(export_model, converted_mglycan[0], input_protein_side_chain_residue, root_chain_id, chainID);
                        clipper::MiniMol clash_free = trim_graft_until_no_clashes_left(export_model, glycan_to_graft, std::make_pair(input_protein_side_chain_residue, root_chain_id), chainID, clashes_after_grafting);
                        export_model = clash_free;
                    }
                        
                    if(enable_user_messages && !debug_output)
                        std::cout << "Glycan has been grafted!" << std::endl;
                    
                    if(debug_output)
                        DBG << "Glycan has been grafted!" << std::endl;

                }
                else
                {
                    clipper::MPolymer converted_mglycan = delete_atom_from_mglycan(converted_mglycan, sugar_vector_point_target);

                    if(enable_user_messages && !debug_output)
                        std::cout << "Grafting input glycan with Chain ID of: " << chainID << " to " << residue_name << "-" << input_protein_side_chain_residue.id().trim() << std::endl;
                    
                    if(debug_output)
                        DBG << "Grafting input glycan with Chain ID of: " << chainID << " to " << residue_name << "-" << input_protein_side_chain_residue.id().trim() << std::endl;
                    
                    this->grafted_glycan = glycan_to_graft;
                    export_model.insert(converted_mglycan);

                    if(trim_donor_when_clashes_detected)
                    {
                        std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > clashes_after_grafting = check_for_clashes(export_model, converted_mglycan[0], input_protein_side_chain_residue, root_chain_id, chainID);
                        clipper::MiniMol clash_free = trim_graft_until_no_clashes_left(export_model, glycan_to_graft, std::make_pair(input_protein_side_chain_residue, root_chain_id), chainID, clashes_after_grafting);
                        export_model = clash_free;
                    }
                    
                    if(enable_user_messages && !debug_output)
                        std::cout << "Glycan has been grafted!" << std::endl;
                    
                    if(debug_output)
                        DBG << "Glycan has been grafted!" << std::endl;
                }
            }  
        }

        // Function adopted from Paul Emsley's Coot software. 
        clipper::Coord_orth Grafter::generate_rotation_matrix_from_rodrigues_rotation_formula(clipper::Coord_orth direction, clipper::Coord_orth position, clipper::Coord_orth origin_shift, double angle)
        {
            clipper::Coord_orth unit_vec = clipper::Coord_orth(direction.unit());
            
            double l = unit_vec[0];
            double m = unit_vec[1];
            double n = unit_vec[2];

            double ll = l*l;
            double mm = m*m;
            double nn = n*n;
            double cosk = cos(angle);
            double sink = sin(angle);
            double I_cosk = 1.0 - cosk;
            
            // The Rotation matrix angle w about vector with direction cosines l,m,n.
            // 
            // ( l**2+(m**2+n**2)cos k     lm(1-cos k)-nsin k        nl(1-cos k)+msin k   )
            // ( lm(1-cos k)+nsin k        m**2+(l**2+n**2)cos k     mn(1-cos k)-lsin k   )
            // ( nl(1-cos k)-msin k        mn(1-cos k)+lsin k        n*2+(l**2+m**2)cos k )
            //
            // (Amore documentation) Thanks for that pointer EJD :).
            
            clipper::Mat33<double> r( ll+(mm+nn)*cosk,    l*m*I_cosk-n*sink,  n*l*I_cosk+m*sink,
                            l*m*I_cosk+n*sink,  mm+(ll+nn)*cosk,    m*n*I_cosk-l*sink,
                            n*l*I_cosk-m*sink,  m*n*I_cosk+l*sink,  nn+(ll+mm)*cosk );
            
            clipper::RTop_orth rtop(r, clipper::Coord_orth(0,0,0));
            return origin_shift + (position-origin_shift).transform(rtop);
        }


        void Grafter::rotate_mglycan_until_torsion_angle_fulfilled(clipper::MPolymer& converted_mglycan, clipper::MMonomer& protein_residue, clipper::Coord_orth direction, clipper::Coord_orth origin_shift, std::vector<std::pair<clipper::MAtom, std::string>>& torsionAtoms, double targetAngle, bool debug_output)
        {
            clipper::MAtom firstTorsionAtom; // always sugar
            clipper::MAtom secondTorsionAtom; 
            clipper::MAtom thirdTorsionAtom; // always protein
            clipper::MAtom fourthTorsionAtom; // always protein

            firstTorsionAtom = converted_mglycan[0].find(torsionAtoms[0].first.id().trim(), clipper::MM::UNIQUE);
            
            if(torsionAtoms[1].second == "sugar")
            {
                clipper::String secondTorsionAtomName = torsionAtoms[1].first.id().trim();
                secondTorsionAtom = converted_mglycan[0].find(secondTorsionAtomName, clipper::MM::UNIQUE);
            }
            else
            {
                clipper::String secondTorsionAtomName = torsionAtoms[1].first.id().trim();
                secondTorsionAtom = protein_residue.find(secondTorsionAtomName, clipper::MM::UNIQUE);
            }
            
            thirdTorsionAtom = protein_residue.find(torsionAtoms[2].first.id().trim(), clipper::MM::UNIQUE);
            fourthTorsionAtom = protein_residue.find(torsionAtoms[3].first.id().trim(), clipper::MM::UNIQUE);

            double currentTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(firstTorsionAtom.coord_orth(), secondTorsionAtom.coord_orth(), thirdTorsionAtom.coord_orth(), fourthTorsionAtom.coord_orth()));
            double rotangle_radian = clipper::Util::d2rad(currentTorsionAngle - targetAngle);

            if(debug_output)
                DBG << "targetTorsion: " << targetAngle << "\t\tCurrent value of torsion: " << currentTorsionAngle << "\trotangle " << clipper::Util::rad2d(rotangle_radian) << "\t\t\t" << firstTorsionAtom.id().trim() << "-" << secondTorsionAtom.id().trim() << "-" << thirdTorsionAtom.id().trim() << "-" << fourthTorsionAtom.id().trim() << std::endl;

            for(int residue = 0; residue < converted_mglycan.size(); residue++)
            {
                clipper::MMonomer currentResidue = converted_mglycan[residue];
                for(int atom = 0; atom < converted_mglycan[residue].size(); atom++)
                {
                    clipper::MAtom currentAtom = converted_mglycan[residue][atom];
                    clipper::Coord_orth old_pos = currentAtom.coord_orth();
                    clipper::Coord_orth new_pos = generate_rotation_matrix_from_rodrigues_rotation_formula(direction, old_pos, origin_shift, rotangle_radian);
                    converted_mglycan[residue][atom].set_coord_orth(new_pos);
                }
            }
        }

        clipper::MPolymer Grafter::rotate_mglycan_until_clashes_are_minimized_parallelized(clipper::MiniMol& export_model, clipper::MPolymer& converted_mglycan, clipper::MMonomer& protein_residue, std::vector<std::pair<clipper::MAtom, std::string>>& phiTorsionAtoms, std::vector<std::pair<clipper::MAtom, std::string>>& psiTorsionAtoms, double phiError, double psiError, clipper::String root_chain_id, clipper::String root_sugar_chain_id, bool debug_output)
        {
            double iteration_step = 1.0;
            clipper::MiniMol single_clash_model = export_model;
            clipper::MPolymer best_performing_glycan = converted_mglycan;
            double bestPerformingPsi, bestPerformingPhi;

            single_clash_model.insert(converted_mglycan);
            std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > reference_clashes = check_for_clashes(single_clash_model, converted_mglycan[0], protein_residue, root_chain_id, root_sugar_chain_id);
            // std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > current_clashes;
            int n_reference_clashes = reference_clashes.size();

            single_clash_model = export_model;

            clipper::MAtom firstInitialTorsionAtomPsi; // always sugar
            clipper::MAtom secondInitialTorsionAtomPsi; 
            clipper::MAtom thirdInitialTorsionAtomPsi; // always protein
            clipper::MAtom fourthInitialTorsionAtomPsi; // always protein

            clipper::MAtom firstInitialTorsionAtomPhi; // always sugar
            clipper::MAtom secondInitialTorsionAtomPhi; 
            clipper::MAtom thirdInitialTorsionAtomPhi; // always protein
            clipper::MAtom fourthInitialTorsionAtomPhi; // always protein

            firstInitialTorsionAtomPsi = converted_mglycan[0].find(psiTorsionAtoms[0].first.id().trim(), clipper::MM::UNIQUE);
            secondInitialTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[1].first.id().trim(), clipper::MM::UNIQUE);
            thirdInitialTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[2].first.id().trim(), clipper::MM::UNIQUE);
            fourthInitialTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[3].first.id().trim(), clipper::MM::UNIQUE);

            firstInitialTorsionAtomPhi = converted_mglycan[0].find(phiTorsionAtoms[0].first.id().trim(), clipper::MM::UNIQUE);
            secondInitialTorsionAtomPhi = converted_mglycan[0].find(phiTorsionAtoms[1].first.id().trim(), clipper::MM::UNIQUE);
            thirdInitialTorsionAtomPhi = protein_residue.find(phiTorsionAtoms[2].first.id().trim(), clipper::MM::UNIQUE);
            fourthInitialTorsionAtomPhi = protein_residue.find(phiTorsionAtoms[3].first.id().trim(), clipper::MM::UNIQUE);

            double referencePhiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(firstInitialTorsionAtomPhi.coord_orth(), secondInitialTorsionAtomPhi.coord_orth(), thirdInitialTorsionAtomPhi.coord_orth(), fourthInitialTorsionAtomPhi.coord_orth()));
            double referencePsiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(firstInitialTorsionAtomPsi.coord_orth(), secondInitialTorsionAtomPsi.coord_orth(), thirdInitialTorsionAtomPsi.coord_orth(), fourthInitialTorsionAtomPsi.coord_orth()));
            
            bestPerformingPsi = referencePhiTorsionAngle;
            bestPerformingPhi = referencePhiTorsionAngle;
            best_performing_glycan = converted_mglycan;

            if(useParallelism && nThreads >= 2)
            {
                std::vector<std::pair<double, double>> iteration_data; // .first = Psi, .second = Phi;
            
                for(double currentPsiIterator = -psiError; currentPsiIterator <= +psiError; currentPsiIterator += iteration_step)
                {
                    for(double currentPhiIterator = -phiError; currentPhiIterator <= +phiError; currentPhiIterator += iteration_step)
                    {
                        auto currentPair = std::make_pair(currentPsiIterator, currentPhiIterator);
                        iteration_data.push_back(currentPair);
                    }
                }

                struct ClashMinimizedGlycan 
                {
                    clipper::MPolymer glycan;
                    std::pair<double, double> torsion_angles;
                    std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > clashes;
                };

                std::vector<std::future<void>> thread_results;
                std::vector<ClashMinimizedGlycan> results(iteration_data.size());
                size_t iterations_per_thread = iteration_data.size() / nThreads + 1;
                size_t start = 0, end;
                for (size_t i = 0; i < nThreads; ++i)
                {
                    end = std::min(start + iterations_per_thread, iteration_data.size());
                    thread_results.push_back(std::async(std::launch::async,
                        [&]( std::vector<std::pair<double, double>>& iteration_data, std::vector<ClashMinimizedGlycan>& results, clipper::MPolymer& converted_mglycan, clipper::MMonomer& protein_residue, 
                            clipper::MiniMol& export_model, std::vector<std::pair<clipper::MAtom, std::string>>& phiTorsionAtoms, std::vector<std::pair<clipper::MAtom, std::string>>& psiTorsionAtoms, 
                            clipper::String root_chain_id, clipper::String root_sugar_chain_id, size_t start, size_t end, bool& debug_output )
                        {
                            for (size_t index = start; index < end; ++index)
                            {    
                                clipper::MPolymer manipulated_glycan = converted_mglycan;
                                clipper::MiniMol tmp_clash_model = export_model;

                                clipper::MAtom firstSubsequentTorsionAtomPsi = manipulated_glycan[0].find(psiTorsionAtoms[0].first.id().trim(), clipper::MM::UNIQUE);
                                clipper::MAtom secondSubsequentTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[1].first.id().trim(), clipper::MM::UNIQUE);
                                clipper::MAtom thirdSubsequentTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[2].first.id().trim(), clipper::MM::UNIQUE);
                                clipper::MAtom fourthSubsequentTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[3].first.id().trim(), clipper::MM::UNIQUE);

                                clipper::MAtom firstSubsequentTorsionAtomPhi = manipulated_glycan[0].find(phiTorsionAtoms[0].first.id().trim(), clipper::MM::UNIQUE);
                                clipper::MAtom secondSubsequentTorsionAtomPhi = manipulated_glycan[0].find(phiTorsionAtoms[1].first.id().trim(), clipper::MM::UNIQUE);
                                clipper::MAtom thirdSubsequentTorsionAtomPhi = protein_residue.find(phiTorsionAtoms[2].first.id().trim(), clipper::MM::UNIQUE);
                                clipper::MAtom fourthSubsequentTorsionAtomPhi = protein_residue.find(phiTorsionAtoms[3].first.id().trim(), clipper::MM::UNIQUE);

                                clipper::Coord_orth psiDirection = thirdSubsequentTorsionAtomPsi.coord_orth() - secondSubsequentTorsionAtomPsi.coord_orth(); // ND2(origin_shift) - C1(base)
                                clipper::Coord_orth phiDirection = thirdSubsequentTorsionAtomPhi.coord_orth() - secondSubsequentTorsionAtomPhi.coord_orth(); // ND2(origin_shift) - C1(base)

                                for(int residue = 0; residue < manipulated_glycan.size(); residue++)
                                {
                                    clipper::MMonomer currentResidue = manipulated_glycan[residue];
                                    for(int atom = 0; atom < manipulated_glycan[residue].size(); atom++)
                                    {
                                        clipper::MAtom currentAtom = manipulated_glycan[residue][atom];
                                        clipper::Coord_orth old_pos = currentAtom.coord_orth();
                                        clipper::Coord_orth new_pos = generate_rotation_matrix_from_rodrigues_rotation_formula(psiDirection, old_pos, thirdSubsequentTorsionAtomPsi.coord_orth(), clipper::Util::d2rad(-iteration_data[index].first));
                                        manipulated_glycan[residue][atom].set_coord_orth(new_pos);
                                    }
                                }

                                firstSubsequentTorsionAtomPsi = manipulated_glycan[0].find(psiTorsionAtoms[0].first.id().trim(), clipper::MM::UNIQUE);
                                secondSubsequentTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[1].first.id().trim(), clipper::MM::UNIQUE);
                                thirdSubsequentTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[2].first.id().trim(), clipper::MM::UNIQUE);
                                fourthSubsequentTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[3].first.id().trim(), clipper::MM::UNIQUE);

                                firstSubsequentTorsionAtomPhi = manipulated_glycan[0].find(phiTorsionAtoms[0].first.id().trim(), clipper::MM::UNIQUE);
                                secondSubsequentTorsionAtomPhi = manipulated_glycan[0].find(phiTorsionAtoms[1].first.id().trim(), clipper::MM::UNIQUE);
                                thirdSubsequentTorsionAtomPhi = protein_residue.find(phiTorsionAtoms[2].first.id().trim(), clipper::MM::UNIQUE);
                                fourthSubsequentTorsionAtomPhi = protein_residue.find(phiTorsionAtoms[3].first.id().trim(), clipper::MM::UNIQUE);


                                for(int residue = 0; residue < manipulated_glycan.size(); residue++)
                                {
                                    clipper::MMonomer currentResidue = manipulated_glycan[residue];
                                    for(int atom = 0; atom < manipulated_glycan[residue].size(); atom++)
                                    {
                                        clipper::MAtom currentAtom = manipulated_glycan[residue][atom];
                                        clipper::Coord_orth old_pos = currentAtom.coord_orth();
                                        clipper::Coord_orth new_pos = generate_rotation_matrix_from_rodrigues_rotation_formula(phiDirection, old_pos, thirdSubsequentTorsionAtomPhi.coord_orth(), clipper::Util::d2rad(-iteration_data[index].second));
                                        manipulated_glycan[residue][atom].set_coord_orth(new_pos);
                                    }
                                }

                                firstSubsequentTorsionAtomPsi = manipulated_glycan[0].find(psiTorsionAtoms[0].first.id().trim(), clipper::MM::UNIQUE);
                                secondSubsequentTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[1].first.id().trim(), clipper::MM::UNIQUE);
                                thirdSubsequentTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[2].first.id().trim(), clipper::MM::UNIQUE);
                                fourthSubsequentTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[3].first.id().trim(), clipper::MM::UNIQUE);

                                firstSubsequentTorsionAtomPhi = manipulated_glycan[0].find(phiTorsionAtoms[0].first.id().trim(), clipper::MM::UNIQUE);
                                secondSubsequentTorsionAtomPhi = manipulated_glycan[0].find(phiTorsionAtoms[1].first.id().trim(), clipper::MM::UNIQUE);
                                thirdSubsequentTorsionAtomPhi = protein_residue.find(phiTorsionAtoms[2].first.id().trim(), clipper::MM::UNIQUE);
                                fourthSubsequentTorsionAtomPhi = protein_residue.find(phiTorsionAtoms[3].first.id().trim(), clipper::MM::UNIQUE);

                                double currentPhiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(firstSubsequentTorsionAtomPhi.coord_orth(), secondSubsequentTorsionAtomPhi.coord_orth(), thirdSubsequentTorsionAtomPhi.coord_orth(), fourthSubsequentTorsionAtomPhi.coord_orth()));
                                double currentPsiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(firstSubsequentTorsionAtomPsi.coord_orth(), secondSubsequentTorsionAtomPsi.coord_orth(), thirdSubsequentTorsionAtomPsi.coord_orth(), fourthSubsequentTorsionAtomPsi.coord_orth()));
                                
                                tmp_clash_model.insert(manipulated_glycan);
                                std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > current_clashes = check_for_clashes(tmp_clash_model, manipulated_glycan[0], protein_residue, root_chain_id, root_sugar_chain_id);

                                ClashMinimizedGlycan output {manipulated_glycan, iteration_data[index], current_clashes };
                                results[index] = output; 
                            }

                        },
                        std::ref(iteration_data), std::ref(results), std::ref(converted_mglycan), std::ref(protein_residue), std::ref(export_model), std::ref(phiTorsionAtoms), std::ref(psiTorsionAtoms), root_chain_id, root_sugar_chain_id, start, end, std::ref(debug_output)
                    ));
                    start += iterations_per_thread;
                }

                for (auto& r: thread_results)
                    r.get();

                thread_results.clear();

                std::sort(results.begin(), results.end(), [](ClashMinimizedGlycan& a, ClashMinimizedGlycan& b) {
                    return a.clashes.size() < b.clashes.size();
                });

                if(!results.empty())
                {
                    if(results[0].clashes.size() == 0)
                    {
                        std::vector<ClashMinimizedGlycan> clashes_eliminated;
                        int index = 0;
                        while(results[index].clashes.size() == 0)
                        {
                            clashes_eliminated.push_back(results[index]);
                            index++;
                        }
                        std::sort(clashes_eliminated.begin(), clashes_eliminated.end(), [](ClashMinimizedGlycan& lhs, ClashMinimizedGlycan& rhs) {
                            auto lhs_torsion_angles = lhs.torsion_angles;
                            auto rhs_torsion_angles = rhs.torsion_angles;

                            double lhs_aggregatedAngle = abs(lhs_torsion_angles.first) + abs(lhs_torsion_angles.second);
                            double rhs_aggregatedAngle = abs(lhs_torsion_angles.first) + abs(lhs_torsion_angles.second);

                            return lhs_aggregatedAngle < rhs_aggregatedAngle;
                        });
                        
                        if(!clashes_eliminated.empty())
                            return clashes_eliminated[0].glycan;
                    }
                    else
                    {
                        int currentNumClashes = results[0].clashes.size();
                        std::vector<ClashMinimizedGlycan> clashes_minimized;
                        int index = 0;
                        while(results[index].clashes.size() == currentNumClashes)
                        {
                            clashes_minimized.push_back(results[index]);
                            index++;
                        }
                        std::sort(clashes_minimized.begin(), clashes_minimized.end(), [&](ClashMinimizedGlycan& lhs, ClashMinimizedGlycan& rhs) {
                            auto lhs_clashes = lhs.clashes;
                            auto rhs_clashes = rhs.clashes;

                            double lhs_averageResidueDistance = get_average_distance_between_clashing_residues(lhs_clashes);
                            double rhs_averageResidueDistance = get_average_distance_between_clashing_residues(rhs_clashes);

                            return lhs_averageResidueDistance > rhs_averageResidueDistance;
                        });
                        
                        if(!clashes_minimized.empty())
                            return clashes_minimized[0].glycan;
                    }
                }

            }

            return best_performing_glycan;
        }        
        
   
        clipper::MPolymer Grafter::rotate_mglycan_until_clashes_are_minimized_singlethreaded(clipper::MiniMol& export_model, clipper::MPolymer& converted_mglycan, clipper::MMonomer& protein_residue, std::vector<std::pair<clipper::MAtom, std::string>>& phiTorsionAtoms, std::vector<std::pair<clipper::MAtom, std::string>>& psiTorsionAtoms, double phiError, double psiError, clipper::String root_chain_id, clipper::String root_sugar_chain_id, bool debug_output)
        {
            clipper::MiniMol tmp_clash_model = export_model;
            clipper::MPolymer manipulated_glycan = converted_mglycan;
            clipper::MPolymer best_performing_glycan = converted_mglycan;
            double bestPerformingPsi, bestPerformingPhi;
            double previousAverageDistanceBetweenResidues = 0.0;

            tmp_clash_model.insert(manipulated_glycan);
            std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > reference_clashes = check_for_clashes(tmp_clash_model, converted_mglycan[0], protein_residue, root_chain_id, root_sugar_chain_id);
            // std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > current_clashes;
            int n_reference_clashes = reference_clashes.size();
            int n_most_recent_clashes = n_reference_clashes;

            tmp_clash_model = export_model;

            clipper::MAtom firstTorsionAtomPsi; // always sugar
            clipper::MAtom secondTorsionAtomPsi; 
            clipper::MAtom thirdTorsionAtomPsi; // always protein
            clipper::MAtom fourthTorsionAtomPsi; // always protein

            clipper::MAtom firstTorsionAtomPhi; // always sugar
            clipper::MAtom secondTorsionAtomPhi; 
            clipper::MAtom thirdTorsionAtomPhi; // always protein
            clipper::MAtom fourthTorsionAtomPhi; // always protein

            firstTorsionAtomPsi = manipulated_glycan[0].find(psiTorsionAtoms[0].first.id().trim(), clipper::MM::UNIQUE);
            secondTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[1].first.id().trim(), clipper::MM::UNIQUE);
            thirdTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[2].first.id().trim(), clipper::MM::UNIQUE);
            fourthTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[3].first.id().trim(), clipper::MM::UNIQUE);

            firstTorsionAtomPhi = manipulated_glycan[0].find(phiTorsionAtoms[0].first.id().trim(), clipper::MM::UNIQUE);
            secondTorsionAtomPhi = manipulated_glycan[0].find(phiTorsionAtoms[1].first.id().trim(), clipper::MM::UNIQUE);
            thirdTorsionAtomPhi = protein_residue.find(phiTorsionAtoms[2].first.id().trim(), clipper::MM::UNIQUE);
            fourthTorsionAtomPhi = protein_residue.find(phiTorsionAtoms[3].first.id().trim(), clipper::MM::UNIQUE);

            double currentPhiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(firstTorsionAtomPhi.coord_orth(), secondTorsionAtomPhi.coord_orth(), thirdTorsionAtomPhi.coord_orth(), fourthTorsionAtomPhi.coord_orth()));
            double currentPsiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(firstTorsionAtomPsi.coord_orth(), secondTorsionAtomPsi.coord_orth(), thirdTorsionAtomPsi.coord_orth(), fourthTorsionAtomPsi.coord_orth()));
            
            bestPerformingPsi = currentPsiTorsionAngle;
            bestPerformingPhi = currentPhiTorsionAngle;
            best_performing_glycan = converted_mglycan;

            for(double currentPsiIterator = -psiError; currentPsiIterator <= +psiError; currentPsiIterator += 5.0)
            {
                for(double currentPhiIterator = -phiError; currentPhiIterator <= +phiError; currentPhiIterator += 5.0)
                {
                    manipulated_glycan = converted_mglycan;
                    tmp_clash_model = export_model;

                    firstTorsionAtomPsi = manipulated_glycan[0].find(psiTorsionAtoms[0].first.id().trim(), clipper::MM::UNIQUE);
                    secondTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[1].first.id().trim(), clipper::MM::UNIQUE);
                    thirdTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[2].first.id().trim(), clipper::MM::UNIQUE);
                    fourthTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[3].first.id().trim(), clipper::MM::UNIQUE);

                    firstTorsionAtomPhi = manipulated_glycan[0].find(phiTorsionAtoms[0].first.id().trim(), clipper::MM::UNIQUE);
                    secondTorsionAtomPhi = manipulated_glycan[0].find(phiTorsionAtoms[1].first.id().trim(), clipper::MM::UNIQUE);
                    thirdTorsionAtomPhi = protein_residue.find(phiTorsionAtoms[2].first.id().trim(), clipper::MM::UNIQUE);
                    fourthTorsionAtomPhi = protein_residue.find(phiTorsionAtoms[3].first.id().trim(), clipper::MM::UNIQUE);

                    clipper::Coord_orth psiDirection = thirdTorsionAtomPsi.coord_orth() - secondTorsionAtomPsi.coord_orth(); // ND2(origin_shift) - C1(base)
                    clipper::Coord_orth phiDirection = thirdTorsionAtomPhi.coord_orth() - secondTorsionAtomPhi.coord_orth(); // ND2(origin_shift) - C1(base)

                    for(int residue = 0; residue < manipulated_glycan.size(); residue++)
                    {
                        clipper::MMonomer currentResidue = manipulated_glycan[residue];
                        for(int atom = 0; atom < manipulated_glycan[residue].size(); atom++)
                        {
                            clipper::MAtom currentAtom = manipulated_glycan[residue][atom];
                            clipper::Coord_orth old_pos = currentAtom.coord_orth();
                            clipper::Coord_orth new_pos = generate_rotation_matrix_from_rodrigues_rotation_formula(psiDirection, old_pos, thirdTorsionAtomPsi.coord_orth(), clipper::Util::d2rad(-currentPsiIterator));
                            manipulated_glycan[residue][atom].set_coord_orth(new_pos);
                        }
                    }

                    firstTorsionAtomPsi = manipulated_glycan[0].find(psiTorsionAtoms[0].first.id().trim(), clipper::MM::UNIQUE);
                    secondTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[1].first.id().trim(), clipper::MM::UNIQUE);
                    thirdTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[2].first.id().trim(), clipper::MM::UNIQUE);
                    fourthTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[3].first.id().trim(), clipper::MM::UNIQUE);

                    firstTorsionAtomPhi = manipulated_glycan[0].find(phiTorsionAtoms[0].first.id().trim(), clipper::MM::UNIQUE);
                    secondTorsionAtomPhi = manipulated_glycan[0].find(phiTorsionAtoms[1].first.id().trim(), clipper::MM::UNIQUE);
                    thirdTorsionAtomPhi = protein_residue.find(phiTorsionAtoms[2].first.id().trim(), clipper::MM::UNIQUE);
                    fourthTorsionAtomPhi = protein_residue.find(phiTorsionAtoms[3].first.id().trim(), clipper::MM::UNIQUE);


                    for(int residue = 0; residue < manipulated_glycan.size(); residue++)
                    {
                        clipper::MMonomer currentResidue = manipulated_glycan[residue];
                        for(int atom = 0; atom < manipulated_glycan[residue].size(); atom++)
                        {
                            clipper::MAtom currentAtom = manipulated_glycan[residue][atom];
                            clipper::Coord_orth old_pos = currentAtom.coord_orth();
                            clipper::Coord_orth new_pos = generate_rotation_matrix_from_rodrigues_rotation_formula(phiDirection, old_pos, thirdTorsionAtomPhi.coord_orth(), clipper::Util::d2rad(-currentPhiIterator));
                            manipulated_glycan[residue][atom].set_coord_orth(new_pos);
                        }
                    }

                    firstTorsionAtomPsi = manipulated_glycan[0].find(psiTorsionAtoms[0].first.id().trim(), clipper::MM::UNIQUE);
                    secondTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[1].first.id().trim(), clipper::MM::UNIQUE);
                    thirdTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[2].first.id().trim(), clipper::MM::UNIQUE);
                    fourthTorsionAtomPsi = protein_residue.find(psiTorsionAtoms[3].first.id().trim(), clipper::MM::UNIQUE);

                    firstTorsionAtomPhi = manipulated_glycan[0].find(phiTorsionAtoms[0].first.id().trim(), clipper::MM::UNIQUE);
                    secondTorsionAtomPhi = manipulated_glycan[0].find(phiTorsionAtoms[1].first.id().trim(), clipper::MM::UNIQUE);
                    thirdTorsionAtomPhi = protein_residue.find(phiTorsionAtoms[2].first.id().trim(), clipper::MM::UNIQUE);
                    fourthTorsionAtomPhi = protein_residue.find(phiTorsionAtoms[3].first.id().trim(), clipper::MM::UNIQUE);

                    currentPhiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(firstTorsionAtomPhi.coord_orth(), secondTorsionAtomPhi.coord_orth(), thirdTorsionAtomPhi.coord_orth(), fourthTorsionAtomPhi.coord_orth()));
                    currentPsiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(firstTorsionAtomPsi.coord_orth(), secondTorsionAtomPsi.coord_orth(), thirdTorsionAtomPsi.coord_orth(), fourthTorsionAtomPsi.coord_orth()));
                    
                    tmp_clash_model.insert(manipulated_glycan);
                    std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > current_clashes = check_for_clashes(tmp_clash_model, manipulated_glycan[0], protein_residue, root_chain_id, root_sugar_chain_id);
                    n_most_recent_clashes = current_clashes.size();

                    if(debug_output)
                        DBG << "currentPsiIterator: " << currentPsiIterator << " currentPhiIterator: " << currentPhiIterator << "\t\tClashes after rotation: " << current_clashes.size() << "\t\tCurrent value of Psi torsion: " << currentPsiTorsionAngle << "\t\tCurrent value of Phi torsion: " << currentPhiTorsionAngle << std::endl;
                    
                    n_most_recent_clashes = current_clashes.size();

                    if(n_most_recent_clashes < n_reference_clashes || n_most_recent_clashes == n_reference_clashes)
                    {
                        if(n_most_recent_clashes < n_reference_clashes)
                        {
                            bestPerformingPsi = currentPsiTorsionAngle;
                            bestPerformingPhi = currentPhiTorsionAngle;
                            best_performing_glycan = manipulated_glycan;
                            n_reference_clashes = n_most_recent_clashes;
                            previousAverageDistanceBetweenResidues = 0.0;
                        }
                        else
                        {
                            double currentAverageDistanceBetweenResidues = 0.0;
                            double totalAveragesOfPerResidueDistances = 0.0;

                            for(int i = 0; i < current_clashes.size(); i++)
                            {
                                int n_protein_side_chain_atoms = current_clashes[i].first.first.size();
                                int n_sugar_atoms = current_clashes[i].second.first.size();
                                double averageDistance = 0.0;
                                double totalDistance = 0.0;
                                int totalIterations = 0;
                                
                                for(int sugarAtom = 0; sugarAtom < n_sugar_atoms; sugarAtom++)
                                {
                                    for(int proteinResidueAtom = 0; proteinResidueAtom < n_protein_side_chain_atoms; proteinResidueAtom++)
                                    {
                                        double currentDistance = clipper::Coord_orth::length(current_clashes[i].first.first[proteinResidueAtom].coord_orth(), current_clashes[i].second.first[sugarAtom].coord_orth());
                                        totalDistance = totalDistance + currentDistance;
                                        totalIterations++;
                                    }
                                }

                                averageDistance = totalDistance / totalIterations;
                                totalAveragesOfPerResidueDistances = totalAveragesOfPerResidueDistances + averageDistance;
                            }
                            
                            currentAverageDistanceBetweenResidues = totalAveragesOfPerResidueDistances / n_most_recent_clashes;

                            if(currentAverageDistanceBetweenResidues > previousAverageDistanceBetweenResidues)
                            {
                                bestPerformingPsi = currentPsiTorsionAngle;
                                bestPerformingPhi = currentPhiTorsionAngle;
                                best_performing_glycan = manipulated_glycan;
                                n_reference_clashes = n_most_recent_clashes;
                                previousAverageDistanceBetweenResidues = currentAverageDistanceBetweenResidues;
                            }
                        }
                    }

                    if(current_clashes.empty())
                        return best_performing_glycan;

                    manipulated_glycan = converted_mglycan;
                    tmp_clash_model = export_model;
                }
            }

            return best_performing_glycan;
        }

        std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > Grafter::check_for_clashes(clipper::MiniMol& input_model, clipper::MMonomer& root_sugar, clipper::MMonomer& input_protein_side_chain_residue, clipper::String root_chain_id, clipper::String root_sugar_chain_id)
        {
            clipper::MAtomNonBond manb = clipper::MAtomNonBond(input_model, 1.0);
            clipper::MGlycology mgl = clipper::MGlycology(input_model, manb, debug_output, "undefined");
            std::vector<clipper::MGlycan> list_of_glycans = mgl.get_list_of_glycans();

            clipper::MGlycan grafted_mglycan;
            clipper::MSugar root_msugar;


            for(int i = 0; i < list_of_glycans.size(); i++)
            {
                clipper::MGlycan currentGlycan = list_of_glycans[i];
                clipper::MSugar currentRootSugar = currentGlycan.get_root().second;

                if(currentGlycan.get_chain().trim() == root_chain_id && currentGlycan.get_root_sugar_chainID().trim() == root_sugar_chain_id && currentRootSugar.type().trim() == root_sugar.type().trim() && currentRootSugar.id().trim() == root_sugar.id().trim())
                {
                    grafted_mglycan = currentGlycan;
                    root_msugar = currentRootSugar;
                }
            }

            std::vector<clipper::MSugar> sugars_in_grafted_mglycan = grafted_mglycan.get_sugars();

            std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > clashing_residues;
            clipper::MAtomNonBond clashmanb( input_model, 3.0 );
            for(int node = 0; node < grafted_mglycan.number_of_nodes(); node++)
            {
                clipper::MGlycan::Node currentNode = grafted_mglycan.get_node(node);
                clipper::MSugar currentSugar = currentNode.get_sugar();
                
                for(int atom = 0; atom < currentSugar.size(); atom++)
                {
                    clipper::Coord_orth currentSugarAtomCoord = currentSugar[atom].coord_orth();
                    std::vector<clipper::MAtomIndexSymmetry> neighbourhood_atoms = clashmanb( currentSugarAtomCoord, 2.5 );


                    for(int i = 0; i < neighbourhood_atoms.size(); i++)
                    {
                        int detected_chain      = neighbourhood_atoms[i].polymer();
                        int detected_monomer    = neighbourhood_atoms[i].monomer();
                        int detected_atom       = neighbourhood_atoms[i].atom();

                        auto search_result = std::find_if(std::begin(sugars_in_grafted_mglycan), std::end(sugars_in_grafted_mglycan),
                        [&input_model, &detected_chain, &detected_monomer](clipper::MSugar& sugar) {
                            return sugar.id().trim() == input_model[detected_chain][detected_monomer].id().trim() && sugar.type().trim() == input_model[detected_chain][detected_monomer].type().trim() &&
                                   sugar.seqnum() == input_model[detected_chain][detected_monomer].seqnum();
                        });

                        if(search_result == std::end(sugars_in_grafted_mglycan) && clipper::Coord_orth::length(currentSugarAtomCoord, input_model[detected_chain][detected_monomer][detected_atom].coord_orth()) < 2.25)
                        {
                            auto previously_identified = std::find_if(std::begin(clashing_residues), std::end(clashing_residues),
                            [&](std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> >& element) {
                                return  element.first.first.type().trim() == input_model[detected_chain][detected_monomer].type().trim() && element.first.first.id().trim() == input_model[detected_chain][detected_monomer].id().trim() && element.first.first.seqnum() == input_model[detected_chain][detected_monomer].seqnum() && 
                                        element.second.first.type().trim() == currentSugar.type().trim() && element.second.first.id().trim() == currentSugar.id().trim() && element.second.first.seqnum() == currentSugar.seqnum() &&
                                        element.first.second == root_chain_id && element.second.second == root_sugar_chain_id;
                            });

                            if(previously_identified == std::end(clashing_residues))
                            {
                                // Make sure that we are not detecting root protein side-chain as a clash.
                                if(root_chain_id.trim() != input_model[detected_chain].id().trim() || input_protein_side_chain_residue.id().trim() != input_model[detected_chain][detected_monomer].id().trim() || input_protein_side_chain_residue.type().trim() != input_model[detected_chain][detected_monomer].type().trim() || input_protein_side_chain_residue.seqnum() != input_model[detected_chain][detected_monomer].seqnum())
                                    clashing_residues.push_back( std::make_pair(std::make_pair(input_model[detected_chain][detected_monomer], root_chain_id), std::make_pair(currentSugar, root_sugar_chain_id)) );
                            }
                        }
                        
                        // DBG << "Detected " << input_model[detected_chain].id().trim() << "/" << input_model[detected_chain][detected_monomer].type().trim() << "-" << input_model[detected_chain][detected_monomer].id().trim()
                            // << " with the distance of: " << clipper::Coord_orth::length(currentSugarAtomCoord, input_model[detected_chain][detected_monomer][detected_atom].coord_orth()) << "\t\t\t" << currentSugar[atom].id().trim() << std::endl;
                    }
                }
            }

            if(debug_output)
            {
                for(int i = 0; i < clashing_residues.size(); i++)
                {
                     DBG << "Detected clash: " << clashing_residues[i].first.second << "/" << clashing_residues[i].first.first.id().trim() << "-" << clashing_residues[i].first.first.type().trim() << "\t\t\tSugar: " << clashing_residues[i].second.second << "/" << clashing_residues[i].second.first.id().trim() << "-" << clashing_residues[i].second.first.type().trim() << std::endl;
                }
            }

            return clashing_residues;
        }

        clipper::MiniMol Grafter::trim_graft_until_no_clashes_left(clipper::MiniMol& current_model, clipper::MGlycan& grafted_glycan, std::pair<clipper::MMonomer, clipper::String> root_residue, clipper::String graft_chain_id, std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > current_clashes)
        {
            clipper::MGlycan lastNode_removed = grafted_glycan, leafNode_removed = grafted_glycan;
            clipper::MiniMol tmp_model = current_model;
            clipper::MiniMol output_model = current_model;
            std::vector<std::pair<clipper::MMonomer, clipper::String>> residues_to_delete_lastNode, residues_to_delete_leafNode;
            int amountOfCurrentClashingResidues = current_clashes.size();
            int totalNodes = grafted_glycan.number_of_nodes();

            int residueDeletions = 0;


            for(int i = totalNodes; i > 1; i--)
            {
                lastNode_removed = remove_last_node(lastNode_removed);
                leafNode_removed = remove_first_leaf_node(leafNode_removed);

                residues_to_delete_lastNode = get_designated_residues_for_deletion(grafted_glycan, lastNode_removed, graft_chain_id);
                residues_to_delete_leafNode = get_designated_residues_for_deletion(grafted_glycan, leafNode_removed, graft_chain_id);


                clipper::MiniMol deleted_lastNode_model = get_model_with_trimmed_glycan(tmp_model, residues_to_delete_lastNode);
                clipper::MiniMol deleted_leafNode_model = get_model_with_trimmed_glycan(tmp_model, residues_to_delete_leafNode);

                clipper::MMonomer root_sugar_lastNode = lastNode_removed.get_node(0).get_sugar();
                clipper::MMonomer root_sugar_leafNode = leafNode_removed.get_node(0).get_sugar();

                std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > lastNode_removed_clashes = check_for_clashes(deleted_lastNode_model, root_sugar_lastNode, root_residue.first, root_residue.second, graft_chain_id);
                std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > leafNode_removed_clashes = check_for_clashes(deleted_leafNode_model, root_sugar_leafNode, root_residue.first, root_residue.second, graft_chain_id);

                if(lastNode_removed_clashes.empty())
                {
                    output_model = deleted_lastNode_model;
                    this->grafted_glycan = lastNode_removed;
                    break;
                }
                else if(leafNode_removed_clashes.empty())
                {
                    output_model = deleted_leafNode_model;
                    this->grafted_glycan = leafNode_removed;
                    break;
                }
                else
                {
                    output_model = deleted_lastNode_model;
                    this->grafted_glycan = lastNode_removed;
                }

                residueDeletions++;
            }
            
            clipper::MMonomer output_model_root_sugar = grafted_glycan.get_node(0).get_sugar();
            std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > output_model_clashes = check_for_clashes(output_model, output_model_root_sugar, root_residue.first, root_residue.second, graft_chain_id);
            this->clashes = output_model_clashes;
            
            if(enable_user_messages)
            {
                std::cout << "Deleted " << residueDeletions << "/" << totalNodes << " sugars in order to remove clashes resulting from glycan grafting. Number of total clashes remaining: " << output_model_clashes.size() << std::endl;
                for(int i = 0; i < output_model_clashes.size(); i++)
                {
                    std::cout << "Detected clash: " << output_model_clashes[i].first.second << "/" << output_model_clashes[i].first.first.id().trim() << "-" << output_model_clashes[i].first.first.type().trim() << "\t\t\tSugar: " << output_model_clashes[i].second.second << "/" << output_model_clashes[i].second.first.id().trim() << "-" << output_model_clashes[i].second.first.type().trim() << std::endl;
                }
            }
            
            if(debug_output)
            {
                DBG << "Deleted " << residueDeletions << "/" << totalNodes << " sugars in order to remove clashes resulting from glycan grafting. Number of total clashes remaining: " << output_model_clashes.size() << std::endl;
                for(int i = 0; i < output_model_clashes.size(); i++)
                {
                     DBG << "Detected clash: " << output_model_clashes[i].first.second << "/" << output_model_clashes[i].first.first.id().trim() << "-" << output_model_clashes[i].first.first.type().trim() << "\t\t\tSugar: " << output_model_clashes[i].second.second << "/" << output_model_clashes[i].second.first.id().trim() << "-" << output_model_clashes[i].second.first.type().trim() << std::endl;
                }
            }
            
            return output_model;
        }

        clipper::MiniMol Grafter::get_model_with_trimmed_glycan(clipper::MiniMol& current_model, std::vector<std::pair<clipper::MMonomer, clipper::String>>& residues_to_delete)
        {
            clipper::MiniMol tmp_model = current_model;
            clipper::String current_chainID = residues_to_delete[0].second.trim();

            clipper::MPolymer chain_to_be_modified = tmp_model.find(current_chainID, clipper::MM::UNIQUE);
            
            for(int monomer = 0; monomer < chain_to_be_modified.size(); monomer++)
            {
                for(int i = 0; i < residues_to_delete.size(); i++)
                {
                    if(chain_to_be_modified[monomer].id().trim() == residues_to_delete[i].first.id().trim())
                        chain_to_be_modified[monomer].set_type("~~~");
                }
            }

            clipper::MiniMol output_model(tmp_model.spacegroup(), tmp_model.cell());

            for(int chain = 0; chain < tmp_model.size(); chain++)
            {
                clipper::MPolymer exportChain;
                exportChain.set_id(tmp_model[chain].id().trim());
                if(current_chainID == exportChain.id().trim())
                {
                    for(int residue = 0; residue < chain_to_be_modified.size(); residue++)
                    {
                        if(chain_to_be_modified[residue].type().trim() != "~~~")
                        {
                            exportChain.insert(chain_to_be_modified[residue]);
                        }
                    }
                }
                else
                {
                    for(int residue = 0; residue < tmp_model[chain].size(); residue++)
                    {
                        exportChain.insert(tmp_model[chain][residue]);
                    }
                }
                output_model.insert(exportChain);
            }

            return output_model;
        }


        std::vector<std::pair<clipper::MMonomer, clipper::String>> Grafter::get_designated_residues_for_deletion(clipper::MGlycan& original_graft, clipper::MGlycan& trimmed_graft, clipper::String graft_chain_id)
        {
            std::vector<std::pair<clipper::MMonomer, clipper::String>> output;
            std::vector<clipper::MGlycan::Node> original_graft_nodes = original_graft.get_node_list_vector();
            std::vector<clipper::MGlycan::Node> trimmed_graft_nodes = trimmed_graft.get_node_list_vector();
            std::vector<clipper::MGlycan::Node> difference;

            std::sort(original_graft_nodes.begin(), original_graft_nodes.end(), compareMGlycanNode);
            std::sort(trimmed_graft_nodes.begin(), trimmed_graft_nodes.end(), compareMGlycanNode);
            std::set_difference(original_graft_nodes.begin(), original_graft_nodes.end(), trimmed_graft_nodes.begin(), trimmed_graft_nodes.end(), std::inserter(difference, difference.end()), compareMGlycanNode);

            for(int i = 0; i < difference.size(); i++)
            {
                clipper::MMonomer convertedMMonomer = difference[i].get_sugar();
                output.push_back(std::make_pair(convertedMMonomer, graft_chain_id));
            }

            return output;
        }

        bool Grafter::check_if_residue_has_hydrogens(clipper::MMonomer residue_to_check)
        {
            for(int atom = 0; atom < residue_to_check.size(); atom++)
            {
                clipper::MAtom currentAtom = residue_to_check[atom];
                if(currentAtom.element().trim() == "H")
                    return true;
            }

            return false;
        }

        clipper::MPolymer Grafter::delete_atom_from_mglycan(clipper::MPolymer& converted_mglycan, clipper::MAtom& atom_to_be_deleted)
        {
            clipper::MPolymer output;
            for(int residue = 0; residue < converted_mglycan.size(); residue++)
            {
                clipper::MMonomer output_residue;
                for(int atom = 0; atom < converted_mglycan[residue].size(); atom++)
                {
                    clipper::MAtom currentAtom = converted_mglycan[residue][atom];
                    if(currentAtom.coord_orth() != atom_to_be_deleted.coord_orth())
                    {
                        output_residue.insert(currentAtom);
                    }
                }
                output_residue.set_id(converted_mglycan[residue].id());
                output_residue.set_type(converted_mglycan[residue].type());
                output_residue.set_seqnum(converted_mglycan[residue].seqnum());
                output.insert(output_residue);
            }
            output.set_id(converted_mglycan.id());
            return output;
        }

        double Grafter::get_average_distance_between_clashing_residues(std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > >& input_clashes)
        {
            double currentAverageDistanceBetweenResidues = 0.0;
            double totalAveragesOfPerResidueDistances = 0.0;

            for(int i = 0; i < input_clashes.size(); i++)
            {
                int n_protein_side_chain_atoms = input_clashes[i].first.first.size();
                int n_sugar_atoms = input_clashes[i].second.first.size();
                double averageDistance = 0.0;
                double totalDistance = 0.0;
                int totalIterations = 0;
                
                for(int sugarAtom = 0; sugarAtom < n_sugar_atoms; sugarAtom++)
                {
                    for(int proteinResidueAtom = 0; proteinResidueAtom < n_protein_side_chain_atoms; proteinResidueAtom++)
                    {
                        double currentDistance = clipper::Coord_orth::length(input_clashes[i].first.first[proteinResidueAtom].coord_orth(), input_clashes[i].second.first[sugarAtom].coord_orth());
                        totalDistance = totalDistance + currentDistance;
                        totalIterations++;
                    }
                }

                averageDistance = totalDistance / totalIterations;
                totalAveragesOfPerResidueDistances = totalAveragesOfPerResidueDistances + averageDistance;
            }
            
            currentAverageDistanceBetweenResidues = totalAveragesOfPerResidueDistances / input_clashes.size();
            return currentAverageDistanceBetweenResidues;
        }

        int Grafter::lookup_protein_backbone_glycosylation_database (clipper::String name)
        {
            for (int i = 0; i < backbone_instructions_size; i++)
            {
                if (name.trim() == backbone_instructions[i].residue_name.trim())
                {
                    return i;
                }
            }
            return -1;
        }

        int Grafter::lookup_glycan_type_glycosylation_database (clipper::String type)
        {
            for (int i = 0; i < sugar_instructions_size; i++)
            {
                if (type.trim() == sugar_instructions[i].glycan_type.trim())
                {
                    return i;
                }
            }

            return -1;
        }
    
    }
}
