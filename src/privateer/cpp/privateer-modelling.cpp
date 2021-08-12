// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#include "privateer-modelling.h"

#define DBG std::cout << "[" << __FUNCTION__ << "] - "

namespace privateer
{
    namespace modelling 
    {
        const protein_sidechain_glycosylation backbone_instructions[] =
        {
            { "ASN", "ND2", "CG", "CB", -100, 180 },
            { "ASP", "OD2", "CG", "CB" -100, 180 }
        };
        const int backbone_instructions_size = sizeof( backbone_instructions ) / sizeof( backbone_instructions[0] );

        const sugar_attachment sugar_instructions[] =
        {
            { "ideal", "C1", "O1", "" },   
            { "non-ideal", "C1", "", "" }
        };
        const int sugar_instructions_size = sizeof( sugar_instructions ) / sizeof( sugar_instructions[0] );
    
        Grafter::Grafter(clipper::MiniMol receiving_model, clipper::MiniMol donor_model, bool enable_user_messages, bool debug_output)
        {
            this->enable_user_messages = enable_user_messages;
            this->debug_output = debug_output;
            this->receiving_model = receiving_model;
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

        clipper::Coord_orth Grafter::get_glycan_target_point(clipper::Coord_orth& connecting_atom, clipper::Coord_orth& vector_origin, clipper::Coord_orth& vector_target, float vectorShiftDistance)
        {
            clipper::Coord_orth coord; 

            clipper::Vec3<clipper::ftype> baseVector((vector_target.x()-vector_origin.x()),(vector_target.y()-vector_origin.y()), (vector_target.z()-vector_origin.z()));
            // Create a 1A unit vector out of baseVector, to be used later in vector shifting
            clipper::Vec3<clipper::ftype> unitVector = baseVector.unit();

            // Obtain coordinates in the middle of suspected glycan density via 5A vector shift. This is the nearest glycan bonded via ND2 atom to ASN residue.
            coord = clipper::Coord_orth( (connecting_atom.x()+(unitVector[0]*vectorShiftDistance)), (connecting_atom.y()+(unitVector[1]*vectorShiftDistance)), (connecting_atom.z()+(unitVector[2]*vectorShiftDistance)) );

            return coord;
        }

        void Grafter::graft_mpolymer_to_receiving_model(clipper::MGlycan& glycan_to_graft, clipper::MMonomer& input_protein_side_chain_residue, bool ANY_search_policy)
        { 

            clipper::Mat33<clipper::ftype> identity_matrix;
            identity_matrix = identity_matrix.identity();

            int receiver_atom_index = lookup_protein_backbone_glycosylation_database(input_protein_side_chain_residue.type().trim());

            clipper::String residue_name;       
            clipper::String connected_atom;          
            clipper::String vector_point_alpha;
            clipper::String vector_point_bravo;
            clipper::ftype targetPhi;        
            clipper::ftype targetPsi; 

            if(receiver_atom_index != -1)
            {
                residue_name = privateer::modelling::backbone_instructions[receiver_atom_index].residue_name;
                connected_atom = privateer::modelling::backbone_instructions[receiver_atom_index].connected_atom;
                vector_point_alpha = privateer::modelling::backbone_instructions[receiver_atom_index].vector_point_alpha;
                vector_point_bravo = privateer::modelling::backbone_instructions[receiver_atom_index].vector_point_bravo;
                targetPhi = privateer::modelling::backbone_instructions[receiver_atom_index].Phi;
                targetPsi = privateer::modelling::backbone_instructions[receiver_atom_index].Psi;

                if(enable_user_messages && !debug_output)
                    std::cout << "Successfully located " << residue_name << " instructions. Will connect glycan to " << connected_atom << " with " << vector_point_alpha << " and " << vector_point_bravo << " used to generate rotation-translation matrix. Targer Phi = " << targetPhi << " target Psi = " << targetPsi << std::endl;

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
            for(int i = 0; i < receiving_model.size(); i++)
            {
                clipper::String currentChainID = receiving_model[i].id().trim();
                used_chain_ids.push_back(currentChainID);
            }

            char new_chain_id = chainids[used_chain_ids.size()];
            std::string tempchainid(1, new_chain_id);
            clipper::String chainID(tempchainid);

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
                clipper::String sugar_connection_atom;

                int glycan_grafting_type = lookup_glycan_type_glycosylation_database(glycan_type);
                if(glycan_grafting_type != -1)
                {
                    sugar_connection_atom = privateer::modelling::sugar_instructions[glycan_grafting_type].connection_atom;
                    

                    if(enable_user_messages && !debug_output)
                        std::cout << "Successfully located " << glycan_type << " glycan grafting instructions. Will connect glycan via " << sugar_connection_atom << " used to generate rotation-translation matrix." << std::endl;

                    if(debug_output)
                        DBG << "Successfully located " << glycan_type << " glycan grafting instructions. Will connect glycan via " << sugar_connection_atom << " used to generate rotation-translation matrix." << std::endl;
                }
                else
                {
                    DBG << "Unable to locate " << glycan_type << " in sugar_attachment_instruction_set from donor model! ...Aborting..." << std::endl;
                    throw std::invalid_argument( "Unable to generate instructions from donor model from input receiver." );
                }

                clipper::Coord_orth protein_connecting_target;
                clipper::Coord_orth protein_vector_point_alpha;
                clipper::Coord_orth protein_vector_point_bravo;
                clipper::Coord_orth sugar_connection_target;

                if(ANY_search_policy)
                {
                    protein_connecting_target = input_protein_side_chain_residue.find(connected_atom, clipper::MM::ANY).coord_orth(); // ND2
                    protein_vector_point_alpha = input_protein_side_chain_residue.find(vector_point_alpha, clipper::MM::ANY).coord_orth(); // CG
                    protein_vector_point_bravo = input_protein_side_chain_residue.find(vector_point_bravo, clipper::MM::ANY).coord_orth(); // CB
                    sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, clipper::MM::ANY).coord_orth(); // C1
                }
                else
                {
                    protein_connecting_target = input_protein_side_chain_residue.find(connected_atom, clipper::MM::UNIQUE).coord_orth(); // ND2
                    protein_vector_point_alpha = input_protein_side_chain_residue.find(vector_point_alpha, clipper::MM::UNIQUE).coord_orth(); // CG
                    protein_vector_point_bravo = input_protein_side_chain_residue.find(vector_point_bravo, clipper::MM::UNIQUE).coord_orth(); // CB
                    sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, clipper::MM::UNIQUE).coord_orth(); // C1
                }


                clipper::Coord_orth potential_C1_position = get_glycan_target_point(protein_connecting_target, protein_vector_point_alpha, protein_vector_point_bravo, 1.5);

                clipper::Coord_orth source(sugar_connection_target);
                std::vector<clipper::Coord_orth> sourceVector;
                sourceVector.push_back(source);

                clipper::Coord_orth target(potential_C1_position);
                std::vector<clipper::Coord_orth> targetVector;
                targetVector.push_back(target);
                
                

                if(enable_user_messages && !debug_output)
                {
                    std::cout << "Coordinates of " << connected_atom << ":\t" << protein_connecting_target.format() << std::endl;
                    std::cout << "Coordinates of " << vector_point_alpha << ":\t" << protein_vector_point_alpha.format() << std::endl;
                    std::cout << "Coordinates of " << vector_point_bravo << ":\t" << protein_vector_point_bravo.format() << std::endl;
                    std::cout << "Coordinates of " << sugar_connection_atom << ":\t" << sugar_connection_target.format() << std::endl;
                    std::cout << "Coordinates of potential C1 position:\t" << potential_C1_position.format() << std::endl;
                    std::cout << "Coordinates of source:\t" << source.format() << std::endl;
                    std::cout << "Coordinates of target:\t" << target.format() << std::endl;
                }

                if(debug_output)
                {
                    DBG << "Coordinates of " << connected_atom << ":\t" << protein_connecting_target.format() << std::endl;
                    DBG << "Coordinates of " << vector_point_alpha << ":\t" << protein_vector_point_alpha.format() << std::endl;
                    DBG << "Coordinates of " << vector_point_bravo << ":\t" << protein_vector_point_bravo.format() << std::endl;
                    DBG << "Coordinates of " << sugar_connection_atom << ":\t" << sugar_connection_target.format() << std::endl;
                    DBG << "Coordinates of potential C1 position:\t" << potential_C1_position.format() << std::endl;
                    DBG << "Coordinates of source:\t" << source.format() << std::endl;
                    DBG << "Coordinates of target:\t" << target.format() << std::endl;
                }

                clipper::RTop_orth relocator(sourceVector, targetVector);
                converted_mglycan.transform(relocator);

                
                if(enable_user_messages && !debug_output)
                    std::cout << "Grafting input glycan with Chain ID of: " << chainID << " to " << residue_name << "-" << input_protein_side_chain_residue.id().trim() << std::endl;
                
                if(debug_output)
                    DBG << "Grafting input glycan with Chain ID of: " << chainID << " to " << residue_name << "-" << input_protein_side_chain_residue.id().trim() << std::endl;
                
                converted_mglycan.set_id(chainID);
                receiving_model.insert(converted_mglycan);

                if(enable_user_messages && !debug_output)
                    std::cout << "Glycan has been grafted!" << std::endl;
                
                if(debug_output)
                    DBG << "Glycan has been grafted!" << std::endl;
            }
            else if(glycan_type == "ideal")
            {
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
                        DBG << "Successfully located " << glycan_type << " glycan grafting instructions. Will connect glycan via " << sugar_connection_atom << " used to generate rotation-translation matrix." << std::endl;
                }
                else
                {
                    DBG << "Unable to locate " << glycan_type << " in sugar_attachment_instruction_set from donor model! ...Aborting..." << std::endl;
                    throw std::invalid_argument( "Unable to generate instructions from donor model from input receiver." );
                }

                clipper::Coord_orth protein_connecting_target;
                clipper::Coord_orth protein_vector_point_alpha;
                clipper::Coord_orth protein_vector_point_bravo;
                clipper::Coord_orth sugar_connection_target;
                clipper::Coord_orth sugar_vector_point_target;

                if(ANY_search_policy)
                {
                    protein_connecting_target = input_protein_side_chain_residue.find(connected_atom, clipper::MM::ANY).coord_orth(); // ND2
                    protein_vector_point_alpha = input_protein_side_chain_residue.find(vector_point_alpha, clipper::MM::ANY).coord_orth(); // CB
                    protein_vector_point_bravo = input_protein_side_chain_residue.find(vector_point_bravo, clipper::MM::ANY).coord_orth(); // CG
                    sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, clipper::MM::ANY).coord_orth(); // C1
                    sugar_vector_point_target = converted_mglycan[0].find(sugar_vector_point, clipper::MM::ANY).coord_orth(); // O1
                }
                else
                {
                    protein_connecting_target = input_protein_side_chain_residue.find(connected_atom, clipper::MM::UNIQUE).coord_orth(); // ND2
                    protein_vector_point_alpha = input_protein_side_chain_residue.find(vector_point_alpha, clipper::MM::UNIQUE).coord_orth(); // CB
                    protein_vector_point_bravo = input_protein_side_chain_residue.find(vector_point_bravo, clipper::MM::UNIQUE).coord_orth(); // CG
                    sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, clipper::MM::UNIQUE).coord_orth(); // C1
                    sugar_vector_point_target = converted_mglycan[0].find(sugar_vector_point, clipper::MM::UNIQUE).coord_orth(); // O1
                }


                clipper::Vec3<clipper::ftype> O1toND2vector(protein_connecting_target.x() - sugar_vector_point_target.x(), 
                                                            protein_connecting_target.y() - sugar_vector_point_target.y(), 
                                                            protein_connecting_target.z() - sugar_vector_point_target.z());
                
    
                clipper::RTop_orth ND2_O1_overlayer(identity_matrix, O1toND2vector);
                converted_mglycan.transform(ND2_O1_overlayer);

                // Update the coordinates of C1 and O1, get O5 coords from the glycan.
                clipper::String ring_oxygen_name = glycan_to_graft.get_sugars()[0].ring_members()[0].id().trim();
                clipper::Coord_orth ring_oxygen_coords;
                if(ANY_search_policy)
                {
                    sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, clipper::MM::ANY).coord_orth(); // C1
                    sugar_vector_point_target = converted_mglycan[0].find(sugar_vector_point, clipper::MM::ANY).coord_orth(); // O1
                    ring_oxygen_coords = converted_mglycan[0].find(ring_oxygen_name, clipper::MM::ANY).coord_orth(); // O5
                }
                else
                {
                    sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, clipper::MM::UNIQUE).coord_orth(); // C1
                    sugar_vector_point_target = converted_mglycan[0].find(sugar_vector_point, clipper::MM::UNIQUE).coord_orth(); // O1
                    ring_oxygen_coords = converted_mglycan[0].find(ring_oxygen_name, clipper::MM::UNIQUE).coord_orth(); // O5
                }

                clipper::Vec3<clipper::ftype> C1toND2vector(protein_connecting_target.x() - sugar_connection_target.x(), 
                                                            protein_connecting_target.y() - sugar_connection_target.y(), 
                                                            protein_connecting_target.z() - sugar_connection_target.z());
                
                clipper::ftype currentPhiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(ring_oxygen_coords, sugar_connection_target, protein_connecting_target, protein_vector_point_alpha));
                
                if(debug_output)
                    DBG << "Phi = " << currentPhiTorsionAngle << "\t\t" << ring_oxygen_name << "-" << sugar_connection_atom << "-" << connected_atom << "-" << vector_point_alpha << std::endl;
                    

                bool isCurrentPhiNegative = std::signbit(currentPhiTorsionAngle);
                clipper::ftype phiDifference;
                if(isCurrentPhiNegative)
                    phiDifference = targetPhi - currentPhiTorsionAngle;
                else
                    phiDifference = targetPhi + currentPhiTorsionAngle;

                if(debug_output)
                {
                    DBG << "Coordinates of " << sugar_connection_atom << ":\t" << sugar_connection_target.format() << std::endl;
                    DBG << "Coordinates of " << sugar_vector_point << ":\t" << sugar_vector_point_target.format() << std::endl;
                    DBG << "Coordinates of " << ring_oxygen_name << ":\t" << ring_oxygen_coords.format() << std::endl;
                }
                
                clipper::Mat33<clipper::ftype> rodrigues_rotation_matrix_for_phi = generate_rotation_matrix_from_rodrigues_rotation_formula(phiDifference, C1toND2vector);

                if(debug_output)
                {
                    clipper::RTop_orth phiRotator(rodrigues_rotation_matrix_for_phi);
                    DBG << "targetPhi: " << targetPhi << "\t\tcurrentPhiTorsionAngle: " << currentPhiTorsionAngle << "\t\tphiDifference: " << phiDifference << "\t\t" << ring_oxygen_name << "-" << sugar_connection_atom << "-" << connected_atom << "-" << vector_point_alpha << std::endl;
                    // DBG << "RTop_orth: \n" << phiRotator.format() << std::endl;
                }
                   
                if(enable_user_messages && !debug_output)
                    std::cout << "Target value of Phi torsion angle: " << targetPhi << "\t\tValue of Phi after translation: " << currentPhiTorsionAngle << "\t\tRotation along will be done by: " << phiDifference << " degrees." << "\t\t" << ring_oxygen_name << "-" << sugar_connection_atom << "-" << connected_atom << "-" << vector_point_alpha << std::endl;
                
                clipper::RTop_orth phiRotator(rodrigues_rotation_matrix_for_phi);
                converted_mglycan.transform(phiRotator);

                if(ANY_search_policy)
                {
                    sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, clipper::MM::ANY).coord_orth(); // C1
                    sugar_vector_point_target = converted_mglycan[0].find(sugar_vector_point, clipper::MM::ANY).coord_orth(); // O1
                    ring_oxygen_coords = converted_mglycan[0].find(ring_oxygen_name, clipper::MM::ANY).coord_orth(); // O5
                }
                else
                {
                    sugar_connection_target = converted_mglycan[0].find(sugar_connection_atom, clipper::MM::UNIQUE).coord_orth(); // C1
                    sugar_vector_point_target = converted_mglycan[0].find(sugar_vector_point, clipper::MM::UNIQUE).coord_orth(); // O1
                    ring_oxygen_coords = converted_mglycan[0].find(ring_oxygen_name, clipper::MM::UNIQUE).coord_orth(); // O5
                }

                if(debug_output)
                {
                    // DBG << "Coordinates of " << connected_atom << ":\t" << protein_connecting_target.format() << std::endl;
                    // DBG << "Coordinates of " << vector_point_alpha << ":\t" << protein_vector_point_alpha.format() << std::endl;
                    // DBG << "Coordinates of " << vector_point_bravo << ":\t" << protein_vector_point_bravo.format() << std::endl;
                    DBG << "Coordinates of " << sugar_connection_atom << ":\t" << sugar_connection_target.format() << std::endl;
                    DBG << "Coordinates of " << sugar_vector_point << ":\t" << sugar_vector_point_target.format() << std::endl;
                    DBG << "Coordinates of " << ring_oxygen_name << ":\t" << ring_oxygen_coords.format() << std::endl;
                    currentPhiTorsionAngle = clipper::Util::rad2d(clipper::Coord_orth::torsion(ring_oxygen_coords, sugar_connection_target, protein_connecting_target, protein_vector_point_alpha));
                    DBG << "targetPhi: " << targetPhi << "\t\tcurrentPhiTorsionAngle: " << currentPhiTorsionAngle << "\t\tphiDifference: " << phiDifference << "\t\t" << ring_oxygen_name << "-" << sugar_connection_atom << "-" << connected_atom << "-" << vector_point_alpha << std::endl;
                }


                if(enable_user_messages && !debug_output)
                {
                    std::cout << "Coordinates of " << connected_atom << ":\t" << protein_connecting_target.format() << std::endl;
                    std::cout << "Coordinates of " << vector_point_alpha << ":\t" << protein_vector_point_alpha.format() << std::endl;
                    std::cout << "Coordinates of " << vector_point_bravo << ":\t" << protein_vector_point_bravo.format() << std::endl;
                    std::cout << "Coordinates of " << sugar_connection_atom << ":\t" << sugar_connection_target.format() << std::endl;
                    std::cout << "Coordinates of " << ring_oxygen_name << ":\t" << ring_oxygen_coords.format() << std::endl;
                }

                if(debug_output)
                {
                    DBG << "Coordinates of " << connected_atom << ":\t" << protein_connecting_target.format() << std::endl;
                    DBG << "Coordinates of " << vector_point_alpha << ":\t" << protein_vector_point_alpha.format() << std::endl;
                    DBG << "Coordinates of " << vector_point_bravo << ":\t" << protein_vector_point_bravo.format() << std::endl;
                    DBG << "Coordinates of " << sugar_connection_atom << ":\t" << sugar_connection_target.format() << std::endl;
                    std::cout << "Coordinates of " << ring_oxygen_name << ":\t" << ring_oxygen_coords.format() << std::endl;
                }

                if(enable_user_messages && !debug_output)
                    std::cout << "Grafting input glycan with Chain ID of: " << chainID << " to " << residue_name << "-" << input_protein_side_chain_residue.id().trim() << std::endl;
                
                if(debug_output)
                    DBG << "Grafting input glycan with Chain ID of: " << chainID << " to " << residue_name << "-" << input_protein_side_chain_residue.id().trim() << std::endl;
                
                converted_mglycan.set_id(chainID);
                receiving_model.insert(converted_mglycan);

                if(enable_user_messages && !debug_output)
                    std::cout << "Glycan has been grafted!" << std::endl;
                
                if(debug_output)
                    DBG << "Glycan has been grafted!" << std::endl;
            }   
        }

        clipper::Mat33<clipper::ftype> Grafter::generate_rotation_matrix_from_rodrigues_rotation_formula(clipper::ftype degrees_to_rotate, clipper::Vec3<clipper::ftype> input_vector)
        {
            clipper::Vec3<clipper::ftype> input_unit_vector = input_vector.unit();
            clipper::ftype input_angle_in_radians = clipper::Util::d2rad(degrees_to_rotate);

            clipper::ftype omega_x = input_unit_vector[0];
            clipper::ftype omega_y = input_unit_vector[1];
            clipper::ftype omega_z = input_unit_vector[2];
            
            /*
            clipper::Mat33<clipper::ftype> matrix(1, 2, 3, 4, 5, 6, 7, 8, 9);
            |         1,         2,         3|
            |         4,         5,         6|
            |         7,         8,         9|
            */
            // clipper::ftype mat00 = ( cos(input_angle_in_radians) + pow(omega_x, 2) * (1 - cos(input_angle_in_radians)) );
            // clipper::ftype mat01 = ( omega_x * omega_y * (1 - cos(input_angle_in_radians)) - omega_z * sin(input_angle_in_radians) );
            // clipper::ftype mat02 = ( omega_y * sin(input_angle_in_radians) + omega_x * omega_z * (1 - cos(input_angle_in_radians)) );
            // clipper::ftype mat10 = ( omega_z * sin(input_angle_in_radians) + omega_x * omega_y * (1 - cos(input_angle_in_radians)) );
            // clipper::ftype mat11 = ( cos(input_angle_in_radians) + pow(omega_y, 2) * (1 - cos(input_angle_in_radians)) );
            // clipper::ftype mat12 = ( (-1 * omega_x) * sin(input_angle_in_radians) + omega_y * omega_z * (1 - cos(input_angle_in_radians)) );
            // clipper::ftype mat20 = ( (-1 * omega_y) * sin(input_angle_in_radians) + omega_x * omega_z * (1 - cos(input_angle_in_radians)) );
            // clipper::ftype mat21 = ( omega_x * sin(input_angle_in_radians ) + omega_y * omega_z * (1 - cos(input_angle_in_radians)) );
            // clipper::ftype mat22 = ( cos(input_angle_in_radians) + pow(omega_z, 2) * (1 - cos(input_angle_in_radians)) );

            clipper::ftype mat00 = ( pow(omega_x, 2) + (1 - pow(omega_x, 2)) * cos(input_angle_in_radians) );
            clipper::ftype mat01 = ( omega_x * omega_y * (1 - cos(input_angle_in_radians)) + omega_z * sin(input_angle_in_radians) );
            clipper::ftype mat02 = ( omega_z * omega_x * (1 - cos(input_angle_in_radians)) + omega_y * sin(input_angle_in_radians) );

            clipper::ftype mat10 = ( omega_x * omega_y * (1 - cos(input_angle_in_radians)) + omega_z * sin(input_angle_in_radians) );
            clipper::ftype mat11 = ( pow(omega_y, 2) + (1 - pow(omega_y, 2)) * cos(input_angle_in_radians) );
            clipper::ftype mat12 = ( omega_y * omega_z * (1 - cos(input_angle_in_radians)) - omega_x * sin(input_angle_in_radians) );

            clipper::ftype mat20 = ( omega_z * omega_x * (1 - cos(input_angle_in_radians)) - omega_y * sin(input_angle_in_radians) );
            clipper::ftype mat21 = ( omega_y * omega_z * (1 - cos(input_angle_in_radians)) + omega_x * sin(input_angle_in_radians) );
            clipper::ftype mat22 = ( pow(omega_z, 2) + (1 - pow(omega_z, 2)) * cos(input_angle_in_radians) );
            
            clipper::Mat33<clipper::ftype> rodrigures_rotation_matrix(mat00, mat01, mat02, mat10, mat11, mat12, mat20, mat21, mat22);

            return rodrigures_rotation_matrix;
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
