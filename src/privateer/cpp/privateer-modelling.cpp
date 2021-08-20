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
    
        Grafter::Grafter(clipper::MiniMol receiving_model, clipper::MiniMol donor_model, bool trim_donor_when_clashes_detected, bool enable_user_messages, bool debug_output)
        {
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

            if(receiver_atom_index != -1)
            {
                residue_name = privateer::modelling::backbone_instructions[receiver_atom_index].residue_name;
                connected_atom = privateer::modelling::backbone_instructions[receiver_atom_index].connected_atom;
                vector_point_alpha = privateer::modelling::backbone_instructions[receiver_atom_index].vector_point_alpha;
                vector_point_bravo = privateer::modelling::backbone_instructions[receiver_atom_index].vector_point_bravo;
                targetPhi = privateer::modelling::backbone_instructions[receiver_atom_index].Phi;
                targetPsi = privateer::modelling::backbone_instructions[receiver_atom_index].Psi;

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

                if(enable_user_messages && !debug_output)
                    std::cout << "Grafting input glycan with Chain ID of: " << chainID << " to " << residue_name << "-" << input_protein_side_chain_residue.id().trim() << std::endl;
                
                if(debug_output)
                    DBG << "Grafting input glycan with Chain ID of: " << chainID << " to " << residue_name << "-" << input_protein_side_chain_residue.id().trim() << std::endl;
                
                converted_mglycan.set_id(chainID);
                export_model.insert(converted_mglycan);

                if(enable_user_messages && !debug_output)
                    std::cout << "Glycan has been grafted!" << std::endl;
                
                if(debug_output)
                    DBG << "Glycan has been grafted!" << std::endl;
            }
            else if(glycan_type == "ideal")
            {
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
                    
                    converted_mglycan.set_id(chainID);
                    export_model.insert(converted_mglycan);

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
                    
                    converted_mglycan.set_id(chainID);
                    export_model.insert(converted_mglycan);

                    if(enable_user_messages && !debug_output)
                        std::cout << "Glycan has been grafted!" << std::endl;
                    
                    if(debug_output)
                        DBG << "Glycan has been grafted!" << std::endl;
                }
            }
            
            std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > clashes_after_grafting = check_for_clashes(trim_donor_when_clashes_detected, export_model, converted_mglycan[0], input_protein_side_chain_residue, root_chain_id, chainID);
            
            // The above std::pair<>std::pair needs to go into above if statements.
            // Ideally we should keep grafting and checking for clashes until we exhaust all Phi, Psi possibilities. Store the best current performing grafted glycan in terms of clashes.

            // Only then try to 
            // if(delete_clashes)
                // trim_graft_until_no_clashes_left
            

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
                        output_residue.insert(currentAtom, -1);
                    }
                }
                output_residue.set_id(converted_mglycan[residue].id());
                output_residue.set_type(converted_mglycan[residue].type());
                output_residue.set_seqnum(converted_mglycan[residue].seqnum());
                output.insert(output_residue, -1);
            }
            output.set_id(converted_mglycan.id());
            return output;
        }

        std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > Grafter::check_for_clashes(bool trim_donor_when_clashes_detected, clipper::MiniMol& export_model, clipper::MMonomer& root_sugar, clipper::MMonomer& input_protein_side_chain_residue, clipper::String root_chain_id, clipper::String root_sugar_chain_id)
        {
            clipper::MAtomNonBond manb = clipper::MAtomNonBond(export_model, 1.0);
            clipper::MGlycology mgl = clipper::MGlycology(export_model, manb, debug_output, "undefined");
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
            clipper::MAtomNonBond clashmanb( export_model, 4.0 );
            for(int node = 0; node < grafted_mglycan.number_of_nodes(); node++)
            {
                clipper::MGlycan::Node currentNode = grafted_mglycan.get_node(node);
                clipper::MSugar currentSugar = currentNode.get_sugar();
                
                for(int atom = 0; atom < currentSugar.size(); atom++)
                {
                    clipper::Coord_orth currentSugarAtomCoord = currentSugar[atom].coord_orth();
                    std::vector<clipper::MAtomIndexSymmetry> neighbourhood_atoms = clashmanb( currentSugarAtomCoord, 3.0 );


                    for(int i = 0; i < neighbourhood_atoms.size(); i++)
                    {
                        int detected_chain      = neighbourhood_atoms[i].polymer();
                        int detected_monomer    = neighbourhood_atoms[i].monomer();
                        int detected_atom       = neighbourhood_atoms[i].atom();

                        auto search_result = std::find_if(std::begin(sugars_in_grafted_mglycan), std::end(sugars_in_grafted_mglycan),
                        [&export_model, &detected_chain, &detected_monomer](clipper::MSugar& sugar) {
                            return sugar.id().trim() == export_model[detected_chain][detected_monomer].id().trim() && sugar.type().trim() == export_model[detected_chain][detected_monomer].type().trim() &&
                                   sugar.seqnum() == export_model[detected_chain][detected_monomer].seqnum();
                        });

                        if(search_result == std::end(sugars_in_grafted_mglycan) && clipper::Coord_orth::length(currentSugarAtomCoord, export_model[detected_chain][detected_monomer][detected_atom].coord_orth()) < 3.0)
                        {
                            auto previously_identified = std::find_if(std::begin(clashing_residues), std::end(clashing_residues),
                            [&](std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> >& element) {
                                return  element.first.first.type().trim() == export_model[detected_chain][detected_monomer].type().trim() && element.first.first.id().trim() == export_model[detected_chain][detected_monomer].id().trim() && element.first.first.seqnum() == export_model[detected_chain][detected_monomer].seqnum() && 
                                        element.second.first.type().trim() == currentSugar.type().trim() && element.second.first.id().trim() == currentSugar.id().trim() && element.second.first.seqnum() == currentSugar.seqnum() &&
                                        element.first.second == root_chain_id && element.second.second == root_sugar_chain_id;
                            });

                            if(previously_identified == std::end(clashing_residues))
                            {
                                // Make sure that we are not detecting root protein side-chain as a clash.
                                if(root_chain_id.trim() != export_model[detected_chain].id().trim() || input_protein_side_chain_residue.id().trim() != export_model[detected_chain][detected_monomer].id().trim() || input_protein_side_chain_residue.type().trim() != export_model[detected_chain][detected_monomer].type().trim() || input_protein_side_chain_residue.seqnum() != export_model[detected_chain][detected_monomer].seqnum())
                                    clashing_residues.push_back( std::make_pair(std::make_pair(export_model[detected_chain][detected_monomer], root_chain_id), std::make_pair(currentSugar, root_sugar_chain_id)) );
                            }
                        }
                        
                        // DBG << "Detected " << export_model[detected_chain].id().trim() << "/" << export_model[detected_chain][detected_monomer].type().trim() << "-" << export_model[detected_chain][detected_monomer].id().trim()
                            // << " with the distance of: " << clipper::Coord_orth::length(currentSugarAtomCoord, export_model[detected_chain][detected_monomer][detected_atom].coord_orth()) << "\t\t\t" << currentSugar[atom].id().trim() << std::endl;
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
