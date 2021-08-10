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
            { "ASN", "ND2", "CB", "CG", -100, 180 },
            { "ASP", "OD2", "CB", "CG" -100, 180 }
        };
        const int backbone_instructions_size = sizeof( backbone_instructions ) / sizeof( backbone_instructions[0] );

        const sugar_attachment sugar_instructions[] =
        {
            { "ideal", "C1", "C1", "O1" },   
            { "non-ideal", "C1", "C5", "O5" }
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

        void Grafter::graft_mpolymer_to_receiving_model(clipper::RTop_orth relocator, clipper::MPolymer input_chain)
        {    
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
            input_chain.transform(relocator);
            
            if(enable_user_messages && !debug_output)
                std::cout << "Grafting input glycan with Chain ID of: " << chainID << std::endl;
            
            if(debug_output)
                DBG << "Grafting input glycan with Chain ID of: " << chainID << std::endl;
            
            input_chain.set_id(chainID);
            receiving_model.insert(input_chain);
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
