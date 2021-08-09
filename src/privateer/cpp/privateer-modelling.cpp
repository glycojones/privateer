// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#include "privateer-modelling.h"



privateer::modelling::Grafter::Grafter(clipper::MiniMol receiving_model, clipper::MiniMol donor_model, bool debug_output)
{
    this->debug_output = debug_output;
    this->receiving_model = receiving_model;
    clipper::MGlycology donor_model_mgl = clipper::MGlycology(donor_model, debug_output, "undefined");
    

    std::vector<clipper::MGlycan> list_of_glycans_from_donor = donor_model_mgl.get_list_of_glycans();
    this->donor_glycans=list_of_glycans_from_donor;

    this->numDonorGlycansDetected=list_of_glycans_from_donor.size();
}

clipper::MPolymer privateer::modelling::Grafter::convert_mglycan_to_mpolymer(clipper::MGlycan input)
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

void privateer::modelling::Grafter::graft_mpolymer_to_receiving_model(clipper::Coord_orth target, clipper::Coord_orth source, std::vector<float> torsions, clipper::MPolymer input_chain)
{
    clipper::Mat33<clipper::ftype> Phi(1, 0, 0, 0, cos(clipper::Util::d2rad(torsions[0])), sin(clipper::Util::d2rad(torsions[0])), 0, -sin(clipper::Util::d2rad(torsions[0])), cos(clipper::Util::d2rad(torsions[0])));
    clipper::Mat33<clipper::ftype> Theta(cos(clipper::Util::d2rad(0)), 0, -sin(clipper::Util::d2rad(0)), 0, 1, 0, sin(clipper::Util::d2rad(0)), 0, cos(clipper::Util::d2rad(0)));
    clipper::Mat33<clipper::ftype> Psi(cos(clipper::Util::d2rad(torsions[1])), sin(clipper::Util::d2rad(torsions[1])), 0, -sin(clipper::Util::d2rad(torsions[1])), cos(clipper::Util::d2rad(torsions[1])), 0, 0, 0, 1);

    clipper::Mat33<clipper::ftype> rotation_matrix = Phi * Theta * Psi;
    clipper::Vec3<clipper::ftype> translation_vector(target.x() - source.x(), target.y() - source.y(), target.z() - source.z());

    clipper::RTop_orth relocator(rotation_matrix, translation_vector);

    const clipper::String chainids = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    std::vector<clipper::String> used_chain_ids;
    for(int i = 0; i < receiving_model.size(); i++)
    {
        clipper::String currentChainID = receiving_model[i].id().trim();
        used_chain_ids.push_back(currentChainID);
    }

    clipper::String new_chain_id = chainids[used_chain_ids.size()];
    input_chain.transform(relocator);
    input_chain.set_id(new_chain_id);
    receiving_model.insert(input_chain);
}

