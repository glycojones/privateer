// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#include "privateer-interactions.h"

#define DBG std::cout << "[" << __FUNCTION__ << "] - "

namespace privateer
{
    namespace interactions
    {

        void hydrogenate_input_model(std::string input_model, std::string output_path)
        {
            std::cout << "Hydrogenating '" << input_model << "' using project-gemmi/gemmi third party library." << std::endl;
            try
            {
                std::string output;
                if (output == "undefined")
                    output = "hydrogenated_input_model.pdb";
                else
                    output = output_path;
                gemmi::HydrogenChange h_change = gemmi::HydrogenChange::ReAdd;
                std::string monomer_dir = privateer::restraints::check_monlib_access();
                if (monomer_dir.empty())
                    throw std::runtime_error("Failed to locate $CLIBD_MON. Have ccp4 env variables been sourced?");

                gemmi::Structure st = gemmi::read_structure_file(input_model);

                if (st.models.empty() || st.models[0].chains.empty())
                    throw std::runtime_error("No atoms in the input file.");

                gemmi::setup_entities(st);
                size_t initial_h = 0;
                initial_h = gemmi::count_hydrogen_sites(st);

                std::vector<std::string> res_names = st.models[0].get_all_residue_names();

                std::printf("Reading %zu monomers and all links from %s\n",
                            res_names.size(), input_model.c_str());

                gemmi::MonLib monlib = gemmi::read_monomer_lib(monomer_dir, res_names,
                                                               gemmi::cif::read_file);

                for (size_t i = 0; i != st.models.size(); ++i)
                    gemmi::prepare_topology(st, monlib, i, h_change, false);

                std::printf("Hydrogen site count: %zu in input, %zu in output.\n",
                            initial_h, gemmi::count_hydrogen_sites(st));

                // Clean up the H atoms that were placed at (0.00, 0.00, 0.00)
                for (gemmi::Model &model : st.models)
                    for (gemmi::Chain &chain : model.chains)
                        for (gemmi::Residue &res : chain.residues)
                        {
                            std::vector<gemmi::Atom> modified_atoms;
                            for (gemmi::Atom &atom : res.atoms)
                            {
                                if (atom.is_hydrogen())
                                {
                                    if (atom.pos.x == 0 && atom.pos.y == 0 && atom.pos.z == 0 && atom.occ == 0)
                                        continue;
                                    else
                                        modified_atoms.push_back(atom);
                                }
                                else
                                {
                                    modified_atoms.push_back(atom);
                                }
                            }
                            res.atoms = modified_atoms;
                        }

                std::printf("Writing coordinates to %s\n", output.c_str());
                gemmi::Ofstream os(output, &std::cout);
                gemmi::write_pdb(st, os.ref());

                std::cout << input_model << " was successfully hydrogenated!" << std::endl;
            }
            catch (std::exception &e)
            {
                std::fprintf(stderr, "ERROR: %s\n", e.what());
                // throw stderr;
            }
        }

        privateer::interactions::CHPiBondsParser::CHPiBondsParser(std::string &input_model, std::string output_path, std::string algorithm)
        {
            this->algorithm = algorithm;
            hydrogenate_input_model(input_model, output_path);

            clipper::MMDBfile mfile;
            clipper::String clipperfied_input_model_path;
            if (output_path == "undefined")
                clipperfied_input_model_path = "hydrogenated_input_model.pdb";
            else
                clipperfied_input_model_path = output_path;
            privateer::util::read_coordinate_file_mtz(mfile, this->hydrogenated_input_model, clipperfied_input_model_path, true);

            this->manb_object = clipper::MAtomNonBond(this->hydrogenated_input_model, 5.0); // 1.2 for sugar initialization, 3.9 for max hbond distance.
            privateer::json::GlobalTorsionZScore torsions_zscore_database;
            this->hydrogenated_mglycology = clipper::MGlycology(this->hydrogenated_input_model, this->manb_object, torsions_zscore_database, false, "undefined");
        }

        std::vector<privateer::interactions::CHPiBond> privateer::interactions::CHPiBondsParser::get_CHPi_interactions(int glycanIndex)
        {
            std::vector<clipper::MGlycan> list_of_glycans = this->hydrogenated_mglycology.get_list_of_glycans();
            if (glycanIndex >= list_of_glycans.size() || glycanIndex < 0)
                throw std::runtime_error("Out of bounds access to std::vector storing MGlycans. Supplied index: " + std::to_string(glycanIndex) + "\tsize of vector: " + std::to_string(list_of_glycans.size()));

            std::vector<privateer::interactions::CHPiBond> result;
            clipper::MGlycan inputGlycan = list_of_glycans[glycanIndex];
            std::vector<clipper::MSugar> sugars_in_glycan = inputGlycan.get_sugars();

            for (int sugar = 0; sugar < sugars_in_glycan.size(); sugar++)
            {
                clipper::MSugar currentSugar = sugars_in_glycan[sugar];
                std::vector<privateer::interactions::CHPiBond> contacts = get_stacked_residues_python(currentSugar, sugar, sugars_in_glycan.size(), this->algorithm);

                if (!contacts.empty())
                {
                    // std::cout << "contacts size" << contacts.size() << std::endl;
                    for(int contact = 0; contact < contacts.size(); contact++)
                    {
                        result.push_back(contacts[contact]);
                    }
                }
            }
            return result;
        }

    clipper::Coord_orth get_ring_atoms (clipper::MMonomer mmon, std::vector<std::string> ring_atoms)
    {
    clipper::Coord_orth coords(0.0, 0.0, 0.0);
    int ring_size = ring_atoms.size();
    
        for (int atom = 0; atom < ring_atoms.size(); atom++)
        {
            coords += mmon.find(ring_atoms[atom], clipper::MM::ANY).coord_orth();
        }
        return clipper::Coord_orth(coords.x() / ring_size, coords.y() / ring_size, coords.z() / ring_size); // get coordinates for centre of π-ring
    }
    
    clipper::Coord_orth get_aromatic_centre(clipper::MMonomer mmon, std::string ring = "A")
        {
            clipper::Coord_orth result(0.0, 0.0, 0.0);

            if (mmon.type().trim() == "TRP")
            {
                if (ring == "A")
                { // smaller pyrrole ring
                    std::vector<std::string> trp_ringa_atoms = {"CE2", "NE1", "CD1", "CG", "CD2"};
                    result = get_ring_atoms(mmon, trp_ringa_atoms);
                }
                else
                {
                    std::vector<std::string> trp_ringb_atoms = {"CE2", "CZ2", "CH2", "CZ3", "CE3", "CD2"};
                    result = get_ring_atoms(mmon, trp_ringb_atoms);
                }
            }
            else if (mmon.type().trim() == "TYR" || mmon.type().trim() == "PHE")
            {
                std::vector<std::string> tyrphe_ring_atoms = {"CE2", "CZ", "CD1", "CE1", "CD2", "CG"};
                result = get_ring_atoms(mmon, tyrphe_ring_atoms);
            }
            else if (mmon.type().trim() == "HIS")
            {
                std::vector<std::string> his_ring_atoms = {"CE1", "ND1", "NE2", "CG", "CD2"};
                result = get_ring_atoms(mmon, his_ring_atoms);
            }
            return result;
        }
    }

        clipper::Vec3<clipper::ftype> get_plane(clipper::MMonomer mmon, std::string atom_1, std::string atom_2, std::string atom_3)
        {    
            clipper::Vec3<clipper::ftype> result(0.0, 0.0, 0.0);
            bool foundAtomAlpha = false, foundAtomBravo = false, foundAtomCharlie = false;
            clipper::Coord_orth coords_a(0.0, 0.0, 0.0), coords_b(0.0, 0.0, 0.0), coords_c(0.0, 0.0, 0.0); // Three atoms in π-ring used to create two vectors to find the normal to the π-ring

            for (int atom = 0; atom < mmon.size(); atom++)
            {
                clipper::MAtom currentAtom = mmon[atom];
                if (currentAtom.id().trim() == atom_1)
                    foundAtomAlpha = true;
                    coords_a += mmon.find(atom_1, clipper::MM::ANY).coord_orth(); 
                if (currentAtom.id().trim() == atom_2)
                    foundAtomBravo = true;
                    coords_b += mmon.find(atom_2, clipper::MM::ANY).coord_orth(); 
                if (currentAtom.id().trim() == atom_3)
                    foundAtomCharlie = true;
                    coords_c += mmon.find(atom_3 , clipper::MM::ANY).coord_orth(); 
            }

            if (foundAtomAlpha == false || foundAtomBravo == false || foundAtomCharlie == false)
                return result;
        
            clipper::Vec3<clipper::ftype> vec1(coords_a - coords_b);
            clipper::Vec3<clipper::ftype> vec2(coords_c - coords_b);
            result = clipper::Vec3<clipper::ftype>::cross(vec1, vec2);
            return result;
        }

    inline clipper::Vec3<clipper::ftype> find_aromatic_plane(clipper::MMonomer mmon)
    {
        clipper::Vec3<clipper::ftype> result(0.0, 0.0, 0.0);
        std::string atom_a, atom_b, atom_c = "";

        if (mmon.type().trim() == "TRP")
        {
            atom_a = "CD1", atom_b = "CD2", atom_c = "CE2";
            result = get_plane(mmon, atom_a, atom_b, atom_c);
        }
        else if (mmon.type().trim() == "TYR" || mmon.type().trim() == "PHE")
        {
            atom_a = "CE1", atom_b = "CG", atom_c = "CE2";
            result = get_plane(mmon, atom_a, atom_b, atom_c);
        }
        else if (mmon.type().trim() == "HIS")
        {
            atom_a = "CE1", atom_b = "CG", atom_c = "NE2";
            result = get_plane(mmon, atom_a, atom_b, atom_c);
        }
            return result.unit();
            return result;
    }

    clipper::ftype get_angle(clipper::Vec3<clipper::ftype> vec1, clipper::Vec3<clipper::ftype> vec2, std::string parameter)
    {
        clipper::ftype angle = acos(clipper::Vec3<clipper::ftype>::dot(vec1, vec2) /
                                    (sqrt(pow(vec1[0], 2) + pow(vec1[1], 2) + pow(vec1[2], 2)) *
                                     sqrt(pow(vec2[0], 2) + pow(vec2[1], 2) + pow(vec2[2], 2))));

        if(parameter == "theta")
            if(clipper::Vec3<clipper::ftype>::dot(vec1, vec2) < 0) // the angle is obtuse
            {
                angle = M_PI - angle;
            }
        return angle;
    }
    
    clipper::Vec3<clipper::ftype> get_sugar_ch_coords(std::pair<clipper::MAtom, clipper::MAtom> pair) // can probably get rid of this function and merge into code in get_stacked_residues
    {
        clipper::Vec3<clipper::ftype> vector(pair.second.coord_orth().x() - pair.first.coord_orth().x(),
                                             pair.second.coord_orth().y() - pair.first.coord_orth().y(),
                                             pair.second.coord_orth().z() - pair.first.coord_orth().z());
        return vector;
    }
    
    bool is_same_side(const clipper::Coord_orth& h_coords, const clipper::Coord_orth& anomeric_O_coords, const clipper::Vec3<clipper::ftype>& sugar_normal)
    { // is h-atom on the same side as anomeric oxygen (beta-face) or not (alpha-face)
        double dot_product = ((h_coords.x() - anomeric_O_coords.x()) * sugar_normal[0] +
                              (h_coords.y() - anomeric_O_coords.y()) * sugar_normal[1] +
                              (h_coords.z() - anomeric_O_coords.z()) * sugar_normal[2]) / 
                              std::sqrt(sugar_normal[0] * sugar_normal[0] +
                                        sugar_normal[1] * sugar_normal[1] +
                                        sugar_normal[2] * sugar_normal[2]);
        
        return dot_product >= 0; // Points on the same side have a non-negative dot product
    }

    std::string get_sugar_face(clipper::MSugar &input_sugar, std::pair<clipper::MAtom, clipper::MAtom> &ch_atoms)
    { // checks whether the h_atom is on the same side as the anomeric oxygen (beta) or the other side (alpha)
        std::string sugar_face = "α";
        const clipper::Vec3<clipper::ftype>& sugar_normal = input_sugar.ring_mean_plane(); // normal vector of the plane
        clipper::Coord_orth anomeric_O_coords = input_sugar.find("O5", clipper::MM::ANY).coord_orth(); // point on the plane
        clipper::Coord_orth h_coords = ch_atoms.second.coord_orth(); // point (h_atom coords) to check
        
        if (is_same_side(h_coords, anomeric_O_coords, sugar_normal)) // if true, h_atom is on same side as anomeric oxygen, then it is beta
        {
            sugar_face = "β";
        }
        return sugar_face;
    }

    std::vector<clipper::ftype> privateer::interactions::CHPiBond::calculate_hudson_parameters(const clipper::MAtomIndexSymmetry &neighbourhood,
                                                                                               const std::pair<clipper::MAtom, clipper::MAtom> &ch_atoms, 
                                                                                               const std::string &trp_ring,
                                                                                               const clipper::MiniMol &hydrogenated_input_model,
                                                                                               const clipper::MMonomer &mmon)
    {
        clipper::ftype distance = 0.0;
        std::vector<clipper::ftype> parameters;
        clipper::Coord_orth aromatic_centre = get_aromatic_centre(mmon, trp_ring);
        
        if (neighbourhood.symmetry() == 0)
        {
            distance = clipper::Coord_orth::length(ch_atoms.first.coord_orth(), aromatic_centre);
        }
        else // this neighbour is actually a symmetry mate
        {
            clipper::Spacegroup spgr = hydrogenated_input_model.spacegroup();
            clipper::Coord_frac f1 = ch_atoms.first.coord_orth().coord_frac(hydrogenated_input_model.cell());
            clipper::Coord_frac f2 = aromatic_centre.coord_frac(hydrogenated_input_model.cell());
            f1 = spgr.symop(neighbourhood.symmetry()) * f1;
            f1 = f1.lattice_copy_near(f2);
            distance = sqrt((f2 - f1).lengthsq(hydrogenated_input_model.cell()));
        }

        if (distance > 4.5) return {};
        clipper::ftype distance_ho = clipper::Coord_orth::length(ch_atoms.second.coord_orth(), aromatic_centre);      
        
        if (distance < distance_ho) return {}; // makes sure CH bond is pointing towards the ring - the C-π-ring distance must be longer than the H-π-ring distance
        clipper::Vec3<clipper::ftype> hx_vector = ch_atoms.second.coord_orth() - ch_atoms.first.coord_orth();
        clipper::Vec3<clipper::ftype> aromatic_vector = find_aromatic_plane(mmon);
        clipper::ftype theta = clipper::Util::rad2d((get_angle(hx_vector, aromatic_vector, "theta")));

        if (theta > 40) return {};
        clipper::Vec3<clipper::ftype> co_vector = aromatic_centre - ch_atoms.first.coord_orth();
        clipper::Vec3<clipper::ftype> theta_aromatic_vector = find_aromatic_plane(mmon);
        clipper::ftype angle1 = clipper::Util::rad2d((get_angle(theta_aromatic_vector, co_vector, "cp_distance")));
        clipper::ftype cp_distance = abs(cos(clipper::Util::d2rad(90) - clipper::Util::d2rad(angle1)) * distance);

        if (mmon.type().trim() == "HIS" || trp_ring == "A")
        {
            if (cp_distance > 1.6) return {}; // different cp restrictions for TrpA & His (<= 1.6)
            parameters.push_back(cp_distance);
        }
        else
        {
            if (cp_distance > 2.0) return {}; // different cp restrictions for TrpB & Tyr & Phe (<= 2.0)
            parameters.push_back(cp_distance);
        }
        parameters.push_back(distance);
        parameters.push_back(theta);
        return {parameters};
    }

    std::vector<clipper::ftype> privateer::interactions::CHPiBond::calculate_plevin_parameters(const clipper::MAtomIndexSymmetry &neighbourhood, 
                                                                                               const std::pair<clipper::MAtom, clipper::MAtom> &ch_atoms, 
                                                                                               const std::string &trp_ring,
                                                                                               const clipper::MiniMol &hydrogenated_input_model,
                                                                                               const clipper::MMonomer &mmon)
    {        
        clipper::ftype distance = 0.0;
        std::vector<clipper::ftype> parameters;
        clipper::Coord_orth aromatic_centre = get_aromatic_centre(mmon, trp_ring);
        
        if (neighbourhood.symmetry() == 0)
        {
            distance = clipper::Coord_orth::length(ch_atoms.first.coord_orth(), aromatic_centre);
        }
        else // this neighbour is actually a symmetry mate
        {
            clipper::Spacegroup spgr = hydrogenated_input_model.spacegroup();
            clipper::Coord_frac f1 = ch_atoms.first.coord_orth().coord_frac(hydrogenated_input_model.cell());
            clipper::Coord_frac f2 = aromatic_centre.coord_frac(hydrogenated_input_model.cell());
            f1 = spgr.symop(neighbourhood.symmetry()) * f1;
            f1 = f1.lattice_copy_near(f2);
            distance = sqrt((f2 - f1).lengthsq(hydrogenated_input_model.cell()));
        }

        if (distance >= 4.3) return {};
        clipper::ftype distance_ho = clipper::Coord_orth::length(ch_atoms.second.coord_orth(), aromatic_centre);             
        
        if (distance < distance_ho) return {};
        clipper::Vec3<clipper::ftype> ox_vector = aromatic_centre - ch_atoms.first.coord_orth(); // centre of ring-x atom
        
        clipper::Vec3<clipper::ftype> aromatic_vector = find_aromatic_plane(mmon);
        clipper::ftype theta = clipper::Util::rad2d(get_angle(ox_vector, aromatic_vector, "theta"));

        if (theta >= 25) return {};
        clipper::Vec3<clipper::ftype> hx_vector = ch_atoms.second.coord_orth() - ch_atoms.first.coord_orth(); // x atom-hydrogen
        clipper::Vec3<clipper::ftype> oh_vector = ch_atoms.second.coord_orth() - aromatic_centre; // centre of ring-hydrogen
        clipper::ftype phi = clipper::Util::rad2d(get_angle(oh_vector, hx_vector, "phi"));

        if (phi >= 120) return {};

        parameters.push_back(distance);
        parameters.push_back(theta);
        parameters.push_back(phi);

        return {parameters};
    }
    
    std::vector<privateer::interactions::CHPiBond> privateer::interactions::CHPiBondsParser::get_stacked_residues_python(   clipper::MSugar &input_sugar,
                                                                                                                            int sugarIndex,
                                                                                                                            int glycanSize,
                                                                                                                            std::string algorithm,
                                                                                                                            float distance,
                                                                                                                            float theta,
                                                                                                                            float phi,
                                                                                                                            std::string sugarFace) const
    {
        std::vector<std::pair<clipper::MAtom, clipper::MAtom>> ch_atoms;
        std::vector<clipper::Vec3<clipper::ftype>> c_to_h_vectors;
        clipper::Vec3<clipper::ftype> vector;
        clipper::MAtom ma;
        clipper::Coord_orth centre_apolar = input_sugar.ring_centre();
        std::vector<privateer::interactions::CHPiBond> results;
        clipper::MAtom h_atom;
        clipper::MAtom x_atom;
        std::vector<clipper::MAtom> bonded_atoms;

        for (int atom = 0; atom < input_sugar.size(); atom++)
        {            
            if (input_sugar[atom].element().trim() == "H") // find all h-atoms in sugar
            {
                h_atom = input_sugar[atom];
                const std::vector<clipper::MAtomIndexSymmetry> neighbourhood = this->manb_object.atoms_near(h_atom.coord_orth(), 2.0); // find neighbourhood atoms of h_atom

                for (int i = 0; i < neighbourhood.size(); i++)
                {   
                    clipper::MAtom sug_atom = this->hydrogenated_input_model[neighbourhood[i].polymer()][neighbourhood[i].monomer()][neighbourhood[i].atom()];
                    clipper::MMonomer sug_res = this->hydrogenated_input_model[neighbourhood[i].polymer()][neighbourhood[i].monomer()];
                    
                    if (input_sugar.id() == sug_res.id() && input_sugar.check_if_bonded(sug_atom, h_atom)) // check if the sug_atom and H-atom are in the same residue and covalently bonded
                    {    
                        if (sug_atom.element().trim() == "C" || 
                            sug_atom.element().trim() == "N" || 
                            sug_atom.element().trim() == "O" ||
                            sug_atom.element().trim() == "S" )
                        {
                            x_atom = sug_atom;
                            std::pair<clipper::MAtom, clipper::MAtom> pair(x_atom, h_atom);
                            vector = get_sugar_ch_coords(pair);
                            ch_atoms.push_back(pair);
                            c_to_h_vectors.push_back(vector);
                        }
                    }
                }
            }    
            continue;
        }

        const std::vector<clipper::MAtomIndexSymmetry> neighbourhood = this->manb_object.atoms_near(centre_apolar, 5.0);
        for (int k = 0; k < neighbourhood.size(); k++)
        {
            clipper::MMonomer mmon = this->hydrogenated_input_model[neighbourhood[k].polymer()][neighbourhood[k].monomer()];

            if ((mmon.type().trim() != "TRP") && // might be worth extending to cover GLU, ASP, GLN, ASN
                (mmon.type().trim() != "TYR") &&
                (mmon.type().trim() != "PHE") &&
                (mmon.type().trim() != "HIS"))
                continue;
            
            for (int j = 0; j < ch_atoms.size(); j++)  // Need to do this for each of the vectors in c_to_h_vectors
            {   
                if (algorithm == "hudson")
                { // Parameters: Theta(CH^normal), CX(C..ring centre), Cp(C..Cprojection)
                    if (mmon.type().trim() == "TRP")
                    {             
                        std::vector<clipper::ftype> parameters_trpA = privateer::interactions::CHPiBond::calculate_hudson_parameters(neighbourhood[k], ch_atoms[j], "A", this->hydrogenated_input_model, mmon);
                        if (!parameters_trpA.empty()) // need this if statement because return an empty function - could turn this into a guard funciton
                        {
                            // std::cout << ch_atoms[j].first.id() << std::endl;
                            privateer::interactions::CHPiBond the_interaction(input_sugar.chain_id(), this->hydrogenated_input_model[neighbourhood[k].polymer()].id(), input_sugar, mmon, theta, "hudson");
                            the_interaction.set_sugar_index(sugarIndex);
                            the_interaction.set_glycan_size(glycanSize);
                            the_interaction.set_distance_cp(parameters_trpA[0]);
                            the_interaction.set_angle_theta_h(parameters_trpA[2]);
                            the_interaction.set_distance_cx(parameters_trpA[1]);
                            the_interaction.set_trp_ring("A");
                            the_interaction.set_xh_pair(ch_atoms[j]);
                            the_interaction.set_sugar_face(get_sugar_face(input_sugar, ch_atoms[j]));

                            int residue;
                            for (residue = 0; residue < results.size(); residue++)
                            {
                                if(results[residue].get_stacked_residue().id() == this->hydrogenated_input_model[neighbourhood[k].polymer()][neighbourhood[k].monomer()].id() && 
                                   results[residue].get_stacked_residue().type().trim() == this->hydrogenated_input_model[neighbourhood[k].polymer()][neighbourhood[k].monomer()].type().trim() && 
                                   results[residue].get_stacked_residue_chainID() == this->hydrogenated_input_model[neighbourhood[k].polymer()].id() &&
                                   results[residue].get_xh_pair().first.id().trim() == ch_atoms[j].first.id().trim() &&
                                   results[residue].get_trp_ring() == "A")
                                break;
                            }
                            if (residue == results.size())
                            results.push_back(the_interaction);
                        }

                        std::vector<clipper::ftype> parameters_trpB = privateer::interactions::CHPiBond::calculate_hudson_parameters(neighbourhood[k], ch_atoms[j], "B", this->hydrogenated_input_model, mmon);                                                
                        if (!parameters_trpB.empty())
                        {
                            // std::cout << ch_atoms[j].first.id() << std::endl;
                            privateer::interactions::CHPiBond the_interaction(input_sugar.chain_id(), this->hydrogenated_input_model[neighbourhood[k].polymer()].id(), input_sugar, mmon, theta, "hudson");
                            the_interaction.set_sugar_index(sugarIndex);
                            the_interaction.set_sugar_face(sugarFace);
                            the_interaction.set_glycan_size(glycanSize);
                            the_interaction.set_distance_cp(parameters_trpB[0]);
                            the_interaction.set_angle_theta_h(parameters_trpB[2]);
                            the_interaction.set_distance_cx(parameters_trpB[1]);
                            the_interaction.set_trp_ring("B");
                            the_interaction.set_xh_pair(ch_atoms[j]);
                            the_interaction.set_sugar_face(get_sugar_face(input_sugar, ch_atoms[j]));

                            int residue;
                            for (residue = 0; residue < results.size(); residue++)
                            {
                                if(results[residue].get_stacked_residue().id() == this->hydrogenated_input_model[neighbourhood[k].polymer()][neighbourhood[k].monomer()].id() && 
                                   results[residue].get_stacked_residue().type().trim() == this->hydrogenated_input_model[neighbourhood[k].polymer()][neighbourhood[k].monomer()].type().trim() && 
                                   results[residue].get_stacked_residue_chainID() == this->hydrogenated_input_model[neighbourhood[k].polymer()].id() &&
                                   results[residue].get_xh_pair().first.id().trim() == ch_atoms[j].first.id().trim() &&
                                   results[residue].get_trp_ring() == "B")
                                break;
                            }
                            if (residue == results.size())
                            results.push_back(the_interaction);
                        }
                    }
                    else if (mmon.type().trim() == "TYR" || mmon.type().trim() == "PHE" || mmon.type().trim() == "HIS")
                    {
                        std::vector<clipper::ftype> parameters = privateer::interactions::CHPiBond::calculate_hudson_parameters(neighbourhood[k], ch_atoms[j], " ", this->hydrogenated_input_model, mmon);
                        if (!parameters.empty()) // only capturing TrpA interactions now so need to see why
                        {   
                            // std::cout << ch_atoms[j].first.id() << std::endl;
                            privateer::interactions::CHPiBond the_interaction(input_sugar.chain_id(), this->hydrogenated_input_model[neighbourhood[k].polymer()].id(), input_sugar, mmon, theta, "hudson");
                            the_interaction.set_sugar_index(sugarIndex);
                            the_interaction.set_glycan_size(glycanSize);
                            the_interaction.set_distance_cp(parameters[0]);
                            the_interaction.set_angle_theta_h(parameters[2]);
                            the_interaction.set_distance_cx(parameters[1]);
                            the_interaction.set_xh_pair(ch_atoms[j]);
                            the_interaction.set_sugar_face(get_sugar_face(input_sugar, ch_atoms[j]));
                            
                            int residue;
                            for (residue = 0; residue < results.size(); residue++)
                            {
                                if(results[residue].get_stacked_residue().id() == this->hydrogenated_input_model[neighbourhood[k].polymer()][neighbourhood[k].monomer()].id() && 
                                   results[residue].get_stacked_residue().type().trim() == this->hydrogenated_input_model[neighbourhood[k].polymer()][neighbourhood[k].monomer()].type().trim() && 
                                   results[residue].get_stacked_residue_chainID() == this->hydrogenated_input_model[neighbourhood[k].polymer()].id() &&
                                   results[residue].get_xh_pair().first.id().trim() == ch_atoms[j].first.id().trim())
                                break;
                            }
                            if (residue == results.size())
                            results.push_back(the_interaction);
                        }
                    }
                }
                else if (algorithm == "plevin")
                { 
                    if (mmon.type().trim() == "TRP")
                    {
                        std::vector<clipper::ftype> parameters_trpA = privateer::interactions::CHPiBond::calculate_plevin_parameters(neighbourhood[k], ch_atoms[j], "A", this->hydrogenated_input_model, mmon);
                        if (!parameters_trpA.empty()) // need this if statement because return an empty function - could turn this into a guard funciton
                        {
                            privateer::interactions::CHPiBond the_interaction(input_sugar.chain_id(), this->hydrogenated_input_model[neighbourhood[k].polymer()].id(), input_sugar, mmon, theta, "plevin");
                            the_interaction.set_sugar_index(sugarIndex);
                            the_interaction.set_glycan_size(glycanSize);
                            the_interaction.set_angle_phi(parameters_trpA[2]);
                            the_interaction.set_angle_theta_p(parameters_trpA[1]);
                            the_interaction.set_distance_cx(parameters_trpA[0]);
                            the_interaction.set_trp_ring("A");
                            the_interaction.set_xh_pair(ch_atoms[j]);
                            
                            int residue;
                                for (residue = 0; residue < results.size(); residue++)
                                {
                                    if(results[residue].get_stacked_residue().id() == this->hydrogenated_input_model[neighbourhood[k].polymer()][neighbourhood[k].monomer()].id() && 
                                        results[residue].get_stacked_residue().type().trim() == this->hydrogenated_input_model[neighbourhood[k].polymer()][neighbourhood[k].monomer()].type().trim() && 
                                        results[residue].get_stacked_residue_chainID() == this->hydrogenated_input_model[neighbourhood[k].polymer()].id() &&
                                        results[residue].get_xh_pair().first.id().trim() == ch_atoms[j].first.id().trim() &&
                                        results[residue].get_trp_ring() == "A")
                                        break;
                                }
                                if (residue == results.size())
                                results.push_back(the_interaction);
                        }

                        std::vector<clipper::ftype> parameters_trpB = privateer::interactions::CHPiBond::calculate_plevin_parameters(neighbourhood[k], ch_atoms[j], "B", this->hydrogenated_input_model, mmon);                                                
                        if (!parameters_trpB.empty())
                        {
                            privateer::interactions::CHPiBond the_interaction(input_sugar.chain_id(), this->hydrogenated_input_model[neighbourhood[k].polymer()].id(), input_sugar, mmon, theta, "plevin");
                            the_interaction.set_sugar_index(sugarIndex);
                            the_interaction.set_glycan_size(glycanSize);
                            the_interaction.set_angle_phi(parameters_trpB[2]);
                            the_interaction.set_angle_theta_p(parameters_trpB[1]);
                            the_interaction.set_distance_cx(parameters_trpB[0]);
                            the_interaction.set_trp_ring("B");
                            the_interaction.set_xh_pair(ch_atoms[j]);
        
                            int residue;
                                for (residue = 0; residue < results.size(); residue++)
                                {
                                    if(results[residue].get_stacked_residue().id() == this->hydrogenated_input_model[neighbourhood[k].polymer()][neighbourhood[k].monomer()].id() && 
                                        results[residue].get_stacked_residue().type().trim() == this->hydrogenated_input_model[neighbourhood[k].polymer()][neighbourhood[k].monomer()].type().trim() && 
                                        results[residue].get_stacked_residue_chainID() == this->hydrogenated_input_model[neighbourhood[k].polymer()].id() &&
                                        results[residue].get_xh_pair().first.id().trim() == ch_atoms[j].first.id().trim() &&
                                        results[residue].get_trp_ring() == "B")
                                        break;
                                }
                            if (residue == results.size())
                            results.push_back(the_interaction);
                        }
                    }
                    else if (mmon.type().trim() == "TYR" || mmon.type().trim() == "PHE" || mmon.type().trim() == "HIS")
                    {
                        std::vector<clipper::ftype> parameters = privateer::interactions::CHPiBond::calculate_plevin_parameters(neighbourhood[k], ch_atoms[j], " ", this->hydrogenated_input_model, mmon);
                        if (!parameters.empty()) // only capturing TrpA interactions now so need to see why
                        { 
                            privateer::interactions::CHPiBond the_interaction(input_sugar.chain_id(), this->hydrogenated_input_model[neighbourhood[k].polymer()].id(), input_sugar, mmon, theta, "plevin");
                            the_interaction.set_sugar_index(sugarIndex);
                            the_interaction.set_glycan_size(glycanSize);
                            the_interaction.set_angle_phi(parameters[2]);
                            the_interaction.set_angle_theta_p(parameters[1]);
                            the_interaction.set_distance_cx(parameters[0]);
                            the_interaction.set_xh_pair(ch_atoms[j]);

                            int residue;
                                for (residue = 0; residue < results.size(); residue++)
                                {
                                    if(results[residue].get_stacked_residue().id() == this->hydrogenated_input_model[neighbourhood[k].polymer()][neighbourhood[k].monomer()].id() && 
                                        results[residue].get_stacked_residue().type().trim() == this->hydrogenated_input_model[neighbourhood[k].polymer()][neighbourhood[k].monomer()].type().trim() && 
                                        results[residue].get_stacked_residue_chainID() == this->hydrogenated_input_model[neighbourhood[k].polymer()].id() &&
                                        results[residue].get_xh_pair().first.id().trim() == ch_atoms[j].first.id().trim())
                                        break;
                                }
                            if (residue == results.size())
                            results.push_back(the_interaction);
                        }
                    }
                }
            }
        }
        return results;
    }

    privateer::interactions::HBondsParser::HBondsParser(std::string &input_model, std::string output_path)
    {
        hydrogenate_input_model(input_model, output_path);

        clipper::MMDBfile mfile;
        clipper::String clipperfied_input_model_path;
        if (output_path == "undefined")
            clipperfied_input_model_path = "hydrogenated_input_model.pdb";
        else
            clipperfied_input_model_path = output_path;
        privateer::util::read_coordinate_file_mtz(mfile, this->hydrogenated_input_model, clipperfied_input_model_path, true);

        this->privateer::interactions::HBondsParser::import_ener_lib();
        this->hydrogenated_input_model = mark_hbond_donors_and_acceptors(this->hydrogenated_input_model);
        this->manb_object = clipper::MAtomNonBond(this->hydrogenated_input_model, 5.0); // 1.2 for sugar initialization, 3.9 for max hbond distance.
        privateer::json::GlobalTorsionZScore torsions_zscore_database;
        this->hydrogenated_mglycology = clipper::MGlycology(this->hydrogenated_input_model, this->manb_object, torsions_zscore_database, false, "undefined");
    }

    clipper::MiniMol privateer::interactions::HBondsParser::mark_hbond_donors_and_acceptors(clipper::MiniMol &input_model)
    {
        clipper::MiniMol output = input_model;
        std::string monomer_dir = privateer::restraints::check_monlib_access();
        if (monomer_dir.empty())
            throw std::runtime_error("Failed to locate $CLIBD_MON. Have ccp4 env variables been sourced?");

        for (int chain = 0; chain < output.size(); chain++)
        {
            for (int residue = 0; residue < output[chain].size(); residue++)
            {
                for (int atom = 0; atom < output[chain][residue].size(); atom++)
                {
                    hb_type h_bond_type = get_h_bond_type(output[chain][residue][atom], output[chain][residue].type().trim());
                    // std::cout << output[chain][residue].type().trim() << " - " << output[chain][residue][atom].id().trim() << "\t" << h_bond_type << std::endl;
                    output[chain][residue][atom].set_property("hb_type", clipper::Property<hb_type>(h_bond_type));
                }
            }
        }

        // for(int chain = 0; chain < output.size(); chain++)
        // {
        //     for(int residue = 0; residue < output[chain].size(); residue++)
        //     {
        //         for(int atom = 0; atom < output[chain][residue].size(); atom++)
        //         {
        //             if(output[chain][residue][atom].exists_property("hb_type"))
        //             {
        //                 const hb_type& atom_hb_potential = dynamic_cast<const clipper::Property<hb_type>& >(output[chain][residue][atom].get_property( "hb_type" )).value(); // kurwa biski ugly, bet kak pap tak
        //                 // std::cout   << "Chain " << output[chain].id().trim() << ": " << output[chain][residue].type().trim() << "/" << output[chain][residue].id().trim()
        //                 //             << " - " << output[chain][residue][atom].id().trim() << "\t\t\t" << atom_hb_potential << std::endl;
        //             }
        //         }
        //     }
        // }

        return output;
    }

    privateer::interactions::HBondsParser::hb_type privateer::interactions::HBondsParser::get_h_bond_type(clipper::MAtom &input_atom, std::string input_residue_type)
    {
        privateer::interactions::HBondsParser::monomer_dictionary dict = check_or_import_monomer_library_chem_comp_for_residue(input_residue_type);
        hb_type hb_t = HB_UNASSIGNED;

        std::string input_atom_name;
        char altconf = privateer::util::get_altconformation(input_atom);
        if (altconf != ' ')
        {
            clipper::String temp_string = input_atom.id().substr(0, 4);
            input_atom_name = temp_string.trim();
        }
        else
            input_atom_name = input_atom.id().trim();

        if (input_residue_type == dict.monomer_name)
        {
            std::vector<residue_monomer_library_chem_comp> dict_atoms = dict.dictionary_of_atoms;
            for (int i = 0; i < dict_atoms.size(); i++)
            {
                if (dict_atoms[i].atom_id == input_atom_name)
                {
                    if (dict_atoms[i].energy_type == "H")
                    {
                        if (is_connected_to_hydrogen_donor(dict_atoms[i].atom_id, dict_atoms))
                        {
                            hb_t = HB_HYDROGEN;
                        }
                    }
                    else
                    {
                        std::string current_atom_energy_type = dict_atoms[i].energy_type;
                        auto search_result_for_hb_type_of_current_atom = std::find_if(std::begin(this->energy_library), std::end(this->energy_library),
                                                                                      [&current_atom_energy_type](energy_library_entry &entry)
                                                                                      {
                                                                                          return entry.type == current_atom_energy_type;
                                                                                      });

                        if (search_result_for_hb_type_of_current_atom != std::end(energy_library))
                        {
                            std::string dict_hb_type = search_result_for_hb_type_of_current_atom->hb_type;
                            if (dict_hb_type == "N")
                                hb_t = HB_NEITHER;
                            else if (dict_hb_type == "D")
                                hb_t = HB_DONOR;
                            else if (dict_hb_type == "A")
                                hb_t = HB_ACCEPTOR;
                            else if (dict_hb_type == "B")
                                hb_t = HB_BOTH;
                            else if (dict_hb_type == "H")
                                hb_t = HB_HYDROGEN;
                        }
                    }
                }
            }
        }
        else
            throw std::runtime_error("get_h_bond_type: input_residue_type " + input_residue_type + " does not match dict.monomer_name " + dict.monomer_name);

        return hb_t;
    }

    bool privateer::interactions::HBondsParser::is_connected_to_hydrogen_donor(std::string atom_name, std::vector<residue_monomer_library_chem_comp> &residue_atoms)
    {
        bool result = false;

        auto search_result_for_input_hydrogen = std::find_if(std::begin(residue_atoms), std::end(residue_atoms),
                                                             [&atom_name](residue_monomer_library_chem_comp &entry)
                                                             {
                                                                 return entry.atom_id == atom_name;
                                                             });

        // std::cout << std::endl;
        // std::cout << "Received " << atom_name << std::endl;
        // for(int j = 0; j < residue_atoms.size(); j++)
        // {
        //     std::cout << "\t\t" << residue_atoms[j].atom_id << "\t"
        //                         << residue_atoms[j].element << "\t"
        //                         << residue_atoms[j].energy_type << "\t"
        //                         << residue_atoms[j].charge << "\t"
        //                         << residue_atoms[j].bonded_to_atom_id << "\t"
        //                         << residue_atoms[j].bond_type << std::endl;
        // }
        // std::cout << std::endl;

        if (search_result_for_input_hydrogen != std::end(residue_atoms))
        {
            residue_monomer_library_chem_comp input_hydrogen_row = *search_result_for_input_hydrogen;
            std::string current_hydrogen_is_bonded_to = input_hydrogen_row.bonded_to_atom_id;
            auto search_result_for_current_hydrogen_is_bonded_to = std::find_if(std::begin(residue_atoms), std::end(residue_atoms),
                                                                                [&current_hydrogen_is_bonded_to](residue_monomer_library_chem_comp &entry)
                                                                                {
                                                                                    return entry.atom_id == current_hydrogen_is_bonded_to;
                                                                                });
            // std::cout << "Bonded_to_atom_id " << current_hydrogen_is_bonded_to << std::endl;
            if (search_result_for_current_hydrogen_is_bonded_to != std::end(residue_atoms))
            {
                residue_monomer_library_chem_comp row_of_hydrogen_bonded_to = *search_result_for_current_hydrogen_is_bonded_to;
                std::string retrieved_energy_type = row_of_hydrogen_bonded_to.energy_type;
                // std::cout << "Received " << atom_name << " and energy is " << retrieved_energy_type << std::endl;
                auto search_result_for_energy_type_of_atom_bonded_to = std::find_if(std::begin(this->energy_library), std::end(this->energy_library),
                                                                                    [&retrieved_energy_type](energy_library_entry &entry)
                                                                                    {
                                                                                        return entry.type == retrieved_energy_type;
                                                                                    });

                if (search_result_for_energy_type_of_atom_bonded_to != std::end(energy_library))
                {
                    energy_library_entry row_of_hydrogen_bonded_to = *search_result_for_energy_type_of_atom_bonded_to;
                    std::string returned_hb_type = row_of_hydrogen_bonded_to.hb_type;
                    // std::cout << returned_hb_type << std::endl;
                    if (returned_hb_type == "D" || returned_hb_type == "B")
                    {
                        result = true;
                    }
                }
            }
        }
        return result;
    }
}

// Need to account for symmetry and distance here also.
std::vector<privateer::interactions::HBond> privateer::interactions::HBondsParser::get_HBonds_via_mcdonald_and_thornton(int glycanIndex, double max_dist)
{
    std::vector<clipper::MGlycan> list_of_hydrogenated_glycans = this->hydrogenated_mglycology.get_list_of_glycans();
    if (glycanIndex >= list_of_hydrogenated_glycans.size() || glycanIndex < 0)
        throw std::runtime_error("Out of bounds access to std::vector storing MGlycans. Supplied index: " + std::to_string(glycanIndex) + "\tsize of vector: " + std::to_string(list_of_hydrogenated_glycans.size()));

    std::vector<privateer::interactions::HBond> output;

    float min_dist = 0.1;
    clipper::MGlycan inputGlycan = list_of_hydrogenated_glycans[glycanIndex];
    // TO DO: This needs to be reworked after mmdb::Atom atom->GetChainID() equivalent is implemented.
    std::vector<clipper::MSugar> inputGlycanSugars = inputGlycan.get_sugars();
    for (int sugar = 0; sugar < inputGlycanSugars.size(); sugar++)
    {
        clipper::MSugar currentSugar = inputGlycanSugars[sugar];
        clipper::String glycanChainID = currentSugar.chain_id().trim();
        for (int atom = 0; atom < currentSugar.size(); atom++)
        {
            clipper::MAtom currentAtom = currentSugar[atom];
            std::vector<clipper::MAtomIndexSymmetry> neighbourhood = manb_object.atoms_near(currentAtom.coord_orth(), max_dist);
            // std::cout << currentAtom.id().trim() << std::endl;

            for (int i = 0; i < neighbourhood.size(); i++)
            {
                clipper::ftype distance = clipper::Coord_orth::length(currentAtom.coord_orth(), hydrogenated_input_model.atom(neighbourhood[i]).coord_orth());
                if (neighbourhood[i].symmetry() == 0)
                {
                    int detected_chain = neighbourhood[i].polymer();
                    int detected_monomer = neighbourhood[i].monomer();
                    int detected_atom = neighbourhood[i].atom();

                    // passing 'this' to lambda function is needed to capture class variable - 'clipper::MiniMol hydrogenated_input_model'
                    auto check_if_neighbour_is_in_currentMGlycan = std::find_if(std::begin(inputGlycanSugars), std::end(inputGlycanSugars),
                                                                                [this, &detected_chain, &detected_monomer, &detected_atom, &glycanChainID](clipper::MSugar &element)
                                                                                {
                                                                                    return element.id().trim() == hydrogenated_input_model[detected_chain][detected_monomer].id().trim() && element.type().trim() == hydrogenated_input_model[detected_chain][detected_monomer].type().trim() &&
                                                                                           element.seqnum() == hydrogenated_input_model[detected_chain][detected_monomer].seqnum() && glycanChainID == hydrogenated_input_model[detected_chain].id().trim() &&
                                                                                           element.size() == hydrogenated_input_model[detected_chain][detected_monomer].size();
                                                                                });

                    if (check_if_neighbour_is_in_currentMGlycan == std::end(inputGlycanSugars) && clipper::Coord_orth::length(currentAtom.coord_orth(), hydrogenated_input_model[detected_chain][detected_monomer][detected_atom].coord_orth()) <= max_dist)
                    {
                        if (currentAtom.exists_property("hb_type") && hydrogenated_input_model[detected_chain][detected_monomer][detected_atom].exists_property("hb_type"))
                        {
                            const hb_type &atom_hb_potential_target = dynamic_cast<const clipper::Property<hb_type> &>(currentAtom.get_property("hb_type")).value();
                            const hb_type &atom_hb_potential_neighbour = dynamic_cast<const clipper::Property<hb_type> &>(hydrogenated_input_model[detected_chain][detected_monomer][detected_atom].get_property("hb_type")).value();

                            if (atom_hb_potential_target == HB_HYDROGEN)
                            {
                                if (atom_hb_potential_neighbour == HB_ACCEPTOR || atom_hb_potential_neighbour == HB_BOTH)
                                {
                                    // Modify this, so that multiple neighbours are returned.
                                    std::vector<std::pair<clipper::MAtom, float>> hb_potential_target_neighbours = get_closest_neighbour_atoms(currentAtom, currentSugar, glycanChainID);
                                    std::vector<std::pair<clipper::MAtom, float>> hb_potential_neighbour_neighbours = get_closest_neighbour_atoms(hydrogenated_input_model[detected_chain][detected_monomer][detected_atom], hydrogenated_input_model[detected_chain][detected_monomer], hydrogenated_input_model[detected_chain].id().trim());
                                    std::pair<bool, HBond> new_h_bond = make_h_bond_from_sugar_hydrogen(glycanChainID, currentSugar, currentAtom, hydrogenated_input_model[detected_chain].id().trim(), hydrogenated_input_model[detected_chain][detected_monomer], hydrogenated_input_model[detected_chain][detected_monomer][detected_atom], hb_potential_target_neighbours, hb_potential_neighbour_neighbours);
                                    new_h_bond.second.sugarIndex = sugar;
                                    new_h_bond.second.glycanSize = inputGlycanSugars.size();
                                    if (new_h_bond.first)
                                        output.push_back(new_h_bond.second);

                                    // std::cout   << "DEBUG:: ===> pushing back new_h_bond from sugar_hydrogen, " << new_h_bond.second.donor_chainID << "/" << new_h_bond.second.donor_residue.type().trim() << "-" << new_h_bond.second.donor_residue.id().trim() << "_H: " << new_h_bond.second.HBonding_hydrogen.id().trim()
                                    //             << " donor: " << new_h_bond.second.donor.id().trim() << " acceptor: " << new_h_bond.second.acceptor.id().trim()
                                    //             << " acceptor_neighbour: " << new_h_bond.second.acceptor_neighbour.id().trim() << "_" << new_h_bond.second.acceptor_chainID << "/" << new_h_bond.second.acceptor_residue.type().trim() << "-" << new_h_bond.second.acceptor_residue.id().trim() << " angles:  " << new_h_bond.second.angle_1 << " " << new_h_bond.second.angle_2 << " " << new_h_bond.second.angle_3
                                    //             << " HBondLength: " << new_h_bond.second.HBondLength << std::boolalpha << " sugar_is_atom_donor: " << new_h_bond.second.sugar_atom_is_donor << " hydrogen_is_sugar_atom " << new_h_bond.second.hydrogen_is_sugar_atom
                                    //             << " bond_has_hydrogen_flag " << new_h_bond.second.bond_has_hydrogen_flag << "\tactually getting pushed in: " << new_h_bond.first << std::endl;
                                }
                            }
                            if (atom_hb_potential_target == HB_ACCEPTOR || atom_hb_potential_target == HB_BOTH)
                            {
                                if (atom_hb_potential_neighbour == HB_HYDROGEN || hydrogenated_input_model[detected_chain][detected_monomer].type().trim() == "HOH")
                                {
                                    std::vector<std::pair<clipper::MAtom, float>> hb_potential_target_neighbours = get_closest_neighbour_atoms(currentAtom, currentSugar, glycanChainID);
                                    std::vector<std::pair<clipper::MAtom, float>> hb_potential_neighbour_neighbours = get_closest_neighbour_atoms(hydrogenated_input_model[detected_chain][detected_monomer][detected_atom], hydrogenated_input_model[detected_chain][detected_monomer], hydrogenated_input_model[detected_chain].id().trim());
                                    std::pair<bool, HBond> new_h_bond = make_h_bond_from_environment_residue_hydrogen(hydrogenated_input_model[detected_chain].id().trim(), hydrogenated_input_model[detected_chain][detected_monomer], hydrogenated_input_model[detected_chain][detected_monomer][detected_atom], glycanChainID, currentSugar, currentAtom, hb_potential_neighbour_neighbours, hb_potential_target_neighbours);
                                    new_h_bond.second.sugarIndex = sugar;
                                    new_h_bond.second.glycanSize = inputGlycanSugars.size();
                                    if (new_h_bond.first)
                                        output.push_back(new_h_bond.second);

                                    // std::cout   << "DEBUG:: ===> pushing back new_h_bond from environment_residue, " << new_h_bond.second.donor_chainID << "/" << new_h_bond.second.donor_residue.type().trim() << "-" << new_h_bond.second.donor_residue.id().trim() << "_H: " << new_h_bond.second.HBonding_hydrogen.id().trim()
                                    //             << " donor: " << new_h_bond.second.donor.id().trim() << " acceptor: " << new_h_bond.second.acceptor.id().trim()
                                    //             << " acceptor_neighbour: " << new_h_bond.second.acceptor_neighbour.id().trim() << "_" << new_h_bond.second.acceptor_chainID << "/" << new_h_bond.second.acceptor_residue.type().trim() << "-" << new_h_bond.second.acceptor_residue.id().trim() << " angles:  " << new_h_bond.second.angle_1 << " " << new_h_bond.second.angle_2 << " " << new_h_bond.second.angle_3
                                    //             << " HBondLength: " << new_h_bond.second.HBondLength << std::boolalpha << " sugar_is_atom_donor: " << new_h_bond.second.sugar_atom_is_donor << " hydrogen_is_sugar_atom " << new_h_bond.second.hydrogen_is_sugar_atom
                                    //             << " bond_has_hydrogen_flag " << new_h_bond.second.bond_has_hydrogen_flag << "\tactually getting pushed in: " << new_h_bond.first << std::endl;
                                }
                            }
                        }
                        else
                        {
                            // std::cout << "hb_type: property does not exist" << std::endl;
                        }
                    }
                    else
                    {
                        // std::cout << hydrogenated_input_model[detected_chain][detected_monomer].id().trim() << "-" << hydrogenated_input_model[detected_chain][detected_monomer].type().trim() << hydrogenated_input_model[detected_chain][detected_monomer][detected_atom].id().trim() << " have been detected to be in the current MGlycan" << std::endl;
                    }
                }
            }
        }
    }
    return output;
}

std::vector<std::pair<clipper::MAtom, float>> privateer::interactions::HBondsParser::get_closest_neighbour_atoms(clipper::MAtom &input_atom, clipper::MMonomer &residue_atom_located_in, clipper::String chainID)
{

    std::vector<clipper::MAtomIndexSymmetry> neighbourhood = manb_object.atoms_near(input_atom.coord_orth(), 1.8);
    // Validate the length via clipper::Coord_orth::length

    std::vector<std::pair<clipper::MAtom, float>> atoms_in_current_residue_only;
    for (int i = 0; i < neighbourhood.size(); i++)
    {
        clipper::ftype distance = clipper::Coord_orth::length(input_atom.coord_orth(), hydrogenated_input_model.atom(neighbourhood[i]).coord_orth());
        if (neighbourhood[i].symmetry() == 0)
        {
            int detected_chain = neighbourhood[i].polymer();
            int detected_monomer = neighbourhood[i].monomer();
            int detected_atom = neighbourhood[i].atom();
            // std::cout   << "Input: " << input_atom.coord_orth().format() << " " << input_atom.id().trim() << "--" << residue_atom_located_in.type().trim() << "-" << residue_atom_located_in.id().trim() << "_" << residue_atom_located_in.seqnum()
            //             << " chain: " << chainID << "\t\tDetected: " << hydrogenated_input_model[detected_chain][detected_monomer][detected_atom].id().trim() << "--" << hydrogenated_input_model[detected_chain][detected_monomer][detected_atom].element() << " " << hydrogenated_input_model[detected_chain][detected_monomer][detected_atom].coord_orth().format()
            //             << " " << hydrogenated_input_model[detected_chain][detected_monomer].type().trim() << "-" << hydrogenated_input_model[detected_chain][detected_monomer].id().trim() << "_" << hydrogenated_input_model[detected_chain][detected_monomer].seqnum() << " chain: " << hydrogenated_input_model[detected_chain].id().trim() << std::endl;

            if (hydrogenated_input_model[detected_chain][detected_monomer][detected_atom].element().trim() != "H" &&
                hydrogenated_input_model[detected_chain][detected_monomer][detected_atom].id().trim() != input_atom.id().trim() &&
                hydrogenated_input_model[detected_chain][detected_monomer][detected_atom].coord_orth() != input_atom.coord_orth() &&
                hydrogenated_input_model[detected_chain][detected_monomer].type().trim() == residue_atom_located_in.type().trim() &&
                hydrogenated_input_model[detected_chain][detected_monomer].id().trim() == residue_atom_located_in.id().trim() &&
                hydrogenated_input_model[detected_chain][detected_monomer].seqnum() == residue_atom_located_in.seqnum() &&
                hydrogenated_input_model[detected_chain].id().trim() == chainID &&
                clipper::Coord_orth::length(input_atom.coord_orth(), hydrogenated_input_model[detected_chain][detected_monomer][detected_atom].coord_orth()) <= 1.8)
            {
                float dist = clipper::Coord_orth::length(input_atom.coord_orth(), hydrogenated_input_model[detected_chain][detected_monomer][detected_atom].coord_orth());
                atoms_in_current_residue_only.push_back(std::make_pair(hydrogenated_input_model[detected_chain][detected_monomer][detected_atom], dist));
            }
        }
    }

    std::sort(atoms_in_current_residue_only.begin(), atoms_in_current_residue_only.end(), [](std::pair<clipper::MAtom, float> &a, std::pair<clipper::MAtom, float> &b)
              { return a.second < b.second; });

    return atoms_in_current_residue_only;
}

std::pair<bool, privateer::interactions::HBond> privateer::interactions::HBondsParser::make_h_bond_from_sugar_hydrogen(std::string input_hydrogen_chainID, clipper::MMonomer &input_hydrogen_residue, clipper::MAtom &hydrogen, std::string input_acceptor_chainID, clipper::MMonomer &input_acceptor_residue, clipper::MAtom &acceptor, std::vector<std::pair<clipper::MAtom, float>> &hydrogen_neighbours, std::vector<std::pair<clipper::MAtom, float>> &acceptor_neighbours)
{
    privateer::interactions::HBond bond(input_hydrogen_chainID, input_hydrogen_residue, hydrogen, input_acceptor_chainID, input_acceptor_residue, acceptor, true);
    bond.HBondLength = clipper::Coord_orth::length(hydrogen.coord_orth(), acceptor.coord_orth());
    bool neighbour_distances_and_angles_are_good = true;
    bool good_donor_acceptor_dist = false;

    // Angle Donor-H-Acceptor
    std::pair<clipper::MAtom, float> hydrogen_neighbour = hydrogen_neighbours.front();
    double D_H_A_angle = clipper::Util::rad2d(clipper::Coord_orth::angle(hydrogen_neighbour.first.coord_orth(), hydrogen.coord_orth(), acceptor.coord_orth()));
    double D_A_dist = clipper::Coord_orth::length(hydrogen_neighbour.first.coord_orth(), acceptor.coord_orth());

    if (D_A_dist < 3.9)
        good_donor_acceptor_dist = true;

    if (bond.donor.is_null())
    {
        bond.donor = hydrogen_neighbour.first;
        bond.angle_1 = D_H_A_angle;
    }
    if (D_H_A_angle < 90)
        neighbour_distances_and_angles_are_good = false;

    // Angle H-Acceptor-AcceptorNeighbour
    for (int i = 0; i < acceptor_neighbours.size(); i++)
    {
        double H_A_AA_angle = clipper::Util::rad2d(clipper::Coord_orth::angle(hydrogen.coord_orth(), acceptor.coord_orth(), acceptor_neighbours[i].first.coord_orth()));

        bond.angle_2 = H_A_AA_angle;
        if (H_A_AA_angle < 90)
            neighbour_distances_and_angles_are_good = false;

        // Angle Donor-Acceptor-AcceptorNeighbour
        double D_A_AA_angle = clipper::Util::rad2d(clipper::Coord_orth::angle(hydrogen_neighbour.first.coord_orth(), acceptor.coord_orth(), acceptor_neighbours[i].first.coord_orth()));

        bond.acceptor_neighbour = acceptor_neighbours[i].first;
        bond.angle_3 = D_A_AA_angle;

        if (D_A_AA_angle < 90)
            neighbour_distances_and_angles_are_good = false;

        if (!neighbour_distances_and_angles_are_good)
            break;
    }

    if (bond.angle_2 < 90 || bond.angle_3 < 90)
        neighbour_distances_and_angles_are_good = false;

    return std::pair<bool, privateer::interactions::HBond>(neighbour_distances_and_angles_are_good && good_donor_acceptor_dist, bond);
}

// std::cout << "from_environment: " << hydrogenated_input_model[detected_chain][detected_monomer].type().trim() << "_" << hb_potential_neighbour_neighbour.first.id().trim() << "-(" << hydrogenated_input_model[detected_chain][detected_monomer][detected_atom].id().trim() << ")--(" << currentAtom.id().trim() << ")-" << hb_potential_target_neighbour.first.id().trim() << "_" << currentSugar.type().trim() << std::endl;
// std::pair<bool, HBond> new_h_bond = make_h_bond_from_environment_residue_hydrogen(currentAtom, hydrogenated_input_model[detected_chain][detected_monomer][detected_atom], "HOH", hb_potential_target_neighbour.first, hb_potential_neighbour_neighbour.first);

std::pair<bool, privateer::interactions::HBond> privateer::interactions::HBondsParser::make_h_bond_from_environment_residue_hydrogen(std::string input_hydrogen_chainID, clipper::MMonomer &input_hydrogen_residue, clipper::MAtom &hydrogen, std::string input_acceptor_on_sugar_chainID, clipper::MMonomer &input_acceptor_on_sugar_residue, clipper::MAtom &acceptor_on_sugar, std::vector<std::pair<clipper::MAtom, float>> &hydrogen_neighbours, std::vector<std::pair<clipper::MAtom, float>> &acceptor_on_sugar_neighbours)
{
    double water_dist_max = 3.25;
    bool ligand_atom_is_H_flag = false;

    // HBond(std::string input_hydrogen_chainID, clipper::MMonomer& input_hydrogen_residue, clipper::MAtom& input_hydrogen, std::string input_acceptor_chainID, clipper::MMonomer& input_acceptor_residue, clipper::MAtom& input_acceptor, bool input_sugar_atom_is_H_flag)
    privateer::interactions::HBond bond(input_hydrogen_chainID, input_hydrogen_residue, hydrogen, input_acceptor_on_sugar_chainID, input_acceptor_on_sugar_residue, acceptor_on_sugar, ligand_atom_is_H_flag);
    bond.HBondLength = clipper::Coord_orth::length(acceptor_on_sugar.coord_orth(), hydrogen.coord_orth());
    bool neighbour_distances_and_angles_are_good = true;
    bool good_donor_acceptor_dist = false;

    // Donor-Acceptor
    std::pair<clipper::MAtom, float> hydrogen_neighbour = std::make_pair(clipper::MAtom().null(), -1.0);
    if (!hydrogen_neighbours.empty())
        hydrogen_neighbour = hydrogen_neighbours.front();
    double D_A_dist = clipper::Coord_orth::length(hydrogen_neighbour.first.coord_orth(), acceptor_on_sugar.coord_orth());

    if (D_A_dist < 3.9)
        good_donor_acceptor_dist = true;

    if (input_hydrogen_residue.type().trim() == "HOH" && hydrogen.element() == "O")
    {
        if (D_A_dist < water_dist_max)
        {
            good_donor_acceptor_dist = true;
            bond.donor = hydrogen;
        }
    }

    // Donor-Hydrogen-Acceptor
    double D_H_A_angle = clipper::Util::rad2d(clipper::Coord_orth::angle(hydrogen_neighbour.first.coord_orth(), hydrogen.coord_orth(), acceptor_on_sugar.coord_orth()));
    if (D_H_A_angle < 90)
        neighbour_distances_and_angles_are_good = false;

    bond.donor = hydrogen_neighbour.first;
    bond.angle_1 = D_H_A_angle;

    // Angle Hydrogen-Acceptor-AcceptorNeighbour
    for (int i = 0; i < acceptor_on_sugar_neighbours.size(); i++)
    {
        double H_A_AA_angle = clipper::Util::rad2d(clipper::Coord_orth::angle(hydrogen.coord_orth(), acceptor_on_sugar.coord_orth(), acceptor_on_sugar_neighbours[i].first.coord_orth()));
        if (H_A_AA_angle < 90)
            neighbour_distances_and_angles_are_good = false;

        bond.acceptor = acceptor_on_sugar;
        bond.angle_2 = H_A_AA_angle;

        // Angle Donor-Acceptor-AcceptorNeighbour
        if (input_hydrogen_residue.type().trim() != "HOH")
        {
            double D_A_AA_angle = clipper::Util::rad2d(clipper::Coord_orth::angle(hydrogen_neighbour.first.coord_orth(), acceptor_on_sugar.coord_orth(), acceptor_on_sugar_neighbours[i].first.coord_orth()));
            if (D_A_AA_angle < 90)
                neighbour_distances_and_angles_are_good = false;

            bond.acceptor_neighbour = acceptor_on_sugar_neighbours[i].first;
            bond.angle_3 = D_A_AA_angle;
        }
        else
        {
            double D_A_AA_angle = clipper::Util::rad2d(clipper::Coord_orth::angle(hydrogen.coord_orth(), acceptor_on_sugar.coord_orth(), acceptor_on_sugar_neighbours[i].first.coord_orth()));
            bond.acceptor_neighbour = acceptor_on_sugar_neighbours[i].first;
            bond.angle_3 = D_A_AA_angle;
        }

        if (!neighbour_distances_and_angles_are_good)
            break;
    }

    if (bond.angle_2 < 90 || bond.angle_3 < 90)
        neighbour_distances_and_angles_are_good = false;

    return std::pair<bool, privateer::interactions::HBond>(neighbour_distances_and_angles_are_good && good_donor_acceptor_dist, bond);
}

// Private methods

void privateer::interactions::HBondsParser::import_ener_lib()
{
    std::string monomer_dir = privateer::restraints::check_monlib_access();
    if (monomer_dir.empty())
        throw std::runtime_error("Failed to locate $CLIBD_MON. Have ccp4 env variables been sourced?");

    std::stringstream path;

    int expected_ncolumns_in_energy_library_table = 9;

    if (!monomer_dir.empty())
    {
        path << monomer_dir << "ener_lib"
             << ".cif";
        std::cout << "Path to 'ener_lib.cif': " << path.str() << std::endl;
        gemmi::cif::Document ener_lib_gemmi_document = gemmi::cif::read_file(path.str());
        for (gemmi::cif::Block &block : ener_lib_gemmi_document.blocks)
        {
            if (!block.name.empty() && block.name == "energy")
            {
                gemmi::cif::Table energy_library_table = block.find("_lib_atom.",
                                                                    {"type", "weight", "hb_type", "vdw_radius", "vdwh_radius", "ion_radius", "element", "valency", "sp"});

                if (energy_library_table.width() != expected_ncolumns_in_energy_library_table)
                    throw std::runtime_error("Number of expected columns doesn't match the number of columns retrieved from 'ener_lib.cif'. Please report this issue to the developers of Privateer.");

                for (size_t i = 0; i != energy_library_table.length(); ++i)
                {

                    gemmi::cif::Table::Row currentGemmiRow = energy_library_table[i];
                    energy_library_entry export_row;
                    export_row.type = currentGemmiRow[0];
                    export_row.weight = currentGemmiRow[1];
                    export_row.hb_type = currentGemmiRow[2];
                    export_row.vdw_radius = currentGemmiRow[3];
                    export_row.vdwh_radius = currentGemmiRow[4];
                    export_row.ion_radius = currentGemmiRow[5];
                    export_row.element = currentGemmiRow[6];
                    export_row.valency = currentGemmiRow[7];
                    export_row.sp = currentGemmiRow[8];

                    if (export_row.type == "." || export_row.type == "n/a" || export_row.type.empty())
                        export_row.type = "null";
                    if (export_row.weight == "." || export_row.weight == "n/a" || export_row.weight.empty())
                        export_row.weight = "null";
                    if (export_row.hb_type == "." || export_row.hb_type == "n/a" || export_row.hb_type.empty())
                        export_row.hb_type = "null";
                    if (export_row.vdw_radius == "." || export_row.vdw_radius == "n/a" || export_row.vdw_radius.empty())
                        export_row.vdw_radius = "null";
                    if (export_row.vdwh_radius == "." || export_row.vdwh_radius == "n/a" || export_row.vdwh_radius.empty())
                        export_row.vdwh_radius = "null";
                    if (export_row.ion_radius == "." || export_row.ion_radius == "n/a" || export_row.ion_radius.empty())
                        export_row.ion_radius = "null";
                    if (export_row.element == "." || export_row.element == "n/a" || export_row.element.empty())
                        export_row.element = "null";
                    if (export_row.valency == "." || export_row.valency == "n/a" || export_row.valency.empty())
                        export_row.valency = "null";
                    if (export_row.sp == "." || export_row.sp == "n/a" || export_row.sp.empty())
                        export_row.sp = "null";

                    this->energy_library.push_back(export_row);
                }
            }
        }
    }
    else
        throw std::runtime_error("Failed to locate $CLIBD_MON. Have ccp4 env variables been sourced?");

    // std::cout << std::endl;
    // for(int i = 0; i < energy_library.size(); i++)
    // {
    //     std::cout   << i << "/" << energy_library.size() << ": "
    //                 << energy_library[i].type << "\t"
    //                 << energy_library[i].weight << "\t"
    //                 << energy_library[i].hb_type << "\t"
    //                 << energy_library[i].vdw_radius << "\t"
    //                 << energy_library[i].vdwh_radius << "\t"
    //                 << energy_library[i].ion_radius << "\t"
    //                 << energy_library[i].element << "\t"
    //                 << energy_library[i].valency << "\t"
    //                 << energy_library[i].sp << std::endl;
    // }
    // std::cout << std::endl;
}

privateer::interactions::HBondsParser::monomer_dictionary privateer::interactions::HBondsParser::check_or_import_monomer_library_chem_comp_for_residue(std::string input_residue_type)
{

    auto search_result = std::find_if(std::begin(monomer_dict), std::end(monomer_dict),
                                      [&input_residue_type](monomer_dictionary &entry)
                                      {
                                          return entry.monomer_name == input_residue_type;
                                      });

    if (search_result != std::end(monomer_dict))
    {
        return *search_result;
    }
    else
    {
        std::string monomer_dir = privateer::restraints::check_monlib_access();
        if (monomer_dir.empty())
            throw std::runtime_error("Failed to locate $CLIBD_MON. Have ccp4 env variables been sourced?");

        std::stringstream path;
        std::locale loc;

        char initial = std::tolower(input_residue_type[0], loc);

        if (!monomer_dir.empty())
        {
            path << monomer_dir << initial << "/" << input_residue_type << ".cif";
            // std::cout << path.str() << std::endl;
            gemmi::cif::Document monomer_dict_gemmi_document = gemmi::cif::read_file(path.str());
            std::string block_name_search_string = "comp_" + input_residue_type;
            std::vector<residue_monomer_library_chem_comp> dict_of_atoms;
            monomer_dictionary new_entry_in_monomer_dict;
            new_entry_in_monomer_dict.monomer_name = input_residue_type;
            std::vector<residue_monomer_library_chem_comp> &export_vector = new_entry_in_monomer_dict.dictionary_of_atoms;
            for (gemmi::cif::Block &block : monomer_dict_gemmi_document.blocks)
            {
                if (!block.name.empty() && block.name == block_name_search_string)
                {
                    gemmi::cif::Table chem_comp_atom_table = block.find("_chem_comp_atom.",
                                                                        {"comp_id", "atom_id", "type_symbol", "type_energy", "?charge", "?partial_charge"});

                    gemmi::cif::Table chem_comp_bond_table = block.find("_chem_comp_bond.",
                                                                        {"comp_id", "atom_id_1", "atom_id_2", "type"});

                    for (size_t i = 0; i != chem_comp_atom_table.length(); ++i)
                    {
                        gemmi::cif::Table::Row currentAtomTableGemmiRow = chem_comp_atom_table[i];
                        residue_monomer_library_chem_comp export_row;

                        if (chem_comp_atom_table.has_column(0))
                            export_row.residue_type = currentAtomTableGemmiRow[0];
                        if (chem_comp_atom_table.has_column(1))
                            export_row.atom_id = currentAtomTableGemmiRow[1];
                        if (chem_comp_atom_table.has_column(2))
                            export_row.element = currentAtomTableGemmiRow[2];
                        if (chem_comp_atom_table.has_column(3))
                            export_row.energy_type = currentAtomTableGemmiRow[3];
                        if (chem_comp_atom_table.has_column(4))
                            export_row.charge = currentAtomTableGemmiRow[4];
                        if (chem_comp_atom_table.has_column(5))
                            export_row.charge = currentAtomTableGemmiRow[5];

                        for (size_t j = 0; j != chem_comp_bond_table.length(); ++j)
                        {
                            gemmi::cif::Table::Row currentBondTableGemmiRow = chem_comp_bond_table[j];
                            if (currentBondTableGemmiRow[1] == export_row.atom_id)
                            {
                                auto bond_table_search_result = std::find_if(std::begin(dict_of_atoms), std::end(dict_of_atoms),
                                                                             [&currentBondTableGemmiRow](residue_monomer_library_chem_comp &element)
                                                                             {
                                                                                 return element.atom_id == currentBondTableGemmiRow[1];
                                                                             });

                                if (bond_table_search_result != std::end(dict_of_atoms))
                                {
                                    residue_monomer_library_chem_comp new_export_row = export_row;
                                    if (chem_comp_bond_table.has_column(2))
                                        new_export_row.bonded_to_atom_id = currentBondTableGemmiRow[2];
                                    if (chem_comp_bond_table.has_column(3))
                                        new_export_row.bond_type = currentBondTableGemmiRow[3];

                                    if (new_export_row.residue_type == "." || new_export_row.residue_type == "n/a" || new_export_row.residue_type.empty())
                                        new_export_row.residue_type = "null";
                                    if (new_export_row.atom_id == "." || new_export_row.atom_id == "n/a" || new_export_row.atom_id.empty())
                                        new_export_row.atom_id = "null";
                                    if (new_export_row.element == "." || new_export_row.element == "n/a" || new_export_row.element.empty())
                                        new_export_row.element = "null";
                                    if (new_export_row.energy_type == "." || new_export_row.energy_type == "n/a" || new_export_row.energy_type.empty())
                                        new_export_row.energy_type = "null";
                                    if (new_export_row.charge == "." || new_export_row.charge == "n/a" || new_export_row.charge.empty())
                                        new_export_row.charge = "null";
                                    if (new_export_row.bonded_to_atom_id == "." || new_export_row.bonded_to_atom_id == "n/a" || new_export_row.bonded_to_atom_id.empty())
                                        new_export_row.bonded_to_atom_id = "null";
                                    if (new_export_row.bond_type == "." || new_export_row.bond_type == "n/a" || new_export_row.bond_type.empty())
                                        new_export_row.bond_type = "null";

                                    dict_of_atoms.push_back(new_export_row);
                                }
                                else
                                {
                                    if (chem_comp_bond_table.has_column(2))
                                        export_row.bonded_to_atom_id = currentBondTableGemmiRow[2];
                                    if (chem_comp_bond_table.has_column(3))
                                        export_row.bond_type = currentBondTableGemmiRow[3];

                                    if (export_row.residue_type == "." || export_row.residue_type == "n/a" || export_row.residue_type.empty())
                                        export_row.residue_type = "null";
                                    if (export_row.atom_id == "." || export_row.atom_id == "n/a" || export_row.atom_id.empty())
                                        export_row.atom_id = "null";
                                    if (export_row.element == "." || export_row.element == "n/a" || export_row.element.empty())
                                        export_row.element = "null";
                                    if (export_row.energy_type == "." || export_row.energy_type == "n/a" || export_row.energy_type.empty())
                                        export_row.energy_type = "null";
                                    if (export_row.charge == "." || export_row.charge == "n/a" || export_row.charge.empty())
                                        export_row.charge = "null";
                                    if (export_row.bonded_to_atom_id == "." || export_row.bonded_to_atom_id == "n/a" || export_row.bonded_to_atom_id.empty())
                                        export_row.bonded_to_atom_id = "null";
                                    if (export_row.bond_type == "." || export_row.bond_type == "n/a" || export_row.bond_type.empty())
                                        export_row.bond_type = "null";

                                    dict_of_atoms.push_back(export_row);
                                }
                            }
                        }

                        auto second_bond_table_search_result = std::find_if(std::begin(dict_of_atoms), std::end(dict_of_atoms),
                                                                            [&export_row](residue_monomer_library_chem_comp &element)
                                                                            {
                                                                                return element.atom_id == export_row.atom_id;
                                                                            });

                        if (second_bond_table_search_result == std::end(dict_of_atoms))
                        {
                            if (export_row.residue_type == "." || export_row.residue_type == "n/a" || export_row.residue_type.empty())
                                export_row.residue_type = "null";
                            if (export_row.atom_id == "." || export_row.atom_id == "n/a" || export_row.atom_id.empty())
                                export_row.atom_id = "null";
                            if (export_row.element == "." || export_row.element == "n/a" || export_row.element.empty())
                                export_row.element = "null";
                            if (export_row.energy_type == "." || export_row.energy_type == "n/a" || export_row.energy_type.empty())
                                export_row.energy_type = "null";
                            if (export_row.charge == "." || export_row.charge == "n/a" || export_row.charge.empty())
                                export_row.charge = "null";
                            if (export_row.bonded_to_atom_id == "." || export_row.bonded_to_atom_id == "n/a" || export_row.bonded_to_atom_id.empty())
                                export_row.bonded_to_atom_id = "null";
                            if (export_row.bond_type == "." || export_row.bond_type == "n/a" || export_row.bond_type.empty())
                                export_row.bond_type = "null";

                            dict_of_atoms.push_back(export_row);
                        }
                    }
                }
            }

            for (int i = 0; i < dict_of_atoms.size(); i++)
            {
                if (dict_of_atoms[i].element == "H" && dict_of_atoms[i].bonded_to_atom_id == "null")
                {
                    for (int j = 0; j < dict_of_atoms.size(); j++)
                    {
                        if (dict_of_atoms[i].atom_id == dict_of_atoms[j].bonded_to_atom_id)
                            dict_of_atoms[i].bonded_to_atom_id = dict_of_atoms[j].atom_id;
                        dict_of_atoms[i].bond_type = "single";
                    }
                }
            }

            export_vector.insert(export_vector.end(), std::make_move_iterator(dict_of_atoms.begin()), std::make_move_iterator(dict_of_atoms.end()));
            this->monomer_dict.push_back(new_entry_in_monomer_dict);

            // std::cout << std::endl;
            // for(int i = 0; i < this->monomer_dict.size(); i++)
            // {
            //     std::cout << i << "/" << monomer_dict.size() << ": " << monomer_dict[i].monomer_name << std::endl;
            //     std::vector<residue_monomer_library_chem_comp> dictionary = monomer_dict[i].dictionary_of_atoms;
            //     for(int j = 0; j < dictionary.size(); j++)
            //     {
            //         std::cout << "\t\t" << dictionary[j].atom_id << "\t"
            //                             << dictionary[j].element << "\t"
            //                             << dictionary[j].energy_type << "\t"
            //                             << dictionary[j].charge << "\t"
            //                             << dictionary[j].bonded_to_atom_id << "\t"
            //                             << dictionary[j].bond_type << std::endl;
            //     }
            // }
            // std::cout << std::endl;

            return new_entry_in_monomer_dict;
        }
        else
            throw std::runtime_error("Failed to locate $CLIBD_MON. Have ccp4 env variables been sourced?");
    }
}
