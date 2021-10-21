// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#include "privateer-interactions.h"

#define DBG std::cout << "[" << __FUNCTION__ << "] - "

namespace privateer {
	namespace interactions {
        
        void hydrogenate_input_model(std::string input_model)
        {
            std::cout << "Hydrogenating '" << input_model << "' using project-gemmi/gemmi third party library." << std::endl;
            try 
            {
                std::string output = "hydrogenated_input_model.pdb";
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
                                                        gemmi::read_cif_gz);
                
                for (size_t i = 0; i != st.models.size(); ++i)
                    gemmi::prepare_topology(st, monlib, i, h_change, false);
                
                std::printf("Hydrogen site count: %zu in input, %zu in output.\n",
                    initial_h, gemmi::count_hydrogen_sites(st));

                // Clean up the H atoms that were placed at (0.00, 0.00, 0.00)
                for (gemmi::Model& model : st.models)
                    for (gemmi::Chain& chain : model.chains)
                        for (gemmi::Residue& res : chain.residues)
                        {
                            std::vector<gemmi::Atom> modified_atoms;
                            for (gemmi::Atom& atom : res.atoms)
                            {
                                if(atom.is_hydrogen())
                                {
                                    if(atom.pos.x == 0 && atom.pos.y == 0 && atom.pos.z == 0 && atom.occ == 0)
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
            catch (std::exception& e) 
            {
                std::fprintf(stderr, "ERROR: %s\n", e.what());
                throw stderr;
            }
        }

        privateer::interactions::HBondsParser::HBondsParser()
        {
			this->privateer::interactions::HBondsParser::import_ener_lib();
            std::vector<atom_info_for_hydrogen_bonding> empty_marked_hbond_donors_and_acceptors;
            this->marked_hbond_donors_and_acceptors = empty_marked_hbond_donors_and_acceptors;
        }

        std::vector<privateer::interactions::HBond> privateer::interactions::HBondsParser::get_HBonds_via_mcdonald_and_thornton(clipper::MGlycan& input_glycan, clipper::MiniMol& input_model, double max_dist)
        {
            std::vector<privateer::interactions::HBond> output;
            std::vector<privateer::interactions::HBondsParser::atom_info_for_hydrogen_bonding> marked_hbond_donors_and_acceptors = mark_hbond_donors_and_acceptors(input_glycan, input_model);

            return output;
        }

        std::vector<privateer::interactions::HBondsParser::atom_info_for_hydrogen_bonding> privateer::interactions::HBondsParser::mark_hbond_donors_and_acceptors(clipper::MGlycan& input_glycan, clipper::MiniMol& input_model)
        {
            std::vector<privateer::interactions::HBondsParser::atom_info_for_hydrogen_bonding> output;
            std::string monomer_dir = privateer::restraints::check_monlib_access();
            if (monomer_dir.empty())
                throw std::runtime_error("Failed to locate $CLIBD_MON. Have ccp4 env variables been sourced?");

            std::vector<clipper::MSugar> input_sugars = input_glycan.get_sugars();
            for(int sugar = 0; sugar < input_sugars.size(); sugar++)
            {
                clipper::MSugar currentSugar = input_sugars[sugar];
                for(int atom = 0; atom < currentSugar.size(); atom++)
                {
                    clipper::MAtom currentAtom = currentSugar[atom];
                    int h_bond_type = get_h_bond_type(currentAtom, currentSugar.type().trim());
                }
            }

            return output;
        }

        int privateer::interactions::HBondsParser::get_h_bond_type(clipper::MAtom& input_atom, std::string input_residue_type)
        {
            check_or_import_monomer_library_chem_comp_for_residue(input_residue_type);

            return 0;
        }
    }
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
        path << monomer_dir << "ener_lib" << ".cif";
        std::cout << "Path to 'ener_lib.cif': " << path.str() << std::endl;
        gemmi::cif::Document ener_lib_gemmi_document = gemmi::cif::read_file( path.str() );
        for (gemmi::cif::Block& block : ener_lib_gemmi_document.blocks)
        {
            if (!block.name.empty() && block.name == "energy")
            {
                gemmi::cif::Table energy_library_table = block.find("_lib_atom.", 
                                                    {"type", "weight", "hb_type", "vdw_radius", "vdwh_radius", "ion_radius", "element", "valency", "sp"});

                if(energy_library_table.width() != expected_ncolumns_in_energy_library_table)
                        throw std::runtime_error("Number of expected columns doesn't match the number of columns retrieved from 'ener_lib.cif'. Please report this issue to the developers of Privateer.");
                
                for(size_t i = 0; i != energy_library_table.length(); ++i)
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

                    if(export_row.type == "." || export_row.type == "n/a")
                        export_row.type = "null";
                    if(export_row.weight == "." || export_row.weight == "n/a")
                        export_row.weight = "null";
                    if(export_row.hb_type == "." || export_row.hb_type == "n/a")
                        export_row.hb_type = "null";
                    if(export_row.vdw_radius == "." || export_row.vdw_radius == "n/a")
                        export_row.vdw_radius = "null";
                    if(export_row.vdwh_radius == "." || export_row.vdwh_radius == "n/a")
                        export_row.vdwh_radius = "null";
                    if(export_row.ion_radius == "." || export_row.ion_radius == "n/a")
                        export_row.ion_radius = "null";
                    if(export_row.element == "." || export_row.element == "n/a")
                        export_row.element = "null";
                    if(export_row.valency == "." || export_row.valency == "n/a")
                        export_row.valency = "null";
                    if(export_row.sp == "." || export_row.sp == "n/a")
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

privateer::interactions::HBondsParser::residue_monomer_library_chem_comp privateer::interactions::HBondsParser::check_or_import_monomer_library_chem_comp_for_residue(std::string input_residue_type)
{

    auto search_result = std::find_if(std::begin(residue_library), std::end(residue_library),
    [&input_residue_type](residue_monomer_library_chem_comp& element) {
        return  element.residue_type == input_residue_type;
    });

    if(search_result != std::end(residue_library))
    {
        std::cout << "returning search result without using gemmi" << std::endl;
        return *search_result;
    }
    else
    {
        std::string monomer_dir = privateer::restraints::check_monlib_access();
        if (monomer_dir.empty())
            throw std::runtime_error("Failed to locate $CLIBD_MON. Have ccp4 env variables been sourced?");

        std::stringstream path;
        std::locale loc;

        char initial = std::tolower(input_residue_type[0],loc);
        
        int expected_ncolumns_in_chem_comp_atom_table = 6;
        int expected_ncolumns_in_chem_comp_bond_table = 4;

        if (!monomer_dir.empty()) 
        {
            path << monomer_dir << initial << "/" << input_residue_type << ".cif";
            std::cout << path.str() << std::endl;
            gemmi::cif::Document monomer_dict_gemmi_document = gemmi::cif::read_file( path.str() );
            std::string block_name_search_string = "comp_" + input_residue_type;
            std::vector<residue_monomer_library_chem_comp> temp_residue_library;
            for (gemmi::cif::Block& block : monomer_dict_gemmi_document.blocks)
            {
                if (!block.name.empty() && block.name == block_name_search_string)
                {
                    gemmi::cif::Table chem_comp_atom_table = block.find("_chem_comp_atom.", 
                                                                        {"comp_id", "atom_id", "type_symbol", "type_energy", "?charge", "?partial_charge"});
                    
                    gemmi::cif::Table chem_comp_bond_table = block.find("_chem_comp_bond.", 
                                                                        {"comp_id", "atom_id_1", "atom_id_2", "type"});
                    
                    std::cout << chem_comp_atom_table.width() << std::endl;
                    std::cout << chem_comp_bond_table.width() << std::endl;
                    if(chem_comp_atom_table.width() != expected_ncolumns_in_chem_comp_atom_table)
                            throw std::runtime_error("chem_comp_atom_table: Number of expected columns doesn't match the number of columns retrieved from '" + input_residue_type + ".cif', Please report this issue to the developers of Privateer.");
                    
                    if(chem_comp_bond_table.width() != expected_ncolumns_in_chem_comp_bond_table)
                            throw std::runtime_error("chem_comp_bond_table: Number of expected columns doesn't match the number of columns retrieved from '" + input_residue_type + ".cif', Please report this issue to the developers of Privateer.");

                    for(size_t i = 0; i != chem_comp_atom_table.length(); ++i)
                    {
                        gemmi::cif::Table::Row currentAtomTableGemmiRow = chem_comp_atom_table[i];
                        residue_monomer_library_chem_comp export_row;

                        export_row.residue_type = currentAtomTableGemmiRow[0];
                        export_row.atom_id = currentAtomTableGemmiRow[1];
                        export_row.element = currentAtomTableGemmiRow[2];
                        export_row.energy_type = currentAtomTableGemmiRow[3];
                        if(chem_comp_atom_table.has_column(4))
                            export_row.charge = currentAtomTableGemmiRow[4];
                        if(chem_comp_atom_table.has_column(5))
                            export_row.charge = currentAtomTableGemmiRow[5];

                        for(size_t j = 0; j != chem_comp_bond_table.length(); ++j)
                        {
                            gemmi::cif::Table::Row currentBondTableGemmiRow = chem_comp_bond_table[j];
                            if(currentBondTableGemmiRow[1] == export_row.atom_id && currentBondTableGemmiRow[0] == export_row.residue_type)
                            {
                                auto bond_table_search_result = std::find_if(std::begin(temp_residue_library), std::end(temp_residue_library),
                                    [&currentBondTableGemmiRow](residue_monomer_library_chem_comp& element) {
                                        return  element.atom_id == currentBondTableGemmiRow[1] && element.bonded_to_atom_id == currentBondTableGemmiRow[2];
                                    });

                                if(bond_table_search_result == std::end(temp_residue_library))
                                {
                                    export_row.bonded_to_atom_id = currentBondTableGemmiRow[2];
                                    export_row.bond_type = currentBondTableGemmiRow[3];

                                    if(export_row.residue_type == "." || export_row.residue_type == "n/a")
                                        export_row.residue_type = "null";
                                    if(export_row.atom_id == "." || export_row.atom_id == "n/a")
                                        export_row.atom_id = "null";
                                    if(export_row.element == "." || export_row.element == "n/a")
                                        export_row.element = "null";
                                    if(export_row.energy_type == "." || export_row.energy_type == "n/a")
                                        export_row.energy_type = "null";
                                    if(export_row.charge == "." || export_row.charge == "n/a")
                                        export_row.charge = "null";
                                    if(export_row.bonded_to_atom_id == "." || export_row.bonded_to_atom_id == "n/a")
                                        export_row.bonded_to_atom_id = "null";
                                    if(export_row.bond_type == "." || export_row.bond_type == "n/a")
                                        export_row.bond_type = "null";

                                    temp_residue_library.push_back(export_row);
                                }   
                            }
                        }
                    }
                }
            }
            this->residue_library.insert(this->residue_library.end(), std::make_move_iterator(temp_residue_library.begin()), std::make_move_iterator(temp_residue_library.end()));

            std::cout << std::endl;
            for(int i = 0; i < residue_library.size(); i++)
            {
                std::cout   << i << "/" << residue_library.size() << ": " 
                            << residue_library[i].residue_type << "\t"
                            << residue_library[i].atom_id << "\t"
                            << residue_library[i].element << "\t"
                            << residue_library[i].energy_type << "\t"
                            << residue_library[i].charge << "\t"
                            << residue_library[i].bonded_to_atom_id << "\t"
                            << residue_library[i].bond_type << std::endl;
            }
            std::cout << std::endl;

            residue_monomer_library_chem_comp most_recent_row = residue_library.back();
            if (most_recent_row.residue_type == input_residue_type)
                return most_recent_row;
            else
            {
                std::cout << "WARNING! Unable to add a monomer library record for " << input_residue_type << " , either it's missing from the monomer library or there is a bug in Privateer." << std::endl;
                return most_recent_row;
            }
        }
        else
            throw std::runtime_error("Failed to locate $CLIBD_MON. Have ccp4 env variables been sourced?");
    }

}