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
            try 
            {
                std::string output = "hydrogenated_input_model.pdb";
                gemmi::HydrogenChange h_change = gemmi::HydrogenChange::ReAdd;
                std::string monomer_dir = privateer::restraints::check_monlib_access();
                if (monomer_dir.empty())
                    throw std::runtime_error("Failed to locate $CLIBD_MON. Have ccp4 env variables been sourced?");;
                
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
            }
            catch (std::exception& e) 
            {
                std::fprintf(stderr, "ERROR: %s\n", e.what());
                throw stderr;
            }
        }
    }
}
