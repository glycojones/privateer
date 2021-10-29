// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#ifndef PRIVATEER_MODELLING_H_INCLUDED
#define PRIVATEER_MODELLING_H_INCLUDED

#include <pybind11/pybind11.h>
#include <fstream>
#include <string>
#include <locale>
#include <algorithm>
#include "clipper-glyco.h"
#include "privateer-restraints.h"
#include "privateer-lib.h"
#include <gemmi/mmread.hpp>
#include <gemmi/monlib.hpp>
#include <gemmi/placeh.hpp>
#include <gemmi/fstream.hpp>
#define GEMMI_READ_CIF_IMPLEMENTATION
#include <gemmi/read_cif.hpp>
#define GEMMI_WRITE_IMPLEMENTATION
#include <gemmi/to_pdb.hpp>



namespace privateer {
	namespace interactions {
		void hydrogenate_input_model(std::string input_model);

		class HBond 
		{
			public: 
				clipper::MAtom HBonding_hydrogen;
				clipper::MAtom donor;
				clipper::MAtom acceptor;
				clipper::MAtom donor_neighbour;
				clipper::MAtom acceptor_neighbour;
				double angle_1;
				double angle_2;
				double angle_3;
				double HBondLength;
				bool sugar_atom_is_donor;
				bool hydrogen_is_sugar_atom;
				bool bond_has_hydrogen_flag;

				HBond() 
				{
					HBonding_hydrogen = clipper::MAtom();
					donor = clipper::MAtom();
					acceptor = clipper::MAtom();
					donor_neighbour = clipper::MAtom();
					acceptor_neighbour = clipper::MAtom();
					sugar_atom_is_donor = true;
					angle_1 = -1;
					angle_2 = -1;
					angle_3 = -1;
					HBondLength = -1;
					
					hydrogen_is_sugar_atom = false;
					bond_has_hydrogen_flag = false;
				}

				HBond(clipper::MAtom& input_donor, clipper::MAtom input_acceptor)
				{
					HBonding_hydrogen = clipper::MAtom();
					donor = input_donor;
					acceptor = input_acceptor;
					donor_neighbour = clipper::MAtom();
					acceptor_neighbour = clipper::MAtom();
					sugar_atom_is_donor = true;
					angle_1 = -1;
					angle_2 = -1;
					angle_3 = -1;
					HBondLength = -1;
					
					hydrogen_is_sugar_atom = false;
					bond_has_hydrogen_flag = false;
				}

				HBond(clipper::MAtom& input_hydrogen, clipper::MAtom input_acceptor, bool input_sugar_atom_is_H_flag)
				{
					HBonding_hydrogen = input_hydrogen;
					bond_has_hydrogen_flag = true;
					acceptor = input_acceptor;
					donor = clipper::MAtom();
					donor_neighbour = clipper::MAtom();
					acceptor_neighbour = clipper::MAtom();
					sugar_atom_is_donor = input_sugar_atom_is_H_flag;
					angle_1 = -1;
					angle_2 = -1;
					angle_3 = -1;
					HBondLength = -1;
					
					hydrogen_is_sugar_atom = input_sugar_atom_is_H_flag;
				}

				bool has_hydrogen() const { return bond_has_hydrogen_flag; }
				bool ligand_atom_is_H() const { return hydrogen_is_sugar_atom; }
		};

		class HBondsParser
		{
			public:
				enum hb_type
				{
					HB_UNASSIGNED = -1,
					HB_NEITHER = 0,
					HB_DONOR = 1,
					HB_ACCEPTOR = 2, 
					HB_BOTH = 3,
					HB_HYDROGEN = 4
				};


				struct energy_library_entry
				{
					std::string type; 
					std::string weight; 
					std::string hb_type;
					std::string vdw_radius;
					std::string vdwh_radius;
					std::string ion_radius;
					std::string element; 
					std::string valency;
					std::string sp;
				};

				struct residue_monomer_library_chem_comp
				{
					std::string residue_type;
					std::string atom_id;
					std::string element;
					std::string energy_type;
					std::string charge;
					std::string bonded_to_atom_id;
					std::string bond_type;
				};

				struct monomer_dictionary
				{
					std::string monomer_name;
					std::vector<residue_monomer_library_chem_comp> dictionary_of_atoms;
				};

				HBondsParser(std::string& input_model);
				std::vector<privateer::interactions::HBond> get_HBonds_via_mcdonald_and_thornton(clipper::MGlycan& input_glycan, clipper::MiniMol& input_model, double max_dist = 3.9);
				clipper::MiniMol mark_hbond_donors_and_acceptors(clipper::MiniMol& input_model);
				hb_type get_h_bond_type(clipper::MAtom& input_atom, std::string input_residue_type);
				bool is_connected_to_hydrogen_donor(std::string atom_name, std::vector<residue_monomer_library_chem_comp>& residue_atoms);
			private:
				std::vector<energy_library_entry> energy_library;
				std::vector<monomer_dictionary> monomer_dict;
				clipper::MiniMol hydrogenated_input_model;

				void import_ener_lib();
				monomer_dictionary check_or_import_monomer_library_chem_comp_for_residue(std::string input_residue_type);
		};

        
  	}
}


#endif