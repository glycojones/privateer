// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#ifndef PRIVATEER_INTERACTIONS_H_INCLUDED
#define PRIVATEER_INTERACTIONS_H_INCLUDED

#include <fstream>
#include <string>
#include <locale>
#include <algorithm>
#include "clipper-glyco.h"
#include "privateer-restraints.h"
#include "privateer-lib.h"

#include "gemmi/mmread.hpp"
#include "gemmi/monlib.hpp"
#include "gemmi/placeh.hpp"
#include "gemmi/fstream.hpp"
#include "gemmi/cif.hpp"
#include "gemmi/to_pdb.hpp"



namespace privateer {
	namespace interactions {
		void hydrogenate_input_model(std::string input_model, std::string output_path);

		class HBond
		{
			public:
				std::string donor_chainID;
				std::string acceptor_chainID;
				clipper::MMonomer donor_residue;
				clipper::MMonomer acceptor_residue;
				clipper::MAtom HBonding_hydrogen;
				clipper::MAtom donor;
				clipper::MAtom acceptor;
				clipper::MAtom acceptor_neighbour;
				int sugarIndex;
				int glycanSize;
				float angle_1;
				float angle_2;
				float angle_3;
				float HBondLength;
				bool sugar_atom_is_donor;
				bool hydrogen_is_sugar_atom;
				bool bond_has_hydrogen_flag;

				HBond()
				{
					donor_chainID = "-1";
					acceptor_chainID = "-1";
					donor_residue = clipper::MMonomer();
					acceptor_residue = clipper::MMonomer();
					HBonding_hydrogen = clipper::MAtom().null();
					donor = clipper::MAtom().null();
					acceptor = clipper::MAtom().null();
					acceptor_neighbour = clipper::MAtom().null();
					sugar_atom_is_donor = true;
					int sugarIndex = -1;
					int glycanSize = -1;
					angle_1 = -1;
					angle_2 = -1;
					angle_3 = -1;
					HBondLength = -1;

					hydrogen_is_sugar_atom = false;
					bond_has_hydrogen_flag = false;
				}

				HBond(std::string input_donor_chainID, clipper::MMonomer& input_donor_residue, clipper::MAtom& input_donor, std::string input_acceptor_chainID, clipper::MMonomer& input_acceptor_residue, clipper::MAtom& input_acceptor)
				{
					donor_chainID = input_donor_chainID;
					acceptor_chainID = input_acceptor_chainID;
					donor_residue = input_donor_residue;
					acceptor_residue = input_acceptor_residue;
					HBonding_hydrogen = clipper::MAtom().null();
					donor = input_donor;
					acceptor = input_acceptor;
					acceptor_neighbour = clipper::MAtom().null();
					sugar_atom_is_donor = false;
					int sugarIndex = -1;
					int glycanSize = -1;
					angle_1 = -1;
					angle_2 = -1;
					angle_3 = -1;
					HBondLength = -1;

					hydrogen_is_sugar_atom = false;
					bond_has_hydrogen_flag = false;
				}

				HBond(std::string input_hydrogen_chainID, clipper::MMonomer& input_hydrogen_residue, clipper::MAtom& input_hydrogen, std::string input_acceptor_chainID, clipper::MMonomer& input_acceptor_residue, clipper::MAtom& input_acceptor, bool input_sugar_atom_is_H_flag)
				{
					donor_chainID = input_hydrogen_chainID;
					acceptor_chainID = input_acceptor_chainID;
					donor_residue = input_hydrogen_residue;
					acceptor_residue = input_acceptor_residue;
					HBonding_hydrogen = input_hydrogen;
					bond_has_hydrogen_flag = true;
					acceptor = input_acceptor;
					donor = clipper::MAtom().null();
					acceptor_neighbour = clipper::MAtom().null();
					sugar_atom_is_donor = input_sugar_atom_is_H_flag;
					int sugarIndex = -1;
					int glycanSize = -1;
					angle_1 = -1;
					angle_2 = -1;
					angle_3 = -1;
					HBondLength = -1;

					hydrogen_is_sugar_atom = input_sugar_atom_is_H_flag;
				}

				bool has_hydrogen() const { return bond_has_hydrogen_flag; }
				bool ligand_atom_is_H() const { return hydrogen_is_sugar_atom; }
		};

		class CHPiBond
		{
			public:
				CHPiBond(std::string input_sugar_chainID, std::string input_stacked_residue_chainID, clipper::MSugar& input_sugar, clipper::MMonomer& input_stacked_residue, float input_angle, std::string algorithm_used = "hudson" )
				{
					sugar_chainID = input_sugar_chainID;
					stacked_residue_chainID = input_stacked_residue_chainID;
					sugar = input_sugar;
					stacked_residue = input_stacked_residue;
					algorithm = algorithm_used;
					angle = input_angle;
					sugarIndex = -1;
					glycanSize = -1;
					angle_theta_h = 0.0f;
					angle_theta_p = 0.0f;
					angle_phi = 0.0f;
					distance_cx = 0.0f;
					distance_cp = 0.0f;
					sugarFace = "β";

				}
				
				struct chpi_hudson_parameters
                {
                    clipper::ftype xo_distance;
                    clipper::ftype theta;
                    clipper::ftype xp_distance;
                
                chpi_hudson_parameters(clipper::ftype xo_distance, clipper::ftype theta, clipper::ftype xp_distance)
                {
                    this->xo_distance = xo_distance;
                    this->theta = theta;
                    this->xp_distance = xp_distance;
                }
                };

				
				std::string get_algorithm ( ) {
					return this->algorithm;
				}

				void set_algorithm ( ) {
					this->algorithm = algorithm;
				}
				
				std::string get_sugar_chainID ( ) {
					return this->sugar_chainID;
				}

				void set_sugar_chainID ( ) {
					this->sugar_chainID = sugar_chainID;
				}

				int get_sugar_index ( ) {
					return this->sugarIndex;
				}

				void set_sugar_index ( int sugar_index ) {
					this->sugarIndex = sugar_index;
				}

				int get_glycan_size ( ) {
					return this->glycanSize;
				}

				void set_glycan_size ( int glycan_size ) {
					this->glycanSize = glycan_size;
				}

				clipper::MSugar get_sugar () {
					return this->sugar;
				}

				void set_sugar ( clipper::MSugar sugar ) {
					this->sugar = sugar;
				}

				clipper::MMonomer get_stacked_residue ( ){
					return this->stacked_residue;
				}

				std::string get_stacked_residue_chainID ( ){
					return this->stacked_residue_chainID;
				}

				float get_angle_theta_h ( ) {
					return this->angle_theta_h;
				}

				void set_angle_theta_h ( float angle_theta_h ) {
					this->angle_theta_h = angle_theta_h;
				}

				float get_angle_theta_p ( ) {
					return this->angle_theta_p;
				}

				void set_angle_theta_p ( float angle_theta_p ) {
					this->angle_theta_p = angle_theta_p;
				}

				float get_angle_phi ( ) {
					return this->angle_phi;
				}

				void set_angle_phi ( float angle_phi ) {
					this->angle_phi = angle_phi;
				}

				float get_distance_cx ( ) {
					return this->distance_cx;
				}

				void set_distance_cx ( float distance_cx ) {
					this->distance_cx = distance_cx;
				}

				float get_distance_cp ( ) {
					return this->distance_cp;
				}

				void set_distance_cp ( float distance_cp ) {
					this->distance_cp = distance_cp;
				}

				std::string get_trp_ring () {
					return trp_ring;
				}

				void set_trp_ring ( std::string trp_ring ) {
					this->trp_ring = trp_ring;
				}

				std::pair<clipper::MAtom, clipper::MAtom> get_xh_pair () {
					return xh_pair;
				}
				
				void set_xh_pair ( std::pair<clipper::MAtom, clipper::MAtom> xh_pair) {
					this->xh_pair = xh_pair;
				}

				std::string get_sugar_face ( ) {
					return this->sugarFace;
				}

				void set_sugar_face ( std::string sugar_face ) {
					this->sugarFace = sugar_face;
				}
				
				static std::vector<clipper::ftype> calculate_hudson_parameters(const clipper::MAtomIndexSymmetry &neighbourhood, 
                                                                  			   const std::pair<clipper::MAtom, clipper::MAtom> &ch_atoms, 
																			   const std::string &trp_ring,
																			   const clipper::MiniMol &hydrogenated_input_model,
																			   const clipper::MMonomer &mmon);

				static std::vector<clipper::ftype> calculate_plevin_parameters(const clipper::MAtomIndexSymmetry &neighbourhood, 
                                                                  			   const std::pair<clipper::MAtom, clipper::MAtom> &ch_atoms, 
																			   const std::string &trp_ring,
																			   const clipper::MiniMol &hydrogenated_input_model,
																			   const clipper::MMonomer &mmon);
				
				private:
				clipper::MiniMol hydrogenated_input_model;
				std::string sugar_chainID;
				std::string stacked_residue_chainID;
				clipper::MSugar sugar;
				clipper::MMonomer stacked_residue;
				float angle;
				std::string algorithm;
				int sugarIndex;
				int glycanSize;
				std::pair<clipper::MAtom, clipper::MAtom> xh_pair;
				float angle_theta_h;
				float angle_theta_p;
				float angle_phi;
				float distance_cx;
				float distance_cp;
				std::string sugarFace;
				std::string trp_ring; // For TRP exclusively, A + B
				// clipper::Coord_orth get_aromatic_centre ( clipper::MMonomer mmon, std::string ring = "A" );
				// clipper::ftype get_angle ( clipper::Vec3<clipper::ftype> vec1, clipper::Vec3<clipper::ftype> vec2 );
				// clipper::Vec3<clipper::ftype> find_aromatic_plane ( clipper::MMonomer mmon );
    	};

		class CHPiBondsParser
		{
			
			typedef std::vector <privateer::interactions::CHPiBond> chpibonds;
	
			enum algorithm_type {HUDSON, PLEVIN};
			
			public:
				CHPiBondsParser() { }
				CHPiBondsParser(std::string& input_model, std::string output_path = "undefined", std::string algorithm = "hudson");
				std::vector<privateer::interactions::CHPiBond> get_CHPi_interactions(int glycanIndex);
				std::vector <privateer::interactions::CHPiBond> get_stacked_residues_python(clipper::MSugar& input_sugar,
																							int sugarIndex,
																							int glycanSize,
																							std::string = "hudson",
																							float = 4.5,
																							float = 40.0,
																							float = 0.0,
																							std::string sugarFace = "β") const ;
			private:
				clipper::MiniMol input_model;
				clipper::MAtomNonBond manb_object;
				clipper::MGlycology mglycology;
				clipper::MiniMol hydrogenated_input_model;
				clipper::MGlycology hydrogenated_mglycology;
				std::string algorithm;

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

				HBondsParser() { }
				HBondsParser(std::string& input_model, std::string output_path = "undefined");
				clipper::MiniMol mark_hbond_donors_and_acceptors(clipper::MiniMol& input_model);
				hb_type get_h_bond_type(clipper::MAtom& input_atom, std::string input_residue_type);
				bool is_connected_to_hydrogen_donor(std::string atom_name, std::vector<residue_monomer_library_chem_comp>& residue_atoms);
				std::vector<privateer::interactions::HBond> get_HBonds_via_mcdonald_and_thornton(int glycanIndex, double max_dist = 3.9);
				std::vector<std::pair<clipper::MAtom, float>> get_closest_neighbour_atoms(clipper::MAtom& input_atom, clipper::MMonomer& residue_atom_located_in, clipper::String chainID);
				std::pair<bool, HBond> make_h_bond_from_sugar_hydrogen(std::string input_hydrogen_chainID, clipper::MMonomer& input_hydrogen_residue, clipper::MAtom& hydrogen, std::string input_acceptor_chainID, clipper::MMonomer& input_acceptor_residue, clipper::MAtom& acceptor, std::vector<std::pair<clipper::MAtom, float>>& hydrogen_neighbours, std::vector<std::pair<clipper::MAtom, float>>& acceptor_neighbours);
				std::pair<bool, HBond> make_h_bond_from_environment_residue_hydrogen(std::string input_hydrogen_chainID, clipper::MMonomer& input_hydrogen_residue, clipper::MAtom& hydrogen, std::string input_acceptor_on_sugar_chainID, clipper::MMonomer& input_acceptor_on_sugar_residue, clipper::MAtom& acceptor_on_sugar, std::vector<std::pair<clipper::MAtom, float>>& hydrogen_neighbours, std::vector<std::pair<clipper::MAtom, float>>& acceptor_on_sugar_neighbours);
			private:
				std::vector<energy_library_entry> energy_library;
				std::vector<monomer_dictionary> monomer_dict;
				clipper::MiniMol hydrogenated_input_model;
				clipper::MAtomNonBond manb_object;
				clipper::MGlycology hydrogenated_mglycology;

				void import_ener_lib();
				monomer_dictionary check_or_import_monomer_library_chem_comp_for_residue(std::string input_residue_type);

				// std::vector < std::pair< clipper::MAtomIndexSymmetry, clipper::ftype > > MSugar::get_stacked_residues
		};


  	}
}


#endif
