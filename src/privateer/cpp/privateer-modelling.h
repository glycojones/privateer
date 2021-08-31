// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#ifndef PRIVATEER_MODELLING_H_INCLUDED
#define PRIVATEER_MODELLING_H_INCLUDED

#include <fstream>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <clipper/clipper.h>
#include <clipper/clipper-cif.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <clipper/contrib/sfcalc_obs.h>
#include <clipper/minimol/minimol_utils.h>
#include "clipper-glyco_data.h"
#include "clipper-glyco.h"
#include "privateer-composition.h"

namespace privateer {
	namespace modelling {

		struct protein_sidechain_glycosylation
        {
            clipper::String residue_name;       
            clipper::String connected_atom;          
            clipper::String vector_point_alpha;
			clipper::String vector_point_bravo;     
            clipper::ftype Phi;        
            clipper::ftype Psi;
			clipper::ftype Phi_error;
			clipper::ftype Psi_error;    
        };


		struct sugar_attachment
        {
            clipper::String glycan_type;       
            clipper::String connection_atom;
			clipper::String vector_point_alpha;        
            clipper::String vector_point_bravo;   
        };

		extern const protein_sidechain_glycosylation backbone_instructions[];
        extern const sugar_attachment sugar_instructions[];
		extern const int backbone_instructions_size;
        extern const int sugar_instructions_size;

		class Grafter
    	{
      		public:
				Grafter() { } //!< null constructor
				Grafter(clipper::MiniMol receiving_model, clipper::MiniMol donor_model, bool trim_donor_when_clashes_detected, bool enable_user_messages, bool debug_output);
				int get_number_of_glycans_detected() { return donor_glycans.size(); };
				clipper::MiniMol& get_final_receiving_model() { return export_model; };
				std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > get_grafted_clashes() { return clashes; };
				clipper::MGlycan& get_grafted_glycan() { return grafted_glycan; };
				clipper::MiniMol& get_donor_model() { return donor_model; };
				std::vector<clipper::MGlycan>& get_donor_glycans() { return donor_glycans; };
				clipper::MPolymer convert_mglycan_to_mpolymer(clipper::MGlycan input);
				clipper::Coord_orth get_glycan_target_point(clipper::Coord_orth connecting_atom, clipper::Coord_orth vector_origin, clipper::Coord_orth vector_target, float vectorShiftDistance);
				void overlay_mglycan_via_atom(clipper::Coord_orth target, clipper::Coord_orth origin, clipper::MPolymer& converted_mglycan);
				void graft_mpolymer_to_receiving_model(clipper::MGlycan& glycan_to_graft, clipper::MMonomer& input_protein_side_chain_residue, clipper::String root_chain_id, bool ANY_search_policy);
				clipper::Coord_orth generate_rotation_matrix_from_rodrigues_rotation_formula(clipper::Coord_orth direction, clipper::Coord_orth position, clipper::Coord_orth origin_shift, double targetAngle);
				void rotate_mglycan_until_torsion_angle_fulfilled(clipper::MPolymer& converted_mglycan, clipper::MMonomer& protein_residue, clipper::Coord_orth direction, clipper::Coord_orth origin_shift, std::vector<std::pair<clipper::MAtom, std::string>>& torsionAtoms, double angle, bool debug_output);
				clipper::MPolymer rotate_mglycan_until_clashes_are_minimized(clipper::MiniMol& export_model, clipper::MPolymer& converted_mglycan, clipper::MMonomer& protein_residue, std::vector<std::pair<clipper::MAtom, std::string>>& phiTorsionAtoms, std::vector<std::pair<clipper::MAtom, std::string>>& psiTorsionAtoms, double phiError, double psiError, clipper::String root_chain_id, clipper::String root_sugar_chain_id, bool debug_output);

				std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > check_for_clashes(clipper::MiniMol& input_model, clipper::MMonomer& root_sugar, clipper::MMonomer& input_protein_side_chain_residue, clipper::String root_chain_id, clipper::String root_sugar_chain_id);
				clipper::MiniMol trim_graft_until_no_clashes_left(clipper::MiniMol& current_model, clipper::MGlycan& grafted_glycan, std::pair<clipper::MMonomer, clipper::String> root_residue, clipper::String graft_chain_id, std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > current_clashes);

				clipper::MiniMol get_model_with_trimmed_glycan(clipper::MiniMol& current_model, std::vector<std::pair<clipper::MMonomer, clipper::String>>& residues_to_delete);
				std::vector<std::pair<clipper::MMonomer, clipper::String>> get_designated_residues_for_deletion(clipper::MGlycan& original_graft, clipper::MGlycan& trimmed_graft, clipper::String graft_chain_id);
				bool check_if_residue_has_hydrogens(clipper::MMonomer residue_to_check);
				clipper::MPolymer delete_atom_from_mglycan(clipper::MPolymer& converted_mglycan, clipper::MAtom& atom_to_be_deleted);

				int lookup_protein_backbone_glycosylation_database( clipper::String name);
				int lookup_glycan_type_glycosylation_database( clipper::String type);
			private:
				bool trim_donor_when_clashes_detected;
				bool enable_user_messages;
				bool debug_output;

				clipper::MiniMol receiving_model;
				clipper::MiniMol donor_model;

				clipper::MiniMol export_model;
				clipper::MGlycan grafted_glycan; 

				std::vector< std::pair< std::pair<clipper::MMonomer, clipper::String>, std::pair<clipper::MSugar, clipper::String> > > clashes;

				std::vector<clipper::MGlycan> donor_glycans;
				int numDonorGlycansDetected;

    	};


  	}
}


#endif