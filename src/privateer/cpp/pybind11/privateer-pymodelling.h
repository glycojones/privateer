// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#ifndef PRIVATEER_PYMODELLING_H_INCLUDED
#define PRIVATEER_PYMODELLING_H_INCLUDED

#include <unordered_map>
#include "privateer-modelling.h"
#include "privateer-pyanalysis.h"
#include <Python.h>
#include <pybind11_json/pybind11_json.hpp> //nlohman::json
#include <nlohmann/json.hpp> //nlohman::json
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace privateer {

	namespace pymodelling {

			class Builder 
			{
				public:
					Builder() { };
					Builder(std::string& path_to_receiving_model_file, std::string& path_to_donor_model, bool trim_donor_when_clashes_detected, bool ANY_search_policy, bool enable_user_messages, bool debug_output);
					
					void read_from_file( std::string& path_to_receiving_model_file, std::string& path_to_donor_model, bool trim_donor_when_clashes_detected, bool ANY_search_policy, bool enable_user_messages, bool debug_output );
					std::string get_path_of_receiving_model_file_used ( ) { return path_to_receiving_model; };
					std::string get_path_of_donor_model_file_used ( ) { return path_to_donor_model; };
					std::string convert_three_letter_code_to_single_letter (std::string three_letter_code);

					pybind11::list get_receiving_model_sequence_info () { return imported_receiving_model_seq_info; };
					pybind11::list get_glycan_summary_from_donor () { return glycan_summary_donor; };

					pybind11::list get_summary_of_grafted_glycans() { return summary_of_grafted_glycans; };
					void graft_glycan_to_receiver(int mglycanindex, int receiver_chain_index, int received_residue_index);
					void export_grafted_model( std::string& output_path );
					


				private:
					bool ANY_search_policy;
					bool enable_user_messages;
					bool debug_output;
					bool trim_donor_when_clashes_detected;
					std::string path_to_receiving_model;
					std::string path_to_donor_model;

					privateer::modelling::Grafter grafter;

					clipper::MiniMol imported_receiving_model;
					clipper::MiniMol imported_donor_model;

					clipper::MiniMol export_model;

					pybind11::list imported_receiving_model_seq_info;
					pybind11::list glycan_summary_donor;
					pybind11::list summary_of_grafted_glycans;

					pybind11::list get_protein_sequence_information (clipper::MiniMol& input);
			};
		

	}
}


#endif