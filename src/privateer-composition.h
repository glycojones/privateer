
// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2018 Haroldas Bagdonas & Kevin Cowtan & Jon Agirre
// York Structural Biology Laboratory
// The University of York
// mailto: hb1115@york.ac.uk
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk


#ifndef COMPOSITION_H_INCLUDED
#define COMPOSITION_H_INCLUDED

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
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-minimol.h>
#include "clipper-glyco.h"
#include "clipper-glyco_data.h"
#include <clipper/minimol/minimol_utils.h>
#include "privateer-dbquery.h"
#include "privateer-lib.h"
#include <nlohmann/json.hpp>


std::vector<std::pair<clipper::MGlycan, std::vector<int>>> generate_closest_matches(clipper::MGlycan& fullglycan, nlohmann::json &jsonObject);
std::vector<int> get_editable_node_list_for_anomer_permutations(std::vector < clipper::MSugar >& sugar_list);
std::vector<int> get_editable_node_list_for_monomer_permutations(std::vector < clipper::MSugar >& sugar_list);
void generate_all_anomer_permutations(std::vector<std::pair<clipper::MGlycan, std::vector<int>>>& result, std::vector<int>& editable_node_list, clipper::MGlycan glycan, nlohmann::json& jsonObject, int residuePermuations, int residueDeletions);
void generate_all_monomer_permutations(std::vector<std::pair<clipper::MGlycan, std::vector<int>>>& result, std::vector<int>& editable_node_list, clipper::MGlycan glycan, nlohmann::json& jsonObject, int residueDeletions);
std::vector<std::vector<int>> generate_all_possible_index_combinations(std::vector<int>& editable_node_list);
void CombinationGenerator(std::vector<int> editable_node_list, int reqLength, int shifter, int currLength, bool checkedEditableNode[], int totalEditableNodes, std::vector<std::vector<int>>& result, std::vector<int>& temp);
clipper::MGlycan get_shorter_matching_fragment(clipper::MGlycan& fullglycan, nlohmann::json &jsonObject);
clipper::MGlycan shorten_fragment_by_removing_last_node(int totalNodes, clipper::MGlycan inputglycan, nlohmann::json& jsonObject);
clipper::MGlycan shorten_fragment_by_removing_first_leaf_node(clipper::MGlycan inputglycan, nlohmann::json& jsonObject);
clipper::MGlycan remove_first_leaf_node(clipper::MGlycan inputglycan);



#endif