//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2018 Haroldas Bagdonas & Kevin Cowtan & Jon Agirre
// York Structural Bioexpy Laboratory
// The University of York
// mailto: hb1115@york.ac.uk
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk

#include "privateer-composition.h"
// #define DUMP 1
// #define DBG std::cout << "[" << __FUNCTION__ << "] - "


std::vector<std::pair<clipper::MGlycan, std::vector<int>>> generate_closest_matches(clipper::MGlycan& fullglycan, nlohmann::json& jsonObject, bool glucose_only, privateer::thread_pool& pool, bool useParallelism)
{
    std::vector<std::pair<clipper::MGlycan, std::vector<int>>> resultAlpha;
    std::vector<std::pair<clipper::MGlycan, std::vector<int>>> resultBravo;
    std::vector<std::pair<clipper::MGlycan, std::vector<int>>> resultFinal; // std::vector.push_back() is not thread safe, need to create temporary vectors.


    int totalNodes = fullglycan.number_of_nodes();
    
    clipper::MGlycan tempGlycanAlpha = fullglycan;
    clipper::MGlycan tempGlycanBravo = fullglycan;

    bool finishedAlpha = false, 
         finishedBravo = false;

    // TO DO: Need to design a seperate function specifically for multithreaded aspect.
    //        In essence, have main thread collate all results from spawned threads. Current approach is not going to work at all.
    if(useParallelism)
    {
        pool.push([&tempGlycanAlpha, &resultAlpha, &jsonObject, glucose_only, totalNodes, &finishedAlpha](int id)
        {
            #if DUMP
                std::cout << std::endl;
                DBG << "Calculating tempGlycanAlpha from Thread ID: " << id << '.' << std::endl;
            #endif

            int residueDeletions = 0;
            
            for(int i = 0; i < totalNodes; i++)
            {
                int residuePermutations = 0;

                std::vector < clipper::MSugar > sugars = tempGlycanAlpha.get_sugars();
                std::vector<int> editable_nodes_for_anomer_permutations = get_editable_node_list_for_anomer_permutations(sugars);

                generate_all_anomer_permutations(resultAlpha, editable_nodes_for_anomer_permutations, tempGlycanAlpha, jsonObject, residuePermutations, residueDeletions);

                std::vector<int> editable_nodes_for_monomer_permutations = get_editable_node_list_for_monomer_permutations(sugars, glucose_only);

                generate_all_monomer_permutations(resultAlpha, editable_nodes_for_monomer_permutations, tempGlycanAlpha, jsonObject, residueDeletions);

                if(tempGlycanAlpha.number_of_nodes() > 1)
                {
                    bool glyConnectTrue = false;
                    
                    tempGlycanAlpha = remove_first_leaf_node(tempGlycanAlpha);
                    residueDeletions++;

                    clipper::String temporaryWURCS = tempGlycanAlpha.generate_wurcs();

                    int valueLocation = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCS);

                    if (valueLocation != -1)
                    {
                        // std::cout << "PERMUTATION: " << temporaryWURCS << std::endl << "valueLocation = " << valueLocation << std::endl;
                        if (jsonObject[valueLocation]["glyconnect"] != "NotFound") glyConnectTrue = true;
                        // std::cout << "GlyConnect status " << std::boolalpha << glyConnectTrue << std::endl;
                    }

                    if (glyConnectTrue)
                        {
                            std::vector<int> mutations(3);
                            mutations[0] = 0;
                            mutations[1] = residuePermutations;
                            mutations[2] = residueDeletions;
                            std::pair <clipper::MGlycan,std::vector<int>> tempPair;
                            tempPair = std::make_pair(tempGlycanAlpha, mutations);
                            resultAlpha.push_back(tempPair);
                        }
                }
            }
            finishedAlpha = true;
        });

        pool.push([&tempGlycanBravo, &resultBravo, &jsonObject, glucose_only, totalNodes, &finishedBravo](int id)
        { 
            #if DUMP
                std::cout << std::endl;
                DBG << "Calculating tempGlycanBravo from Thread ID: " << id << '.' << std::endl;
            #endif

            int residueDeletions = 0;
            
            for(int i = 0; i < totalNodes; i++)
            {
                int residuePermutations = 0;

                std::vector < clipper::MSugar > sugars = tempGlycanBravo.get_sugars();
                std::vector<int> editable_nodes_for_anomer_permutations = get_editable_node_list_for_anomer_permutations(sugars);

                generate_all_anomer_permutations(resultBravo, editable_nodes_for_anomer_permutations, tempGlycanBravo, jsonObject, residuePermutations, residueDeletions);

                std::vector<int> editable_nodes_for_monomer_permutations = get_editable_node_list_for_monomer_permutations(sugars, glucose_only);

                generate_all_monomer_permutations(resultBravo, editable_nodes_for_monomer_permutations, tempGlycanBravo, jsonObject, residueDeletions);

                if(tempGlycanBravo.number_of_nodes() > 1)
                {
                    bool glyConnectTrue = false;
                    
                    tempGlycanBravo = remove_last_node(tempGlycanBravo);
                    residueDeletions++;

                    clipper::String temporaryWURCS = tempGlycanBravo.generate_wurcs();

                    int valueLocation = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCS);

                    if (valueLocation != -1)
                    {
                        // std::cout << "PERMUTATION: " << temporaryWURCS << std::endl << "valueLocation = " << valueLocation << std::endl;
                        if (jsonObject[valueLocation]["glyconnect"] != "NotFound") glyConnectTrue = true;
                        // std::cout << "GlyConnect status " << std::boolalpha << glyConnectTrue << std::endl;
                    }

                    if (glyConnectTrue)
                        {
                            std::vector<int> mutations(3);
                            mutations[0] = 0;
                            mutations[1] = residuePermutations;
                            mutations[2] = residueDeletions;
                            std::pair <clipper::MGlycan,std::vector<int>> tempPair;
                            tempPair = std::make_pair(tempGlycanBravo, mutations);
                            resultBravo.push_back(tempPair);
                        }
                }
            }
            finishedBravo = true;
        });

    
        #if DUMP
            std::cout << std::endl;
            DBG << "Number of jobs in the queue: " << pool.n_remaining_jobs() << std::endl;
        #endif
        
        while(!finishedAlpha && !finishedBravo)
        {
             #if DUMP
                std::cout << std::endl;
                DBG << "pool.n_idle() = " << pool.n_idle() << "\t pool.size() = " << pool.size() << std::boolalpha << "\tfinishedAlpha = " << finishedAlpha << "\tfinishedBravo = " << finishedBravo << std::endl;
            #endif
            pool.sync();
        }

        #if DUMP
            std::cout << std::endl;
            DBG << std::boolalpha << "finishedAlpha = " << finishedAlpha << "\tfinishedBravo = " << finishedBravo << std::endl;
        #endif
            
        #if DUMP
            DBG << "Number of jobs in the queue: " << pool.n_remaining_jobs() << " after sync operation!" << std::endl;
        #endif
    
    
        resultFinal.reserve( resultAlpha.size() + resultBravo.size() ); // preallocate memory
        resultFinal.insert( resultFinal.end(), resultAlpha.begin(), resultAlpha.end() );
        resultFinal.insert( resultFinal.end(), resultBravo.begin(), resultBravo.end() );
    }
    else
    {
        int residueDeletions = 0;

        for(int i = 0; i < totalNodes; i++)
        {
            int residuePermutations = 0;

            std::vector < clipper::MSugar > sugars = tempGlycanAlpha.get_sugars();
            std::vector<int> editable_nodes_for_anomer_permutations = get_editable_node_list_for_anomer_permutations(sugars);

            generate_all_anomer_permutations(resultFinal, editable_nodes_for_anomer_permutations, tempGlycanAlpha, jsonObject, residuePermutations, residueDeletions);

            std::vector<int> editable_nodes_for_monomer_permutations = get_editable_node_list_for_monomer_permutations(sugars, glucose_only);

            generate_all_monomer_permutations(resultFinal, editable_nodes_for_monomer_permutations, tempGlycanAlpha, jsonObject, residueDeletions);

            if(tempGlycanAlpha.number_of_nodes() > 1)
            {
                bool glyConnectTrue = false;
                
                tempGlycanAlpha = remove_first_leaf_node(tempGlycanAlpha);
                residueDeletions++;

                clipper::String temporaryWURCS = tempGlycanAlpha.generate_wurcs();

                int valueLocation = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCS);

                if (valueLocation != -1)
                {
                    // std::cout << "PERMUTATION: " << temporaryWURCS << std::endl << "valueLocation = " << valueLocation << std::endl;
                    if (jsonObject[valueLocation]["glyconnect"] != "NotFound") glyConnectTrue = true;
                    // std::cout << "GlyConnect status " << std::boolalpha << glyConnectTrue << std::endl;
                }

                if (glyConnectTrue)
                {
                    std::vector<int> mutations(3);
                    mutations[0] = 0;
                    mutations[1] = residuePermutations;
                    mutations[2] = residueDeletions;
                    std::pair <clipper::MGlycan,std::vector<int>> tempPair;
                    tempPair = std::make_pair(tempGlycanAlpha, mutations);
                    resultFinal.push_back(tempPair);
                }
            }
        }
    }
        
    // do the score calculations here, return the final vector already containing those float values. 

    return resultFinal;
}

std::vector<int> get_editable_node_list_for_anomer_permutations(std::vector < clipper::MSugar >& sugar_list)
{
    std::vector<int> list;
    for(int i = 0; i < sugar_list.size(); i++)
        {
            if(clipper::data::residue_has_alternate_anomer(sugar_list[i].type().trim())) 
                list.push_back(i);
        }
    return list; 
}

std::vector<int> get_editable_node_list_for_monomer_permutations(std::vector < clipper::MSugar >& sugar_list, bool glucose_only)
{
    std::vector<int> list;
    for(int i = 0; i < sugar_list.size(); i++)
        {
            if(clipper::data::residue_has_alternate_monomer(sugar_list[i].type().trim(), glucose_only)) 
                list.push_back(i);
        }
    return list; std::cout << std::endl;
}



void generate_all_anomer_permutations(std::vector<std::pair<clipper::MGlycan, std::vector<int>>>& result, std::vector<int>& editable_node_list, clipper::MGlycan glycan, nlohmann::json& jsonObject, int residuePermuations, int residueDeletions)
{
    std::vector<std::vector<int>> totalCombinations = generate_all_possible_index_combinations(editable_node_list);

    #pragma omp parallel for
    for(int i = 0; i < totalCombinations.size(); i++)
    {
        clipper::MGlycan tempGlycan = glycan;
        bool glyConnectTrue = false;
        int anomerPermutations = 0;
        
        for(int j = 0; j < totalCombinations[i].size(); j++)
        {
            int nodeID = totalCombinations[i][j];

            clipper::MSugar msug = tempGlycan.get_node(nodeID).get_sugar();
            msug.set_type(clipper::data::alternative_anomer(msug.type().trim()));

            tempGlycan.replace_sugar_at_index(nodeID, msug);

            if(nodeID == 0) tempGlycan.update_msugar_in_root(msug);

        }     
        clipper::String temporaryWURCS = tempGlycan.generate_wurcs();
        int valueLocation = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCS);

        // std::cout << "ANOMER PERMUTATION: " << temporaryWURCS << std::endl << "valueLocation = " << valueLocation << std::endl;

        if (valueLocation != -1)
            {
                // std::cout << "PERMUTATION: " << temporaryWURCS << std::endl << "valueLocation = " << valueLocation << std::endl;
                if (jsonObject[valueLocation]["glyconnect"] != "NotFound") glyConnectTrue = true;
                // std::cout << "GlyConnect status " << std::boolalpha << glyConnectTrue << std::endl;
            }

        if (glyConnectTrue)
            {
                anomerPermutations = totalCombinations[i].size();
                std::vector<int> mutations(3);
                mutations[0] = anomerPermutations;
                mutations[1] = residuePermuations;
                mutations[2] = residueDeletions;
                std::pair <clipper::MGlycan,std::vector<int>> tempPair;
                tempPair = std::make_pair(tempGlycan, mutations);
                result.push_back(tempPair);
            }
    }
}

// need to rework this.
void generate_all_monomer_permutations(std::vector<std::pair<clipper::MGlycan, std::vector<int>>>& result, std::vector<int>& editable_node_list, clipper::MGlycan glycan, nlohmann::json& jsonObject, int residueDeletions)
{
    std::vector<std::vector<int>> totalCombinations = generate_all_possible_index_combinations(editable_node_list);


    for(int i = 0; i < totalCombinations.size(); i++)
    {

        clipper::MGlycan tempGlycanAlpha = glycan, 
                         tempGlycanBravo = glycan;

        bool glyConnectTrueAlpha = false,
             glyConnectTrueBravo = false;

        int residuePermutationsAlpha = 0;
        int residuePermutationsBravo = 0;
        
        int anomerPermutationsAlpha = 0; 
        int anomerPermutationsBravo = 0;

        for(int j = 0; j < totalCombinations[i].size(); j++)
        {
            int nodeID = totalCombinations[i][j];

            clipper::MSugar msugAlpha = tempGlycanAlpha.get_node(nodeID).get_sugar(), 
                            msugBravo = tempGlycanBravo.get_node(nodeID).get_sugar();

            std::vector<std::string> alternative_monomers = clipper::data::alternative_monomer(msugAlpha.type().trim());
            
            msugAlpha.set_type(alternative_monomers[0]);
            msugBravo.set_type(alternative_monomers[1]);

            tempGlycanAlpha.replace_sugar_at_index(nodeID, msugAlpha);
            tempGlycanBravo.replace_sugar_at_index(nodeID, msugBravo);
        }
        clipper::String temporaryWURCSAlpha = tempGlycanAlpha.generate_wurcs();
        clipper::String temporaryWURCSBravo = tempGlycanBravo.generate_wurcs();
        
        int valueLocationAlpha = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCSAlpha);
        int valueLocationBravo = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCSBravo);

        // std::cout << "MONOMER PERMUTATION ALPHA: " << temporaryWURCSAlpha << std::endl << "valueLocation = " << valueLocationAlpha << std::endl;
        // std::cout << "MONOMER PERMUTATION BRAVO: " << temporaryWURCSBravo << std::endl << "valueLocation = " << valueLocationBravo << std::endl;

        if (valueLocationAlpha != -1)
            {
                // std::cout << "MONOMER PERMUTATION ALPHA: " << temporaryWURCSAlpha << std::endl << "valueLocation = " << valueLocationAlpha << std::endl;
                if (jsonObject[valueLocationAlpha]["glyconnect"] != "NotFound") glyConnectTrueAlpha = true;
                // std::cout << "GlyConnect status " << std::boolalpha << glyConnectTrueAlpha << std::endl;
            }
        
        if (valueLocationBravo != -1)
            {
                // std::cout << "MONOMER PERMUTATION BRAVO: " << temporaryWURCSBravo << std::endl << "valueLocation = " << valueLocationBravo << std::endl;
                if (jsonObject[valueLocationBravo]["glyconnect"] != "NotFound") glyConnectTrueBravo = true;
                // std::cout << "GlyConnect status " << std::boolalpha << glyConnectTrueBravo << std::endl;
            }

        if (glyConnectTrueAlpha)
            {
                residuePermutationsAlpha = totalCombinations[i].size();
                std::vector<int> mutations(3);
                mutations[0] = anomerPermutationsAlpha;
                mutations[1] = residuePermutationsAlpha;
                mutations[2] = residueDeletions;
                std::pair <clipper::MGlycan,std::vector<int>> tempPair;
                tempPair = std::make_pair(tempGlycanAlpha, mutations);
                result.push_back(tempPair);
            } 

        if (glyConnectTrueBravo)
            {
                residuePermutationsBravo = totalCombinations[i].size();
                std::vector<int> mutations(3);
                mutations[0] = anomerPermutationsBravo;
                mutations[1] = residuePermutationsBravo;
                mutations[2] = residueDeletions;
                std::pair <clipper::MGlycan,std::vector<int>> tempPair;
                tempPair = std::make_pair(tempGlycanBravo, mutations);
                result.push_back(tempPair);
            }
        
        std::vector < clipper::MSugar > sugarsAlpha = tempGlycanAlpha.get_sugars();
        std::vector<int> editable_nodes_for_anomer_permutations_Alpha = get_editable_node_list_for_anomer_permutations(sugarsAlpha);

        std::vector < clipper::MSugar > sugarsBravo = tempGlycanBravo.get_sugars();
        std::vector<int> editable_nodes_for_anomer_permutations_Bravo = get_editable_node_list_for_anomer_permutations(sugarsBravo);

        generate_all_anomer_permutations(result, editable_nodes_for_anomer_permutations_Alpha, tempGlycanAlpha, jsonObject, residuePermutationsAlpha, residueDeletions);
        generate_all_anomer_permutations(result, editable_nodes_for_anomer_permutations_Bravo, tempGlycanBravo, jsonObject, residuePermutationsBravo, residueDeletions);
    }
}


std::vector<std::vector<int>> generate_all_possible_index_combinations(std::vector<int>& editable_node_list)
{

    std::vector<std::vector<int>> result;
    std::vector<int> temp;

    int totalEditableNodes = editable_node_list.size();
    bool checkedEditableNode[totalEditableNodes];
    memset( checkedEditableNode, false, totalEditableNodes*sizeof(bool) );

    
    for(int i = 1; i <= totalEditableNodes; i++)
    {
        CombinationGenerator(editable_node_list, i, 0, 0, checkedEditableNode, totalEditableNodes, result, temp);
    }
    
    return result;
}

void CombinationGenerator(std::vector<int> editable_node_list, int reqLength, int shifter, int currLength, bool checkedEditableNode[], int totalEditableNodes, std::vector<std::vector<int>>& result, std::vector<int>& temp)
{
   if(currLength > reqLength)
   return;
   else if (currLength == reqLength) {
      for (int i = 0; i < totalEditableNodes; i++) {
         if (checkedEditableNode[i] == true) {
         	temp.push_back(editable_node_list[i]);
         }
      }
      result.push_back(temp);
      temp.clear();
      return;
   }
   if (shifter == totalEditableNodes) {
      return;
   }
   checkedEditableNode[shifter] = true;
   CombinationGenerator(editable_node_list, reqLength, shifter + 1, currLength + 1, checkedEditableNode, totalEditableNodes, result, temp);
   checkedEditableNode[shifter] = false;
   CombinationGenerator(editable_node_list, reqLength, shifter + 1, currLength, checkedEditableNode, totalEditableNodes, result, temp);
}

clipper::MGlycan remove_first_leaf_node( clipper::MGlycan inputglycan )
{
    for(int i = 0; i < inputglycan.number_of_nodes(); i++)
        {
            int numberOfConnectionsFromSugar = inputglycan.get_number_of_connections_at_index(i);
            if(numberOfConnectionsFromSugar == 0)
            {
                inputglycan.remove_node_at_index(i);
                return inputglycan;
            }
        }
        return inputglycan;
}

clipper::MGlycan remove_last_node( clipper::MGlycan inputglycan )
{
    int numberOfNodes = inputglycan.number_of_nodes();
    int lastNode = numberOfNodes - 1;

    inputglycan.remove_node_at_index(lastNode);

    return inputglycan;
}


// clipper::MGlycan get_shorter_matching_fragment(clipper::MGlycan& fullglycan, nlohmann::json& jsonObject)
// {
//     clipper::MGlycan tempGlycan = fullglycan;
//     int totalNodes = tempGlycan.number_of_nodes();

//     std::vector<clipper::MGlycan> truncatedGlycanFragments;

//     clipper::MGlycan deletedLastNodeGlycan = shorten_fragment_by_removing_last_node(totalNodes, tempGlycan, jsonObject);
//     truncatedGlycanFragments.push_back(deletedLastNodeGlycan);

//     clipper::MGlycan deletedLeafNodeGlycan = shorten_fragment_by_removing_first_leaf_node(tempGlycan, jsonObject);
//     truncatedGlycanFragments.push_back(deletedLeafNodeGlycan);

//     std::sort(truncatedGlycanFragments.begin(), truncatedGlycanFragments.end(), [](clipper::MGlycan & one, clipper::MGlycan & two) { return one.number_of_nodes() > two.number_of_nodes(); } ); 

//     if (truncatedGlycanFragments[0].number_of_nodes() > 0) 
//     {
//         for(int i = 0; i < truncatedGlycanFragments.size(); i++)
//         {
//             if(truncatedGlycanFragments[i].number_of_nodes() == truncatedGlycanFragments[0].number_of_nodes())
//             {
//                 clipper::String temporaryWURCS = truncatedGlycanFragments[i].generate_wurcs();
//                 int valueLocation = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCS);
//                 std::cout << "Found a closest match after deleting some units:" << std::endl;
//                 std::cout << temporaryWURCS << std::endl;
//                 print_output_from_database(jsonObject, valueLocation, truncatedGlycanFragments[i]);
//             }
//         }
//     }
//     else std::cout << "ERROR: Still unable to find a matching GlyTouCanID for WURCS sequence from this Glycan sequence!" << std::endl;

//     return truncatedGlycanFragments[0];
// }

// clipper::MGlycan shorten_fragment_by_removing_last_node(int totalNodes, clipper::MGlycan inputglycan, nlohmann::json& jsonObject)
// {
//         clipper::MGlycan empty;
//         clipper::String temporaryWURCS;
//         int valueLocation = -1;
//         bool glyConnectTrue = false;

//         for(int i = (totalNodes - 1); i > 0; i--)
//         {
//             inputglycan.remove_node_at_index(i);
            
//             temporaryWURCS = inputglycan.generate_wurcs();

//             valueLocation = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCS);
//             if (valueLocation != -1)
//             {
//                 if (jsonObject[valueLocation]["glyconnect"] != "NotFound") glyConnectTrue = true;
//             }

//             // std::cout << "REMOVING_LAST_NODE: " << temporaryWURCS << std::endl << "valueLocation = " << valueLocation << std::endl;

//             // if (valueLocation != -1) return inputglycan;
//             if (glyConnectTrue) return inputglycan;
//         }
//         if (!glyConnectTrue) return empty;
// }

// clipper::MGlycan shorten_fragment_by_removing_first_leaf_node(clipper::MGlycan inputglycan, nlohmann::json& jsonObject)
// {
//         clipper::MGlycan empty;
//         clipper::String temporaryWURCS;
//         int valueLocation = -1;
//         bool glyConnectTrue = false;
//         int totalNodes = inputglycan.number_of_nodes();

//         for(int i = 0; i < totalNodes; i++)
//         {
//             inputglycan = remove_first_leaf_node(inputglycan);
//             if(inputglycan.number_of_nodes() > 1)
//             {

//                     temporaryWURCS = inputglycan.generate_wurcs();

//                     valueLocation = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCS);

//                     if (valueLocation != -1)
//                     {
//                         if (jsonObject[valueLocation]["glyconnect"] != "NotFound") glyConnectTrue = true;
//                     }

//                     // std::cout << "REMOVING_LEAF_NODE: " << temporaryWURCS << std::endl << "valueLocation = " << valueLocation << std::endl;
                    
//                     // if (valueLocation != -1) return inputglycan;
//                     if (glyConnectTrue) return inputglycan;
//             }
//             else return empty;
//         }
//         if (!glyConnectTrue) return empty;
// }

