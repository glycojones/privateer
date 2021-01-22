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
#define DUMP 1
#define DBG std::cout << "[" << __FUNCTION__ << "] - "


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

// TO DO: Need to redesign this function and its child functions. 
// TO DO: Change the paradigm, so that first all the various glycan chains are generated, vectors populated with potential chains. In this stage, do not query the database for glycan matches.
// TO DO: So, truncated glycans, monomer and anomer permutations done in first stage
// TO DO: In second stage, after the above is done, take the resultant vector and check for entries within the vector itself in the database(in main thread) 
// TO DO: Resize the vector to only contain matches of glycan chains that are found on databases and return that as the final result. 
// TO DO: In essence, glycan mutations first, database checking(on single thread) afterwards.

std::vector<std::pair<clipper::MGlycan, std::vector<int>>> generate_closest_matches_parallel(clipper::MGlycan& fullglycan, nlohmann::json &jsonObject, bool glucose_only, privateer::thread_pool& pool, bool useParallelism)
{
    std::vector<std::pair<clipper::MGlycan, std::vector<int>>> result; // std::vector.push_back() is not thread safe, need to create temporary vectors.

    int totalNodes = fullglycan.number_of_nodes();
    clipper::MGlycan permutatedGlycanLeafNode = fullglycan,
                     permutatedGlycanLastNode = fullglycan;

    int residueDeletions = 0;

    // Return to main thread all the time basically - figure this shit out. 
    for(int i = totalNodes; i > 0; i--)
    {
        if(i == totalNodes)
        {
            int residuePermutations = 0;
        
            std::vector < clipper::MSugar > sugarsLeafNode = permutatedGlycanLeafNode.get_sugars();
            
            std::vector<int> editable_nodes_for_anomer_permutations = get_editable_node_list_for_anomer_permutations(sugarsLeafNode);

            generate_all_anomer_permutations(result, editable_nodes_for_anomer_permutations, permutatedGlycanLeafNode, jsonObject, residuePermutations, residueDeletions);

            std::vector<int> editable_nodes_for_monomer_permutations = get_editable_node_list_for_monomer_permutations(sugarsLeafNode, glucose_only);

            bool statusControl = false;
            std::vector<std::pair<clipper::MGlycan, std::vector<int>>> tempResultMonomerPermutations;
            generate_all_monomer_permutations_parallel_initial(tempResultMonomerPermutations, editable_nodes_for_monomer_permutations, permutatedGlycanLeafNode, jsonObject, residueDeletions, pool, useParallelism, statusControl);

            std::cout << "tempResultMonomerPermutations.size() " << tempResultMonomerPermutations.size() << std::endl;
            std::cout << "Back in the else loop after initiating monomer permutations!" << "\tStatus of statusControl " << std::boolalpha << statusControl << "." << std::endl;
            
            while(!statusControl)
            {
                pool.sync();
            }

            std::cout << "Finished monomer permutations!" << "\tStatus of statusControl " << std::boolalpha << statusControl << "." << std::endl;

            result.insert( result.end(), tempResultMonomerPermutations.begin(), tempResultMonomerPermutations.end() );

            remove_first_leaf_node_and_check_db(result, permutatedGlycanLeafNode, jsonObject, residueDeletions);
            remove_last_node_and_check_db(result, permutatedGlycanLastNode, jsonObject, residueDeletions);

            residueDeletions++;
        }
        else
        {
            int residuePermutations = 0;
        
            std::vector<std::pair<clipper::MGlycan, std::vector<int>>> tempResultMonomerAlpha, tempResultMonomerBravo, tempResultMonomerCharlie, tempResultMonomerDelta;
            std::vector < clipper::MSugar > sugarsLeafNode = permutatedGlycanLeafNode.get_sugars();
            std::vector < clipper::MSugar > sugarsLastNode = permutatedGlycanLastNode.get_sugars();
            
            std::vector<int> editable_leaf_nodes_for_anomer_permutations = get_editable_node_list_for_anomer_permutations(sugarsLeafNode);
            std::vector<int> editable_last_nodes_for_anomer_permutations = get_editable_node_list_for_anomer_permutations(sugarsLastNode);

            // Might even initiate pool launch here. 
            generate_all_anomer_permutations(result, editable_leaf_nodes_for_anomer_permutations, permutatedGlycanLeafNode, jsonObject, residuePermutations, residueDeletions);
            generate_all_anomer_permutations(result, editable_last_nodes_for_anomer_permutations, permutatedGlycanLastNode, jsonObject, residuePermutations, residueDeletions);

            std::vector<int> editable_leaf_nodes_for_monomer_permutations = get_editable_node_list_for_monomer_permutations(sugarsLeafNode, glucose_only);
            std::vector<int> editable_last_nodes_for_monomer_permutations = get_editable_node_list_for_monomer_permutations(sugarsLastNode, glucose_only);

            bool statusControl = false;
            std::vector<std::pair<clipper::MGlycan, std::vector<int>>> tempResultMonomerPermutations;
            generate_all_monomer_permutations_parallel_subsequent(tempResultMonomerPermutations, editable_leaf_nodes_for_monomer_permutations, editable_last_nodes_for_monomer_permutations, permutatedGlycanLeafNode, permutatedGlycanLastNode, jsonObject, residueDeletions, pool, useParallelism, statusControl);
            // (std::vector<std::pair<clipper::MGlycan, std::vector<int>>>& result, std::vector<int>& editable_node_list_leaf, std::vector<int>& editable_node_list_last, clipper::MGlycan glycan_node_removed_leaf, clipper::MGlycan glycan_node_removed_last, nlohmann::json& jsonObject, int residueDeletions, privateer::thread_pool& pool, bool useParallelism, bool& statusControl)
            
            std::cout << "Back in the else loop after initiating monomer permutations!" << "\tStatus of statusControl " << std::boolalpha << statusControl << "." << std::endl;
            
            while(!statusControl)
            {
                pool.sync();
            }

            std::cout << "Finished monomer permutations!" << "\tStatus of statusControl " << std::boolalpha << statusControl << "." << std::endl;

            result.insert( result.end(), tempResultMonomerPermutations.begin(), tempResultMonomerPermutations.end() );

            std::cout << "Size of result = " << result.size() << std::endl;
            
            remove_first_leaf_node_and_check_db(result, permutatedGlycanLeafNode, jsonObject, residueDeletions);
            remove_last_node_and_check_db(result, permutatedGlycanLastNode, jsonObject, residueDeletions);

            residueDeletions++;
            //functions to remove glycans here...And that's how I am going to be using 4 threads max...
        }
    }
    std::cout << "Returning result vector!" << std::endl;
    return result;
}

//last idea: after initiating jobs, check whether indices that were supposed to have shit pushed through are empty or not. after they not empty, only then proceed.
//another idea: revert to multiple containers. Somehow do a long pause after the very last one. 
//Desperate idea: revert to singlethreaded mode when combination numbers become very low. 
//Even more depsperate idea: std::this_thread::sleep_for(5 seconds), but might be best solution. Combine with multiple container idea at the very end. 
void generate_all_monomer_permutations_parallel_initial(std::vector<std::pair<clipper::MGlycan, std::vector<int>>>& result, std::vector<int>& editable_node_list, clipper::MGlycan glycan, nlohmann::json& jsonObject, int residueDeletions, privateer::thread_pool& pool, bool useParallelism, bool& statusControl)
{
    std::vector<std::vector<int>> totalCombinations = generate_all_possible_index_combinations(editable_node_list);

    int totalThreads = pool.size();

    int sleepTimer = 1;

    if(useParallelism)
    {
        std::vector<std::pair<clipper::MGlycan, std::vector<int>>>  threadUnsafeContainerAlpha,
                                                                    threadUnsafeContainerBravo;

        std::vector<std::vector<std::pair<clipper::MGlycan, std::vector<int>>>> threadSafeContainerAlpha(totalCombinations.size()),
                                                                                threadSafeContainerBravo(totalCombinations.size());

        int localIterator = totalThreads;

        for(int globalOffset = 0; globalOffset < totalCombinations.size(); globalOffset+=localIterator)
        {
            while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                pool.sync();

            for(int job = 0; job < totalThreads; job++)
            {
                if( (job + globalOffset) < totalCombinations.size() )
                {
                    pool.push([job, globalOffset, &totalCombinations, &glycan, &jsonObject, residueDeletions, &threadSafeContainerAlpha](int id)
                    {

                        int index = job + globalOffset;

                        #if DUMP
                            std::cout << std::endl;
                            DBG << "Calculating monomer permutations for tempGlycanAlpha from Thread ID: " << id << '.' << std::endl;
                            DBG << "MONOMER PERMUTATION tempGlycanAlpha combination: " << index << std::endl << "/" << totalCombinations.size() << "." << std::endl;
                        #endif

                        std::vector<std::pair<clipper::MGlycan, std::vector<int>>> tempContainer; 
                        
                        clipper::MGlycan tempGlycanAlpha = glycan;

                        bool glyConnectTrue = false;

                        int residuePermutations = 0;
                        
                        int anomerPermutations = 0;

                    
                        for(int j = 0; j < totalCombinations[index].size(); j++)
                        {
                            int nodeID = totalCombinations[index][j];

                            clipper::MSugar msug = tempGlycanAlpha.get_node(nodeID).get_sugar();

                            std::vector<std::string> alternative_monomers = clipper::data::alternative_monomer(msug.type().trim());
                            
                            msug.set_type(alternative_monomers[0]);

                            tempGlycanAlpha.replace_sugar_at_index(nodeID, msug);
                        }
                        
                        clipper::String temporaryWURCS = tempGlycanAlpha.generate_wurcs();
                        
                        int valueLocation = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCS);

                        // std::cout << "MONOMER PERMUTATION ALPHA: " << temporaryWURCS << std::endl << "valueLocation = " << valueLocation << std::endl;

                        if (valueLocation != -1)
                            {
                                if (jsonObject[valueLocation]["glyconnect"] != "NotFound") glyConnectTrue = true;
                            }
                        

                        if (glyConnectTrue)
                            {
                                residuePermutations = totalCombinations[index].size();
                                std::vector<int> mutations(3);
                                mutations[0] = anomerPermutations;
                                mutations[1] = residuePermutations;
                                mutations[2] = residueDeletions;
                                std::pair <clipper::MGlycan,std::vector<int>> tempPair;
                                tempPair = std::make_pair(tempGlycanAlpha, mutations);

                                tempContainer.push_back(tempPair);
                            } 

                        std::vector < clipper::MSugar > sugars = tempGlycanAlpha.get_sugars();
                        std::vector<int> editable_nodes_for_anomer_permutations_Alpha = get_editable_node_list_for_anomer_permutations(sugars);
                    
                        generate_all_anomer_permutations(tempContainer, editable_nodes_for_anomer_permutations_Alpha, tempGlycanAlpha, jsonObject, residuePermutations, residueDeletions);      
                        threadSafeContainerAlpha[index] = tempContainer;
                    });
                }
            }

            while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                pool.sync();

            if((globalOffset + localIterator) >= totalCombinations.size())
            {
                localIterator = totalCombinations.size() - globalOffset;
            }

            if((globalOffset + localIterator) >= totalCombinations.size() || localIterator < 0)
            {
                if(totalCombinations.size() < 65)
                    std::this_thread::sleep_for(std::chrono::seconds(sleepTimer));
                
                while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                    pool.sync();
                
                for(int i = 0; i < threadSafeContainerAlpha.size(); i++)
                {
                    std::vector<std::pair<clipper::MGlycan, std::vector<int>>> tempVector = threadSafeContainerAlpha[i];
                    threadUnsafeContainerAlpha.insert(threadUnsafeContainerAlpha.end(), tempVector.begin(), tempVector.end());
                }
                break;
            }
        }

        
        while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
            pool.sync();

        // while(threadSafeContainerAlpha[(totalCombinations.size() - 1)].empty())
        //     if(!threadSafeContainerAlpha[(totalCombinations.size() - 1)].empty())
        //         result.insert( result.end(), threadUnsafeContainerAlpha.begin(), threadUnsafeContainerAlpha.end());

        // result.insert( result.end(), threadUnsafeContainer.begin(), threadUnsafeContainer.end());

        localIterator = totalThreads;


        for(int globalOffset = 0; globalOffset < totalCombinations.size(); globalOffset+=localIterator)
        {
            while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                pool.sync();

            for(int job = 0; job < totalThreads; job++)
            {
                if( (job + globalOffset) < totalCombinations.size() )
                {
                    pool.push([job, globalOffset, &totalCombinations, &glycan, &jsonObject, residueDeletions, &threadSafeContainerBravo](int id)
                    {

                        int index = job + globalOffset;

                        #if DUMP
                            std::cout << std::endl;
                            DBG << "Calculating monomer permutations for tempGlycanBravo from Thread ID: " << id << '.' << std::endl;
                            DBG << "MONOMER PERMUTATION tempGlycanBravo combination: " << index << std::endl << "/" << totalCombinations.size() << "." << std::endl;
                        #endif

                        std::vector<std::pair<clipper::MGlycan, std::vector<int>>> tempContainer; 
                        
                        clipper::MGlycan tempGlycanBravo = glycan;

                        bool glyConnectTrue = false;

                        int residuePermutations = 0;
                        
                        int anomerPermutations = 0;

                    
                        for(int j = 0; j < totalCombinations[index].size(); j++)
                        {
                            int nodeID = totalCombinations[index][j];

                            clipper::MSugar msug = tempGlycanBravo.get_node(nodeID).get_sugar();

                            std::vector<std::string> alternative_monomers = clipper::data::alternative_monomer(msug.type().trim());
                            
                            msug.set_type(alternative_monomers[1]);

                            tempGlycanBravo.replace_sugar_at_index(nodeID, msug);
                        }
                        
                        clipper::String temporaryWURCS = tempGlycanBravo.generate_wurcs();
                        
                        int valueLocation = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCS);

                        // std::cout << "MONOMER PERMUTATION BRAVO: " << temporaryWURCS << std::endl << "valueLocation = " << valueLocation << std::endl;

                        if (valueLocation != -1)
                            {
                                if (jsonObject[valueLocation]["glyconnect"] != "NotFound") glyConnectTrue = true;
                            }
                        

                        if (glyConnectTrue)
                            {
                                residuePermutations = totalCombinations[index].size();
                                std::vector<int> mutations(3);
                                mutations[0] = anomerPermutations;
                                mutations[1] = residuePermutations;
                                mutations[2] = residueDeletions;
                                std::pair <clipper::MGlycan,std::vector<int>> tempPair;
                                tempPair = std::make_pair(tempGlycanBravo, mutations);

                                tempContainer.push_back(tempPair);
                            } 

                        std::vector < clipper::MSugar > sugars = tempGlycanBravo.get_sugars();
                        std::vector<int> editable_nodes_for_anomer_permutations_Bravo = get_editable_node_list_for_anomer_permutations(sugars);
                    
                        generate_all_anomer_permutations(tempContainer, editable_nodes_for_anomer_permutations_Bravo, tempGlycanBravo, jsonObject, residuePermutations, residueDeletions);      
                        threadSafeContainerBravo[index] = tempContainer;
                    });
                }
            }

            while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                pool.sync();

            if((globalOffset + localIterator) >= totalCombinations.size())
            {
                localIterator = totalCombinations.size() - globalOffset;
            }

            if((globalOffset + localIterator) >= totalCombinations.size() || localIterator < 0)
            {
                if(totalCombinations.size() < 65)
                    std::this_thread::sleep_for(std::chrono::seconds(sleepTimer));
                
                while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                    pool.sync();
                for(int i = 0; i < threadSafeContainerBravo.size(); i++)
                {
                    std::vector<std::pair<clipper::MGlycan, std::vector<int>>> tempVector = threadSafeContainerBravo[i];
                    threadUnsafeContainerBravo.insert(threadUnsafeContainerBravo.end(), tempVector.begin(), tempVector.end());
                }
                break;
            }
        }

        while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
            pool.sync();

        // while(threadSafeContainerBravo[(totalCombinations.size() - 1)].empty())
        //     if(!threadSafeContainerBravo[(totalCombinations.size() - 1)].empty())
        //         result.insert( result.end(), threadUnsafeContainerBravo.begin(), threadUnsafeContainerBravo.end());

        result.insert( result.end(), threadUnsafeContainerAlpha.begin(), threadUnsafeContainerAlpha.end());
        result.insert( result.end(), threadUnsafeContainerBravo.begin(), threadUnsafeContainerBravo.end());

    }
    
    while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
        {
            pool.sync();
            // statusControl = true;
        }
    statusControl = true;
}

void generate_all_monomer_permutations_parallel_subsequent(std::vector<std::pair<clipper::MGlycan, std::vector<int>>>& result, std::vector<int>& editable_node_list_leaf, std::vector<int>& editable_node_list_last, clipper::MGlycan glycan_node_removed_leaf, clipper::MGlycan glycan_node_removed_last, nlohmann::json& jsonObject, int residueDeletions, privateer::thread_pool& pool, bool useParallelism, bool& statusControl)
{

    std::vector<std::vector<int>> totalCombinationsNodeRemovedLeaf = generate_all_possible_index_combinations(editable_node_list_leaf);
    std::vector<std::vector<int>> totalCombinationsNodeRemovedLast = generate_all_possible_index_combinations(editable_node_list_last);

    int totalThreads = pool.size();

    int sleepTimer = 1;
    if(useParallelism)
    {
        std::vector<std::pair<clipper::MGlycan, std::vector<int>>>  threadUnsafeContainerLeafAlpha,
                                                                    threadUnsafeContainerLeafBravo,
                                                                    threadUnsafeContainerLastAlpha,
                                                                    threadUnsafeContainerLastBravo;
        
        std::vector<std::vector<std::pair<clipper::MGlycan, std::vector<int>>>> threadSafeContainerLeafAlpha(totalCombinationsNodeRemovedLeaf.size()),
                                                                                threadSafeContainerLeafBravo(totalCombinationsNodeRemovedLeaf.size()),
                                                                                threadSafeContainerLastAlpha(totalCombinationsNodeRemovedLast.size()),
                                                                                threadSafeContainerLastBravo(totalCombinationsNodeRemovedLast.size());

        if(totalThreads > totalCombinationsNodeRemovedLeaf.size())
        {
            totalThreads = totalCombinationsNodeRemovedLeaf.size();
        }
        
        int localIterator = totalThreads;
        for(int globalOffset = 0; globalOffset < totalCombinationsNodeRemovedLeaf.size(); globalOffset+=localIterator)
        {
            while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                pool.sync();

            for(int job = 0; job < totalThreads; job++)
            {
                if( (job + globalOffset) < totalCombinationsNodeRemovedLeaf.size() )
                {
                    pool.push([job, globalOffset, &totalCombinationsNodeRemovedLeaf, &glycan_node_removed_leaf, &jsonObject, residueDeletions, &threadSafeContainerLeafAlpha](int id)
                    {

                        int index = job + globalOffset;

                        #if DUMP
                            std::cout << std::endl;
                            DBG << "Calculating monomer permutations for tempGlycanAlpha_LEAF from Thread ID: " << id << '.' << std::endl;
                            DBG << "MONOMER PERMUTATION tempGlycanAlpha_LEAF combination: " << index << std::endl << "/" << totalCombinationsNodeRemovedLeaf.size() << "." << std::endl;
                        #endif

                        std::vector<std::pair<clipper::MGlycan, std::vector<int>>> tempContainer; 
                        
                        clipper::MGlycan tempGlycanAlpha = glycan_node_removed_leaf;

                        bool glyConnectTrue = false;

                        int residuePermutations = 0;
                        
                        int anomerPermutations = 0;

                    
                        for(int j = 0; j < totalCombinationsNodeRemovedLeaf[index].size(); j++)
                        {
                            int nodeID = totalCombinationsNodeRemovedLeaf[index][j];

                            clipper::MSugar msug = tempGlycanAlpha.get_node(nodeID).get_sugar();

                            std::vector<std::string> alternative_monomers = clipper::data::alternative_monomer(msug.type().trim());
                            
                            msug.set_type(alternative_monomers[0]);

                            tempGlycanAlpha.replace_sugar_at_index(nodeID, msug);
                        }
                        
                        clipper::String temporaryWURCS = tempGlycanAlpha.generate_wurcs();
                        
                        int valueLocation = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCS);

                        // std::cout << "MONOMER PERMUTATION ALPHA: " << temporaryWURCS << std::endl << "valueLocation = " << valueLocation << std::endl;

                        if (valueLocation != -1)
                            {
                                if (jsonObject[valueLocation]["glyconnect"] != "NotFound") glyConnectTrue = true;
                            }
                        

                        if (glyConnectTrue)
                            {
                                residuePermutations = totalCombinationsNodeRemovedLeaf[index].size();
                                std::vector<int> mutations(3);
                                mutations[0] = anomerPermutations;
                                mutations[1] = residuePermutations;
                                mutations[2] = residueDeletions;
                                std::pair <clipper::MGlycan,std::vector<int>> tempPair;
                                tempPair = std::make_pair(tempGlycanAlpha, mutations);

                                tempContainer.push_back(tempPair);
                            } 

                        std::vector < clipper::MSugar > sugars = tempGlycanAlpha.get_sugars();
                        std::vector<int> editable_nodes_for_anomer_permutations_Alpha = get_editable_node_list_for_anomer_permutations(sugars);
                    
                        generate_all_anomer_permutations(tempContainer, editable_nodes_for_anomer_permutations_Alpha, tempGlycanAlpha, jsonObject, residuePermutations, residueDeletions);      
                        threadSafeContainerLeafAlpha[index] = tempContainer;
                    });
                }
            }

            while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                pool.sync();

            if((globalOffset + localIterator) >= totalCombinationsNodeRemovedLeaf.size())
            {
                localIterator = totalCombinationsNodeRemovedLeaf.size() - globalOffset;
            }

            if((globalOffset + localIterator) >= totalCombinationsNodeRemovedLeaf.size() || localIterator < 0)
            {
                if(totalCombinationsNodeRemovedLeaf.size() < 65)
                    std::this_thread::sleep_for(std::chrono::seconds(sleepTimer));
                
                while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                    pool.sync();

                for(int i = 0; i < threadSafeContainerLeafAlpha.size(); i++)
                {
                    std::vector<std::pair<clipper::MGlycan, std::vector<int>>> tempVector = threadSafeContainerLeafAlpha[i];
                    threadUnsafeContainerLeafAlpha.insert(threadUnsafeContainerLeafAlpha.end(), tempVector.begin(), tempVector.end());
                }
                break;
            }
        }

        
        while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
            pool.sync();
        
        
        // while(threadSafeContainerLeafAlpha[(totalCombinationsNodeRemovedLeaf.size() - 1)].empty())
        //     if(!threadSafeContainerLeafAlpha[(totalCombinationsNodeRemovedLeaf.size() - 1)].empty())
        //         result.insert( result.end(), threadUnsafeContainerLeaf.begin(), threadUnsafeContainerLeaf.end());

        localIterator = totalThreads;


        for(int globalOffset = 0; globalOffset < totalCombinationsNodeRemovedLeaf.size(); globalOffset+=localIterator)
        {
            while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                pool.sync();

            for(int job = 0; job < totalThreads; job++)
            {
                if( (job + globalOffset) < totalCombinationsNodeRemovedLeaf.size() )
                {
                    pool.push([job, globalOffset, &totalCombinationsNodeRemovedLeaf, &glycan_node_removed_leaf, &jsonObject, residueDeletions, &threadSafeContainerLeafBravo](int id)
                    {

                        int index = job + globalOffset;

                        #if DUMP
                            std::cout << std::endl;
                            DBG << "Calculating monomer permutations for tempGlycanBravo_LEAF from Thread ID: " << id << '.' << std::endl;
                            DBG << "MONOMER PERMUTATION tempGlycanBravo_LEAF combination: " << index << std::endl << "/" << totalCombinationsNodeRemovedLeaf.size() << "." << std::endl;
                        #endif

                        std::vector<std::pair<clipper::MGlycan, std::vector<int>>> tempContainer; 
                        
                        clipper::MGlycan tempGlycanBravo = glycan_node_removed_leaf;

                        bool glyConnectTrue = false;

                        int residuePermutations = 0;
                        
                        int anomerPermutations = 0;

                    
                        for(int j = 0; j < totalCombinationsNodeRemovedLeaf[index].size(); j++)
                        {
                            int nodeID = totalCombinationsNodeRemovedLeaf[index][j];

                            clipper::MSugar msug = tempGlycanBravo.get_node(nodeID).get_sugar();

                            std::vector<std::string> alternative_monomers = clipper::data::alternative_monomer(msug.type().trim());
                            
                            msug.set_type(alternative_monomers[1]);

                            tempGlycanBravo.replace_sugar_at_index(nodeID, msug);
                        }
                        
                        clipper::String temporaryWURCS = tempGlycanBravo.generate_wurcs();
                        
                        int valueLocation = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCS);

                        // std::cout << "MONOMER PERMUTATION Bravo: " << temporaryWURCS << std::endl << "valueLocation = " << valueLocation << std::endl;

                        if (valueLocation != -1)
                            {
                                if (jsonObject[valueLocation]["glyconnect"] != "NotFound") glyConnectTrue = true;
                            }
                        

                        if (glyConnectTrue)
                            {
                                residuePermutations = totalCombinationsNodeRemovedLeaf[index].size();
                                std::vector<int> mutations(3);
                                mutations[0] = anomerPermutations;
                                mutations[1] = residuePermutations;
                                mutations[2] = residueDeletions;
                                std::pair <clipper::MGlycan,std::vector<int>> tempPair;
                                tempPair = std::make_pair(tempGlycanBravo, mutations);

                                tempContainer.push_back(tempPair);
                            } 

                        std::vector < clipper::MSugar > sugars = tempGlycanBravo.get_sugars();
                        std::vector<int> editable_nodes_for_anomer_permutations_Bravo = get_editable_node_list_for_anomer_permutations(sugars);
                    
                        generate_all_anomer_permutations(tempContainer, editable_nodes_for_anomer_permutations_Bravo, tempGlycanBravo, jsonObject, residuePermutations, residueDeletions);      
                        threadSafeContainerLeafBravo[index] = tempContainer;
                    });
                }
            }

            while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                pool.sync();

            if((globalOffset + localIterator) >= totalCombinationsNodeRemovedLeaf.size())
            {
                localIterator = totalCombinationsNodeRemovedLeaf.size() - globalOffset;
            }

            if((globalOffset + localIterator) >= totalCombinationsNodeRemovedLeaf.size() || localIterator < 0)
            {
                if(totalCombinationsNodeRemovedLeaf.size() < 65)
                    std::this_thread::sleep_for(std::chrono::seconds(sleepTimer));
                
                while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                    pool.sync();
                
                for(int i = 0; i < threadSafeContainerLeafBravo.size(); i++)
                {
                    std::vector<std::pair<clipper::MGlycan, std::vector<int>>> tempVector = threadSafeContainerLeafBravo[i];
                    threadUnsafeContainerLeafBravo.insert(threadUnsafeContainerLeafBravo.end(), tempVector.begin(), tempVector.end());
                }
                break;
            }
        }

        
        while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
            pool.sync();
        
        // while(threadSafeContainerLeafBravo[(totalCombinationsNodeRemovedLeaf.size() - 1)].empty())
        //     if(!threadSafeContainerLeafBravo[(totalCombinationsNodeRemovedLeaf.size() - 1)].empty())
        //         result.insert( result.end(), threadUnsafeContainerLeafBravo.begin(), threadUnsafeContainerLeafBravo.end());
        
        // result.insert( result.end(), threadUnsafeContainerLeaf.begin(), threadUnsafeContainerLeaf.end());

        localIterator = totalThreads;

        if(pool.size() > totalCombinationsNodeRemovedLeaf.size())
        {
            totalThreads = totalCombinationsNodeRemovedLeaf.size();
        }
        
        for(int globalOffset = 0; globalOffset < totalCombinationsNodeRemovedLast.size(); globalOffset+=localIterator)
        {
            while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                pool.sync();

            for(int job = 0; job < totalThreads; job++)
            {
                if( (job + globalOffset) < totalCombinationsNodeRemovedLast.size() )
                {
                    pool.push([job, globalOffset, &totalCombinationsNodeRemovedLast, &glycan_node_removed_last, &jsonObject, residueDeletions, &threadSafeContainerLastAlpha](int id)
                    {

                        int index = job + globalOffset;

                        #if DUMP
                            std::cout << std::endl;
                            DBG << "Calculating monomer permutations for tempGlycanAlpha_LAST from Thread ID: " << id << '.' << std::endl;
                            DBG << "MONOMER PERMUTATION tempGlycanAlpha_LAST combination: " << index << std::endl << "/" << totalCombinationsNodeRemovedLast.size() << "." << std::endl;
                        #endif

                        std::vector<std::pair<clipper::MGlycan, std::vector<int>>> tempContainer; 
                        
                        clipper::MGlycan tempGlycanAlpha = glycan_node_removed_last;

                        bool glyConnectTrue = false;

                        int residuePermutations = 0;
                        
                        int anomerPermutations = 0;

                    
                        for(int j = 0; j < totalCombinationsNodeRemovedLast[index].size(); j++)
                        {
                            int nodeID = totalCombinationsNodeRemovedLast[index][j];

                            clipper::MSugar msug = tempGlycanAlpha.get_node(nodeID).get_sugar();

                            std::vector<std::string> alternative_monomers = clipper::data::alternative_monomer(msug.type().trim());
                            
                            msug.set_type(alternative_monomers[0]);

                            tempGlycanAlpha.replace_sugar_at_index(nodeID, msug);
                        }
                        
                        clipper::String temporaryWURCS = tempGlycanAlpha.generate_wurcs();
                        
                        int valueLocation = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCS);

                        // std::cout << "MONOMER PERMUTATION ALPHA: " << temporaryWURCS << std::endl << "valueLocation = " << valueLocation << std::endl;

                        if (valueLocation != -1)
                            {
                                if (jsonObject[valueLocation]["glyconnect"] != "NotFound") glyConnectTrue = true;
                            }
                        

                        if (glyConnectTrue)
                            {
                                residuePermutations = totalCombinationsNodeRemovedLast[index].size();
                                std::vector<int> mutations(3);
                                mutations[0] = anomerPermutations;
                                mutations[1] = residuePermutations;
                                mutations[2] = residueDeletions;
                                std::pair <clipper::MGlycan,std::vector<int>> tempPair;
                                tempPair = std::make_pair(tempGlycanAlpha, mutations);

                                tempContainer.push_back(tempPair);
                            } 

                        std::vector < clipper::MSugar > sugars = tempGlycanAlpha.get_sugars();
                        std::vector<int> editable_nodes_for_anomer_permutations_Alpha = get_editable_node_list_for_anomer_permutations(sugars);
                    
                        generate_all_anomer_permutations(tempContainer, editable_nodes_for_anomer_permutations_Alpha, tempGlycanAlpha, jsonObject, residuePermutations, residueDeletions);      
                        threadSafeContainerLastAlpha[index] = tempContainer;
                    });
                }
            }

            while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                pool.sync();

            if((globalOffset + localIterator) >= totalCombinationsNodeRemovedLast.size())
            {
                localIterator = totalCombinationsNodeRemovedLast.size() - globalOffset;
            }

            if((globalOffset + localIterator) >= totalCombinationsNodeRemovedLast.size() || localIterator < 0)
            {
                if(totalCombinationsNodeRemovedLast.size() < 65)
                    std::this_thread::sleep_for(std::chrono::seconds(sleepTimer));

                while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                    pool.sync();

                for(int i = 0; i < threadSafeContainerLastAlpha.size(); i++)
                {
                    std::vector<std::pair<clipper::MGlycan, std::vector<int>>> tempVector = threadSafeContainerLastAlpha[i];
                    threadUnsafeContainerLastAlpha.insert(threadUnsafeContainerLastAlpha.end(), tempVector.begin(), tempVector.end());
                }
                break;
            }
        }

        
        while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
            pool.sync();

        // while(threadSafeContainerLastAlpha[(totalCombinationsNodeRemovedLast.size() - 1)].empty())
        //     if(!threadSafeContainerLastAlpha[(totalCombinationsNodeRemovedLast.size() - 1)].empty())
        //         result.insert( result.end(), threadUnsafeContainerLastAlpha.begin(), threadUnsafeContainerLastAlpha.end());

        // result.insert( result.end(), threadUnsafeContainerLast.begin(), threadUnsafeContainerLast.end());

        localIterator = totalThreads;


        for(int globalOffset = 0; globalOffset < totalCombinationsNodeRemovedLast.size(); globalOffset+=localIterator)
        {
            while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                pool.sync();

            for(int job = 0; job < totalThreads; job++)
            {
                if( (job + globalOffset) < totalCombinationsNodeRemovedLast.size() )
                {
                    pool.push([job, globalOffset, &totalCombinationsNodeRemovedLast, &glycan_node_removed_last, &jsonObject, residueDeletions, &threadSafeContainerLastBravo](int id)
                    {

                        int index = job + globalOffset;

                        #if DUMP
                            std::cout << std::endl;
                            DBG << "Calculating monomer permutations for tempGlycanBravo_LAST from Thread ID: " << id << '.' << std::endl;
                            DBG << "MONOMER PERMUTATION tempGlycanBravo_LAST combination: " << index << std::endl << "/" << totalCombinationsNodeRemovedLast.size() << "." << std::endl;
                        #endif

                        std::vector<std::pair<clipper::MGlycan, std::vector<int>>> tempContainer; 
                        
                        clipper::MGlycan tempGlycanBravo = glycan_node_removed_last;

                        bool glyConnectTrue = false;

                        int residuePermutations = 0;
                        
                        int anomerPermutations = 0;

                    
                        for(int j = 0; j < totalCombinationsNodeRemovedLast[index].size(); j++)
                        {
                            int nodeID = totalCombinationsNodeRemovedLast[index][j];

                            clipper::MSugar msug = tempGlycanBravo.get_node(nodeID).get_sugar();

                            std::vector<std::string> alternative_monomers = clipper::data::alternative_monomer(msug.type().trim());
                            
                            msug.set_type(alternative_monomers[1]);

                            tempGlycanBravo.replace_sugar_at_index(nodeID, msug);
                        }
                        
                        clipper::String temporaryWURCS = tempGlycanBravo.generate_wurcs();
                        
                        int valueLocation = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCS);

                        // std::cout << "MONOMER PERMUTATION Bravo: " << temporaryWURCS << std::endl << "valueLocation = " << valueLocation << std::endl;

                        if (valueLocation != -1)
                            {
                                if (jsonObject[valueLocation]["glyconnect"] != "NotFound") glyConnectTrue = true;
                            }
                        

                        if (glyConnectTrue)
                            {
                                residuePermutations = totalCombinationsNodeRemovedLast[index].size();
                                std::vector<int> mutations(3);
                                mutations[0] = anomerPermutations;
                                mutations[1] = residuePermutations;
                                mutations[2] = residueDeletions;
                                std::pair <clipper::MGlycan,std::vector<int>> tempPair;
                                tempPair = std::make_pair(tempGlycanBravo, mutations);

                                tempContainer.push_back(tempPair);
                            } 

                        std::vector < clipper::MSugar > sugars = tempGlycanBravo.get_sugars();
                        std::vector<int> editable_nodes_for_anomer_permutations_Bravo = get_editable_node_list_for_anomer_permutations(sugars);
                    
                        generate_all_anomer_permutations(tempContainer, editable_nodes_for_anomer_permutations_Bravo, tempGlycanBravo, jsonObject, residuePermutations, residueDeletions);      
                        threadSafeContainerLastBravo[index] = tempContainer;
                    });
                }
            }

            while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                pool.sync();

            if((globalOffset + localIterator) >= totalCombinationsNodeRemovedLast.size())
            {
                localIterator = totalCombinationsNodeRemovedLast.size() - globalOffset;
            }

            if((globalOffset + localIterator) >= totalCombinationsNodeRemovedLast.size() || localIterator < 0)
            {
                if(totalCombinationsNodeRemovedLast.size() < 65)
                    std::this_thread::sleep_for(std::chrono::seconds(sleepTimer));

                while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
                    pool.sync();
                
                for(int i = 0; i < threadSafeContainerLastBravo.size(); i++)
                {
                    std::vector<std::pair<clipper::MGlycan, std::vector<int>>> tempVector = threadSafeContainerLastBravo[i];
                    threadUnsafeContainerLastBravo.insert(threadUnsafeContainerLastBravo.end(), tempVector.begin(), tempVector.end());
                }
                break;
            }
        }

        
        while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
            pool.sync();
        
        // while(threadSafeContainerLastBravo[(totalCombinationsNodeRemovedLast.size() - 1)].empty())
        //     if(!threadSafeContainerLastBravo[(totalCombinationsNodeRemovedLast.size() - 1)].empty())
        //         result.insert( result.end(), threadUnsafeContainerLastBravo.begin(), threadUnsafeContainerLastBravo.end());

        result.insert( result.end(), threadUnsafeContainerLeafAlpha.begin(), threadUnsafeContainerLeafAlpha.end());
        result.insert( result.end(), threadUnsafeContainerLeafBravo.begin(), threadUnsafeContainerLeafBravo.end());
        result.insert( result.end(), threadUnsafeContainerLastAlpha.begin(), threadUnsafeContainerLastAlpha.end());
        result.insert( result.end(), threadUnsafeContainerLastBravo.begin(), threadUnsafeContainerLastBravo.end());
    }
    
    while(pool.n_idle() != pool.size() || pool.n_remaining_jobs() > 0)
        {
            pool.sync();
            // statusControl = true;
        }
        statusControl = true;
}

std::vector<std::pair<clipper::MGlycan, std::vector<int>>> generate_closest_matches_singlethreaded(clipper::MGlycan& fullglycan, nlohmann::json& jsonObject, bool glucose_only)
{
    std::vector<std::pair<clipper::MGlycan, std::vector<int>>> result;


    int totalNodes = fullglycan.number_of_nodes();
    clipper::MGlycan tempGlycan = fullglycan;

    int residueDeletions = 0;
    for(int i = 0; i < totalNodes; i++)
    {
        int residuePermutations = 0;

        std::vector < clipper::MSugar > sugars = tempGlycan.get_sugars();
        std::vector<int> editable_nodes_for_anomer_permutations = get_editable_node_list_for_anomer_permutations(sugars);

        generate_all_anomer_permutations(result, editable_nodes_for_anomer_permutations, tempGlycan, jsonObject, residuePermutations, residueDeletions);

        std::vector<int> editable_nodes_for_monomer_permutations = get_editable_node_list_for_monomer_permutations(sugars, glucose_only);

        generate_all_monomer_permutations_singlethreaded(result, editable_nodes_for_monomer_permutations, tempGlycan, jsonObject, residueDeletions);

        if(tempGlycan.number_of_nodes() > 1)
        {
            bool glyConnectTrue = false;
            
            tempGlycan = remove_first_leaf_node(tempGlycan);
            residueDeletions++;

            clipper::String temporaryWURCS = tempGlycan.generate_wurcs();

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
                    tempPair = std::make_pair(tempGlycan, mutations);
                    result.push_back(tempPair);
                }
        }
    }

    // do the score calculations here, return the final vector already containing those float values. 

    return result;
}

void generate_all_monomer_permutations_singlethreaded(std::vector<std::pair<clipper::MGlycan, std::vector<int>>>& result, std::vector<int>& editable_node_list, clipper::MGlycan glycan, nlohmann::json& jsonObject, int residueDeletions)
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


void generate_all_anomer_permutations(std::vector<std::pair<clipper::MGlycan, std::vector<int>>>& result, std::vector<int>& editable_node_list, clipper::MGlycan glycan, nlohmann::json& jsonObject, int residuePermuations, int residueDeletions)
{
    std::vector<std::vector<int>> totalCombinations = generate_all_possible_index_combinations(editable_node_list);

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

void remove_first_leaf_node_and_check_db(std::vector<std::pair<clipper::MGlycan, std::vector<int>>>& result, clipper::MGlycan& inputglycan, nlohmann::json& jsonObject, int residueDeletions)
{
    bool glyConnectTrue = false;
    inputglycan = remove_first_leaf_node(inputglycan);
    
    clipper::String temporaryWURCS = inputglycan.generate_wurcs();
    
    int valueLocation = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCS);

    // std::cout << "MONOMER PERMUTATION ALPHA: " << temporaryWURCS << std::endl << "valueLocation = " << valueLocation << std::endl;

    if (valueLocation != -1)
        {
            if (jsonObject[valueLocation]["glyconnect"] != "NotFound") glyConnectTrue = true;
        }
    

    if (glyConnectTrue)
        {
            std::vector<int> mutations(3);
            mutations[0] = 0;
            mutations[1] = 0;
            mutations[2] = residueDeletions;
            std::pair <clipper::MGlycan,std::vector<int>> tempPair;
            tempPair = std::make_pair(inputglycan, mutations);
            result.push_back(tempPair);
        } 
}

void remove_last_node_and_check_db(std::vector<std::pair<clipper::MGlycan, std::vector<int>>>& result, clipper::MGlycan& inputglycan, nlohmann::json& jsonObject, int residueDeletions)
{
    bool glyConnectTrue = false;
    inputglycan = remove_last_node(inputglycan);
    
    clipper::String temporaryWURCS = inputglycan.generate_wurcs();
    
    int valueLocation = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCS);

    // std::cout << "MONOMER PERMUTATION ALPHA: " << temporaryWURCS << std::endl << "valueLocation = " << valueLocation << std::endl;

    if (valueLocation != -1)
        {
            if (jsonObject[valueLocation]["glyconnect"] != "NotFound") glyConnectTrue = true;
        }
    

    if (glyConnectTrue)
        {
            std::vector<int> mutations(3);
            mutations[0] = 0;
            mutations[1] = 0;
            mutations[2] = residueDeletions;
            std::pair <clipper::MGlycan,std::vector<int>> tempPair;
            tempPair = std::make_pair(inputglycan, mutations);
            result.push_back(tempPair);
        }
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

