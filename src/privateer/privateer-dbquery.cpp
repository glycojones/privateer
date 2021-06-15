// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York



#include "privateer-dbquery.h"

#define DBG std::cout << "[" << __FUNCTION__ << "] - "

void output_dbquery(std::vector<privateer::json::Database>& glycomics_database, clipper::String glycanWURCS, clipper::MGlycan &currentGlycan, bool closest_match_disable, std::vector<std::pair<std::pair<clipper::MGlycan, std::vector<int>>,float>>& finalGlycanPermutationContainer, bool glucose_only, bool debug_output, int sleepTimer, privateer::thread_pool& pool, bool useParallelism)
{
    int valueLocation;
    valueLocation = privateer::util::find_index_of_value_from_wurcs(glycomics_database, glycanWURCS);
    if (valueLocation != -1 && glycomics_database[valueLocation].GlyConnectID != "NotFound") print_output_from_database(glycomics_database, valueLocation, currentGlycan);
    else
    {
        if (valueLocation != -1 && currentGlycan.number_of_nodes() > 1)
        {
            print_output_from_database(glycomics_database, valueLocation, currentGlycan);
            if(closest_match_disable) std::cout << "\nWARNING: Unable to find a matching GlyConnectID for WURCS sequence from this Glycan sequence! If you would like to find the closest match on GlyConnect database please delete the -closest_match_disable option!" << std::endl;
            else
            {
                std::cout << "\nWARNING: Unable to find a matching GlyConnectID for WURCS sequence from this Glycan sequence! Attempting to find the closest matches on GlyConnect database by carrying out permutations" << std::endl;

                std::vector<std::pair<clipper::MGlycan, std::vector<int>>> alternativeGlycans;
                
                if(useParallelism) alternativeGlycans = generate_closest_matches_parallel(currentGlycan, glycomics_database, glucose_only, debug_output, sleepTimer, pool, useParallelism);
                else               alternativeGlycans = generate_closest_matches_singlethreaded(currentGlycan, glycomics_database, glucose_only, debug_output);   
                
                if (!alternativeGlycans.empty()) push_data_to_final_permutation_container(glycomics_database, currentGlycan, alternativeGlycans, finalGlycanPermutationContainer);    
                else std::cout << "ERROR: Unable to generate permutations that would be found in GlyConnect database!" << std::endl;
            }
        }
        else if (valueLocation != -1 && currentGlycan.number_of_nodes() <= 1)
        {
            print_output_from_database(glycomics_database, valueLocation, currentGlycan);
            if(!closest_match_disable) std::cout << "ERROR: Glycan is too short for permutations, therefore unable to find the closest match on GlyConnect." << std::endl;
        }
        else
        {
            if(closest_match_disable) std::cout << "\nWARNING: Unable to find a matching GlyTouCanID for WURCS sequence from this Glycan sequence! If you would like to find the closest match on GlyConnect database please delete the -closest_match_disable option!" << std::endl;
            else
            {
                if ( currentGlycan.number_of_nodes() > 1)
                {
                    std::vector<std::pair<clipper::MGlycan, std::vector<int>>> alternativeGlycans;

                    if(useParallelism) alternativeGlycans = generate_closest_matches_parallel(currentGlycan, glycomics_database, glucose_only, debug_output, sleepTimer, pool, useParallelism);
                    else               alternativeGlycans = generate_closest_matches_singlethreaded(currentGlycan, glycomics_database, glucose_only, debug_output); 

                    if (useParallelism)
                    {
                        
                        if(debug_output)
                        {
                            std::cout << std::endl;
                            DBG << "Number of jobs in the queue: " << pool.n_remaining_jobs() << std::endl;
                        }
                        
                        while(pool.n_remaining_jobs() > 0)
                            pool.sync();
                        
                        if(debug_output)
                        {
                            DBG << "Number of jobs in the queue: " << pool.n_remaining_jobs() << " after sync operation!" << std::endl;
                        }
                    }
                
                    if (!alternativeGlycans.empty()) push_data_to_final_permutation_container(glycomics_database, currentGlycan, alternativeGlycans, finalGlycanPermutationContainer);    
                    else std::cout << "ERROR: Unable to generate permutations that would be found on GlyConnect database!" << std::endl;
                }
                else std::cout << "ERROR: Glycan is too short for permutations, therefore unable to generate alternative GlyTouCan and GlyConnect IDs!" << std::endl;
            }
         }
    }
}


void push_data_to_final_permutation_container(std::vector<privateer::json::Database>& glycomics_database, clipper::MGlycan &currentGlycan, std::vector<std::pair<clipper::MGlycan, std::vector<int>>>& alternativeGlycans, std::vector<std::pair<std::pair<clipper::MGlycan, std::vector<int>>,float>>& finalGlycanPermutationContainer)
{

    std::vector<std::pair<std::pair<clipper::MGlycan, std::vector<int>>,float>> tempGlycanPermutationContainer;

    for (int j = 0; j < alternativeGlycans.size(); j++)
    {
        int originalGlycanLength = currentGlycan.number_of_nodes(),
            currentGlycanLength = alternativeGlycans[j].first.number_of_nodes();

        int anomerPermutations = alternativeGlycans[j].second[0],
            residuePermutations = alternativeGlycans[j].second[1],
            residueDeletions = alternativeGlycans[j].second[2];
        
        float finalScore, maxPermutationScore, currentPermutationScore;

        maxPermutationScore = ( ( (currentGlycanLength * 5) + (currentGlycanLength * 25) + ((originalGlycanLength - 1) * 100) ) / originalGlycanLength );
        currentPermutationScore = ( ( (anomerPermutations * 5) + (residuePermutations * 25) + (residueDeletions * 100) ) / originalGlycanLength );
        finalScore = (currentPermutationScore / maxPermutationScore) * 100;

        auto tempObject = std::make_pair(alternativeGlycans[j], finalScore);
        tempGlycanPermutationContainer.push_back(tempObject);
    }
    
    // sort by permutation score for sensical output to the user.
    std::sort(tempGlycanPermutationContainer.begin(), tempGlycanPermutationContainer.end(), [](std::pair<std::pair<clipper::MGlycan, std::vector<int>>,float> a, std::pair<std::pair<clipper::MGlycan, std::vector<int>>,float> b){
        return a.second < b.second;
    });

    //Make sure that only unique permutations are outputted to the user.
    for(int i = 0; i < tempGlycanPermutationContainer.size(); i++)
    {
        auto tempObject = tempGlycanPermutationContainer[i];
        clipper::MGlycan currentGlycan = tempGlycanPermutationContainer[i].first.first;
        clipper::String currentGlycanWURCS = currentGlycan.generate_wurcs();
        
        if(finalGlycanPermutationContainer.empty())
        {
            finalGlycanPermutationContainer.push_back(tempObject);
            continue;
        }
            
        for(int j = 0; j < finalGlycanPermutationContainer.size(); j++)
        {
            clipper::MGlycan currentGlycanInFinalContainer = finalGlycanPermutationContainer[j].first.first;
            clipper::String currentGlycanWURCSInFinalContainer = currentGlycanInFinalContainer.generate_wurcs();

            if(currentGlycanWURCS == currentGlycanWURCSInFinalContainer)
                break;

            if(j == (finalGlycanPermutationContainer.size() - 1))
                finalGlycanPermutationContainer.push_back(tempObject);
        }
    }
    
    // output the permutation summary to the console for the user.
    for(int i = 0; i < finalGlycanPermutationContainer.size(); i++)
    {
        clipper::String temporaryWURCS = finalGlycanPermutationContainer[i].first.first.generate_wurcs();
        int valueLocation = privateer::util::find_index_of_value_from_wurcs(glycomics_database, temporaryWURCS);

        std::cout << "\tGenerated WURCS Sequence: " << temporaryWURCS << std::endl;
        std::cout << "\tAnomer Permutations = " << finalGlycanPermutationContainer[i].first.second[0] << "\t\tResidue Permutations = " << finalGlycanPermutationContainer[i].first.second[1] << "\tResidue Deletions = " << finalGlycanPermutationContainer[i].first.second[2] << std::endl; 
        std::cout << std::fixed << std::setprecision(2) << "\tPermutation Score(out of 100): " << finalGlycanPermutationContainer[i].second << std::endl;

        std::string glytoucanID, glyconnectID;
        glytoucanID = glycomics_database[valueLocation].GlyTouCanID;
        if (glytoucanID.front() == '"' && glytoucanID.front() == '"')
        {
            glytoucanID.erase(0, 1);
            glytoucanID.pop_back();
        }
        
        glyconnectID = glycomics_database[valueLocation].GlyConnectID;
        if (glyconnectID.front() == '"' && glyconnectID.front() == '"')
        {
            glyconnectID.erase(0, 1);
            glyconnectID.pop_back();
        }
        std::cout << "\tGlyTouCan Accession ID: " << glytoucanID << std::endl;
        std::cout << "\tGlyConnect ID: " << glyconnectID << std::endl;

        std::cout << std::endl;
    }     
}


void print_output_from_database(std::vector<privateer::json::Database>& glycomics_database, int valueLocation, clipper::MGlycan &currentGlycan)
{
    std::string glytoucanID;
    glytoucanID = glycomics_database[valueLocation].GlyTouCanID;
    if (glytoucanID.front() == '"' && glytoucanID.front() == '"')
    {
        glytoucanID.erase(0, 1);
        glytoucanID.pop_back();
    }
    std::cout << "\tManaged to find a matching GlyTouCan ID for WURCS sequence for this Glycan sequence!" << std::endl;
    std::cout << "\tGlyTouCan Accession ID: " << glytoucanID << std::endl;
    std::cout << "\tGlyTouCan link: "
              << "https://glytoucan.org/Structures/Glycans/" << glytoucanID << std::endl;

    if (glycomics_database[valueLocation].GlyConnectID != "NotFound")
    {
        std::string glyconnectID;
        glyconnectID = glycomics_database[valueLocation].GlyConnectID;
        if (glyconnectID.front() == '"' && glyconnectID.front() == '"')
        {
            glyconnectID.erase(0, 1);
            glyconnectID.pop_back();
        }
        std::cout << "\t\tFound a GlyConnect entry for this GlyTouCan ID!" << std::endl;
        std::cout << "\t\tGlyConnect ID: " << glyconnectID << std::endl;
        std::cout << "\t\tGlyConnect link: "
                    << "https://glyconnect.expasy.org/browser/structures/" << glyconnectID << std::endl;
    }
    else
    {
        std::cout << "\t\tThis GlyTouCan ID is not deposited on GlyConnect." << std::endl;
    }
}

