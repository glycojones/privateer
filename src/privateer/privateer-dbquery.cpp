// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2018 Haroldas Bagdonas & Kevin Cowtan & Jon Agirre
// York Structural Biology Laboratory
// The University of York
// mailto: hb1115@york.ac.uk
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk


#include "privateer-dbquery.h"

void output_dbquery(nlohmann::json &jsonObject, clipper::String glycanWURCS, clipper::MGlycan &currentGlycan, std::vector<std::pair<std::pair<clipper::MGlycan, std::vector<int>>,float>>& finalGlycanPermutationContainer, bool glucose_only, bool clean_output, int sleepTimer, privateer::thread_pool& pool, bool useParallelism)
{
    int valueLocation;
    valueLocation = privateer::util::find_index_of_value(jsonObject, "Sequence", glycanWURCS);
    if (valueLocation != -1 && jsonObject[valueLocation]["glyconnect"] != "NotFound") print_output_from_database(jsonObject, valueLocation, currentGlycan);
    else
    {
        if (valueLocation != -1 && currentGlycan.number_of_nodes() > 1)
        {
            print_output_from_database(jsonObject, valueLocation, currentGlycan);
            std::cout << "\nWARNING: Unable to find a matching GlyTouCanID for WURCS sequence from this Glycan sequence! Attempting to find the closest matches by carrying out permutations" << std::endl;

            std::vector<std::pair<clipper::MGlycan, std::vector<int>>> alternativeGlycans;
            
            if(useParallelism) alternativeGlycans = generate_closest_matches_parallel(currentGlycan, jsonObject, glucose_only, clean_output, sleepTimer, pool, useParallelism);
            else               alternativeGlycans = generate_closest_matches_singlethreaded(currentGlycan, jsonObject, clean_output, glucose_only);   
            

            if (!alternativeGlycans.empty()) push_data_to_final_permutation_container(jsonObject, currentGlycan, alternativeGlycans, finalGlycanPermutationContainer);    
            else std::cout << "ERROR: Unable to generate permutations that would be found in GlyConnect database!" << std::endl;
        }
        else if (valueLocation != -1 && currentGlycan.number_of_nodes() <= 1)
        {
            print_output_from_database(jsonObject, valueLocation, currentGlycan);
            std::cout << "ERROR: Glycan is too short for permutations, therefore unable to find the closest match on GlyConnect." << std::endl;
        }
        else
        {
            if ( currentGlycan.number_of_nodes() > 1)
            {
                std::vector<std::pair<clipper::MGlycan, std::vector<int>>> alternativeGlycans;

                if(useParallelism) alternativeGlycans = generate_closest_matches_parallel(currentGlycan, jsonObject, glucose_only, clean_output, sleepTimer, pool, useParallelism);
                else               alternativeGlycans = generate_closest_matches_singlethreaded(currentGlycan, jsonObject, clean_output, glucose_only); 

                if (useParallelism)
                {
                    
                    #if DUMP
                        std::cout << std::endl;
                        DBG << "Number of jobs in the queue: " << pool.n_remaining_jobs() << std::endl;
                    #endif
                    
                    while(pool.n_remaining_jobs() > 0)
                        pool.sync();
                    
                    #if DUMP
                        DBG << "Number of jobs in the queue: " << pool.n_remaining_jobs() << " after sync operation!" << std::endl;
                    #endif
                }
            
                if (!alternativeGlycans.empty()) push_data_to_final_permutation_container(jsonObject, currentGlycan, alternativeGlycans, finalGlycanPermutationContainer);    
                else std::cout << "ERROR: Unable to generate permutations that would be found in GlyConnect database!" << std::endl;
            }
            else std::cout << "ERROR: Glycan is too short for permutations, therefore unable to generate alternative GlyTouCan and GlyConnect IDs!" << std::endl;
         }
    }
}


void push_data_to_final_permutation_container(nlohmann::json &jsonObject, clipper::MGlycan &currentGlycan, std::vector<std::pair<clipper::MGlycan, std::vector<int>>>& alternativeGlycans, std::vector<std::pair<std::pair<clipper::MGlycan, std::vector<int>>,float>>& finalGlycanPermutationContainer)
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
        int valueLocation = privateer::util::find_index_of_value(jsonObject, "Sequence", temporaryWURCS);

        std::cout << "\tGenerated WURCS Sequence: " << temporaryWURCS << std::endl;
        std::cout << "\tAnomer Permutations = " << finalGlycanPermutationContainer[i].first.second[0] << "\t\tResidue Permutations = " << finalGlycanPermutationContainer[i].first.second[1] << "\tResidue Deletions = " << finalGlycanPermutationContainer[i].first.second[2] << std::endl; 
        std::cout << std::fixed << std::setprecision(2) << "\tPermutation Score(out of 100): " << finalGlycanPermutationContainer[i].second << std::endl;

        std::string glytoucanID;
        glytoucanID = jsonObject[valueLocation]["AccessionNumber"];
        if (glytoucanID.front() == '"' && glytoucanID.front() == '"')
        {
            glytoucanID.erase(0, 1);
            glytoucanID.pop_back();
        }
        std::cout << "\tGlyTouCan Accession ID: " << glytoucanID << std::endl;
        std::cout << "\tGlyConnect ID: " << jsonObject[valueLocation]["glyconnect"]["id"] << std::endl;

        std::cout << std::endl;
    }     
}


void print_output_from_database(nlohmann::json& jsonObject, int valueLocation, clipper::MGlycan &currentGlycan)
{
    std::string glytoucanID;
    glytoucanID = jsonObject[valueLocation]["AccessionNumber"];
    if (glytoucanID.front() == '"' && glytoucanID.front() == '"')
    {
        glytoucanID.erase(0, 1);
        glytoucanID.pop_back();
    }
    std::cout << "\tManaged to find a matching GlyTouCan ID for WURCS sequence for this Glycan sequence!" << std::endl;
    std::cout << "\tGlyTouCan Accession ID: " << glytoucanID << std::endl;
    std::cout << "\tGlyTouCan link: "
              << "https://glytoucan.org/Structures/Glycans/" << glytoucanID << std::endl;

    if (jsonObject[valueLocation]["glyconnect"] != "NotFound")
    {
        if (!jsonObject[valueLocation]["glyconnect"]["comment"].is_null())
        {
            std::cout << "\t\tFound a GlyConnect entry for this GlyTouCan ID!" << std::endl;
            std::cout << "\t\tGlyConnect ID: " << jsonObject[valueLocation]["glyconnect"]["id"] << std::endl;
            std::cout << "\t\tGlycan Type: " << jsonObject[valueLocation]["glyconnect"]["type"] << std::endl;
            std::cout << "\t\tGlycan Core: " << jsonObject[valueLocation]["glyconnect"]["core"] << std::endl;
            std::cout << "\t\t" << jsonObject[valueLocation]["glyconnect"]["comment"] << std::endl;

            std::cout << std::endl
                      << "\t\tPrivateer checks: " << std::endl;
            std::cout << "\t\tGlycosylation type detected in the model: " << currentGlycan.get_type() << "\tGlycosylation type deposited on GlyConnect: " << jsonObject[valueLocation]["glyconnect"]["type"] << std::endl;
        }
        else
        {
            bool typeMatch = false;

            auto sourcesArray = jsonObject[valueLocation]["glyconnect"]["sources"];

            std::cout << "\t\tFound a GlyConnect entry for this GlyTouCan ID!" << std::endl;
            std::cout << "\t\tGlyConnect ID: " << jsonObject[valueLocation]["glyconnect"]["id"] << std::endl;
            std::cout << "\t\tGlycan Type: " << jsonObject[valueLocation]["glyconnect"]["type"] << std::endl;
            std::cout << "\t\tGlycan Core: " << jsonObject[valueLocation]["glyconnect"]["core"] << std::endl;
            std::cout << "\t\tGlycomics composition: " << jsonObject[valueLocation]["glyconnect"]["composition_string"] << std::endl;
            std::cout << "\t\tExpression system(s): " << std::endl;
            for (auto &element : sourcesArray)
            {
                std::cout << "\t\t\t\t\t" << element["species"] << std::endl;
            }
            std::cout << "\t\tTissue(s): " << std::endl;
            for (auto &element : sourcesArray)
            {
                std::cout << "\t\t\t\t\t" << element["system"] << std::endl;
            }
            std::cout << "\t\tProtein(s): " << std::endl;
            for (auto &element : sourcesArray)
            {
                std::cout << "\t\t\t\t\t" << element["protein"]["name"] << std::endl;
            }
            std::cout << "\t\tReviewed by GlyConnect: " << jsonObject[valueLocation]["glyconnect"]["reviewed"] << std::endl;
            std::cout << "\t\tGlyConnect link: "
                      << "https://glyconnect.expasy.org/browser/structures/" << jsonObject[valueLocation]["glyconnect"]["id"] << std::endl;

            std::cout << std::endl
                      << "\t\tPrivateer checks: " << std::endl;
            std::cout << "\t\tGlycosylation type detected in the model: " << currentGlycan.get_type() << "\tGlycosylation type deposited on GlyConnect: " << jsonObject[valueLocation]["glyconnect"]["type"] << std::endl;
        }
    }
    else
    {
        std::cout << "\t\tThis GlyTouCan ID is not deposited on GlyConnect." << std::endl;
    }
}

