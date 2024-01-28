// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York


#ifndef PRIVATEER_JSON_H_INCLUDED
#define PRIVATEER_JSON_H_INCLUDED

#include <iostream>
#include "privateer-error.h"
// #include "third-party/sajson.h" // Have to include it straight from gemmi, as otherwise it leads to "error: multiple definition of ‘enum sajson::type’" errors.
#include "gemmi/third_party/sajson.h" 
#include <vector>
#include <unordered_map>
#include <cstdio>    // for FILE, fopen, fclose
#include <memory>    // for unique_ptr

#include <algorithm>
#include <iterator>

#if defined(_WIN32)
#include "utf.hpp"
#endif

#include "simdjson.h"


namespace privateer
{
    namespace json
    {

        inline std::string json_type_as_string(sajson::type t) 
        {
            switch (t) 
            {
                case sajson::TYPE_INTEGER: return "<integer>";
                case sajson::TYPE_DOUBLE:  return "<double>";
                case sajson::TYPE_NULL:    return "<null>";
                case sajson::TYPE_FALSE:   return "<false>";
                case sajson::TYPE_TRUE:    return "<true>";
                case sajson::TYPE_STRING:  return "<string>";
                case sajson::TYPE_ARRAY:   return "<array>";
                case sajson::TYPE_OBJECT:  return "<object>";
                default:           return "<unknown type>";
            }
        }

        // file operations
        typedef std::unique_ptr<std::FILE, decltype(&std::fclose)> fileptr_t;

        inline fileptr_t file_open(const char* path, const char* mode) 
        {
            std::FILE* file;
            #if defined(_WIN32)
            std::wstring wpath = UTF8_to_wchar(path);
            std::wstring wmode = UTF8_to_wchar(mode);
            if ((file = ::_wfopen(wpath.c_str(), wmode.c_str())) == nullptr)
            #else
            if ((file = std::fopen(path, mode)) == nullptr)
            #endif
                fail(std::string("Failed to open file") +
                    (*mode == 'w' ? " for writing: " : ": ") + path);
            return fileptr_t(file, &std::fclose);
        }

        inline std::size_t file_size(std::FILE* f, const std::string& path) 
        {
            if (std::fseek(f, 0, SEEK_END) != 0)
                fail(path + ": fseek failed");
            long length = std::ftell(f);
            if (length < 0)
                fail(path + ": ftell failed");
            if (std::fseek(f, 0, SEEK_SET) != 0)
                fail(path + ": fseek failed");
            return length;
        }

        // ____________________________________SPECIFIC JSON FILES BEGIN________________________________________________ //
        
        // ____________________________________GLYCOMICS DATABASE BEGIN__________________________________________________ //

        struct GlycomicsDatabase 
        {
            std::string GlyTouCanID;
            std::string WURCS;
            std::string GlyConnectID;
        };

        inline void generate_glycomics_database(std::vector<GlycomicsDatabase>& glycomics_database, const sajson::document& jsonObject) 
        {
            sajson::value root = jsonObject.get_root();
            if (root.get_type() != sajson::TYPE_OBJECT)
                fail("not Privateer Glycomics Database JSON file");
            std::string database_name_key = root.get_object_key(1).as_string();
            std::string database_name_value = root.get_object_value(1).as_string();

            if (database_name_key != "database_name" || database_name_value != "glycomics_database")
                fail("not Glycomics Database JSON file - should be \"database_name\": \"glycomics_database\", instead received \"" + database_name_key + "\": \"" + database_name_value + "\"\n");
            
            std::string database_last_update = root.get_object_value(2).as_string();
            
            sajson::value entries_array = root.get_object_value(0);
            std::string entries_name = root.get_object_key(0).as_string();
            if (entries_array.get_type() != sajson::TYPE_ARRAY || entries_name != "data")
                fail("Expected TYPE_ARRAY, got " + json_type_as_string(entries_array.get_type()) + " with key value equal to \"data\", instead got: \"" + entries_name + "\"");
            
            for (size_t i = 0; i != entries_array.get_length(); ++i) 
            {
                sajson::value entry = entries_array.get_array_element(i);
                if (entry.get_type() != sajson::TYPE_OBJECT)
                fail("Expected TYPE_OBJECT, got " + json_type_as_string(entry.get_type()));
                
                GlycomicsDatabase temp;
                for(size_t k = 0; k != entry.get_length(); k++) 
                {
                    if(entry.get_object_key(k).as_string() == "AccessionNumber") 
                    {
                        if (entry.get_object_value(k).get_type() != sajson::TYPE_STRING)
                            fail("Expected TYPE_STRING, got " + json_type_as_string(entry.get_object_value(k).get_type()));
                        temp.GlyTouCanID = entry.get_object_value(k).as_string();
                    }
                    else if(entry.get_object_key(k).as_string() == "Sequence") 
                    {
                        if (entry.get_object_value(k).get_type() != sajson::TYPE_STRING)
                            fail("Expected TYPE_STRING, got " + json_type_as_string(entry.get_object_value(k).get_type()));
                        temp.WURCS = entry.get_object_value(k).as_string();
                    }
                    else if(entry.get_object_key(k).as_string() == "glyconnect") 
                    {
                        if (entry.get_object_value(k).get_type() == sajson::TYPE_STRING || entry.get_object_value(k).get_type() == sajson::TYPE_ARRAY)
                        {
                            if(entry.get_object_value(k).get_type() == sajson::TYPE_STRING)
                                temp.GlyConnectID = entry.get_object_value(k).as_string();
                            else
                            {
                                sajson::value glyconnect_id_array = entry.get_object_value(k);
                                if(glyconnect_id_array.get_type() != sajson::TYPE_ARRAY)
                                    fail("Expected TYPE_ARRAY, got " + json_type_as_string(entry.get_object_value(k).get_type()));
                                else
                                {
                                    std::string GlyConnectID_list = "";
                                    for (size_t j = 0; j != glyconnect_id_array.get_length(); ++j) 
                                    {
                                        sajson::value glyconnect_id_entry = glyconnect_id_array.get_array_element(j);
                                        if (glyconnect_id_entry.get_type() != sajson::TYPE_STRING)
                                            fail("Expected TYPE_STRING, got " + json_type_as_string(glyconnect_id_entry.get_type()));
                                        GlyConnectID_list.append(glyconnect_id_entry.as_string());
                                        GlyConnectID_list.append(", ");
                                    }
                                    GlyConnectID_list.pop_back(); // removes empty space
                                    GlyConnectID_list.pop_back(); // removes the comma. Yes, I am this lazy.
                                    temp.GlyConnectID = GlyConnectID_list;
                                }
                            }
                        }
                        else
                            fail("Expected TYPE_STRING or TYPE_ARRAY, got " + json_type_as_string(entry.get_object_value(k).get_type()));
                    }
                }
                glycomics_database.push_back(temp);
            }

            std::cout << std::endl << "Successfully imported Privateer's Glycomics database for cross-checking with GlyTouCan and GlyConnect databases.\nLast database update: " << database_last_update << std::endl << std::endl; 
                          
        }

        inline std::vector<GlycomicsDatabase> read_json_glycomics_insitu(char* buffer, size_t size, const std::string& name) 
        {
            std::vector<GlycomicsDatabase> glycomics_database;
            sajson::document json = sajson::parse(sajson::dynamic_allocation(),
                                                sajson::mutable_string_view(size, buffer));
            if (!json.is_valid())
                fail(name, ":", std::to_string(json.get_error_line()), " error: ",
                    json.get_error_message_as_string());
            generate_glycomics_database(glycomics_database, json);
            return glycomics_database;
        }

        inline std::vector<GlycomicsDatabase> read_json_file_for_glycomics_database(const std::string& path) 
        {
            std::string path_copy = path;
            if(path_copy == "nopath" || path_copy.empty()) 
            {
                std::string env;
                if(std::getenv("PRIVATEERDATA"))
                {
                    std::string env( std::getenv("PRIVATEERDATA") );
                    path_copy = env + "/glycomics/privateer_glycomics_database.json";
                }
                else
                {
                    env = std::getenv("CLIBD");
                    path_copy = env + "/privateer_data/glycomics/privateer_glycomics_database.json";
                }
            }

            std::cout << "Reading " << path_copy << " for Glycomics database" << std::endl;

            fileptr_t f = file_open(path_copy.c_str(), "rb");
            size_t buf_size = file_size(f.get(), path_copy);
            std::vector<char> buffer(buf_size);
            if (std::fread(buffer.data(), buffer.size(), 1, f.get()) != 1)
                fail(path_copy + ": fread failed");
            
            return read_json_glycomics_insitu(buffer.data(), buffer.size(), path_copy);
        }
        // ____________________________________GLYCOMICS DATABASE END__________________________________________________ //
        
        // ____________________________________TORSIONS DATABASE BEGIN__________________________________________________ //
        
        struct TorsionsDatabase 
        {
            std::string type;
            std::string first_residue;
            std::string second_residue;
            std::string donor_position;
            std::string acceptor_position;
            std::vector<std::pair<float, float>> torsions; // .first = Phi, .second = Psi
        };


        inline void generate_torsions_database(std::vector<TorsionsDatabase>& torsions_database, const sajson::document& jsonObject)
        {
            sajson::value root = jsonObject.get_root();
            if (root.get_type() != sajson::TYPE_OBJECT)
                fail("not Privateer Torsions Database JSON file");
            std::string database_name_key = root.get_object_key(1).as_string();
            std::string database_name_value = root.get_object_value(1).as_string();

            if (database_name_key != "database_name" || database_name_value != "torsion_database")
                fail("not Torsions Database JSON file - should be \"database_name\": \"torsion_database\", instead received \"" + database_name_key + "\": \"" + database_name_value + "\"\n");
            
            std::string database_last_update = root.get_object_value(2).as_string();
            
            sajson::value entries_array = root.get_object_value(0);
            std::string entries_name = root.get_object_key(0).as_string();
            if (entries_array.get_type() != sajson::TYPE_ARRAY || entries_name != "data")
                fail("Expected TYPE_ARRAY, got " + json_type_as_string(entries_array.get_type()) + " with key value equal to \"data\", instead got: \"" + entries_name + "\"");
            
            for (size_t i = 0; i != entries_array.get_length(); ++i) 
            {
                sajson::value entry = entries_array.get_array_element(i);
                if (entry.get_type() != sajson::TYPE_OBJECT)
                fail("Expected TYPE_OBJECT, got " + json_type_as_string(entry.get_type()));
                
                TorsionsDatabase temp;
                for(size_t j = 0; j != entry.get_length(); j++) 
                {
                    if(entry.get_object_key(j).as_string() == "first") 
                    {
                        if (entry.get_object_value(j).get_type() != sajson::TYPE_STRING)
                            fail("Expected TYPE_STRING, got " + json_type_as_string(entry.get_object_value(j).get_type()));
                        temp.first_residue = entry.get_object_value(j).as_string();
                    }
                    else if(entry.get_object_key(j).as_string() == "type") 
                    {
                        if (entry.get_object_value(j).get_type() != sajson::TYPE_STRING)
                            fail("Expected TYPE_STRING, got " + json_type_as_string(entry.get_object_value(j).get_type()));
                        temp.type = entry.get_object_value(j).as_string();
                    }
                    else if(entry.get_object_key(j).as_string() == "second") 
                    {
                        sajson::value second_array = entry.get_object_value(j);
                        if (second_array.get_type() != sajson::TYPE_ARRAY)
                            fail("Expected TYPE_ARRAY, got " + json_type_as_string(second_array.get_type()));
                        
                        for(size_t k = 0; k != second_array.get_length(); k++)
                        {
                            sajson::value second_object = second_array.get_array_element(k);

                            if (second_object.get_type() != sajson::TYPE_OBJECT)
                                fail("Expected TYPE_OBJECT, got " + json_type_as_string(second_object.get_type()));

                            for(size_t l = 0; l != second_object.get_length(); l++)
                            {
                                if(second_object.get_object_key(l).as_string() == "sugar")
                                {
                                    if (second_object.get_object_value(l).get_type() != sajson::TYPE_STRING)
                                        fail("Expected TYPE_STRING, got " + json_type_as_string(second_object.get_object_value(l).get_type()));
                                    
                                    temp.second_residue = second_object.get_object_value(l).as_string();
                                }

                                if(second_object.get_object_key(l).as_string() == "donor_position")
                                {
                                    if (second_object.get_object_value(l).get_type() != sajson::TYPE_STRING)
                                        fail("Expected TYPE_STRING, got " + json_type_as_string(second_object.get_object_value(l).get_type()));
                                    
                                    temp.donor_position = second_object.get_object_value(l).as_string();
                                }
                                
                                if(second_object.get_object_key(l).as_string() == "acceptor_position")
                                {
                                    if (second_object.get_object_value(l).get_type() != sajson::TYPE_STRING)
                                        fail("Expected TYPE_STRING, got " + json_type_as_string(second_object.get_object_value(l).get_type()));
                                    
                                    temp.acceptor_position = second_object.get_object_value(l).as_string();
                                }

                                else if(second_object.get_object_key(l).as_string() == "torsions")
                                {
                                    if (second_object.get_object_value(l).get_type() != sajson::TYPE_ARRAY)
                                        fail("Expected TYPE_ARRAY, got " + json_type_as_string(second_object.get_object_value(l).get_type()));

                                    std::vector<std::pair<float, float>> torsion_container;

                                    sajson::value torsion_array = second_object.get_object_value(l);

                                    // if (second_object.get_object_value(l).get_type() != sajson::TYPE_ARRAY)
                                        // fail("Expected TYPE_ARRAY, got " + json_type_as_string(second_object.get_object_value(l).get_type()));

                                    for(size_t m = 0; m != torsion_array.get_length(); m++)
                                    {
                                        sajson::value individual_torsions = torsion_array.get_array_element(m);

                                        if(individual_torsions.get_type() != sajson::TYPE_OBJECT)
                                            fail("Expected TYPE_OBJECT, got " + json_type_as_string(individual_torsions.get_type()));
                                        
                                        float Psi = -42069, Phi = -42069;
                                        for(size_t n = 0; n != individual_torsions.get_length(); n++)
                                        {   
                                            if(individual_torsions.get_object_key(n).as_string() == "Phi") 
                                            {
                                                if (individual_torsions.get_object_value(n).get_type() != sajson::TYPE_DOUBLE)
                                                    fail("Expected TYPE_STRING, got " + json_type_as_string(individual_torsions.get_object_value(n).get_type()));
                                                
                                                Phi = individual_torsions.get_object_value(n).get_number_value();
                                            }
                                            else if(individual_torsions.get_object_key(n).as_string() == "Psi") 
                                            {
                                                if (individual_torsions.get_object_value(n).get_type() != sajson::TYPE_DOUBLE)
                                                    fail("Expected TYPE_STRING, got " + json_type_as_string(individual_torsions.get_object_value(j).get_type()));
                                                
                                                Psi = individual_torsions.get_object_value(n).get_number_value();
                                            }
                                        }

                                        torsion_container.push_back(std::make_pair(Phi, Psi));
                                    }

                                    temp.torsions = torsion_container;
                                }
                            }
                            torsions_database.push_back(temp);                    
                        }
                    }
                }
            }

            std::cout << std::endl << "Successfully imported Privateer's Torsions database for validating sugar linkages Ramachandran style.\nLast database update: " << database_last_update << std::endl << std::endl;                   
        }

        inline std::vector<TorsionsDatabase> read_json_torsions_insitu(char* buffer, size_t size, const std::string& name) 
        {
            std::vector<TorsionsDatabase> torsions_database;
            sajson::document json = sajson::parse(sajson::dynamic_allocation(),
                                                sajson::mutable_string_view(size, buffer));
            if (!json.is_valid())
                fail(name, ":", std::to_string(json.get_error_line()), " error: ",
                    json.get_error_message_as_string());
            generate_torsions_database(torsions_database, json);
            return torsions_database;
        }

        inline std::vector<TorsionsDatabase> read_json_file_for_torsions_database(const std::string& path) 
        {
            std::string path_copy = path;
            if(path_copy == "nopath" || path_copy.empty()) 
            {
                std::string env;
                if(std::getenv("PRIVATEERDATA"))
                {
                    std::string env( std::getenv("PRIVATEERDATA") );
                    path_copy = env + "/linkage_torsions/privateer_torsion_database.json";
                }
                else
                {
                    env = std::getenv("CLIBD");
                    path_copy = env + "/privateer_data/linkage_torsions/privateer_torsion_database.json";
                }
            }

            std::cout << "Reading " << path_copy << " for Torsions database" << std::endl;

            fileptr_t f = file_open(path_copy.c_str(), "rb");
            size_t buf_size = file_size(f.get(), path_copy);
            std::vector<char> buffer(buf_size);
            if (std::fread(buffer.data(), buffer.size(), 1, f.get()) != 1)
                fail(path_copy + ": fread failed");
            
            return read_json_torsions_insitu(buffer.data(), buffer.size(), path_copy);
        }

        // ____________________________________TORSIONS DATABASE END__________________________________________________ //

        // ____________________________________TORSIONS Z SCORE DATABASE BEGIN________________________________________ //
        
        struct TorsionsZScoreStatistics {
            float average_z_score_for_database; 
            float stddev_z_score_for_database;
        };

        struct TorsionsZScoreDatabase 
        {
            std::string donor_sugar;
            std::string acceptor_sugar;
            std::string donor_end;
            std::string acceptor_end;
            std::pair<float, float> summary; // .first = count_mean, .second = count_stdev
            std::vector<std::unordered_map<std::string, int>> bin_data;
            std::vector<std::string> pdb_list;
        };

        struct GlobalTorsionZScore { 
            std::vector<TorsionsZScoreDatabase> database_array; 
            TorsionsZScoreStatistics statistics;
        };

        inline void generate_torsions_zscore_database_simd(GlobalTorsionZScore& torsions_database, const simdjson::ondemand::document) {

        }

        inline void generate_torsions_zscore_database(GlobalTorsionZScore& global_torsion_z_score, const sajson::document& jsonObject) 
        {

            std::vector<TorsionsZScoreDatabase> torsions_zscore_database;
            TorsionsZScoreStatistics temp_statistics;

            sajson::value root = jsonObject.get_root();
            if (root.get_type() != sajson::TYPE_OBJECT)
                fail("not Privateer Torsions Z-Score Database JSON file");
            std::string database_name_key = root.get_object_key(1).as_string();
            std::string database_name_value = root.get_object_value(1).as_string();
            if (database_name_key != "database_name" || database_name_value != "torsions_z_score_database")
                fail("not Torsions Database JSON file - should be \"database_name\": \"torsions_z_score_database\", instead received \"" + database_name_key + "\": \"" + database_name_value + "\"\n");
            
            std::string database_last_update = root.get_object_value(2).as_string();

            sajson::value entries_array = root.get_object_value(0);
            std::string entries_name = root.get_object_key(0).as_string();
            if (entries_array.get_type() != sajson::TYPE_OBJECT || entries_name != "data")
                fail("Expected TYPE_OBJECT, got " + json_type_as_string(entries_array.get_type()) + " with key value equal to \"data\", instead !!!got: \"" + entries_name + "\"");
            
            for (size_t data_index = 0; data_index != entries_array.get_length(); data_index++) { 
                
                sajson::value data_array = entries_array.get_object_value(data_index);
                std::string data_key = entries_array.get_object_key(data_index).as_string();
             
                if (data_key == "linkage_data") {

                    for (size_t linkage_data_index = 0; linkage_data_index != data_array.get_length(); ++linkage_data_index) 
                    {
                        sajson::value entry = data_array.get_array_element(linkage_data_index);
                        if (entry.get_type() != sajson::TYPE_OBJECT)
                        fail("Expected TYPE_OBJECT, got " + json_type_as_string(entry.get_type()));
                        
                        TorsionsZScoreDatabase temp;
                        for(size_t donor_index = 0; donor_index != entry.get_length(); donor_index++) 
                        {
                            if(entry.get_object_key(donor_index).as_string() == "donor") 
                            {
                                if (entry.get_object_value(donor_index).get_type() != sajson::TYPE_STRING)
                                    fail("Expected TYPE_STRING, got " + json_type_as_string(entry.get_object_value(donor_index).get_type()));
                                temp.donor_sugar = entry.get_object_value(donor_index).as_string();
                            }
                            else if(entry.get_object_key(donor_index).as_string() == "acceptor") 
                            {
                                sajson::value acceptor_array = entry.get_object_value(donor_index);
                                if (acceptor_array.get_type() != sajson::TYPE_ARRAY)
                                    fail("Expected TYPE_ARRAY, got " + json_type_as_string(acceptor_array.get_type()));
                                
                                for(size_t acceptor_index = 0; acceptor_index != acceptor_array.get_length(); acceptor_index++)
                                {
                                    sajson::value acceptor_object = acceptor_array.get_array_element(acceptor_index);

                                    if (acceptor_object.get_type() != sajson::TYPE_OBJECT)
                                        fail("Expected TYPE_OBJECT, got " + json_type_as_string(acceptor_object.get_type()));

                                    for(size_t acceptor_object_index = 0; acceptor_object_index != acceptor_object.get_length(); acceptor_object_index++)
                                    {
                                        if(acceptor_object.get_object_key(acceptor_object_index).as_string() == "sugar")
                                        {
                                            if (acceptor_object.get_object_value(acceptor_object_index).get_type() != sajson::TYPE_STRING)
                                                fail("Expected TYPE_STRING, got " + json_type_as_string(acceptor_object.get_object_value(acceptor_object_index).get_type()));
                                            
                                            temp.acceptor_sugar = acceptor_object.get_object_value(acceptor_object_index).as_string();
                                        }
                                        else if(acceptor_object.get_object_key(acceptor_object_index).as_string() == "Linkage")
                                        {
                                            if (acceptor_object.get_object_value(acceptor_object_index).get_type() != sajson::TYPE_ARRAY)
                                                fail("Expected TYPE_ARRAY, got " + json_type_as_string(acceptor_object.get_object_value(acceptor_object_index).get_type()));

                                            sajson::value Linkage_array = acceptor_object.get_object_value(acceptor_object_index);

                                            for(size_t Linkage_index = 0; Linkage_index != Linkage_array.get_length(); Linkage_index++)
                                            {
                                                sajson::value Linkage_object = Linkage_array.get_array_element(Linkage_index);

                                                if(Linkage_object.get_type() != sajson::TYPE_OBJECT)
                                                    fail("Expected TYPE_OBJECT, got " + json_type_as_string(Linkage_object.get_type()));
                                                
                                                for(size_t Linkage_object_index = 0; Linkage_object_index != Linkage_object.get_length(); Linkage_object_index++)
                                                {
                                                    if(Linkage_object.get_object_key(Linkage_object_index).as_string() == "donor_end")
                                                    {
                                                        if (Linkage_object.get_object_value(Linkage_object_index).get_type() != sajson::TYPE_STRING)
                                                            fail("Expected TYPE_STRING, got " + json_type_as_string(Linkage_object.get_object_value(Linkage_object_index).get_type()));
                                                        
                                                        temp.donor_end = Linkage_object.get_object_value(Linkage_object_index).as_string();
                                                    }
                                                    
                                                    else if(Linkage_object.get_object_key(Linkage_object_index).as_string() == "acceptor_end")
                                                    {
                                                        if (Linkage_object.get_object_value(Linkage_object_index).get_type() != sajson::TYPE_STRING)
                                                            fail("Expected TYPE_STRING, got " + json_type_as_string(Linkage_object.get_object_value(Linkage_object_index).get_type()));
                                                        
                                                        temp.acceptor_end = Linkage_object.get_object_value(Linkage_object_index).as_string();
                                                    }

                                                    else if(Linkage_object.get_object_key(Linkage_object_index).as_string() == "Linkage_data")
                                                    {
                                                        sajson::value Linkage_data_object = Linkage_object.get_object_value(Linkage_object_index);
                                                        if(Linkage_data_object.get_type() != sajson::TYPE_OBJECT)
                                                            fail("Expected TYPE_OBJECT, got " + json_type_as_string(Linkage_data_object.get_type()));

                                                        
                                                        std::vector<std::unordered_map<std::string, int>> tmp_bins_for_linkage; 
                                                        for(size_t Linkage_data_object_index = 0; Linkage_data_object_index != Linkage_data_object.get_length(); Linkage_data_object_index++)
                                                        {
                                                        // std::cout << "KEY : " << Linkage_data_object.get_object_key(Linkage_data_object_index).as_string() << std::endl;
                                                            if (Linkage_data_object.get_object_key(Linkage_data_object_index).as_string() == "summary")
                                                            {
                                                                sajson::value Linkage_data_summary_object = Linkage_data_object.get_object_value(Linkage_data_object_index);
                                                                if(Linkage_data_summary_object.get_type() != sajson::TYPE_OBJECT)
                                                                    fail("Expected TYPE_OBJECT, got " + json_type_as_string(Linkage_data_summary_object.get_type()));

                                                                for(size_t Linkage_data_object_summary_index = 0; Linkage_data_object_summary_index != Linkage_data_summary_object.get_length(); Linkage_data_object_summary_index++)
                                                                {
                                                                    if(Linkage_data_summary_object.get_object_key(Linkage_data_object_summary_index).as_string() == "count_mean")
                                                                    {
                                                                        if (Linkage_data_summary_object.get_object_value(Linkage_data_object_summary_index).get_type() != sajson::TYPE_DOUBLE)
                                                                            fail("Expected TYPE_DOUBLE, got " + json_type_as_string(Linkage_data_summary_object.get_object_value(Linkage_data_object_summary_index).get_type()));
                                                                        
                                                                        temp.summary.first = Linkage_data_summary_object.get_object_value(Linkage_data_object_summary_index).get_double_value();
                                                                    }
                                                                    else if(Linkage_data_summary_object.get_object_key(Linkage_data_object_summary_index).as_string() == "count_stdev")
                                                                    {
                                                                        if (Linkage_data_summary_object.get_object_value(Linkage_data_object_summary_index).get_type() != sajson::TYPE_DOUBLE)
                                                                            fail("Expected TYPE_DOUBLE, got " + json_type_as_string(Linkage_data_summary_object.get_object_value(Linkage_data_object_summary_index).get_type()));
                                                                        
                                                                        temp.summary.second = Linkage_data_summary_object.get_object_value(Linkage_data_object_summary_index).get_double_value();
                                                                    }
                                                                }
                                                            }
                                                            else if(Linkage_data_object.get_object_key(Linkage_data_object_index).as_string() == "bin_data")
                                                            {
                                                                sajson::value bin_data_array = Linkage_data_object.get_object_value(Linkage_data_object_index);
                                                                if(bin_data_array.get_type() != sajson::TYPE_ARRAY)
                                                                    fail("Expected TYPE_ARRAY, got " + json_type_as_string(bin_data_array.get_type()));
                                                                
                                                                for(size_t bin_data_index = 0; bin_data_index != bin_data_array.get_length(); bin_data_index++)
                                                                {
                                                                    std::unordered_map<std::string, int> current_bin;

                                                                    sajson::value current_bin_data_entry = bin_data_array.get_array_element(bin_data_index);
                                                                    if(current_bin_data_entry.get_type() != sajson::TYPE_OBJECT)
                                                                        fail("Expected TYPE_OBJECT, got " + json_type_as_string(current_bin_data_entry.get_type()));
                                                                    for(size_t current_bin_data_entry_index = 0; current_bin_data_entry_index != current_bin_data_entry.get_length(); current_bin_data_entry_index++)
                                                                    {
                                                                        if(current_bin_data_entry.get_object_key(current_bin_data_entry_index).as_string() == "lower_phi")
                                                                        {
                                                                            if (current_bin_data_entry.get_object_value(current_bin_data_entry_index).get_type() != sajson::TYPE_DOUBLE)
                                                                                fail("Expected TYPE_DOUBLE, got " + json_type_as_string(current_bin_data_entry.get_object_value(current_bin_data_entry_index).get_type()));
                                                                            
                                                                            current_bin["lower_phi"] = int(current_bin_data_entry.get_object_value(current_bin_data_entry_index).get_double_value());
                                                                        }
                                                                        else if(current_bin_data_entry.get_object_key(current_bin_data_entry_index).as_string() == "higher_phi")
                                                                        {
                                                                            if (current_bin_data_entry.get_object_value(current_bin_data_entry_index).get_type() != sajson::TYPE_DOUBLE)
                                                                                fail("Expected TYPE_DOUBLE, got " + json_type_as_string(current_bin_data_entry.get_object_value(current_bin_data_entry_index).get_type()));
                                                                            
                                                                            current_bin["higher_phi"] = int(current_bin_data_entry.get_object_value(current_bin_data_entry_index).get_double_value());
                                                                        }
                                                                        else if(current_bin_data_entry.get_object_key(current_bin_data_entry_index).as_string() == "lower_psi")
                                                                        {
                                                                            if (current_bin_data_entry.get_object_value(current_bin_data_entry_index).get_type() != sajson::TYPE_DOUBLE)
                                                                                fail("Expected TYPE_DOUBLE, got " + json_type_as_string(current_bin_data_entry.get_object_value(current_bin_data_entry_index).get_type()));
                                                                            
                                                                            current_bin["lower_psi"] = int(current_bin_data_entry.get_object_value(current_bin_data_entry_index).get_double_value());
                                                                        }
                                                                        else if(current_bin_data_entry.get_object_key(current_bin_data_entry_index).as_string() == "higher_psi")
                                                                        {
                                                                            if (current_bin_data_entry.get_object_value(current_bin_data_entry_index).get_type() != sajson::TYPE_DOUBLE)
                                                                                fail("Expected TYPE_DOUBLE, got " + json_type_as_string(current_bin_data_entry.get_object_value(current_bin_data_entry_index).get_type()));
                                                                            
                                                                            current_bin["higher_psi"] = int(current_bin_data_entry.get_object_value(current_bin_data_entry_index).get_double_value());
                                                                        }
                                                                        else if(current_bin_data_entry.get_object_key(current_bin_data_entry_index).as_string() == "count")
                                                                        {
                                                                            if (current_bin_data_entry.get_object_value(current_bin_data_entry_index).get_type() != sajson::TYPE_DOUBLE)
                                                                                fail("Expected TYPE_DOUBLE, got " + json_type_as_string(current_bin_data_entry.get_object_value(current_bin_data_entry_index).get_type()));
                                                                            
                                                                            current_bin["count"] = int(current_bin_data_entry.get_object_value(current_bin_data_entry_index).get_double_value());
                                                                        }
                                                                    }
                                                                    tmp_bins_for_linkage.push_back(current_bin);
                                                                }
                                                            }
                                                            else if(Linkage_data_object.get_object_key(Linkage_data_object_index).as_string() == "pdb_list") { 

                                                                sajson::value pdb_list_array = Linkage_data_object.get_object_value(Linkage_data_object_index);

                                                                if(pdb_list_array.get_type() != sajson::TYPE_ARRAY)
                                                                    fail("Expected TYPE_ARRAY, got " + json_type_as_string(pdb_list_array.get_type()));

                                                                std::vector<std::string> temp_pdb_list;
        
                                                                for(size_t pdb_list_array_index = 0; pdb_list_array_index != pdb_list_array.get_length(); pdb_list_array_index++) {
                                                                    std::string pdb_code = pdb_list_array.get_array_element(pdb_list_array_index).as_string();
                                                                    temp_pdb_list.push_back(pdb_code);
                                                                    
                                                                }
                                                                temp.pdb_list = temp_pdb_list;
                                                            }
                                                        }
                                                        temp.bin_data = tmp_bins_for_linkage;
                                                    }
                                                }
                                                torsions_zscore_database.push_back(temp);
                                            }
                                        }
                                    }                  
                                }
                            }
                        }
                    }
                }

                if (data_key == "statistics") { 

                    for (size_t statistics_data_index = 0; statistics_data_index != data_array.get_length(); statistics_data_index++) { 

                        std::string statistics_key = data_array.get_object_key(statistics_data_index).as_string();
                        float statistics_value = static_cast<float>(data_array.get_object_value(statistics_data_index).get_double_value());


                        if (statistics_key != "Database Z Score Average" && statistics_key != "Database Z Score StdDev") {
                            fail("Key not defined in statistics obeject");
                        }

                        if (statistics_key == "Database Z Score Average") {
                            temp_statistics.average_z_score_for_database = statistics_value;
                        }

                        if (statistics_key == "Database Z Score StdDev") { 
                            temp_statistics.stddev_z_score_for_database = statistics_value;
                        }

                    }

                }
            }

            global_torsion_z_score.database_array = torsions_zscore_database;
            global_torsion_z_score.statistics = temp_statistics;
            // std::cout << std::endl << "Successfully imported Privateer's Torsions Z-score database.\nLast database update: " << database_last_update << std::endl << std::endl;                   
        }

        inline GlobalTorsionZScore read_json_zscore_torsions_insitu(char* buffer, size_t size, const std::string& name) 
        {
            GlobalTorsionZScore global_torsion_z_score;
            sajson::document json = sajson::parse(sajson::dynamic_allocation(),
                                                sajson::mutable_string_view(size, buffer));
            if (!json.is_valid())
                fail(name, ":", std::to_string(json.get_error_line()), " error: ",
                    json.get_error_message_as_string());
            generate_torsions_zscore_database(global_torsion_z_score, json);
            // std::cout << "___Length of torsions_zscore_database = " << global_torsion_z_score.database_array.size() << std::endl;
            return global_torsion_z_score;
        }

        inline GlobalTorsionZScore read_json_zscore_db(const std::string& path) {
            simdjson::ondemand::parser parser;
            auto json = simdjson::padded_string::load(path);
            simdjson::ondemand::document doc = parser.iterate(json);

            std::vector<TorsionsZScoreDatabase> torsions_zscore_database;
            TorsionsZScoreStatistics temp_statistics;
            auto data = doc["data"];

            auto stats = data["statistics"];
            temp_statistics.average_z_score_for_database = stats["Database Z Score Average"].get_double();
            temp_statistics.stddev_z_score_for_database = stats["Database Z Score StdDev"].get_double();

            auto linkage_data = data["linkage_data"];
            for (auto l: linkage_data) {
                auto donor = l["donor"];
                auto acceptors = l["acceptor"];
                for (auto a: acceptors) {
                    auto sugar = a["sugar"];
                    auto linkage = a["Linkage"];
                    for (auto lk: linkage) { 

                        TorsionsZScoreDatabase tmp_db; 
                        tmp_db.acceptor_sugar = std::string(sugar.get_string().value()); 
                        tmp_db.donor_sugar = std::string(donor.get_string().value());
                        tmp_db.acceptor_end = std::string(lk["acceptor_end"].get_string().value());
                        tmp_db.donor_end = std::string(lk["donor_end"].get_string().value());

                        auto linkage_data = lk["Linkage_data"];
                        auto summary = linkage_data["summary"];
                        tmp_db.summary = std::make_pair(
                            summary["count_mean"].get_double().value(),
                            summary["count_stdev"].get_double().value()
                        );

                        std::vector<std::unordered_map<std::string, int>> tmp_bins_for_linkage; 

                        auto bin_data = linkage_data["bin_data"];
                        for (auto bin: bin_data) { 
                            std::unordered_map<std::string, int> bin_map; 
                            bin_map["lower_phi"] = int(bin["lower_phi"].get_double().value());
                            bin_map["higher_phi"] = int(bin["higher_phi"].get_double().value());
                            bin_map["lower_psi"] = int(bin["lower_psi"].get_double().value());
                            bin_map["higher_psi"] = int(bin["higher_psi"].get_double().value());
                            bin_map["count"] = int(bin["count"].get_double());
                            tmp_bins_for_linkage.push_back(bin_map);
                        }
                        tmp_db.bin_data = tmp_bins_for_linkage;

                        auto pdb_list = linkage_data["pdb_list"];
                        std::vector<std::string> temp_pdb_list;

                        for (auto pdb: pdb_list) { 
                            temp_pdb_list.push_back(std::string(pdb.get_string().value()));  
                        }
                        tmp_db.pdb_list = temp_pdb_list;
                        torsions_zscore_database.push_back(tmp_db);

                    }
                }
             }

            GlobalTorsionZScore global_torsion_z_score;
            global_torsion_z_score.statistics = temp_statistics;
            global_torsion_z_score.database_array = torsions_zscore_database;       

            return global_torsion_z_score;
        }

        inline GlobalTorsionZScore read_json_file_for_torsions_zscore_database(const std::string& path) 
        {
            std::string path_copy = path;
            if(path_copy == "nopath" || path_copy.empty())
            {
                std::string env;
                if(std::getenv("PRIVATEERDATA"))
                {
                    std::string env( std::getenv("PRIVATEERDATA") );
                    path_copy = env + "/linkage_torsions/privateer_torsions_z_score_database.json";
                }
                else
                {
                    env = std::getenv("CLIBD");
                    path_copy = env + "/privateer_data/linkage_torsions/privateer_torsions_z_score_database.json";
                }
            }
                std::cout << path_copy << std::endl;
            return read_json_zscore_db(path_copy);

            // fileptr_t f = file_open(path_copy.c_str(), "rb");
            // size_t buf_size = file_size(f.get(), path_copy);
            // std::vector<char> buffer(buf_size);
            // if (std::fread(buffer.data(), buffer.size(), 1, f.get()) != 1)
            //     fail(path_copy + ": fread failed");
            
            // return read_json_zscore_torsions_insitu(buffer.data(), buffer.size(), path_copy);
        }

        // ____________________________________TORSIONS Z SCORE DATABASE END_________________________________________ //

    }
}



#endif