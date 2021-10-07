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
#include "third-party/sajson.h"

#include <vector>
#include <cstdio>    // for FILE, fopen, fclose
#include <memory>    // for unique_ptr

#include <algorithm>
#include <iterator>

#if defined(_WIN32)
#include "utf.hpp"
#endif

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
                        if (entry.get_object_value(k).get_type() != sajson::TYPE_STRING)
                            fail("Expected TYPE_STRING, got " + json_type_as_string(entry.get_object_value(k).get_type()));
                        temp.GlyConnectID = entry.get_object_value(k).as_string();
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
                    path_copy = env + "/privateer_glycomics_database.json";
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

                                else if(second_object.get_object_key(l).as_string() == "torsions")
                                {
                                    if (second_object.get_object_value(l).get_type() != sajson::TYPE_ARRAY)
                                        fail("Expected TYPE_ARRAY, got " + json_type_as_string(second_object.get_object_value(l).get_type()));

                                    std::vector<std::pair<float, float>> torsion_container;

                                    sajson::value torsion_array = second_object.get_object_value(l);

                                    if (second_object.get_object_value(l).get_type() != sajson::TYPE_ARRAY)
                                        fail("Expected TYPE_ARRAY, got " + json_type_as_string(second_object.get_object_value(l).get_type()));

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
                    path_copy = env + "/privateer_torsion_database.json";
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


        // ____________________________________SPECIFIC JSON FILES END________________________________________________ //
        

    }
}



#endif