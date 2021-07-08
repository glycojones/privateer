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

#if defined(_WIN32)
    #include "third-party/utf.hpp"
#endif

#include <vector>
#include <cstdio>    // for FILE, fopen, fclose
#include <memory>    // for unique_ptr

#include <algorithm>
#include <iterator>

// Solution obtained and adopted from Marcin Wojdyr's Gemmi project - https://github.com/project-gemmi/gemmi 
namespace privateer
{
    namespace json
    {
        struct Database 
        {
            std::string GlyTouCanID;
            std::string WURCS;
            std::string GlyConnectID;
        };

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

        inline void generate_database(std::vector<Database>& database, const sajson::document& jsonObject) 
        {
            sajson::value root = jsonObject.get_root();
            if (root.get_type() != sajson::TYPE_OBJECT || root.get_length() != 1)
                fail("not Privateer Database JSON file");
            std::string block_name = root.get_object_key(0).as_string();
            std::cout << block_name << std::endl;
            if (block_name != "privateer_database")
                fail("not Privateer Database JSON file - top level key should equal to privateer_database\n");
            sajson::value entries = root.get_object_value(0);
            std::string entries_name = entries.get_object_key(0).as_string();
            std::cout << entries_name << std::endl;
            if (entries.get_type() != sajson::TYPE_OBJECT)
                fail("Expected TYPE_OBJECT, got " + json_type_as_string(entries.get_type()));
            
            sajson::value entries_array = entries.get_object_value(0);


            if (entries_array.get_type() != sajson::TYPE_ARRAY)
                fail("Expected TYPE_ARRAY, got " + json_type_as_string(entries_array.get_type()));
            
            for (size_t i = 0; i != entries_array.get_length(); ++i) 
            {
                sajson::value entry = entries_array.get_array_element(i);
                if (entry.get_type() != sajson::TYPE_OBJECT)
                fail("Expected TYPE_OBJECT, got " + json_type_as_string(entry.get_type()));
                
                Database temp;
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
                database.push_back(temp);
            }
        }

        inline std::vector<Database> read_json_insitu(char* buffer, size_t size, const std::string& name) 
        {
            std::vector<Database> database;
            sajson::document json = sajson::parse(sajson::dynamic_allocation(),
                                                sajson::mutable_string_view(size, buffer));
            if (!json.is_valid())
                fail(name, ":", std::to_string(json.get_error_line()), " error: ",
                    json.get_error_message_as_string());
            generate_database(database, json);
            return database;
        }

        inline std::vector<Database> read_json_file(const std::string& path) 
        {
            std::string path_copy = path;
            if(path_copy == "nopath" || path_copy.empty()) 
                {
                    std::string env(std::getenv ( "CLIBD" ));

                    path_copy = env + "/privateer_database.json";
                }

            std::cout << "Reading " << path_copy << " for Glycomics database" << std::endl;

            fileptr_t f = file_open(path_copy.c_str(), "rb");
            size_t buf_size = file_size(f.get(), path_copy);
            std::vector<char> buffer(buf_size);
            if (std::fread(buffer.data(), buffer.size(), 1, f.get()) != 1)
                fail(path_copy + ": fread failed");
            
            return read_json_insitu(buffer.data(), buffer.size(), path_copy);
        }
    }
}



#endif