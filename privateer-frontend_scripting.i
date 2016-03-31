
/* privateer-frontend_scripting.i */

%module privateer
%include "std_string.i"
%include "std_vector.i"

%{

#include "clipper-glyco_data.h"
#include "privateer-lib.h"
bool clipper::data::found_in_database ( std::string name );

%}

namespace clipper 
{
    namespace data
    {
        bool found_in_database ( std::string name );
    }
}

namespace privateer
{
    namespace scripting
    {
        std::string get_annotated_glycans ( std::string pdb_filename );
    }
}

