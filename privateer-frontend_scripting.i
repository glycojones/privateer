


/* privateer-python_interface.i */
%module privateer
%include "std_string.i"
%include "std_vector.i"

%{

#include "clipper-glyco_data.h"
bool clipper::data::found_in_database ( std::string name );

%}

namespace clipper 
{
    namespace data
    {
        bool found_in_database ( std::string name );
    }
}

