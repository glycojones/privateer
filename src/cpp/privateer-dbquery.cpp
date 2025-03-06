// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York



#include "privateer-dbquery.h"

std::pair<std::string, std::string> output_dbquery(std::vector<privateer::json::GlycomicsDatabase>& glycomics_database, clipper::String glycanWURCS, clipper::MGlycan &currentGlycan)
{
    int valueLocation = privateer::util::find_index_of_value_from_wurcs(glycomics_database, glycanWURCS);
    std::string glytoucanID = glycomics_database[valueLocation].GlyTouCanID;
    std::string glyconnectID = glycomics_database[valueLocation].GlyConnectID;
    return std::make_pair(glytoucanID, glyconnectID);
}


