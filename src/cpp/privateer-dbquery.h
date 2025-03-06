// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York



#ifndef DBQUERY_H_INCLUDED
#define DBQUERY_H_INCLUDED

#include <fstream>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <future>
#include <clipper/clipper.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-minimol.h>
#include "clipper-glyco.h"
#include "clipper-glyco_data.h"
#include "privateer-composition.h"
#include "privateer-lib.h"
#include <clipper/minimol/minimol_utils.h>
#include "privateer-json.h"

namespace privateer
{
    namespace dbquery
    {
        inline std::pair<std::string, std::string> output_dbquery(std::vector<privateer::json::GlycomicsDatabase>& glycomics_database, clipper::String glycanWURCS, clipper::MGlycan &currentGlycan)
        {
            int valueLocation = privateer::util::find_index_of_value_from_wurcs(glycomics_database, glycanWURCS);
            std::string glytoucanID = glycomics_database[valueLocation].GlyTouCanID;
            std::string glyconnectID = glycomics_database[valueLocation].GlyConnectID;
            return std::make_pair(glytoucanID, glyconnectID);
        }
    }
}



#endif