// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#ifndef PRIVATEER_MODELLING_H_INCLUDED
#define PRIVATEER_MODELLING_H_INCLUDED

#include <pybind11/pybind11.h>
#include <fstream>
#include <string>
#include <locale>
#include <algorithm>
#include "clipper-glyco.h"
#include "privateer-restraints.h"
#include <gemmi/mmread.hpp>
#include <gemmi/monlib.hpp>
#include <gemmi/placeh.hpp>
#include <gemmi/fstream.hpp>
#define GEMMI_READ_CIF_IMPLEMENTATION
#include <gemmi/read_cif.hpp>
#define GEMMI_WRITE_IMPLEMENTATION
#include <gemmi/to_pdb.hpp>



namespace privateer {
	namespace interactions {

        void hydrogenate_input_model(std::string input_model);

  	}
}


#endif