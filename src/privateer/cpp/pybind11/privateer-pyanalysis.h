// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York



#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <exception>
#include "clipper-glyco.h"
#include "privateer-lib.h"


#ifndef PRIVATEER_PYANALYSIS_H_INCLUDED
#define PRIVATEER_PYANALYSIS_H_INCLUDED

namespace privateer {

  namespace pyanalysis {

    class GlycosylationStructure {
      public:
        GlycosylationStructure() { };
        GlycosylationStructure(std::string& path_to_model_file, std::string expression_system) {
          this->read_from_file ( path_to_model_file, expression_system );
        };
        ~GlycosylationStructure() { };
        void read_from_file( std::string path_to_model_file, std::string expression_system );
        void initialize_mglycan_summary( clipper::MGlycology& mglObject );

        std::string get_path_of_model_file_used ( ) { return path_to_model_file; };
        std::string get_expression_system_used ( ) { return expression_system; };
        int get_number_of_glycan_chains_detected ( ) { return numberOfGlycanChains; };
        pybind11::list get_summary_of_detected_glycans () { return mglycanSummary; };

      private:
        std::string path_to_model_file;
        std::string expression_system;
        int numberOfGlycanChains;
        pybind11::list mglycanSummary;
        clipper::MGlycology mgl;
    };
  }
}

#endif
