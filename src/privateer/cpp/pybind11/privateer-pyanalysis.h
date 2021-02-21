// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York



#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include "clipper-glyco.h"

#ifndef PRIVATEER_PYANALYSIS_H_INCLUDED
#define PRIVATEER_PYANALYSIS_H_INCLUDED

namespace privateer {

  namespace pyanalysis {

    std::string check_monlib_access ();

    class GlycosylationComposition {
      public:
        GlycosylationComposition() { };
      private:
    };

    class GlycanStructure {
      public:
        GlycanStructure() { };
        
      private:
    };

    class CarbohydrateStructure {
    public:
        CarbohydrateStructure() { };
      private:
    };


  }
}

#endif
