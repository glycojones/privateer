
// Pybind11 Python bindings for the YSBL program Privateer
// (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2019 Jon Agirre
// York Structural Biology Laboratory
// The University of York
// mailto: jon.agirre@york.ac.uk
//

#include <pybind11/pybind11.h>
#include "privateer-lib.h"

// pybind11 module definition
//

PYBIND11_MODULE(privateer, m)
{
  m.doc() = "Privateer's Python interface.\nVersion history:\n- 2016-2018 MKIII (SWIG)\n2019- MKIV (pybind11)"; // docstring
  //m.def("carbname_of", &::privateer::scripting::carbname_of, "Returns IUPAC name for a monosaccharide's three-letter code (CCD)");
  m.def("svg_graphics_demo", &privateer::scripting::svg_graphics_demo, "SVG graphics demo");
}
