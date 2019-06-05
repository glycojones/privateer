
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
#include <pybind11/stl.h>
#include "privateer-lib.h"
#include "privateer-restraints.h"

using namespace pybind11::literals;

// pybind11 module definition
//
PYBIND11_MODULE(privateer_core, m)
{
  m.doc() = "Privateer's Python interface.\nVersion history:\n- 2016-2018 MKIII (SWIG)\n2019-present MKIV (pybind11)"; // docstring

  m.def("carbname_of",
        &privateer::scripting::carbname_of,
        "Returns IUPAC name for a monosaccharide's three-letter code (CCD)" );

  m.def("found_in_database",
        &privateer::scripting::found_in_database,
        "Returns boolean if found in Privateer's database" );

  m.def("svg_graphics_demo",
        &privateer::scripting::svg_graphics_demo,
        "Creates an SVG file with a SNFG demo",
        "original_colour_scheme"_a = false,
        "inverted_background"_a = false );

  m.def("get_annotated_glycans",
        &privateer::scripting::get_annotated_glycans,
        "Produces XML with validation info for protein glycosylation",
        "pdb_filename"_a,
        "original_colour_scheme"_a = true,
        "expression_system"_a = "undefined");

  m.def("get_annotated_glycans_hierarchical",
        &privateer::scripting::get_annotated_glycans_hierarchical,
        "Produces XML with hierarchical validation info for protein glycosylation",
        "pdb_filename"_a,
        "original_colour_scheme"_a = true,
        "expression_system"_a = "undefined" );

  m.def("write_refmac_keywords",
        &privateer::scripting::write_refmac_keywords,
        "Writes refmac5 keywords",
        "code_list"_a );

  m.def("write_libraries",
        &privateer::scripting::write_libraries,
        "Writes refmac5 libraries",
        "code_list"_a,
        "esd"_a );

  m.def("test_monlib_access",
        &privateer::restraints::test_monlib_access,
        "Checks if the CCP4 monomer library is accessible via environment, returns boolean");

  pybind11::enum_<privateer::glycoplot::Colour>(m, "Colour")
            .value("blue",    privateer::glycoplot::blue)
            .value("red" ,    privateer::glycoplot::red )
            .value("yellow" , privateer::glycoplot::yellow )
            .value("orange" , privateer::glycoplot::orange )
            .value("green" ,  privateer::glycoplot::green )
            .value("purple" , privateer::glycoplot::purple )
            .value("cyan" ,   privateer::glycoplot::cyan )
            .value("tan" ,    privateer::glycoplot::tan )
            .value("black" ,  privateer::glycoplot::black )
            .value("white" ,  privateer::glycoplot::white )
            .export_values();

  m.def("get_colour",
        &privateer::glycoplot::get_colour,
        "Returns RGB colour values in two keys (SNFG or Privateer)",
        "colour"_a,
        "original_style"_a,
        "inverted"_a );

}
