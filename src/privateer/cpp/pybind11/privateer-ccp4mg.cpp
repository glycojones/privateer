// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "privateer-lib.h"

using namespace pybind11::literals;
namespace py=pybind11;

// pybind11 module definition
//
void init_ccp4mg(py::module& m)
{
//   m.doc() = "Privateer's Python interface for ccp4mg.\nVersion history:\n- 2016-2018 MKIII (SWIG)\n2019-present MKIV (pybind11)"; // docstring

  m.def("carbname_of",
        &privateer::scripting::carbname_of,
        "Returns IUPAC name for a monosaccharide's three-letter code (CCD)" );

  m.def("found_in_database",
        &privateer::scripting::found_in_database,
        "Returns boolean if found in Privateer's database" );

  m.def("svg_graphics_demo",
        &privateer::scripting::svg_graphics_demo,
        "Creates an SVG file with a SNFG demo",
        "original_colour_scheme"_a = true,
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


  pybind11::enum_<privateer::glycanbuilderplot::Colour>(m, "Colour")
            .value("blue",    privateer::glycanbuilderplot::blue)
            .value("red" ,    privateer::glycanbuilderplot::red )
            .value("yellow" , privateer::glycanbuilderplot::yellow )
            .value("orange" , privateer::glycanbuilderplot::orange )
            .value("green" ,  privateer::glycanbuilderplot::green )
            .value("purple" , privateer::glycanbuilderplot::purple )
            .value("cyan" ,   privateer::glycanbuilderplot::cyan )
            .value("tan" ,    privateer::glycanbuilderplot::tan )
            .value("black" ,  privateer::glycanbuilderplot::black )
            .value("white" ,  privateer::glycanbuilderplot::white )
            .export_values();

  pybind11::enum_<privateer::Conformation>(m, "Conformation")
            .value("pyranose_4C1", privateer::pyranose_4C1 )
            .value("pyranose_1C4", privateer::pyranose_1C4 )
            .value("pyranose_3OB", privateer::pyranose_3OB )
            .value("pyranose_B25", privateer::pyranose_B25 )
            .value("pyranose_14B", privateer::pyranose_14B )
            .value("pyranose_B3O", privateer::pyranose_B3O )
            .value("pyranose_25B", privateer::pyranose_25B )
            .value("pyranose_B14", privateer::pyranose_B14 )
            .value("pyranose_OE",  privateer::pyranose_OE  )
            .value("pyranose_E5",  privateer::pyranose_E5  )
            .value("pyranose_4E",  privateer::pyranose_4E  )
            .value("pyranose_E3",  privateer::pyranose_E3  )
            .value("pyranose_2E",  privateer::pyranose_2E  )
            .value("pyranose_E1",  privateer::pyranose_E1  )
            .value("pyranose_3E",  privateer::pyranose_3E  )
            .value("pyranose_E2",  privateer::pyranose_E2  )
            .value("pyranose_1E",  privateer::pyranose_1E  )
            .value("pyranose_EO",  privateer::pyranose_EO  )
            .value("pyranose_5E",  privateer::pyranose_5E  )
            .value("pyranose_E4",  privateer::pyranose_E4  )
            .value("pyranose_OH5", privateer::pyranose_OH5 )
            .value("pyranose_4H5", privateer::pyranose_4H5 )
            .value("pyranose_4H3", privateer::pyranose_4H3 )
            .value("pyranose_2H3", privateer::pyranose_2H3 )
            .value("pyranose_2H1", privateer::pyranose_2H1 )
            .value("pyranose_OH1", privateer::pyranose_OH1 )
            .value("pyranose_3H2", privateer::pyranose_3H2 )
            .value("pyranose_1H2", privateer::pyranose_1H2 )
            .value("pyranose_1HO", privateer::pyranose_1HO )
            .value("pyranose_5HO", privateer::pyranose_5HO )
            .value("pyranose_5H4", privateer::pyranose_5H4 )
            .value("pyranose_3H4", privateer::pyranose_3H4 )
            .value("pyranose_OS2", privateer::pyranose_OS2 )
            .value("pyranose_1S5", privateer::pyranose_1S5 )
            .value("pyranose_1S3", privateer::pyranose_1S3 )
            .value("pyranose_2SO", privateer::pyranose_2SO )
            .value("pyranose_5S1", privateer::pyranose_5S1 )
            .value("pyranose_3S1", privateer::pyranose_3S1 )
            .value("furanose_3T2", privateer::furanose_3T2 )
            .value("furanose_3EV", privateer::furanose_3EV )
            .value("furanose_3T4", privateer::furanose_3T4 )
            .value("furanose_EV4", privateer::furanose_EV4 )
            .value("furanose_OT4", privateer::furanose_OT4 )
            .value("furanose_OEV", privateer::furanose_OEV )
            .value("furanose_OT1", privateer::furanose_OT1 )
            .value("furanose_EV1", privateer::furanose_EV1 )
            .value("furanose_2T1", privateer::furanose_2T1 )
            .value("furanose_2EV", privateer::furanose_2EV )
            .value("furanose_2T3", privateer::furanose_2T3 )
            .value("furanose_EV3", privateer::furanose_EV3 )
            .value("furanose_4T3", privateer::furanose_4T3 )
            .value("furanose_4EV", privateer::furanose_4EV )
            .value("furanose_4TO", privateer::furanose_4TO )
            .value("furanose_EVO", privateer::furanose_EVO )
            .value("furanose_1TO", privateer::furanose_1TO )
            .value("furanose_1EV", privateer::furanose_1EV )
            .value("furanose_1T2", privateer::furanose_1T2 )
            .value("furanose_EV2", privateer::furanose_EV2 )
            .export_values();

  m.def("get_colour",
        &privateer::glycanbuilderplot::get_colour,
        "Returns RGB colour values in two keys (SNFG or Privateer)",
        "colour"_a,
        "original_style"_a,
        "inverted"_a);

  m.def("print_wurcs",
        &privateer::scripting::print_wurcs,
        "Returns a WURCS string of all glycans in the glycoprotein model",
        "pdb_filename"_a,
        "expression_system"_a = "undefined");
}
