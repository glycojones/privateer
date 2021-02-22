// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#include "privateer-pyanalysis.h"

using namespace pybind11::literals;


// Class GlycanStructure

void privateer::pyanalysis::GlycosylationStructure::read_from_file( std::string path_to_model_file, std::string expression_system ) {

    if(path_to_model_file == "undefined")
    {
        throw std::invalid_argument( "Error: No path was provided for model file input! Aborting." );

    }

    this->path_to_model_file = path_to_model_file;
    this->expression_system = expression_system;
    
    clipper::MiniMol mmol;
    clipper::MMDBfile mfile;
    clipper::String path_to_model_file_clipper = path_to_model_file;

    privateer::util::read_coordinate_file_mtz(mfile, mmol, path_to_model_file_clipper, true);

    this->mgl = clipper::MGlycology(mmol, expression_system);

    initialize_mglycan_summary(mgl);
}

void privateer::pyanalysis::GlycosylationStructure::initialize_mglycan_summary( clipper::MGlycology& mglObject )
{
    std::vector<clipper::MGlycan> list_of_glycans = mglObject.get_list_of_glycans();

    this->numberOfGlycanChains = list_of_glycans.size();

    auto list = pybind11::list();
    for(int i = 0; i < list_of_glycans.size(); i++)
    {
        std::string wurcsNotation = list_of_glycans[i].generate_wurcs();
        std::string kindOfGlycan = list_of_glycans[i].get_type();
        std::string rootSummary = list_of_glycans[i].get_root_by_name();
        std::string torsions = list_of_glycans[i].print_torsions();
        
        auto dict = pybind11::dict ("id"_a=i, "wurcs"_a=wurcsNotation, "glycosylation_type"_a=kindOfGlycan, "root_info"_a=rootSummary, "linkage_torsion"_a=torsions);
        list.append(dict);
    }

    this->mglycanSummary = list;
}

namespace py = pybind11;
namespace pa = privateer::pyanalysis;

void init_pyanalysis(py::module& m)
{
    py::class_<pa::GlycosylationStructure>(m, "GlycosylationStructure")
        .def(py::init<>())
        .def(py::init<std::string&, std::string>(), py::arg("path_to_model_file")="undefined", py::arg("expression_system")="undefined")
        .def("get_path_of_model_file_used",  &pa::GlycosylationStructure::get_path_of_model_file_used)
        .def("get_expression_system_used",  &pa::GlycosylationStructure::get_expression_system_used)
        .def("get_number_of_glycan_chains_detected",  &pa::GlycosylationStructure::get_number_of_glycan_chains_detected)
        .def("get_summary_of_detected_glycans",  &pa::GlycosylationStructure::get_summary_of_detected_glycans);
}