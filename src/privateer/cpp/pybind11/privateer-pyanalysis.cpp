// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#include "privateer-pyanalysis.h"

using namespace pybind11::literals;


///////////////////////////////////////////////// Class GlycosylationComposition ////////////////////////////////////////////////////////////////////

void privateer::pyanalysis::GlycosylationComposition::read_from_file( std::string path_to_model_file, std::string expression_system ) {

    if(path_to_model_file == "undefined")
    {
        throw std::invalid_argument( "No path was provided for model file input! Aborting." );
    }

    this->path_to_model_file = path_to_model_file;
    this->expression_system = expression_system;
    
    clipper::MiniMol mmol;
    clipper::MMDBfile mfile;
    clipper::String path_to_model_file_clipper = path_to_model_file;

    privateer::util::read_coordinate_file_mtz(mfile, mmol, path_to_model_file_clipper, true);

    this->mgl = clipper::MGlycology(mmol, expression_system);

    initialize_summary_of_detected_glycans(mgl);
}

void privateer::pyanalysis::GlycosylationComposition::initialize_summary_of_detected_glycans( clipper::MGlycology& mglObject )
{
    std::vector<clipper::MGlycan> list_of_glycans = mglObject.get_list_of_glycans();

    this->numberOfGlycanChains = list_of_glycans.size();

    auto list = pybind11::list();
    for(int i = 0; i < list_of_glycans.size(); i++)
    {
        std::string wurcsNotation = list_of_glycans[i].generate_wurcs();
        std::string kindOfGlycan = list_of_glycans[i].get_type();
        
        list_of_glycans[i].get_root_by_name();
        std::string proteinResidue = list_of_glycans[i].get_root().first.type().trim();
        std::string proteinResidueID = list_of_glycans[i].get_root().first.id().trim();
        std::string proteinChainID = list_of_glycans[i].get_chain().substr(0,1);

        auto rootSummary = pybind11::dict ("ProteinResidueType"_a=proteinResidue, "ProteinResidueID"_a=std::stoi(proteinResidueID), "ProteinChainID"_a=proteinChainID);
        
        std::vector<float> torsions = list_of_glycans[i].get_glycosylation_torsions();
        auto protein_glycan_linkage_torsion = pybind11::dict ("Phi"_a=torsions[0], "Psi"_a=torsions[1]);
        
        auto dict = pybind11::dict ("GlycanID"_a=i, "WURCS"_a=wurcsNotation, "GlycosylationType"_a=kindOfGlycan, "RootInfo"_a=rootSummary, "ProteinGlycanLinkageTorsion"_a=protein_glycan_linkage_torsion);
        list.append(dict);
    }
    this->glycosylationSummary = list;
}


privateer::pyanalysis::GlycanStructure privateer::pyanalysis::GlycosylationComposition::get_glycan(const int glycanID)
{
    if(glycanID >= numberOfGlycanChains || glycanID < 0)
    {
        throw std::invalid_argument( "Provided ID is out of bounds and exceeds/inceeds number of glycans detected in the model. \nInput: " + std::to_string(glycanID) + "\tPermitted Range: [0-" + std::to_string(numberOfGlycanChains - 1) + "]");
    }

    auto glycanObject = privateer::pyanalysis::GlycanStructure(mgl, glycanID);
    
    return glycanObject;
}
///////////////////////////////////////////////// Class GlycosylationComposition END ////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////// Class GlycanStructure ////////////////////////////////////////////////////////////////////
void privateer::pyanalysis::GlycanStructure::pyinit( const clipper::MGlycology& mgl, const int id )
{
    std::vector<clipper::MGlycan> list_of_glycans = mgl.get_list_of_glycans();

    clipper::MGlycan inputGlycan = list_of_glycans[id];

    this->glycan = inputGlycan;
    this->sugars_in_glycan = inputGlycan.get_sugars();
    
    this->id = id;
    this->numberOfSugars = inputGlycan.number_of_nodes();
    this->glycanWURCS = inputGlycan.generate_wurcs();

    auto uniqueMonosaccharidesPyList = pybind11::list();
    std::vector<std::string> uniqueMonosaccharidesVector = inputGlycan.obtain_unique_residue_codes();
    for(int i = 0; i < uniqueMonosaccharidesVector.size(); i++)
    {
        uniqueMonosaccharidesPyList.append(uniqueMonosaccharidesVector[i]);
    }
    this->uniqueMonosaccharides = uniqueMonosaccharidesPyList;
    this->numberOfGlycosidicBonds = inputGlycan.obtain_total_number_of_glycosidic_bonds();
    this->glycosylationType = inputGlycan.get_type();
    
    std::string proteinResidue = inputGlycan.get_root().first.type().trim();
    std::string proteinResidueID = inputGlycan.get_root().first.id().trim();
    std::string proteinChainID = inputGlycan.get_chain().substr(0,1);

    auto rootSummary = pybind11::dict("ProteinResidueType"_a=proteinResidue, "ProteinResidueID"_a=std::stoi(proteinResidueID), "ProteinChainID"_a=proteinChainID);
    this->rootSummary = rootSummary;

    std::vector<float> torsions = inputGlycan.get_glycosylation_torsions();
    auto protein_glycan_linkage_torsion = pybind11::dict("Phi"_a=torsions[0], "Psi"_a=torsions[1]);
    this->protein_glycan_linkage_torsion = protein_glycan_linkage_torsion;

    initialize_summary_of_glycan();
}

void privateer::pyanalysis::GlycanStructure::initialize_summary_of_glycan( )
{
    auto dict = pybind11::dict ("GlycanID"_a=id, "WURCS"_a=glycanWURCS, "UniqueMonosaccharides"_a=uniqueMonosaccharides, "TotalSugars"_a=numberOfSugars, "NumberOfGlycosidicBonds"_a=numberOfGlycosidicBonds, "GlycosylationType"_a=glycosylationType, "RootInfo"_a=rootSummary, "ProteinGlycanLinkageTorsion"_a=protein_glycan_linkage_torsion);

    this->glycanSummary = dict;
}

privateer::pyanalysis::CarbohydrateStructure privateer::pyanalysis::GlycanStructure::get_monosaccharide(const int sugarID)
{
    if(sugarID >= numberOfSugars || sugarID < 0)
    {
        throw std::invalid_argument( "Provided ID is out of bounds and exceeds/inceeds number of sugars detected in the glycan. \nInput: " + std::to_string(sugarID) + "\tPermitted Range: [0-" + std::to_string(numberOfSugars - 1) + "]");
    }

    auto sugarObject = privateer::pyanalysis::CarbohydrateStructure(glycan, sugarID);
    
    return sugarObject;
}
///////////////////////////////////////////////// Class GlycanStructure END ////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////// Class CarbohydrateStructure  ////////////////////////////////////////////////////////////////////
void privateer::pyanalysis::CarbohydrateStructure::pyinit(clipper::MGlycan& mglycan, const int id )
{
    // std::vector<clipper::MGlycan> list_of_glycans = mgl.get_list_of_glycans();
    std::vector<clipper::MSugar> list_of_sugars = mglycan.get_sugars();

    clipper::MSugar inputSugar = list_of_sugars[id];

    this->sugar_type = inputSugar.full_type();

    initialize_summary_of_sugar();
}


void privateer::pyanalysis::GlycanStructure::initialize_summary_of_glycan( )
{
    auto dict = pybind11::dict ("GlycanID"_a=id, "WURCS"_a=glycanWURCS, "UniqueMonosaccharides"_a=uniqueMonosaccharides, "TotalSugars"_a=numberOfSugars, "NumberOfGlycosidicBonds"_a=numberOfGlycosidicBonds, "GlycosylationType"_a=glycosylationType, "RootInfo"_a=rootSummary, "ProteinGlycanLinkageTorsion"_a=protein_glycan_linkage_torsion);

    this->glycanSummary = dict;
}




///////////////////////////////////////////////// Class CarbohydrateStructure END ////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////// PYBIND11 BINDING DEFINITIONS ////////////////////////////////////////////////////////////////////
namespace py = pybind11;
namespace pa = privateer::pyanalysis;

void init_pyanalysis(py::module& m)
{
    py::class_<pa::GlycosylationComposition>(m, "GlycosylationComposition")
        .def(py::init<>())
        .def(py::init<std::string&, std::string>(), py::arg("path_to_model_file")="undefined", py::arg("expression_system")="undefined")
        .def("get_path_of_model_file_used",  &pa::GlycosylationComposition::get_path_of_model_file_used)
        .def("get_expression_system_used",  &pa::GlycosylationComposition::get_expression_system_used)
        .def("get_number_of_glycan_chains_detected",  &pa::GlycosylationComposition::get_number_of_glycan_chains_detected)
        .def("get_summary_of_detected_glycans",  &pa::GlycosylationComposition::get_summary_of_detected_glycans)
        .def("get_glycan",  &pa::GlycosylationComposition::get_glycan);

    py::class_<pa::GlycanStructure>(m, "GlycanStructure")
        .def(py::init<>())
        .def(py::init<const clipper::MGlycology&, const int>())
        .def("get_glycan_id", &pa::GlycanStructure::get_glycan_id)
        .def("get_total_number_of_sugars", &pa::GlycanStructure::get_total_number_of_sugars)
        .def("get_wurcs_notation", &pa::GlycanStructure::get_wurcs_notation)
        .def("get_unique_monosaccharides", &pa::GlycanStructure::get_unique_monosaccharides)
        .def("get_total_of_glycosidic_bonds", &pa::GlycanStructure::get_total_of_glycosidic_bonds)
        .def("get_glycosylation_type", &pa::GlycanStructure::get_glycosylation_type)
        .def("get_root_info", &pa::GlycanStructure::get_root_info)
        .def("get_protein_glycan_linkage_torsions", &pa::GlycanStructure::get_protein_glycan_linkage_torsions)
        .def("get_glycan_summary", &pa::GlycanStructure::get_glycan_summary);
}

///////////////////////////////////////////////// PYBIND11 BINDING DEFINITIONS END////////////////////////////////////////////////////////////////////