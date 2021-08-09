// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#include "privateer-pymodelling.h"

privateer::pymodelling::Builder::Builder(std::string& path_to_receiving_model_file, std::string& path_to_donor_model_file, bool debug_output) 
{
    this->debug_output = debug_output;
    this->read_from_file ( path_to_receiving_model_file, path_to_donor_model_file, debug_output );
};

void privateer::pymodelling::Builder::read_from_file( std::string& path_to_receiving_model_file, std::string& path_to_donor_model_file, bool debug_output ) {

    if(path_to_receiving_model == "undefined")
    {
        throw std::invalid_argument( "No path was provided for model file input! Aborting." );
    }

    if(path_to_donor_model == "undefined")
    {
        throw std::invalid_argument( "No path was provided for model file input! Aborting." );
    }

    this->path_to_receiving_model = path_to_receiving_model_file;
    this->path_to_donor_model = path_to_donor_model_file;
    
    clipper::MiniMol receiver_mmol;
    clipper::MiniMol donor_mmol;
    clipper::MMDBfile receiver_mfile;
    clipper::MMDBfile donor_mfile;
    clipper::String path_to_model_file_clipper_receiver = path_to_receiving_model_file;
    clipper::String path_to_model_file_clipper_donor = path_to_donor_model_file;

    privateer::util::read_coordinate_file_mtz(receiver_mfile, receiver_mmol, path_to_model_file_clipper_receiver, true);
    privateer::util::read_coordinate_file_mtz(donor_mfile, donor_mmol, path_to_model_file_clipper_donor, true);

    this->imported_receiving_model = receiver_mmol;
    this->imported_donor_model = donor_mmol;
    this->imported_receiving_model_seq_info = get_protein_sequence_information(imported_receiving_model);

    privateer::pyanalysis::GlycosylationComposition glycosylation(path_to_donor_model, "undefined", false);
    this->glycan_summary_donor = glycosylation.get_summary_of_detected_glycans();
}

pybind11::list privateer::pymodelling::Builder::get_protein_sequence_information( clipper::MiniMol& input ) 
{
    pybind11::list output;
    for(int i = 0; i < input.size(); i++)
    {
        pybind11::list residues;
        for(int j = 0; j < input[i].size(); j++)
        {
            std::string residueType = input[i][j].type().trim();
            int residueSeqnum = input[i][j].seqnum();
            auto residue_dict = pybind11::dict("index"_a=j, "residueType"_a=residueType, "residueSeqnum"_a=residueSeqnum);
            residues.append(residue_dict);
        }
        std::string chainID = input[i].id().trim();
        auto mpolymer_dict = pybind11::dict("index"_a = i, "ChainID"_a=chainID, "Residues"_a=residues);
        output.append(mpolymer_dict);
    }
    return output;
}

void privateer::pymodelling::Builder::graft_glycan_to_receiver(int mglycanindex, int receiver_chain_index, int received_residue_index)
{
    privateer::modelling::Grafter Grafter(imported_receiving_model, imported_donor_model, debug_output);

    std::vector<clipper::MGlycan> donor_glycans = Grafter.get_donor_glycans();
    clipper::MGlycan glycan_to_graft = donor_glycans[mglycanindex];
    clipper::MPolymer converted_mglycan = Grafter.convert_mglycan_to_mpolymer(glycan_to_graft);

    clipper::MMonomer receiver_monomer = imported_receiving_model[receiver_chain_index][received_residue_index];
    clipper::MAtom receiver_atom = receiver_monomer.find("ND2", clipper::MM::ANY);
    clipper::Coord_orth target = receiver_atom.coord_orth();

    clipper::MAtom donor_atom = converted_mglycan[0].find("O1", clipper::MM::ANY);
    clipper::Coord_orth source = donor_atom.coord_orth();

    std::vector<float> torsions;
    torsions.push_back(-110.0);
    torsions.push_back(-175.0);

    Grafter.graft_mpolymer_to_receiving_model(target, source, torsions, converted_mglycan);
    this->imported_receiving_model = Grafter.get_receiving_model();
}

void privateer::pymodelling::Builder::export_grafted_model(std::string& output_path)
{
    clipper::MMDBfile pdbfile;
    pdbfile.export_minimol( imported_receiving_model );
    pdbfile.write_file( output_path );
}



///////////////////////////////////////////////// PYBIND11 BINDING DEFINITIONS ////////////////////////////////////////////////////////////////////
namespace py = pybind11;
namespace pm = privateer::pymodelling;

void init_pymodelling(py::module& m)
{
    py::class_<pm::Builder>(m, "Builder")
        .def(py::init<>())
        .def(py::init<std::string&, std::string&, bool>(), py::arg("path_to_receiving_model_file")="undefined", py::arg("path_to_donor_model")="undefined", py::arg("debug_output")=false)
        .def("get_path_of_receiving_model_file_used", &pm::Builder::get_path_of_receiving_model_file_used)
        .def("get_path_of_donor_model_file_used", &pm::Builder::get_path_of_donor_model_file_used)
        .def("get_receiving_model_sequence_info", &pm::Builder::get_receiving_model_sequence_info)
        .def("get_glycan_summary_from_donor", &pm::Builder::get_glycan_summary_from_donor)
        .def("graft_glycan_to_receiver", &pm::Builder::graft_glycan_to_receiver)
        .def("export_grafted_model", &pm::Builder::export_grafted_model);
}

///////////////////////////////////////////////// PYBIND11 BINDING DEFINITIONS END////////////////////////////////////////////////////////////////////

