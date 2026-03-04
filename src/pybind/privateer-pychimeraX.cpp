#include <pybind11/pybind11.h>
#include <privateer-lib.h>
#include <privateer-json.h>
#include <privateer-dbquery.h>
#include "gemmi/mmread.hpp"
#include <clippergemmi/clipper-minimol.h>
#include <clippergemmi/clipper-gemmi.h>

using namespace pybind11::literals;

pybind11::list validate(std::string& path_to_model_file, std::string& path_to_zscores, std::string& path_to_glycomics)
{
    char *c_data = (char *)path_to_model_file.c_str();
    size_t size = path_to_model_file.length();

    gemmi::Structure structure = gemmi::read_structure_from_char_array(c_data, size, path_to_model_file);
    clipper::GEMMIfile gemmi_file;
    // clipper::GemmiStructure *gemmi_structure = &gemmi_file;
    gemmi_file.set_gemmi_structure(structure);

    clipper::MiniMol mmol;
    //clipper::MMDBfile mfile;
    //clipper::String path_to_model_file_clipper = path_to_model_file;
    gemmi_file.import_minimol(mmol);
    privateer::json::GlobalTorsionZScore torsions_zscore_database = privateer::json::read_json_file_for_torsions_zscore_database(path_to_zscores);
    std::vector<privateer::json::GlycomicsDatabase> glycomics_database = privateer::json::read_json_file_for_glycomics_database(path_to_glycomics);
    
    //privateer::util::read_coordinate_file_mtz(mfile, mmol, path_to_model_file_clipper, true);

    const clipper::MAtomNonBond &manb = clipper::MAtomNonBond(mmol, 1.0); // was 1.0
    clipper::MGlycology mgl = clipper::MGlycology(mmol, manb, torsions_zscore_database, false);


    std::vector<clipper::MGlycan> list_of_glycans = mgl.get_list_of_glycans();

    auto resultslist = pybind11::list();
    if (list_of_glycans.size() > 0)
    {
        std::string current_chain = "";
        for (int i = 0; i < list_of_glycans.size(); i++)
        {
            if (current_chain != list_of_glycans[i].get_chain())
            {
                current_chain = list_of_glycans[i].get_chain();
            }
            clipper::String wurcs = list_of_glycans[i].generate_wurcs();
            std::string wurcs_string = list_of_glycans[i].generate_wurcs();
            std::string kindOfGlycan = list_of_glycans[i].get_type();
            std::pair<std::string, std::string> GlycanIDs = privateer::dbquery::output_dbquery(glycomics_database, wurcs, list_of_glycans[i]);
            std::string GlyConnectID = GlycanIDs.second;
            std::string GlyToucanID = GlycanIDs.first;
            

            privateer::glycanbuilderplot::GlycanErrorCount* err = new privateer::glycanbuilderplot::GlycanErrorCount;

            privateer::glycanbuilderplot::Plot plot(true, true, list_of_glycans[i].get_root_by_name());
            plot.plot_glycan(list_of_glycans[i], err);

            std::vector<clipper::MGlycan::MGlycanTorsion> torsion_list = list_of_glycans[i].return_torsion_collection(); //FLAG: Edit this loop based on the new torsion collection structure I added.
            std::vector<clipper::MGlycan::MGlycanTorsionSummary> torsion_summary = list_of_glycans[i].return_torsion_summary_within_glycan();
            auto torsionlist = pybind11::list();
            for(int j = 0; j < torsion_list.size(); j++) {
                std::string type = torsion_list[j].type;
                std::string sugar_1 = torsion_list[j].first_residue_name;
                std::string sugar_2 = torsion_list[j].second_residue_name;
                std::string firstchainID = torsion_list[j].firstsugchainID.substr(0,1);
                std::string secondchainID = torsion_list[j].secondsugchainID.substr(0,1);
                int firstresID = torsion_list[j].first_seqid;
                int secondresID = torsion_list[j].second_seqid;
                float phi = torsion_list[j].phi;
                float psi = torsion_list[j].psi;
                std::string atom_number_1 = torsion_list[j].donor_atom;
                std::string atom_number_2 = torsion_list[j].acceptor_atom;
                std::vector < clipper::MSugar > list_of_sugars = list_of_glycans[i].get_sugars();
                clipper::MAtom pos1, pos2;
                float pos1_x, pos1_y, pos1_z, pos2_x, pos2_y, pos2_z;
                if (type == "protein-sugar")
                {
                    for(int m = 0; m < list_of_sugars.size(); m++) {
                        clipper::MSugar sugar = list_of_sugars[m];
                        pos1 = torsion_summary[j].atoms[0].first;
                        if (sugar.chain_id().trim().substr(0,1) == secondchainID && sugar.get_seqnum() == secondresID && sugar.type().trim() == sugar_2){
                            if (atom_number_2 == "6"){
                                std::string atom_name_2 = "C5";
                                pos2 = sugar[sugar.lookup(atom_name_2,clipper::MM::ANY)];
                            }
                            else{
                                std::string atom_name_2 = "C" + atom_number_2;
                                pos2 = sugar[sugar.lookup(atom_name_2,clipper::MM::ANY)];   
                            }
                        }
                    }
                }
                else{
                    for(int m = 0; m < list_of_sugars.size(); m++) {
                        clipper::MSugar sugar = list_of_sugars[m];
                        if (sugar.chain_id().trim().substr(0,1) == firstchainID && sugar.get_seqnum() == firstresID && sugar.type().trim() == sugar_1){
                            if (atom_number_1 == "6"){
                                std::string atom_name_1 = "C5";
                                pos1 = sugar[sugar.lookup(atom_name_1,clipper::MM::ANY)];
                            }
                            else{
                                std::string atom_name_1 = "C" + atom_number_1;
                                pos1 = sugar[sugar.lookup(atom_name_1,clipper::MM::ANY)];   
                            }
                        }
                        if (sugar.chain_id().trim().substr(0,1) == secondchainID && sugar.get_seqnum() == secondresID && sugar.type().trim() == sugar_2){
                            if (atom_number_2 == "6"){
                                std::string atom_name_2 = "C5";
                                pos2 = sugar[sugar.lookup(atom_name_2,clipper::MM::ANY)];
                            }
                            else{
                                std::string atom_name_2 = "C" + atom_number_2;
                                pos2 = sugar[sugar.lookup(atom_name_2,clipper::MM::ANY)];   
                            }
                        }
                    }
                }
                pos1_x = pos1.coord_orth().x();
                pos1_y = pos1.coord_orth().y();
                pos1_z = pos1.coord_orth().z();
                pos2_x = pos2.coord_orth().x();
                pos2_y = pos2.coord_orth().y();
                pos2_z = pos2.coord_orth().z();
                auto torsiondict = pybind11::dict("chainID"_a = secondchainID, "sugar_1_resID"_a = firstresID, "sugar_2_resID"_a = secondresID, "sugar_1"_a=sugar_1,"sugar_2"_a=sugar_2,
                                                    "atom_number_1"_a=atom_number_1, "atom_number_2"_a=atom_number_2, "phi"_a=phi, "psi"_a=psi,
                                                    "x1"_a= pos1_x, "y1"_a= pos1_y, "z1"_a= pos1_z, "x2"_a= pos2_x, "y2"_a= pos2_y, "z2"_a= pos2_z);
                torsionlist.append(torsiondict);
            }

/*
            std::vector<clipper::MGlycan::MGlycanTorsionSummary> torsion_list = list_of_glycans[i].return_torsion_summary_within_glycan();
            auto torsionlist = pybind11::list();
            for(int j = 0; j < torsion_list.size(); j++) {
                std::string secondchainID = torsion_list[j].secondsugchainID.substr(0,1);
                std::string firstchainID = torsion_list[j].firstsugchainID.substr(0,1);
                int secondresID = torsion_list[j].secondsugresID;
                int firstresID = torsion_list[j].firstsugresID;
                std::string sugar_1 = torsion_list[j].first_residue_name;
                std::string sugar_2 = torsion_list[j].second_residue_name;
                std::string type = torsion_list[j].type;
                for(int k = 0; k < torsion_list[j].combined_torsions.size(); k++)
                {
                    std::pair<std::pair<std::string, std::string>, std::vector<std::pair<float,float>>> torsion = torsion_list[j].combined_torsions[k];
                    std::string atom_number_1 = torsion.first.first;
                    std::string atom_number_2 = torsion.first.second;
                    for (int l = 0; l < torsion.second.size(); l++) {
                        float phi = torsion.second[l].first;
                        float psi = torsion.second[l].second;
                        std::vector < clipper::MSugar > list_of_sugars = list_of_glycans[i].get_sugars();
                        clipper::MAtom pos1, pos2;
                        float pos1_x, pos1_y, pos1_z, pos2_x, pos2_y, pos2_z;
                        if (type == "protein-sugar")
                        {
                            for(int m = 0; m < list_of_sugars.size(); m++) {
                                clipper::MSugar sugar = list_of_sugars[m];
                                pos1 = torsion_list[j].atoms[0].first;
                                if (sugar.chain_id().trim().substr(0,1) == secondchainID && sugar.get_seqnum() == secondresID && sugar.type().trim() == sugar_2){
                                    if (atom_number_2 == "6"){
                                        std::string atom_name_2 = "C5";
                                        pos2 = sugar[sugar.lookup(atom_name_2,clipper::MM::ANY)];
                                    }
                                    else{
                                        std::string atom_name_2 = "C" + atom_number_2;
                                        pos2 = sugar[sugar.lookup(atom_name_2,clipper::MM::ANY)];   
                                    }
                                }
                            }
                        }
                        else{
                            for(int m = 0; m < list_of_sugars.size(); m++) {
                                clipper::MSugar sugar = list_of_sugars[m];
                                if (sugar.chain_id().trim().substr(0,1) == firstchainID && sugar.get_seqnum() == firstresID && sugar.type().trim() == sugar_1){
                                    if (atom_number_1 == "6"){
                                        std::string atom_name_1 = "C5";
                                        pos1 = sugar[sugar.lookup(atom_name_1,clipper::MM::ANY)];
                                    }
                                    else{
                                        std::string atom_name_1 = "C" + atom_number_1;
                                        pos1 = sugar[sugar.lookup(atom_name_1,clipper::MM::ANY)];   
                                    }
                                }
                                if (sugar.chain_id().trim().substr(0,1) == secondchainID && sugar.get_seqnum() == secondresID && sugar.type().trim() == sugar_2){
                                    if (atom_number_2 == "6"){
                                        std::string atom_name_2 = "C5";
                                        pos2 = sugar[sugar.lookup(atom_name_2,clipper::MM::ANY)];
                                    }
                                    else{
                                        std::string atom_name_2 = "C" + atom_number_2;
                                        pos2 = sugar[sugar.lookup(atom_name_2,clipper::MM::ANY)];   
                                    }
                                }
                            }
                        }
                        pos1_x = pos1.coord_orth().x();
                        pos1_y = pos1.coord_orth().y();
                        pos1_z = pos1.coord_orth().z();
                        pos2_x = pos2.coord_orth().x();
                        pos2_y = pos2.coord_orth().y();
                        pos2_z = pos2.coord_orth().z();
                        
                        auto torsiondict = pybind11::dict("chainID"_a = secondchainID, "sugar_1_resID"_a = firstresID, "sugar_2_resID"_a = secondresID, "sugar_1"_a=sugar_1,"sugar_2"_a=sugar_2,
                                                            "atom_number_1"_a=atom_number_1, "atom_number_2"_a=atom_number_2, "phi"_a=phi, "psi"_a=psi, 
                                                            "x1"_a= pos1_x, "y1"_a= pos1_y, "z1"_a= pos1_z, "x2"_a= pos2_x, "y2"_a= pos2_y, "z2"_a= pos2_z);
                        torsionlist.append(torsiondict);
                    }
                }
            }
*/
            auto sugarcoordlist = pybind11::list();
            std::vector < clipper::MSugar > list_of_sugars = list_of_glycans[i].get_sugars();
            for(int m = 0; m < list_of_sugars.size(); m++) {
                clipper::MSugar sugar = list_of_sugars[m];
                clipper::Coord_orth sugarcentre = sugar.ring_centre();
                float sugarcentre_x = sugarcentre.x();
                float sugarcentre_y = sugarcentre.y();
                float sugarcentre_z = sugarcentre.z();
                clipper::Vec3<clipper::ftype> sugarplane = sugar.ring_mean_plane();
                float sugarplane_x = sugarplane[0];
                float sugarplane_y = sugarplane[1];
                float sugarplane_z = sugarplane[2];
                std::string sugarname = sugar.type().trim();
                std::string sugar_chain_ID = sugar.chain_id().trim().substr(0,1);
                int sugar_res_ID = sugar.get_seqnum();
                clipper::MAtom pos1 = sugar[sugar.lookup("C1",clipper::MM::ANY)];
                clipper::MAtom pos2 = sugar[sugar.lookup("C2",clipper::MM::ANY)];
                clipper::MAtom pos3 = sugar[sugar.lookup("C3",clipper::MM::ANY)];
                clipper::MAtom pos4 = sugar[sugar.lookup("C4",clipper::MM::ANY)];
                clipper::MAtom pos5 = sugar[sugar.lookup("C5",clipper::MM::ANY)];
                float C1_x = pos1.coord_orth().x();
                float C1_y = pos1.coord_orth().y();
                float C1_z = pos1.coord_orth().z();
                float C2_x = pos2.coord_orth().x();
                float C2_y = pos2.coord_orth().y();
                float C2_z = pos2.coord_orth().z();
                float C3_x = pos3.coord_orth().x();
                float C3_y = pos3.coord_orth().y();
                float C3_z = pos3.coord_orth().z();
                float C4_x = pos4.coord_orth().x();
                float C4_y = pos4.coord_orth().y();
                float C4_z = pos4.coord_orth().z();
                float C5_x = pos5.coord_orth().x();
                float C5_y = pos5.coord_orth().y();
                float C5_z = pos5.coord_orth().z();
                int count = 0;
                auto linkageatomslist = pybind11::list();
                for(int j = 0; j < torsion_list.size(); j++) {
                    std::string sugar_1 = torsion_list[j].first_residue_name;
                    std::string sugar_2 = torsion_list[j].second_residue_name;
                    std::string firstchainID = torsion_list[j].firstsugchainID.substr(0,1);
                    std::string secondchainID = torsion_list[j].secondsugchainID.substr(0,1);
                    int firstresID = torsion_list[j].first_seqid;
                    int secondresID = torsion_list[j].second_seqid;
                    if (sugar.chain_id().trim().substr(0,1) == firstchainID && sugar.get_seqnum() == firstresID && sugar.type().trim() == sugar_1){
                        count += 1;
                        linkageatomslist.append(torsion_list[j].donor_atom);
                        // FLAG: Add in here a list of atoms so that I know if both C2 and C3 are in use on square and diamond sugars
                    }
                    if (sugar.chain_id().trim().substr(0,1) == secondchainID && sugar.get_seqnum() == secondresID && sugar.type().trim() == sugar_2){
                        count +=1;
                        linkageatomslist.append(torsion_list[j].acceptor_atom);
                    }
                }
                auto coorddict = pybind11::dict("sugarname"_a = sugarname, "chainID"_a = sugar_chain_ID, "resID"_a = sugar_res_ID,
                                                "sugar_centre_x"_a = sugarcentre_x,"sugar_centre_y"_a = sugarcentre_y,"sugar_centre_z"_a = sugarcentre_z,
                                                "sugar_plane_i"_a = sugarplane_x, "sugar_plane_j"_a = sugarplane_y, "sugar_plane_k"_a = sugarplane_z, 
                                                "C1_x"_a=C1_x, "C1_y"_a=C1_y, "C1_z"_a=C1_z, "C2_x"_a=C2_x, "C2_y"_a=C2_y, "C2_z"_a=C2_z,
                                                "C3_x"_a=C3_x, "C3_y"_a=C3_y, "C3_z"_a=C3_z, "C4_x"_a=C4_x, "C4_y"_a=C4_y, "C4_z"_a=C4_z,
                                                "C5_x"_a=C5_x, "C5_y"_a=C5_y, "C5_z"_a=C5_z, "num_bonds"_a=count, "link_atoms"_a=linkageatomslist);
                sugarcoordlist.append(coorddict);

            }
            auto resultsdict = pybind11::dict ("GlycanNum"_a=i, "WURCS"_a=wurcs_string, "GlycosylationType"_a=kindOfGlycan, "RootID"_a=list_of_glycans[i].get_root_by_name(), "glycanChainID"_a=current_chain,
                                                "TorsionErr"_a=err->torsion_err, "ConformationErr"_a=err->conformation_err, "AnomerErr"_a=err->anomer_err, "PuckeringErr"_a=err->puckering_err, "ChiralityErr"_a=err->chirality_err,
                                                "Torsions"_a=torsionlist, "svg"_a=plot.write_to_string(),"GlyToucanID"_a = GlyToucanID, "GlyConnectID"_a = GlyConnectID, "Sugars"_a=sugarcoordlist);
            resultslist.append(resultsdict);

            delete err;
        }
    }
    return resultslist;
}

namespace py=pybind11;

PYBIND11_MODULE(privateer_core, m) {
    m.doc() = "Python wrapper for privateer_core(C++) exposed via pybind11.";

    py::register_exception_translator([](std::exception_ptr p) {
        try {
            if (p) std::rethrow_exception(p);
        } catch (const std::exception& e) {
            PyErr_SetString(PyExc_RuntimeError, e.what());
        }
    });
   m.def("validate",&validate, "A function that produces a validation report on the glycans in a model",
   py::arg("path_to_model_file"),py::arg("path_to_zscores"),py::arg("path_to_glycomics"));
}