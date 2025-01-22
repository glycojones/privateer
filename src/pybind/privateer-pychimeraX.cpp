#include <pybind11/pybind11.h>
#include <privateer-lib.h>
#include <privateer-json.h>

using namespace pybind11::literals;

pybind11::list validate(std::string& path_to_model_file, std::string& path_to_zscores)
{
    clipper::MiniMol mmol;
    clipper::MMDBfile mfile;
    clipper::String path_to_model_file_clipper = path_to_model_file;

    privateer::util::read_coordinate_file_mtz(mfile, mmol, path_to_model_file_clipper, true);

//   std::cout << "[Privateer] Molecule generated" << std::endl;

    privateer::json::GlobalTorsionZScore torsions_zscore_database = privateer::json::read_json_file_for_torsions_zscore_database(path_to_zscores);
    //privateer::json::GlobalTorsionZScore torsions_zscore_database = privateer::json::read_json_file_for_torsions_zscore_database("/Users/lah583/Development/privateer_chimeraX_bundle/data/linkage_torsions/privateer_torsions_z_score_database.json");
    const clipper::MAtomNonBond &manb = clipper::MAtomNonBond(mmol, 1.0); // was 1.0

//   clipper::MGlycology mgl = clipper::MGlycology(mol, false, ""); <- use this constructor if you do not want to use torsions DB
    clipper::MGlycology mgl = clipper::MGlycology(mmol, manb, torsions_zscore_database, false);

    std::vector<clipper::MGlycan> list_of_glycans = mgl.get_list_of_glycans();

    auto resultslist = pybind11::list();
    if (list_of_glycans.size() > 0)
    {
        std::string current_chain = "";
        for (int i = 0; i < list_of_glycans.size(); i++)
        {
            std::string wurcs_string;
            if (current_chain != list_of_glycans[i].get_chain())
            {
                current_chain = list_of_glycans[i].get_chain();
            }
            wurcs_string = list_of_glycans[i].generate_wurcs();
            std::string kindOfGlycan = list_of_glycans[i].get_type();

            privateer::glycanbuilderplot::GlycanErrorCount* err = new privateer::glycanbuilderplot::GlycanErrorCount;

            privateer::glycanbuilderplot::Plot plot(true, true, list_of_glycans[i].get_root_by_name());
            plot.plot_glycan(list_of_glycans[i], err);

            std::ostringstream os;
            os << list_of_glycans[i].get_root_for_filename() << ".svg";




            std::vector<clipper::MGlycan::MGlycanTorsionSummary> torsion_list = list_of_glycans[i].return_torsion_summary_within_glycan();
            auto torsionlist = pybind11::list();
            for(int i = 0; i < torsion_list.size(); i++) {
                for(int j = 0; j < torsion_list[i].combined_torsions.size(); j++)
                {
                    std::pair<std::pair<std::string, std::string>, std::vector<std::pair<float,float>>> torsion = torsion_list[i].combined_torsions[j];
                    for (int k = 0; k < torsion.second.size(); k++) {
                        std::string sugar_1 = torsion_list[i].first_residue_name;
                        std::string sugar_2 = torsion_list[i].second_residue_name;
                        std::string atom_number_1 = torsion.first.first;
                        std::string atom_number_2 = torsion.first.second;
                        float phi = torsion.second[k].first;
                        float psi = torsion.second[k].second;
                        auto torsiondict = pybind11::dict("sugar_1"_a=sugar_1,"sugar_2"_a=sugar_2,
                                                            "atom_number_1"_a=atom_number_1, "atom_number_2"_a=atom_number_2,
                                                            "phi"_a=phi, "psi"_a=psi);
                        torsionlist.append(torsiondict);
                    }
                }
            }

            // table_entry.description = list_of_glycans[i].get_description();
            auto resultsdict = pybind11::dict ("GlycanNum"_a=i, "WURCS"_a=wurcs_string, "GlycosylationType"_a=kindOfGlycan, "RootID"_a=list_of_glycans[i].get_root_by_name(), "glycanChainID"_a=current_chain,
                                                "TorsionErr"_a=err->torsion_err, "ConformationErr"_a=err->conformation_err, "AnomerErr"_a=err->anomer_err, "PuckeringErr"_a=err->puckering_err, "ChiralityErr"_a=err->chirality_err,
                                                "Torsions"_a=torsionlist, "svg"_a=plot.write_to_string());
            resultslist.append(resultsdict);

            //table_list.emplace_back(table_entry);
            delete err;
            // svg_list.emplace_back(plot.write_to_string());
        }

        //return table_list;
    }
//   std::cout << "[Privateer] No Glycans Found" << std::endl;
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
   py::arg("path_to_model_file"),py::arg("path_to_zscores"));
}