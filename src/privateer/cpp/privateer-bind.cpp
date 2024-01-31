#include <emscripten/bind.h>
#include <privateer-lib.h>
#include <privateer-json.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>
#include <clipper/clipper-contrib.h>

#include <gemmi/mmread.hpp>
#include "clipper-glyco.h"
#include "clipper-glyco.cpp"

using namespace emscripten;

struct TorsionEntry
{
  std::string sugar_1;
  std::string sugar_2;
  std::string atom_number_1;
  std::string atom_number_2;
  float phi;
  float psi;
};

struct ResultsEntry
{
  std::string svg;
  std::string wurcs;
  std::string chain;
  std::string glyconnect_id = "NotFound";
  std::string glytoucan_id = "NotFound";
  std::string id;
  int torsion_err = 0;
  int conformation_err = 0;
  int anomer_err = 0;
  int puckering_err = 0;
  int chirality_err = 0;

  std::vector<TorsionEntry> torsions;
};

clipper::MiniMol read_molecule(const std::string &file, const std::string &name)
{
  char *c_data = (char *)file.c_str();
  size_t size = file.length();

  if (size == 0)
  {
    std::cout << "[Privateer] The supplied file has no content, returning with nothing." << std::endl;
    return {};
  }
  ::gemmi::Structure structure = ::gemmi::read_structure_from_char_array(c_data, size, name);
  std::cout << "structure read" << std::endl;
  clipper::GEMMIFile gemmi_file;
  clipper::GemmiStructure *gemmi_structure = &gemmi_file;
  gemmi_structure->structure_ = structure;
  clipper::MiniMol mol;
  gemmi_file.import_minimol(mol);
  return mol;
}

std::vector<clipper::MGlycan> calculate_validation(clipper::MiniMol &mol)
{
  privateer::json::GlobalTorsionZScore torsions_zscore_database = privateer::json::read_json_file_for_torsions_zscore_database("privateer_torsions_z_score_database.json");
  const clipper::MAtomNonBond &manb = clipper::MAtomNonBond(mol, 1.0);
  clipper::MGlycology mgl = clipper::MGlycology(mol, manb, torsions_zscore_database, false);
  return mgl.get_list_of_glycans();
}

std::vector<ResultsEntry> format_validation(std::vector<clipper::MGlycan> &glycans)
{

  if (glycans.size() == 0) return {};

  std::vector<ResultsEntry> results_list = {};

  clipper::String current_chain = "";

  for (int i = 0; i < glycans.size(); i++)
  {
    clipper::String wurcs_string;
    if (current_chain != glycans[i].get_chain())
    {
      current_chain = glycans[i].get_chain();
    }
    wurcs_string = glycans[i].generate_wurcs();

    privateer::glycanbuilderplot::GlycanErrorCount *err = new privateer::glycanbuilderplot::GlycanErrorCount;

    privateer::glycanbuilderplot::Plot plot(true, true, glycans[i].get_root_by_name());
    plot.plot_glycan(glycans[i], err);

    std::ostringstream os;
    os << glycans[i].get_root_for_filename() << ".svg";

    ResultsEntry results_entry;
    results_entry.svg = plot.write_to_string();
    results_entry.wurcs = wurcs_string;
    results_entry.chain = current_chain;
    results_entry.id = glycans[i].get_root_by_name();
    results_entry.torsion_err = err->torsion_err;
    results_entry.conformation_err = err->conformation_err;
    results_entry.anomer_err = err->anomer_err;
    results_entry.puckering_err = err->puckering_err;
    results_entry.chirality_err = err->chirality_err;

    std::vector<clipper::MGlycan::MGlycanTorsionSummary> torsion_list = glycans[i].return_torsion_summary_within_glycan();
    for (int i = 0; i < torsion_list.size(); i++)
    {
      for (int j = 0; j < torsion_list[i].combined_torsions.size(); j++)
      {
        std::pair<std::pair<std::string, std::string>, std::vector<std::pair<float, float>>> torsion = torsion_list[i].combined_torsions[j];
        for (int k = 0; k < torsion.second.size(); k++)
        {
          TorsionEntry te;
          te.sugar_1 = torsion_list[i].first_residue_name;
          te.sugar_2 = torsion_list[i].second_residue_name;
          te.atom_number_1 = torsion.first.first;
          te.atom_number_2 = torsion.first.second;
          te.phi = torsion.second[k].first;
          te.psi = torsion.second[k].second;
          results_entry.torsions.emplace_back(te);
        }
      }
    }

    results_list.emplace_back(results_entry);
    delete err;
  }
  return results_list;
}

std::vector<ResultsEntry> validate(const std::string &file, const std::string &name)
{

  clipper::MiniMol mol = read_molecule(file, name);
  std::vector<clipper::MGlycan> glycans = calculate_validation(mol);
  std::vector<ResultsEntry> results = format_validation(glycans);
  return results;
}

EMSCRIPTEN_BINDINGS(privateer_module)
{
  value_object<TorsionEntry>("TorsionEntry")
      .field("sugar_1", &TorsionEntry::sugar_1)
      .field("sugar_2", &TorsionEntry::sugar_2)
      .field("atom_number_1", &TorsionEntry::atom_number_1)
      .field("atom_number_2", &TorsionEntry::atom_number_2)
      .field("phi", &TorsionEntry::phi)
      .field("psi", &TorsionEntry::psi);

  register_vector<TorsionEntry>("vector<TorsionEntry>");

  value_object<ResultsEntry>("ResultsEntry")
      .field("svg", &ResultsEntry::svg)
      .field("wurcs", &ResultsEntry::wurcs)
      .field("chain", &ResultsEntry::chain)
      .field("glyconnect_id", &ResultsEntry::glyconnect_id)
      .field("glytoucan_id", &ResultsEntry::glytoucan_id)
      .field("id", &ResultsEntry::id)
      .field("torsion_err", &ResultsEntry::torsion_err)
      .field("conformation_err", &ResultsEntry::conformation_err)
      .field("anomer_err", &ResultsEntry::anomer_err)
      .field("puckering_err", &ResultsEntry::puckering_err)
      .field("chirality_err", &ResultsEntry::chirality_err)
      .field("torsions", &ResultsEntry::torsions);

  function("validate", &validate);
  register_vector<ResultsEntry>("Results");
}
