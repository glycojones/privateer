#include <privateer-lib.h>
#include <privateer-json.h>

struct TorsionEntry {
  std::string sugar_1;
  std::string sugar_2;
  std::string atom_number_1;
  std::string atom_number_2;
  float phi;
  float psi;
};

struct TableEntry
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


std::vector<TableEntry> validate(const std::string &file, const std::string &name)
{

  char *c_data = (char *)file.c_str();
  size_t size = file.length();

  if (size == 0) {
    return {};
  }

  ::gemmi::Structure structure = ::gemmi::read_structure_from_char_array(c_data, size, name);
  clipper::GEMMIfile gemmi_file;
  // clipper::GemmiStructure *gemmi_structure = &gemmi_file;
  gemmi_file.set_gemmi_structure(structure);
//   std::cout << "[Privateer] Structure successfully imported" << std::endl;

  // CHECK FOR AN PROVIDED MTZ
  std::string filename = "/input.mtz";
  std::ifstream infile(filename);
  if (infile.good()) {
    // clipper::CCP4MTZfile mtzin;
    // mtzin.open_read( filename );
    // clipper::HKL_info myhkl();
    // mtzin.import_hkl_info( myhkl );
    // clipper::HKL_data<clipper::data32::F_phi> fphidata( myhkl );
    // mtzin.import_hkl_data( fphidata, "/*/*/[FWT,PHWT]" );
//     std::cout << "Found the input mtz, ready to use..." << std::endl;
  } else {
//     std::cout << "[Privateer] No MTZ supplied, continuing..." << std::endl;
  }

  clipper::MiniMol mol;
  gemmi_file.import_minimol(mol);

//   std::cout << "[Privateer] Molecule generated" << std::endl;

  privateer::json::GlobalTorsionZScore torsions_zscore_database = privateer::json::read_json_file_for_torsions_zscore_database("/data/privateer_data/privateer_torsions_z_score_database.json");

  const clipper::MAtomNonBond &manb = clipper::MAtomNonBond(mol, 1.0); // was 1.0

//   clipper::MGlycology mgl = clipper::MGlycology(mol, false, ""); <- use this constructor if you do not want to use torsions DB
    clipper::MGlycology mgl = clipper::MGlycology(mol, manb, torsions_zscore_database, false);

  std::vector<clipper::MGlycan> list_of_glycans = mgl.get_list_of_glycans();

  std::vector<TableEntry> table_list  = {};

  if (list_of_glycans.size() > 0)
  {
    clipper::String current_chain = "";

    for (int i = 0; i < list_of_glycans.size(); i++)
    {
      clipper::String wurcs_string;
      if (current_chain != list_of_glycans[i].get_chain())
      {
        current_chain = list_of_glycans[i].get_chain();
      }
      wurcs_string = list_of_glycans[i].generate_wurcs();

      privateer::glycanbuilderplot::GlycanErrorCount* err = new privateer::glycanbuilderplot::GlycanErrorCount;

      privateer::glycanbuilderplot::Plot plot(true, true, list_of_glycans[i].get_root_by_name());
      plot.plot_glycan(list_of_glycans[i], err);

      std::ostringstream os;
      os << list_of_glycans[i].get_root_for_filename() << ".svg";

      TableEntry table_entry;
      table_entry.svg = plot.write_to_string();

      // GlycanData glycan_data = query_glycomics_database(list_of_glycans[i], wurcs_string, importedDatabase);

      // table_entry.glyconnect_id = glycan_data.glyconnect_id;
      // table_entry.glytoucan_id = glycan_data.glytoucan_id;
      table_entry.wurcs = wurcs_string;
      table_entry.chain = current_chain;
      table_entry.id = list_of_glycans[i].get_root_by_name();

      table_entry.torsion_err = err->torsion_err;
      table_entry.conformation_err = err->conformation_err;
      table_entry.anomer_err = err->anomer_err;
      table_entry.puckering_err = err->puckering_err;
      table_entry.chirality_err = err->chirality_err;


      std::vector<clipper::MGlycan::MGlycanTorsionSummary> torsion_list = list_of_glycans[i].return_torsion_summary_within_glycan();
       for(int i = 0; i < torsion_list.size(); i++) {
        for(int j = 0; j < torsion_list[i].combined_torsions.size(); j++)
        {
            std::pair<std::pair<std::string, std::string>, std::vector<std::pair<float,float>>> torsion = torsion_list[i].combined_torsions[j];
            for (int k = 0; k < torsion.second.size(); k++) {
              TorsionEntry te;
              te.sugar_1 = torsion_list[i].first_residue_name;
              te.sugar_2 = torsion_list[i].second_residue_name;
              te.atom_number_1 = torsion.first.first;
              te.atom_number_2 = torsion.first.second;
              te.phi = torsion.second[k].first;
              te.psi = torsion.second[k].second;
              table_entry.torsions.emplace_back(te);
            }
        }
       }

      // table_entry.description = list_of_glycans[i].get_description();
      table_list.emplace_back(table_entry);
      delete err;
      // svg_list.emplace_back(plot.write_to_string());
    }

    return table_list;
  }
//   std::cout << "[Privateer] No Glycans Found" << std::endl;
  return {};
}