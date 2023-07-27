#include <emscripten/bind.h>
#include <privateer-lib.h>
#include <clipper/clipper-minimol.h>
#include <gemmi/mmread.hpp>
#include "clipper-glyco.h"
#include "clipper-glyco.cpp"

clipper::MMDBfile mfile;
//   clipper::MiniMol mmol;
//   mfile.read_file( "input.pdb" );
//   mfile.import_minimol( mmol );
using namespace emscripten;

extern "C" std::vector<std::string> read_file(const std::string &file, const std::string &name)
{

  char *c_data = (char *)file.c_str();
  size_t size = file.length();
  ::gemmi::Structure structure = ::gemmi::read_structure_from_char_array(c_data, size, name);

  clipper::GEMMIFile gemmi_file; 
  clipper::GemmiStructure* gemmi_structure = &gemmi_file; 
  gemmi_structure->structure_ = structure; 

  clipper::MiniMol mol; 
  gemmi_file.import_minimol(mol);

  // clipper::Spacegroup spacegroup = clipper::Spacegroup(clipper::Spacegroup::TYPE::P1);
  // clipper::Cell_descr cell_desc = {structure.cell.a, structure.cell.b, structure.cell.c, structure.cell.alpha, structure.cell.beta, structure.cell.gamma};
  // clipper::Cell cell = clipper::Cell(cell_desc);
  // clipper::MiniMol mol(spacegroup, cell);

  // std::cout << "Created minimol" << std::endl;

  // for (int m = 0; m < structure.children().size(); m++)
  // {
  //   clipper::MModel model;
  //   gemmi::Model gemmi_model = structure.children()[m];

  //   for (int c = 0; c < gemmi_model.children().size(); c++)
  //   {
  //     gemmi::Chain gemmi_chain = gemmi_model.children()[c];
  //     clipper::MPolymer polymer;
  //     polymer.set_id(c);
  //     for (int r = 0; r < gemmi_chain.children().size(); r++)
  //     {
  //       gemmi::Residue gemmi_residue = gemmi_chain.children()[r];
  //       clipper::MMonomer monomer;
  //       monomer.set_id(r);
  //       // monomer.set_type(gemmi_residue.entity_type);
  //       for (int a = 0; a < gemmi_residue.children().size(); a++)
  //       {
  //         gemmi::Atom gemmi_atom = gemmi_residue.children()[a];
  //         clipper::MAtom atom;
  //         clipper::Coord_orth coordinate = {gemmi_atom.pos.x, gemmi_atom.pos.y, gemmi_atom.pos.z};
  //         atom.set_coord_orth(coordinate);

  //         atom.set_id(a);
  //         atom.set_name(gemmi_atom.name);
  //         atom.set_occupancy(gemmi_atom.occ);

  //         atom.set_element(gemmi_atom.element.name());
  //         monomer.insert(atom);
  //       }
  //       polymer.insert(monomer);
  //     }
  //     model.insert(polymer);
  //   }
  //   mol.model() = model;
  // }

  // for (int p = 0; p < mol.model().size(); p++) { 
  //   for (int m = 0; m < mol.model()[p].size(); m++) { 
  //     for (int a = 0; a < mol.model()[p][m].size(); a++) { 
  //       std::cout << mol.model()[p].id() << "\t" << mol.model()[p][m].id() << "\t" << mol.model()[p][m][a].id() << std::endl;
  //     }
  //   }
  // }

  privateer::json::GlobalTorsionZScore torsions_zscore_database;
  // = privateer::json::read_json_file_for_torsions_zscore_database("/home/jordan/dev/privateer_wasm/privateer/data/linkage_torsions/privateer_torsions_z_score_database.json");
  const clipper::MAtomNonBond& manb = clipper::MAtomNonBond( mol, 1.0 ); // was 1.0

    clipper::MGlycology mgl = clipper::MGlycology(mol, manb, torsions_zscore_database, false);
  // std::cout << mgl.get_list_of_glycans().size() << std::endl;

  std::vector<clipper::MGlycan> list_of_glycans = mgl.get_list_of_glycans();

  std::vector<std::string> svg_list; 

  if (list_of_glycans.size() > 0)
  {
    int glycansPermutated = 0;
    clipper::String current_chain = "";
    float z_score_total_for_protein = 0;
    int number_of_indiviual_glycans = 0;

    for (int i = 0; i < list_of_glycans.size(); i++)
    {
      clipper::String wurcs_string;
      if (current_chain != list_of_glycans[i].get_chain())
      {
        current_chain = list_of_glycans[i].get_chain();
        // std::cout << std::endl
        //           << std::endl
        //           << "Chain " << current_chain[0] << std::endl
        //           << "-------" << std::endl;
      }
      // std::cout << std::endl
      //           << list_of_glycans[i].print_linear(true, false, true) << std::endl;

      wurcs_string = list_of_glycans[i].generate_wurcs();
      // std::cout << wurcs_string << std::endl;

      privateer::glycanbuilderplot::Plot plot(false, false, list_of_glycans[i].get_root_by_name());
      plot.plot_glycan(list_of_glycans[i]);
      std::ostringstream os;
      os << list_of_glycans[i].get_root_for_filename() << ".svg";
      svg_list.emplace_back(plot.write_to_string());
    }

    return svg_list; 
  }
}

EMSCRIPTEN_BINDINGS(privateer_module)
{

  function("read_structure", &read_file);

  register_vector<std::string>("vector<string>");


}
