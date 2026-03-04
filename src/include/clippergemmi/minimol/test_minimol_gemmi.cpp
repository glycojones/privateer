#include "test_minimol_gemmi.h"
#include <gemmi/model.hpp>
#include <gemmi/symmetry.hpp>
//#include "../gemmi/clipper_gemmi_model.h"

namespace clipper{
  bool Test_minimol_gemmi::test(const String &id, const Spacegroup &sg1, const gemmi::SpaceGroup &sg2)
  {
    int pass = 0;
    //Spacegroup st_spg = GEMMI::spacegroup(*st.find_spacegroup());
    if (sg1.symbol_hall().compare("Unknown") != 0){
      if (sg1.symbol_hall().compare(sg2.hall) != 0) {
        std::cout << id << " (symbol_hall) - Result : ";
        std::cout << sg1.symbol_hall() << " != " << sg2.hall << std::endl;
        pass++;
        error_count++;
      }
    } else { std::cout << "Unknown spacegroup!\n"; }
    return (pass == 0);
  }

  bool Test_minimol_gemmi::test(const String &id, const MiniMol &mmol,
                                const GemmiModel &gmol) {
    int err = 0;
    for (int c = 0; c < mmol.size(); c++) {
      if (!mmol[c].id_match(mmol[c].id(), String(gmol.chains[c].name),
                            MM::UNIQUE)) {
        std::cout << "Chain : " << mmol[c].id();
        std::cout << " != " << gmol.chains[c].name << std::endl;
        err++;
      }
      for (int r = 0; r < mmol[c].size(); r++) {
        if (mmol[c][r].type().compare(gmol.chains[c].residues[r].name) != 0) {
          std::cout << "Res : " << mmol[c][r].type() << " ; ";
          std::cout << gmol.chains[c].residues[r].name << std::endl;
          err++;
        }
        if (mmol[c][r].seqnum() != gmol.chains[c].residues[r].seqid.num.value) {
          std::cout << "SeqNum : " << mmol[c][r].seqnum() << " ; ";
          std::cout << gmol.chains[c].residues[r].seqid.num.value << std::endl;
          err++;
        }
        for (int a = 0; a < mmol[c][r].size(); a++) {
          auto gatm = gmol.chains[c].residues[r].atoms[a];
          if (mmol[c][r][a].name().trim().compare(gatm.name) != 0) {
            std::cout << "Atom : " << mmol[c][r][a].name().trim() << " ; ";
            std::cout << gatm.name << std::endl;
            err++;
          }
          if (mmol[c][r][a].exists_property("AltConf")){
            String altconf = dynamic_cast<const Property<String> &>(mmol[c][r][a].get_property("AltConf")).value();
            if(!gatm.altloc_matches(altconf[0])){
              std::cout << " altloc : " << altconf << " ; " << gatm.altloc << std::endl;
              err++;
            }
          }
          if (mmol[c][r][a].occupancy() != gatm.occ) {
            std::cout << " occupancy : " << mmol[c][r][a].occupancy() << " ; " << gatm.occ << std::endl;
              err++;
          }
	  if (mmol[c][r][a].element().compare(gatm.element.uname()) != 0){
            std::cout << " element : " << mmol[c][r][a].element() << " ; " << gatm.element.uname() << std::endl;
              err++;
          }
          if (mmol[c][r][a].u_iso() != Util::b2u(gatm.b_iso)) {
            std::cout << " u_iso : " << mmol[c][r][a].u_iso() << " ; " << Util::b2u(gatm.b_iso) << std::endl;
              err++;
          }
	  if (!mmol[c][r][a].coord_orth().equals(Vec3<ftype>(gatm.pos.x,gatm.pos.y,gatm.pos.z), 1e-6)){
            std::cout << " coord_orth : " << mmol[c][r][a].coord_orth().format() << " ; ";
            std::cout << gatm.pos.x << ", " << gatm.pos.y << ", " << gatm.pos.z << std::endl;
            err++;
          } 
        }
      }
    }
    error_count += err;
    return (err == 0);
  }

  bool Test_minimol_gemmi::run(const String &filename)
  {
    gemmi::Structure st;
    GEMMIfile gfile;
    MiniMol mmol;
    try { // ground truth read structure using Gemmi function
      st = gemmi::read_structure_file(filename, gemmi::coor_format_from_ext(filename));
    } catch (const std::exception &e) {
      String mess = "Test_minimol_gemmi: run - \n ";
      mess += e.what();
      mess += " \n";
      Message::message(Message_fatal(mess));
    }
    // cast gemmi::Structure to Clipper class
    GemmiStructure gst(st);
    GemmiAtom_list gst_atmlist(gst.models[0].all());
    gfile.read_file(filename);
    gfile.import_minimol(mmol);
    Atom_list catmlist = mmol.atom_list();
    std::cout << "Coordinate file read : " << filename << std::endl;
    std::cout << "Comparing contents from GEMMIfile::read_file with those read "
                 "directly from Gemmi: \n";
    test("Spacegroup", gfile.spacegroup(), *st.find_spacegroup());
    test("MiniMol", mmol, gst.models[0]);
 
    return ( error_count == 0);
  }
}
