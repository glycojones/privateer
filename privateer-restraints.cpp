
// Restraint-handling code for Privateer
// (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2019 Jon Agirre
// York Structural Biology Laboratory
// The University of York
// mailto: jon.agirre@york.ac.uk
//

#include "privateer-restraints.h"


bool privateer::restraints::check_monlib_access () {

  // try complete path first
  std::string mlib_path = std::getenv ( "CLIBD_MON" );
  std::cout << "this is mlib path " << mlib_path << std::endl;
  if (mlib_path.length() == 0) {
    // try general CCP4 environment variable
    std::string ccp4_path = std::getenv ( "CCP4" );
    if (ccp4_path.length() == 0)
      return false;
    else
      return true;
  }
  else
    return true;
}

void privateer::restraints::create_library () {

}

void privateer::restraints::add_to_library () {

}

void privateer::restraints::sign_library_header() {

}

void privateer::restraints::read_dictionary() {

}

void privateer::restraints::write_dictionary() {

}

void privateer::restraints::read_library ( gemmi::cif::Document &doc,
                                           std::vector<gemmi::ChemComp>& list_of_chemicals,
                                           std::string filename ) {

  doc = gemmi::cif::read_file( filename );
  for (gemmi::cif::Block& block : doc.blocks)
    if (!block.name.empty() && block.name != "comp_list") {
      gemmi::ChemComp chem_comp = gemmi::make_chemcomp_from_block(block);
      list_of_chemicals.push_back(chem_comp);
    }
}

void privateer::restraints::write_library( gemmi::cif::Document &doc,
                    std::string filename ) {

  std::ofstream of(filename);
  of << "# " << filename << '\n';
  of << "# modified by Privateer\n";
  gemmi::cif::write_cif_to_stream(of, doc);
  of.close();
}

void privateer::restraints::add_torsion_set (float phi) {

}

void privateer::restraints::add_torsion_set (float phi, float theta) {

}

void privateer::restraints::add_torsion_set ( gemmi::ChemComp &chem_comp,
                                              privateer::Conformation id) {
// TODO: everything
  for (gemmi::Restraints::Torsion& tor : chem_comp.rt.torsions) {
    std::printf("[%s] torsion %3s - %3s - %3s - %3s  %f +/- %f\n",
                chem_comp.name.c_str(),
                tor.id1.atom.c_str(), tor.id2.atom.c_str(),
                tor.id3.atom.c_str(), tor.id4.atom.c_str(),
                tor.value, tor.esd);
    tor.value += 3.5;
    tor.esd = 0.3;
  }
}

void privateer::restraints::restrain_conformation (privateer::Conformation) {

}

privateer::Conformation privateer::restraints::get_conformation ( clipper::MMonomer sugar ) {

  return privateer::pyranose_4C1;
}

void privateer::restraints::replace_conformer () {

}

void privateer::restraints::read_conformer( clipper::MMonomer &sugar ) {

}

void privateer::restraints::calculate_conformer () {

}

void privateer::restraints::refine_conformer () {

}
