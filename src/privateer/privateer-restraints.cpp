
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
using namespace pybind11::literals;

std::string privateer::restraints::check_monlib_access ( ) {

  std::string pathname("");

  // try complete path first
  char *tmp = std::getenv ( "CLIBD_MON" );

  if (tmp == NULL) {
    // try general CCP4 environment variable
    tmp = std::getenv ( "CCP4" );
    if (tmp == NULL)
      // empty pathname means no access
      return pathname;
    else {
      pathname.assign(tmp);
      pathname.append("/lib/data/monomers");
      return pathname;
    }
  }
  else {
    pathname.assign(tmp);
    return pathname;
  }
}


// Class CarbohydrateDictionary

void privateer::restraints::CarbohydrateDictionary::read_from_file( std::string filename ) {

  this->path_to_cif_file = filename;
  this->cif_document = gemmi::cif::read_file( filename );
  for (gemmi::cif::Block& block : cif_document.blocks)
    if (!block.name.empty() && block.name != "comp_list") {
      this->chemical_component = gemmi::make_chemcomp_from_block(block);
    }
}

void privateer::restraints::CarbohydrateDictionary::read_from_monlib ( std::string ccd_id ) {

  this->from_monlib = true;
  std::string path_to_lib = privateer::restraints::check_monlib_access();
  std::stringstream str;
  std::locale loc;

  char initial = std::tolower(ccd_id[0],loc);

  if (!path_to_lib.empty()) {
    str << path_to_lib << initial << "/" << ccd_id << ".cif";
    std::cout << str.str() << std::endl;
    this->cif_document = gemmi::cif::read_file( str.str() );
    for (gemmi::cif::Block& block : cif_document.blocks)
      if (!block.name.empty() && block.name != "comp_list") {
        this->chemical_component = gemmi::make_chemcomp_from_block(block);
      }
  }
}

void privateer::restraints::CarbohydrateDictionary::write_to_file( std::string filename ) {
  std::ofstream of;
  of.open(filename);

  of << "# " << filename << '\n';
  of << "# modified by Privateer\n";
  gemmi::cif::write_cif_to_stream(of, this->cif_document);
  of.close();
}

void privateer::restraints::CarbohydrateDictionary::restrain_rings_unimodal () {
  for (gemmi::cif::Block& block : cif_document.blocks)
    if (!block.name.empty() && block.name != "comp_list") {
      gemmi::cif::Table chem_comp_tor = block.find("_chem_comp_tor.",
                               {"id", "value_angle", "value_angle_esd", "period"});
      assert(chemical_component.rt.torsions.size() == chem_comp_tor.length());
      for (size_t j = 0; j != chemical_component.rt.torsions.size(); ++j) {
        gemmi::Restraints::Torsion& tor = chemical_component.rt.torsions[j];
        auto ring = chemical_component.rt.find_shortest_path(tor.id4, tor.id1, {tor.id2, tor.id3});
        if (!ring.empty()) {
          auto row = chem_comp_tor[j];
          row[0] = "Privateer_"+tor.label;
          row[2] = "3.0";
          row[3] = "1"; // unimodal
        }
      }
    }
}

void privateer::restraints::CarbohydrateDictionary::restrain_rings_unimodal_from_conformer () {
  for (gemmi::cif::Block& block : cif_document.blocks)
    if (!block.name.empty() && block.name != "comp_list") {
      gemmi::cif::Table chem_comp_tor = block.find("_chem_comp_tor.",
                               {"id", "value_angle", "value_angle_esd", "period"});
      assert(chemical_component.rt.torsions.size() == chem_comp_tor.length());
      for (size_t j = 0; j != chemical_component.rt.torsions.size(); ++j) {
        // need to find a way to access ideal cooords
        gemmi::Restraints::Torsion& tor = chemical_component.rt.torsions[j];
        auto ring = chemical_component.rt.find_shortest_path(tor.id4, tor.id1, {tor.id2, tor.id3});
        if (!ring.empty()) {
          auto row = chem_comp_tor[j];
          row[0] = "Privateer_" + tor.label;
          row[2] = "3.0";
          row[3] = "1"; // unimodal
        }
      }
    }
}



void privateer::restraints::CarbohydrateDictionary::add_inverted_torsions () {
  for (gemmi::cif::Block& block : cif_document.blocks)
    if (!block.name.empty() && block.name != "comp_list") {

      gemmi::cif::Table chem_comp_tor = block.find_or_add("_chem_comp_tor.",{"comp_id", "id", "atom_id_1",
                                                          "atom_id_2", "atom_id_3", "atom_id_4", "value_angle",
                                                          "value_angle_esd", "period"});
      assert(chemical_component.rt.torsions.size() == chem_comp_tor.length());
      for (size_t j = 0; j != chemical_component.rt.torsions.size(); ++j) {
        gemmi::Restraints::Torsion& tor = chemical_component.rt.torsions[j];
        auto ring = chemical_component.rt.find_shortest_path(tor.id4, tor.id1, {tor.id2, tor.id3});
        if (!ring.empty()) {
          chem_comp_tor.append_row({ this->chemical_component.name, "1C4_"+tor.label,
                                     tor.id1.atom, tor.id2.atom, tor.id3.atom, tor.id4.atom, std::to_string(-tor.value),
                                     std::to_string(tor.esd), std::to_string(tor.period) } );
        }
      }
    }
}

void privateer::restraints::CarbohydrateDictionary::print_torsion_restraints () {
  int i = 0;
  for (gemmi::Restraints::Torsion& tor : chemical_component.rt.torsions) {
    i++;
    std::cout << "Torsion " << i << "\t" << tor.id1.atom.c_str() << "\t"
              << tor.id2.atom.c_str() << "\t" << tor.id3.atom.c_str() << "\t" << tor.id4.atom.c_str() << "\t"
              << tor.value << "\t" << tor.esd << "\t" << tor.period << std::endl;
  }
}

// This function is similar to gemmi's but returns a standard python object
pybind11::dict privateer::restraints::CarbohydrateDictionary::get_bond (std::string atom_1, std::string atom_2) {

  for ( auto bond : this->chemical_component.rt.bonds ) {
    if (( bond.id1 == atom_1 ) && (bond.id2 == atom_2)) {
      auto result = pybind11::dict ("length"_a=bond.value, "esd"_a=bond.esd);
      return result;
    }
    else if (( bond.id1 == atom_2 ) && (bond.id2 == atom_1)) {
      auto result = pybind11::dict ("length"_a=bond.value, "esd"_a=bond.esd);
      return result;
    }
  }
  auto result = pybind11::dict ("length"_a="", "esd"_a="");
  return result;
}
// End CarbohydrateDictionary class


// Class CarbohydrateLibrary

void privateer::restraints::CarbohydrateLibrary::read_from_file ( std::string filename ) {

  this->path_to_cif_file = filename;
  this->cif_document = gemmi::cif::read_file( filename );
  for (gemmi::cif::Block& block : cif_document.blocks)
    if (!block.name.empty() && block.name != "comp_list") {
      gemmi::ChemComp chem_comp = gemmi::make_chemcomp_from_block(block);
      this->list_of_chemicals.push_back(privateer::restraints::CarbohydrateDictionary(chem_comp));
    }
}

void privateer::restraints::CarbohydrateLibrary::write_to_file( std::string filename = "" ) {

  std::ofstream of;
  if (filename.empty()) {
    of.open(this->path_to_cif_file); // write to input cif file
  }
  else {
    of.open(filename);
  }

  of << "# " << filename << '\n';
  of << "# modified by Privateer\n";
  gemmi::cif::write_cif_to_stream(of, cif_document);
  of.close();
}

void privateer::restraints::CarbohydrateLibrary::add_to_library () {

}

// End CarbohydrateLibrary class


void privateer::restraints::create_library () {

}


void privateer::restraints::sign_library_header() {

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
