// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

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


float privateer::restraints::CarbohydrateDictionary::get_torsion_from_conformer (std::string a1, std::string a2, std::string a3, std::string a4) {
  std::vector<size_t> atom_indices;

  for (gemmi::cif::Block& block : cif_document.blocks) {
    if (!block.name.empty() && block.name != "comp_list") {

      gemmi::cif::Table chem_comp_atom = block.find("_chem_comp_atom.",
                                                    {"comp_id", "atom_id", "x", "y", "z"});

      for (size_t j = 0; j != chem_comp_atom.length(); j++ ) {
        auto atom_in_conformer = chem_comp_atom[j];
        if ( atom_in_conformer[1] == a1 )
          atom_indices.push_back(j);
      }

      for (size_t j = 0; j != chem_comp_atom.length(); j++ ) {
        auto atom_in_conformer = chem_comp_atom[j];
        if ( atom_in_conformer[1] == a2 )
          atom_indices.push_back(j);
      }

      for (size_t j = 0; j != chem_comp_atom.length(); j++ ) {
        auto atom_in_conformer = chem_comp_atom[j];
        if ( atom_in_conformer[1] == a3 )
          atom_indices.push_back(j);
      }

      for (size_t j = 0; j != chem_comp_atom.length(); j++ ) {
        auto atom_in_conformer = chem_comp_atom[j];
        if ( atom_in_conformer[1] == a4 )
          atom_indices.push_back(j);
      }

      if ( atom_indices.size() != 4 )
        return 361.0; // Return out-of-range value to indicate that atoms were not found
      else { // We've got our four atoms identified
        auto atom1_row = chem_comp_atom[atom_indices[0]];
        clipper::Coord_orth atom1_coords ( std::stof(atom1_row[2]),
                                           std::stof(atom1_row[3]),
                                           std::stof(atom1_row[4]) );

        auto atom2_row = chem_comp_atom[atom_indices[1]];
        clipper::Coord_orth atom2_coords ( std::stof(atom2_row[2]),
                                           std::stof(atom2_row[3]),
                                           std::stof(atom2_row[4]) );

        auto atom3_row = chem_comp_atom[atom_indices[2]];
        clipper::Coord_orth atom3_coords ( std::stof(atom3_row[2]),
                                           std::stof(atom3_row[3]),
                                           std::stof(atom3_row[4]) );

        auto atom4_row = chem_comp_atom[atom_indices[3]];
        clipper::Coord_orth atom4_coords ( std::stof(atom4_row[2]),
                                           std::stof(atom4_row[3]),
                                           std::stof(atom4_row[4]) );

        float torsion_value = clipper::Coord_orth::torsion ( atom1_coords,
                                                             atom2_coords,
                                                             atom3_coords,
                                                             atom4_coords );

        return clipper::Util::rad2d(torsion_value);
      }
    }
  }
  return 361.0; // Error: CIF block not found
}


void privateer::restraints::CarbohydrateDictionary::restrain_rings_unimodal_from_conformer () {
  int restraint_index = 1;
  for (gemmi::cif::Block& block : cif_document.blocks)
    if (!block.name.empty() && block.name != "comp_list") {
      gemmi::cif::Table chem_comp_tor  = block.find("_chem_comp_tor.",
                               {"comp_id", "id", "atom_id_1", "atom_id_2", "atom_id_3", "atom_id_4", "value_angle", "value_angle_esd", "period"});

      assert(chemical_component.rt.torsions.size() == chem_comp_tor.length());
      for (int j = chem_comp_tor.length() -1; j >=0 ; j--) {

        gemmi::Restraints::Torsion tor = chemical_component.rt.torsions[j];
        auto ring = chemical_component.rt.find_shortest_path(tor.id4, tor.id1, {tor.id2, tor.id3});
        if (!ring.empty()) {
          chem_comp_tor.remove_row(j); // we want to produce our own torsional set
        }
        else { // regular, non-ring torsion restraint
          auto row = chem_comp_tor[j];
          row[1] = "Privateer_tors_" + std::to_string (restraint_index++);
          float torsion_value = this->get_torsion_from_conformer(row[2], row[3], row[4], row[5]);
          if ( torsion_value != 361.0 )
          {
            std::cout << "Non-ring torsion ";
            for (int i = 2; i < 6; i++)
              std::cout << row[i] << " ";
            std::cout << "of value "<< row[6] << " to be patched with " << torsion_value << " " << std::endl;
            row[6] = std::to_string(torsion_value);
          }
          else {
            std::cout << "Warning: the atoms named in " << row[1] << " do not match the coordinates!" << std::endl;
          }
        }
      }
      // End of the input torsions, need to generate new ring torsions now
      restraint_index = 1;
      auto row = chem_comp_tor[0]; // to get basic info from existing torsions
      gemmi::Restraints::AtomId a_1, a_2, a_3, a_4;
      a_1.comp = a_2.comp = a_3.comp = a_4.comp = 1;

      a_1.atom = "C5";
      a_2.atom = "O5";
      a_3.atom = "C1";
      a_4.atom = "C2";

      if (chemical_component.rt.find_bond (a_2, a_3) != chemical_component.rt.bonds.end()) // general aldopyranose case
      {
        double torsion_value = this->get_torsion_from_conformer(a_1.atom, a_2.atom, a_3.atom, a_4.atom);
        chem_comp_tor.append_row({row[0],
                                  "Privateer_ring_" + std::to_string(restraint_index++),
                                  a_1.atom,
                                  a_2.atom,
                                  a_3.atom,
                                  a_4.atom,
                                  std::to_string(torsion_value),
                                  "3.0",
                                  "1"}); // unimodal

        a_1.atom = "O5";
        a_2.atom = "C1";
        a_3.atom = "C2";
        a_4.atom = "C3";
        auto ring_1 = chemical_component.rt.find_shortest_path(a_1, a_4, {a_2, a_3});
        //if (!ring_1.empty())
          chem_comp_tor.append_row({row[0],
                                    "Privateer_ring_" + std::to_string(restraint_index++),
                                    a_1.atom,
                                    a_2.atom,
                                    a_3.atom,
                                    a_4.atom,
                                    std::to_string(this->get_torsion_from_conformer(a_1.atom, a_2.atom, a_3.atom, a_4.atom)),
                                    "3.0",
                                    "1"}); // unimodal
        a_1.atom = "C1";
        a_2.atom = "C2";
        a_3.atom = "C3";
        a_4.atom = "C4";
        auto ring_2 = chemical_component.rt.find_shortest_path(a_1, a_4, {a_2, a_3});
        //if (!ring_2.empty())
          chem_comp_tor.append_row({row[0],
                                    "Privateer_ring_" + std::to_string(restraint_index++),
                                    a_1.atom,
                                    a_2.atom,
                                    a_3.atom,
                                    a_4.atom,
                                    std::to_string(this->get_torsion_from_conformer(a_1.atom, a_2.atom, a_3.atom, a_4.atom)),
                                    "3.0",
                                    "1"}); // unimodal
        a_1.atom = "C2";
        a_2.atom = "C3";
        a_3.atom = "C4";
        a_4.atom = "C5";
        auto ring_3 = chemical_component.rt.find_shortest_path(a_1, a_4, {a_2, a_3});
        //if (!ring_3.empty())
          chem_comp_tor.append_row({row[0],
                                    "Privateer_ring_" + std::to_string(restraint_index++),
                                    a_1.atom,
                                    a_2.atom,
                                    a_3.atom,
                                    a_4.atom,
                                    std::to_string(this->get_torsion_from_conformer(a_1.atom, a_2.atom, a_3.atom, a_4.atom)),
                                    "3.0",
                                    "1"}); // unimodal

        a_1.atom = "C3";
        a_2.atom = "C4";
        a_3.atom = "C5";
        a_4.atom = "O5";
        auto ring_4 = chemical_component.rt.find_shortest_path(a_1, a_4, {a_2, a_3});
        //if (!ring_4.empty())
          chem_comp_tor.append_row({row[0],
                                    "Privateer_ring_" + std::to_string(restraint_index++),
                                    a_1.atom,
                                    a_2.atom,
                                    a_3.atom,
                                    a_4.atom,
                                    std::to_string(this->get_torsion_from_conformer(a_1.atom, a_2.atom, a_3.atom, a_4.atom)),
                                    "3.0",
                                    "1"}); // unimodal

        a_1.atom = "C4";
        a_2.atom = "C5";
        a_3.atom = "O5";
        a_4.atom = "C1";
        auto ring_5 = chemical_component.rt.find_shortest_path(a_1, a_4, {a_2, a_3});
        //if (!ring_5.empty())
          chem_comp_tor.append_row({row[0],
                                    "Privateer_ring_" + std::to_string(restraint_index++),
                                    a_1.atom,
                                    a_2.atom,
                                    a_3.atom,
                                    a_4.atom,
                                    std::to_string(this->get_torsion_from_conformer(a_1.atom, a_2.atom, a_3.atom, a_4.atom)),
                                    "3.0",
                                    "1"}); // unimodal

      }
      else {
        a_1.atom = "C6";
        a_2.atom = "O6";
        a_3.atom = "C2";
        a_4.atom = "C3";

        if (chemical_component.rt.find_bond (a_2, a_3) != chemical_component.rt.bonds.end()) // general ketopyranose case
        {
          double torsion_value = this->get_torsion_from_conformer(a_1.atom, a_2.atom, a_3.atom, a_4.atom);
          chem_comp_tor.append_row({row[0],
                                    "Privateer_ring_" + std::to_string(restraint_index++),
                                    a_1.atom,
                                    a_2.atom,
                                    a_3.atom,
                                    a_4.atom,
                                    std::to_string(torsion_value),
                                    "3.0",
                                    "1"}); // unimodal

          a_1.atom = "O6";
          a_2.atom = "C2";
          a_3.atom = "C3";
          a_4.atom = "C4";
          auto ring_1 = chemical_component.rt.find_shortest_path(a_1, a_4, {a_2, a_3});
          //if (!ring_1.empty())
            chem_comp_tor.append_row({row[0],
                                      "Privateer_ring_" + std::to_string(restraint_index++),
                                      a_1.atom,
                                      a_2.atom,
                                      a_3.atom,
                                      a_4.atom,
                                      std::to_string(this->get_torsion_from_conformer(a_1.atom, a_2.atom, a_3.atom, a_4.atom)),
                                      "3.0",
                                      "1"}); // unimodal
          a_1.atom = "C2";
          a_2.atom = "C3";
          a_3.atom = "C4";
          a_4.atom = "C5";
          auto ring_2 = chemical_component.rt.find_shortest_path(a_1, a_4, {a_2, a_3});
          //if (!ring_2.empty())
            chem_comp_tor.append_row({row[0],
                                      "Privateer_ring_" + std::to_string(restraint_index++),
                                      a_1.atom,
                                      a_2.atom,
                                      a_3.atom,
                                      a_4.atom,
                                      std::to_string(this->get_torsion_from_conformer(a_1.atom, a_2.atom, a_3.atom, a_4.atom)),
                                      "3.0",
                                      "1"}); // unimodal
          a_1.atom = "C3";
          a_2.atom = "C4";
          a_3.atom = "C5";
          a_4.atom = "C6";
          auto ring_3 = chemical_component.rt.find_shortest_path(a_1, a_4, {a_2, a_3});
          if (!ring_3.empty())
            chem_comp_tor.append_row({row[0],
                                      "Privateer_ring_" + std::to_string(restraint_index++),
                                      a_1.atom,
                                      a_2.atom,
                                      a_3.atom,
                                      a_4.atom,
                                      std::to_string(this->get_torsion_from_conformer(a_1.atom, a_2.atom, a_3.atom, a_4.atom)),
                                      "3.0",
                                      "1"}); // unimodal

          a_1.atom = "C4";
          a_2.atom = "C5";
          a_3.atom = "C6";
          a_4.atom = "O6";
          auto ring_4 = chemical_component.rt.find_shortest_path(a_1, a_4, {a_2, a_3});
          //if (!ring_4.empty())
            chem_comp_tor.append_row({row[0],
                                      "Privateer_ring_" + std::to_string(restraint_index++),
                                      a_1.atom,
                                      a_2.atom,
                                      a_3.atom,
                                      a_4.atom,
                                      std::to_string(this->get_torsion_from_conformer(a_1.atom, a_2.atom, a_3.atom, a_4.atom)),
                                      "3.0",
                                      "1"}); // unimodal

          a_1.atom = "C5";
          a_2.atom = "C6";
          a_3.atom = "O6";
          a_4.atom = "C2";
          auto ring_5 = chemical_component.rt.find_shortest_path(a_1, a_4, {a_2, a_3});
          //if (!ring_5.empty())
            chem_comp_tor.append_row({row[0],
                                      "Privateer_ring_" + std::to_string(restraint_index++),
                                      a_1.atom,
                                      a_2.atom,
                                      a_3.atom,
                                      a_4.atom,
                                      std::to_string(this->get_torsion_from_conformer(a_1.atom, a_2.atom, a_3.atom, a_4.atom)),
                                      "3.0",
                                      "1"}); // unimodal

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
