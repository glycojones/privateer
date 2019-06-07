
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

void privateer::restraints::CarbohydrateDictionary::write_dictionary( std::string filename ) {

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
