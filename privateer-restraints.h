
// Restraint-handling code for Privateer
// (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2019 Jon Agirre
// York Structural Biology Laboratory
// The University of York
// mailto: jon.agirre@york.ac.uk
//

#include <gemmi/chemcomp.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/to_cif.hpp>  // for write_cif_to_stream
#include <string>
#include <locale>
#include "clipper-glyco.h"

namespace privateer {

  enum Conformation {
    pyranose_4C1 = 1, pyranose_1C4, pyranose_3OB, pyranose_B25, pyranose_14B,
    pyranose_B3O, pyranose_25B, pyranose_B14, pyranose_OE, pyranose_E5,
    pyranose_4E, pyranose_E3, pyranose_2E, pyranose_E1, pyranose_3E, pyranose_E2,
    pyranose_1E, pyranose_EO, pyranose_5E, pyranose_E4, pyranose_OH5,
    pyranose_4H5, pyranose_4H3, pyranose_2H3, pyranose_2H1, pyranose_OH1,
    pyranose_3H2, pyranose_1H2, pyranose_1HO, pyranose_5HO, pyranose_5H4,
    pyranose_3H4, pyranose_OS2, pyranose_1S5, pyranose_1S3, pyranose_2SO,
    pyranose_5S1, pyranose_3S1, furanose_3T2, furanose_3EV, furanose_3T4,
    furanose_EV4, furanose_OT4, furanose_OEV, furanose_OT1, furanose_EV1,
    furanose_2T1, furanose_2EV, furanose_2T3, furanose_EV3, furanose_4T3,
    furanose_4EV, furanose_4TO, furanose_EVO, furanose_1TO, furanose_1EV,
    furanose_1T2, furanose_EV2
  };

  namespace restraints {

    std::string check_monlib_access ();

    class Ring {
      public:
        Ring() { };
        Ring(std::vector<gemmi::Restraints::AtomId> list_of_atoms) {
          this->list_of_atoms = list_of_atoms;
        }
        ~Ring() { };
        std::vector<gemmi::Restraints::AtomId> get_list_of_atoms () {
          return list_of_atoms;
        }
        void set_list_of_atoms ( std::vector<gemmi::Restraints::AtomId> list_of_atoms ) {
          this->list_of_atoms = list_of_atoms;
        }
      private:
        std::vector<gemmi::Restraints::AtomId> list_of_atoms;
    };

    class CarbohydrateDictionary {
      public:
        CarbohydrateDictionary() { };
        CarbohydrateDictionary(std::string& path_to_cif_file) {
          this->read_from_file ( path_to_cif_file );
        };
        CarbohydrateDictionary(gemmi::ChemComp chem_comp) {
          this->chemical_component = chem_comp;
          this->path_to_cif_file = "";
          this->from_monlib = false;
        };
        ~CarbohydrateDictionary() { };
        std::string get_chemcomp_id () {
          return this->chemical_component.name;
        }
        void read_from_file( std::string filename );
        void read_from_monlib ( std::string ccd_id );
        void write_to_file( std::string filename );
        void restrain_rings_unimodal ();
        void add_inverted_torsions ();
        void print_torsion_restraints ();

      private:
        gemmi::ChemComp chemical_component;
        gemmi::cif::Document cif_document;
        std::string path_to_cif_file;
        bool from_monlib; // we don't want to write to mon_lib, right?
        std::vector<Ring> list_of_rings;
    };

    class CarbohydrateLibrary {
    public:
        CarbohydrateLibrary() { };
        CarbohydrateLibrary( std::string filename ) {
          this->read_from_file(filename);
        };
        ~CarbohydrateLibrary() { };
        void read_from_file   ( std::string filename );
        void write_to_file    ( std::string filename );
        int number_of_entries ( ) {
          return list_of_chemicals.size();
        };
        void add_dictionary ( privateer::restraints::CarbohydrateDictionary dict) {
          list_of_chemicals.push_back(dict);
        }

      private:
        std::vector<CarbohydrateDictionary> list_of_chemicals;
        gemmi::cif::Document cif_document;
        std::string path_to_cif_file = "";

        void add_to_library ();
    };


    void create_library ();
    void sign_library_header();

    void add_torsion_set (float phi);
    void add_torsion_set (float phi, float theta);
    void add_torsion_set ( gemmi::ChemComp &cc, privateer::Conformation id);
    void restrain_conformation (privateer::Conformation);
    privateer::Conformation get_conformation ( clipper::MMonomer sugar );
    void replace_conformer ();
    void read_conformer(clipper::MMonomer &sugar);
    void calculate_conformer ();
    void refine_conformer ();

  }
}
