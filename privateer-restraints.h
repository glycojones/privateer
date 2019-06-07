
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
    furanose_4EV, furanose_4TO, furanose_EVO, furanose, furanose_1EV,
    furanose_1T2, furanose_EV2
  };

  namespace restraints {

    std::string check_monlib_access ();

    class CarbohydrateDictionary {
      public:
        CarbohydrateDictionary();
        CarbohydrateDictionary(std::string& path_to_cif_file);
        ~CarbohydrateDictionary();
      private:
        gemmi::ChemComp chemical_component;
        std::string path_to_cif_file;
        bool from_monlib; // we don't want to write to mon_lib, right?
    };

    class CarbohydrateLibrary {
    public:
        CarbohydrateLibrary();
        ~CarbohydrateLibrary();
      private:
        std::vector<CarbohydrateDictionary> list_of_entries;
    };


    void create_library ();
    void add_to_library ();
    void sign_library_header();
    void read_dictionary();
    void write_dictionary();
    void read_library ( gemmi::cif::Document &doc,
                        std::vector<gemmi::ChemComp>& list_of_chemicals,
                        std::string filename );
    void write_library ( gemmi::cif::Document &doc, std::string filename );
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
