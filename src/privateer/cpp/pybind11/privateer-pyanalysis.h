// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York



#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <string>
#include <exception>
#include "clipper-glyco.h"
#include "privateer-lib.h"


#ifndef PRIVATEER_PYANALYSIS_H_INCLUDED
#define PRIVATEER_PYANALYSIS_H_INCLUDED

namespace privateer {

  namespace pyanalysis {

    class GlycanStructure;
    class CarbohydrateStructure;

    class GlycosylationComposition 
    {
      public:
        GlycosylationComposition() { };
        GlycosylationComposition(std::string& path_to_model_file, std::string expression_system) {
          this->read_from_file ( path_to_model_file, expression_system );
        };
        ~GlycosylationComposition() { };
        void read_from_file( std::string path_to_model_file, std::string expression_system );
        void initialize_summary_of_detected_glycans( clipper::MGlycology& mglObject );

        std::string get_path_of_model_file_used ( ) { return path_to_model_file; };
        std::string get_expression_system_used ( ) { return expression_system; };
        int get_number_of_glycan_chains_detected ( ) { return numberOfGlycanChains; };
        
        pybind11::list get_summary_of_detected_glycans () { return glycosylationSummary; };
        
        GlycanStructure get_glycan(const int id);
      private:
        clipper::MGlycology mgl;
        std::string path_to_model_file;
        std::string expression_system;
        int numberOfGlycanChains;
        pybind11::list glycosylationSummary;
    };

    class GlycanStructure 
    {
      public:
        GlycanStructure() { };
        GlycanStructure(const clipper::MGlycology& mgl, const int glycanID){
          this->pyinit ( mgl, glycanID );
        };
        ~GlycanStructure() { };
        void pyinit ( const clipper::MGlycology& mgl, const int glycanID);
        void initialize_summary_of_glycan();

        int get_glycan_id( ) const { return glycanID; };
        int get_total_number_of_sugars( ) { return numberOfSugars; };
        std::string get_wurcs_notation( ) { return glycanWURCS; };
        pybind11::list get_unique_monosaccharides( ) { return uniqueMonosaccharides; };
        int get_total_of_glycosidic_bonds( ) { return numberOfGlycosidicBonds; };
        std::string get_glycosylation_type( ) { return glycosylationType; };
        pybind11::dict get_root_info( ) { return rootSummary; };
        pybind11::dict get_protein_glycan_linkage_torsions( ) { return protein_glycan_linkage_torsion; };

        pybind11::dict get_glycan_summary( ) { return glycanSummary; };

        CarbohydrateStructure get_monosaccharide(const int glycanID);
        // pybind11::list get_all_monosaccharides();
      private:
        clipper::MGlycan glycan;
        std::vector<clipper::MSugar> sugars_in_glycan;
        
        int glycanID;
        int numberOfSugars;
        std::string glycanWURCS;
        pybind11::list uniqueMonosaccharides;
        int numberOfGlycosidicBonds;
        std::string glycosylationType;
        pybind11::dict rootSummary;
        pybind11::dict protein_glycan_linkage_torsion;
        pybind11::dict glycanSummary;
    };
  
    class CarbohydrateStructure 
    {
      public:
        CarbohydrateStructure() { };
        CarbohydrateStructure(clipper::MGlycan& mglycan, const int sugarID, const int glycanID){
          this->pyinit ( mglycan, sugarID, glycanID );
        };
        ~CarbohydrateStructure() { };
        void pyinit (clipper::MGlycan& mglycan, const int sugarID, const int glycanID);
        void initialize_summary_of_sugar();

        bool operator==(const CarbohydrateStructure& inputSugar) const { return sugarID == inputSugar.get_sugar_id(); }

        int get_sugar_id( ) const { return sugarID; };
        int get_glycan_id( ) const { return glycanID; };

        std::string get_anomer() { return sugar_anomer; };
        std::string get_type() { return sugar_type; };

        void set_context ( std::string context ) { this->sugar_context = context; }
        
      private:
        clipper::MSugar sugar;
        clipper::MGlycan parentGlycan;

        int sugarID;
        int glycanID;

        std::string sugar_conformation_name;
        std::string sugar_conformation_name_iupac;
        float sugar_puckering_amplitude;
        int sugar_conformation;
        std::string sugar_anomer;
        std::string sugar_handedness;
        std::string sugar_denomination;
        int         sugar_ring_cardinality;
        // const clipper::Coord_orth& sugar_centre;
        // const clipper::Vec3<ftype>& sugar_mean_plane;
        pybind11::list sugar_cremer_pople_params; // std::vector<ftype>
        bool sugar_sane;
        std::string sugar_name;
        std::string sugar_type;
        pybind11::list sugar_ring_angles; // std::vector<ftype>
        pybind11::list sugar_ring_bonds; // std::vector<ftype>
        pybind11::list sugar_ring_torsion; // std::vector<ftype>
        float sugar_ring_bond_rmsd;
        float sugar_ring_angle_rmsd;
        float sugar_bfactor;
        bool sugar_found_db;
        bool sugar_supported;
        bool sugar_diag_ring;
        bool sugar_diag_bonds_rmsd;
        bool sugar_diag_angles_rmsd;
        bool sugar_diag_anomer;
        bool sugar_diag_chirality;
        bool sugar_diag_conformation;
        bool sugar_diag_puckering;
        float sugar_rscc;
        float sugar_accum;
        bool sugar_occupancy_check;
        std::string sugar_context;
    };
  }
}

#endif
