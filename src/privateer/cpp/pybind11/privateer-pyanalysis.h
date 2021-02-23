// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York



#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
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
        GlycanStructure(const clipper::MGlycology& mgl, const int id){
          this->pyinit ( mgl, id );
        };
        ~GlycanStructure() { };
        void pyinit ( const clipper::MGlycology& mgl, const int id);
        void initialize_summary_of_glycan();

        int get_glycan_id( ) { return id; };
        int get_total_number_of_sugars( ) { return numberOfSugars; };
        std::string get_wurcs_notation( ) { return glycanWURCS; };
        pybind11::list get_unique_monosaccharides( ) { return uniqueMonosaccharides; };
        int get_total_of_glycosidic_bonds( ) { return numberOfGlycosidicBonds; };
        std::string get_glycosylation_type( ) { return glycosylationType; };
        pybind11::dict get_root_info( ) { return rootSummary; };
        pybind11::dict get_protein_glycan_linkage_torsions( ) { return protein_glycan_linkage_torsion; };

        pybind11::dict get_glycan_summary( ) { return glycanSummary; };

        CarbohydrateStructure get_monosaccharide(const int id);
        pybind11::list get_all_monosaccharides();
      private:
        clipper::MGlycan glycan;
        std::vector<clipper::MSugar> sugars_in_glycan;
        
        int id;
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
        CarbohydrateStructure(clipper::MGlycan& mglycan, const int id){
          this->pyinit ( mglycan, id );
        };
        ~CarbohydrateStructure() { };
        void pyinit (clipper::MGlycan& mglycan, const int id);
        void initialize_summary_of_sugar();
      private:
        std::string sugar_type;
    };
  }
}

#endif
