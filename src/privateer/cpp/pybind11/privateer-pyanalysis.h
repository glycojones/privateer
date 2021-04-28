// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York


#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <string>
#include <exception>
#include "clipper-glyco.h"
#include "privateer-lib.h"
#include "privateer-xray.h"


#ifndef PRIVATEER_PYANALYSIS_H_INCLUDED
#define PRIVATEER_PYANALYSIS_H_INCLUDED

namespace privateer {

    // class GilManager
    // {
    //     public:
    //         GilManager()
    //         {
    //             mThreadState = PyEval_SaveThread();
    //         }

    //         ~GilManager()
    //         {
    //             if (mThreadState)
    //             PyEval_RestoreThread(mThreadState);
    //         }

    //         GilManager(const GilManager&) = delete;
    //         GilManager& operator=(const GilManager&) = delete;
    //     private:
    //         PyThreadState* mThreadState;
    // };


  namespace pyanalysis {

    class GlycanStructure;
    class CarbohydrateStructure;
    class XRayData;
    class CryoEMData;

    class GlycosylationComposition 
    {
      public:
        GlycosylationComposition() { };
        GlycosylationComposition(std::string& path_to_model_file, std::string expression_system) {
          this->read_from_file ( path_to_model_file, expression_system );
        };
        // GlycosylationComposition(std::string& path_to_model_file, std::string& path_to_model_file); // figure out how to differentiate between the over overloaded operator. 
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
        pybind11::list get_unique_monosaccharide_codes( ) { return uniqueMonosaccharides; };
        int get_total_of_glycosidic_bonds( ) { return numberOfGlycosidicBonds; };
        std::string get_glycosylation_type( ) { return glycosylationType; };
        pybind11::dict get_root_info( ) { return rootSummary; };
        pybind11::dict get_protein_glycan_linkage_torsions( ) { return protein_glycan_linkage_torsion; };

        pybind11::dict get_glycan_summary( ) { return glycanSummary; };

        CarbohydrateStructure get_monosaccharide(const int glycanID);
        pybind11::list get_all_monosaccharides( ) { return sugars; };
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
        pybind11::list sugars;
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

        bool operator==(const CarbohydrateStructure& inputSugar) const { return (sugar_pdb_id == inputSugar.get_sugar_pdb_id() && sugar_pdb_chain == inputSugar.get_sugar_pdb_chain()); }

        pybind11::dict get_sugar_summary( ) { return sugarSummary; };
        
        int get_sugar_id( ) const { return sugarID; };
        int get_glycan_id( ) const { return glycanID; };
        int get_sugar_pdb_id() const { return sugar_pdb_id; };
        std::string get_sugar_pdb_chain() const { return sugar_pdb_chain; };

        std::string get_conformation_name() { return sugar_conformation_name; };
        std::string get_conformation_name_iupac() { return sugar_conformation_name_iupac; };
        float get_puckering_amplitude() { return sugar_puckering_amplitude; };
        std::string get_anomer() { return sugar_anomer; };
        std::string get_handedness() { return sugar_handedness; };
        std::string get_denomination() { return sugar_denomination; };
        int get_ring_cardinality() { return sugar_ring_cardinality; };
        pybind11::list get_cremer_pople_params() { return sugar_cremer_pople_params; };
        bool is_sane() { return sugar_sane; };
        std::string get_name_full() { return sugar_name_full; };
        std::string get_name_short() { return sugar_name_short; };
        std::string get_type() { return sugar_type; };
        pybind11::list get_ring_angles() { return sugar_ring_angles; };
        pybind11::list get_ring_bonds() { return sugar_ring_bonds; };
        pybind11::list get_ring_torsion() { return sugar_ring_torsion; };
        float get_ring_bond_rmsd() { return sugar_ring_bond_rmsd; };
        float get_ring_angle_rmsd() { return sugar_ring_angle_rmsd; };
        float get_bfactor() { return sugar_bfactor; };
        bool is_supported() { return sugar_supported; };
        bool ok_with_ring() { return sugar_diag_ring; };
        bool ok_with_bonds_rmsd() { return sugar_diag_bonds_rmsd; };
        bool ok_with_angles_rmsd() { return sugar_diag_angles_rmsd; };
        bool ok_with_anomer() { return sugar_diag_anomer; };
        bool ok_with_chirality() { return sugar_diag_chirality; };
        bool ok_with_conformation() { return sugar_diag_conformation; };
        bool ok_with_puckering() { return sugar_diag_puckering; };
        std::string get_glycosylation_context() { return sugar_context; };
      
      private:
        clipper::MSugar sugar;
        clipper::MGlycan parentGlycan;

        pybind11::dict sugarSummary;

        int sugarID;
        int glycanID;
        int sugar_pdb_id;
        std::string sugar_pdb_chain;

        int sugar_conformation_code;
        std::string sugar_conformation_name;
        std::string sugar_conformation_name_iupac;
        float sugar_puckering_amplitude;
        std::string sugar_anomer;
        std::string sugar_handedness;
        std::string sugar_denomination;
        int         sugar_ring_cardinality;
        // const clipper::Coord_orth& sugar_centre; // leave this for future implementation.
        // const clipper::Vec3<ftype>& sugar_mean_plane; // leave this for future implementation.
        pybind11::list sugar_cremer_pople_params; // std::vector<ftype>
        bool sugar_sane;
        std::string sugar_name_full;
        std::string sugar_name_short;
        std::string sugar_type;
        pybind11::list sugar_ring_angles; // std::vector<ftype>
        pybind11::list sugar_ring_bonds; // std::vector<ftype>
        pybind11::list sugar_ring_torsion; // std::vector<ftype>
        float sugar_ring_bond_rmsd;
        float sugar_ring_angle_rmsd;
        float sugar_bfactor;
        bool sugar_supported;
        bool sugar_diag_ring;
        bool sugar_diag_bonds_rmsd;
        bool sugar_diag_angles_rmsd;
        bool sugar_diag_anomer;
        bool sugar_diag_chirality;
        bool sugar_diag_conformation;
        bool sugar_diag_puckering;
        float sugar_rscc; // need to develop setter method as in privateer.cpp
        float sugar_accum; // need to develop setter method as in privateer.cpp
        bool sugar_occupancy_check; // need to develop setter method as in privateer.cpp
        std::string sugar_context;
    };

    class XRayData 
    {
      public:
        XRayData() { };
        XRayData(std::string& path_to_mtz_file, std::string& path_to_model_file, std::string& input_column_fobs_user) {
          this->read_from_file ( path_to_mtz_file, path_to_model_file, input_column_fobs_user );
        };
        ~XRayData() { };
        void read_from_file( std::string& path_to_mtz_file, std::string& path_to_model_file, std::string& input_column_fobs_user);
        
      private:
        clipper::MiniMol mmol;
        clipper::HKL_info hklinfo;
        clipper::CCP4MTZfile mtzin;
        clipper::String input_column_fobs; // need to convert user input std::string to clipper::string for internal functions not visible to user. 
        clipper::Xmap<float> sigmaa_omit_fd;
        clipper::Xmap<float> ligandmap; // equivalent to lignadmap in xray implementation
        clipper::Xmap<float> mask;
        clipper::Grid_sampling mygrid;
        clipper::Coord_orth origin;
        clipper::Coord_orth destination;
    };
    
  }
}

#endif
