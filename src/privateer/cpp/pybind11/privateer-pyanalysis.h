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
#include <algorithm>
#include "clipper-glyco.h"
#include "privateer-lib.h"
#include "privateer-xray.h"
#include "privateer-cryo_em.h"
#include "privateer-dbquery.h"
#include "privateer-interactions.h"


#ifndef PRIVATEER_PYANALYSIS_H_INCLUDED
#define PRIVATEER_PYANALYSIS_H_INCLUDED

using namespace pybind11::literals;

namespace privateer {

  namespace pyanalysis {

    // TODO: add [] operator for these classes.

    // class Privateer/Session/Core/Settings; // for stuff like setting output directory, calling a pipeline, user settings, generating coot files etc etc.
    // Maybe even that's the place to force users to default nCores arguments and stuff like that.
    // Leaning the closest to class Session or class Core;
    class CrystallographicData;
    class GlycosylationInteractions;
    class GlycosylationComposition;
    class GlycanStructure;
    class CarbohydrateStructure;
    class XRayData;
    class CryoEMData;
    class OfflineGlycomicsDatabase;
    class OfflineTorsionsDatabase;

    class CrystallographicData
    {
      public:
        CrystallographicData() { }
        CrystallographicData(std::string path_to_model_file, std::string path_to_mtz_file)
        {
          if(path_to_model_file != "undefined")
            this->parse_model_file(path_to_model_file);
          
          if(path_to_mtz_file != "undefined")
            this->parse_mtz_data_file(path_to_mtz_file);
        }

        pybind11::dict get_model_data() { return model_data; };
        pybind11::dict get_mtz_data() { return mtz_data; };
      
      private:
        pybind11::dict model_data;
        pybind11::dict mtz_data;

        void parse_model_file(std::string& path_to_model_file);
        void parse_mtz_data_file(std::string& path_to_mtz_file);

    };

    class GlycosylationInteractions
    {
      public:
        GlycosylationInteractions() { }
        GlycosylationInteractions(std::string& path_to_model_file) {
          this->input_path = path_to_model_file;
          this->read_from_file(path_to_model_file);
        }
        void read_from_file (std::string& path_to_model_file);
        std::string get_path_of_model_file_used ( ) { return input_path; };
        
        pybind11::list get_all_detected_interactions();
        pybind11::list get_all_detected_hbonds();
        pybind11::list get_all_detected_chpibonds();
        pybind11::dict get_all_interactions_for_specific_glycan(int glycanIndex);
        pybind11::dict get_hbonds_for_specific_glycan(int glycanIndex);
        pybind11::dict get_chpibonds_for_specific_glycan(int glycanIndex);

      private:
        std::string input_path;
        clipper::MiniMol input_model;
				clipper::MAtomNonBond manb_object;
				clipper::MGlycology mglycology;
        privateer::interactions::HBondsParser hbonds;
        privateer::interactions::CHPiBondsParser chpibonds;
    };

    class GlycosylationComposition 
    {
      public:
      // need to add constructor with both model and experimental data too. 
        GlycosylationComposition() { };
        GlycosylationComposition(std::string& path_to_model_file, std::string expression_system, bool debug_output) {
          this->read_from_file ( path_to_model_file, expression_system, debug_output );
        };
        GlycosylationComposition(std::string& path_to_model_file, std::string& path_to_mtz_file, std::string& input_column_fobs_user, int nThreads, float ipradius, std::string expression_system, bool debug_output);
        GlycosylationComposition(std::string& path_to_model_file, std::string& path_to_mrc_file, float resolution, int nThreads, float ipradius, std::string expression_system, bool debug_output);
        ~GlycosylationComposition() { };
        void read_from_file( std::string path_to_model_file, std::string expression_system, bool debug_output );
        void initialize_summary_of_detected_glycans();

        std::string get_path_of_model_file_used ( ) { return path_to_model_file; };
        std::string get_expression_system_used ( ) { return expression_system; };
        int get_number_of_glycan_chains_detected ( ) { return numberOfGlycanChains; };
        
        pybind11::list get_summary_of_detected_glycans () { return glycosylationSummary; };
        
        GlycanStructure get_glycan(const int id);

        pybind11::list get_ligands() { return ligands; };
        
        pybind11::list get_torsions_summary(OfflineTorsionsDatabase& importedDatabase);

        void update_with_experimental_data(privateer::pyanalysis::XRayData& xray_data);
        void update_with_experimental_data(privateer::pyanalysis::CryoEMData& cryoem_data);
        bool check_if_updated_with_experimental_data() { return updatedWithExperimentalData; };

      private:
        bool debug_output;

        clipper::MGlycology mgl;
        std::vector<std::pair<clipper::String, clipper::MSugar>> ligandList;
        std::vector<std::pair<clipper::String, clipper::MSugar>> ligandsOnly;
        std::string path_to_model_file;
        std::string expression_system;
        int numberOfGlycanChains;
        pybind11::list glycosylationSummary;
        pybind11::list glycans;
        pybind11::list ligands;
        pybind11::list torsions;

        bool updatedWithExperimentalData;
    };

    class GlycosylationComposition_memsafe
    {
      public:
      // need to add constructor with both model and experimental data too. 
        GlycosylationComposition_memsafe() { };
        GlycosylationComposition_memsafe(std::string& path_to_model_file, std::string expression_system, bool debug_output) {
          this->read_from_file ( path_to_model_file, expression_system, debug_output );
        };
        ~GlycosylationComposition_memsafe() { };
        void read_from_file( std::string path_to_model_file, std::string expression_system, bool debug_output );
        void initialize_summary_of_detected_glycans( clipper::MGlycology& mglObject );

        std::string get_path_of_model_file_used ( ) { return path_to_model_file; };
        std::string get_expression_system_used ( ) { return expression_system; };
        int get_number_of_glycan_chains_detected ( ) { return numberOfGlycanChains; };
        
        pybind11::list get_summary_of_detected_glycans () { return glycosylationSummary; };
        
        GlycanStructure get_glycan(const int id);

        pybind11::list get_torsions_summary(OfflineTorsionsDatabase& importedDatabase);

      private:
        bool debug_output;

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
          this->updatedWithExperimentalData = false;
          this->pyinit_memsafe ( mgl, glycanID );
        };
        GlycanStructure(const clipper::MGlycology& mgl, const int glycanID, privateer::pyanalysis::GlycosylationComposition& parentGlycosylationComposition){
          this->updatedWithExperimentalData = false;
          this->pyinit ( mgl, glycanID, parentGlycosylationComposition );
        };
        GlycanStructure(const clipper::MGlycology& mgl, const int glycanID, privateer::pyanalysis::GlycosylationComposition& parentGlycosylationComposition, std::vector<std::pair< clipper::String , clipper::MSugar> >& finalLigandList){
          this->updatedWithExperimentalData = true;
          this->pyinitWithExperimentalData ( mgl, glycanID, parentGlycosylationComposition, finalLigandList );
        };
        ~GlycanStructure() { };
        void pyinit ( const clipper::MGlycology& mgl, const int glycanID, privateer::pyanalysis::GlycosylationComposition& parentGlycosylationComposition);
        void pyinit_memsafe ( const clipper::MGlycology& mgl, const int glycanID);
        void pyinitWithExperimentalData (const clipper::MGlycology& mgl, const int glycanID, privateer::pyanalysis::GlycosylationComposition& parentGlycosylationComposition, std::vector<std::pair< clipper::String , clipper::MSugar> >& finalLigandList);
        void initialize_summary_of_glycan();
         
        
        int get_glycan_id( ) const { return glycanID; };
        int get_total_number_of_sugars( ) { return numberOfSugars; };
        std::string get_wurcs_notation( ) { return glycanWURCS; };
        pybind11::list get_unique_monosaccharide_codes( ) { return uniqueMonosaccharides; };
        int get_total_of_glycosidic_bonds( ) { return numberOfGlycosidicBonds; };
        std::string get_glycosylation_type( ) { return glycosylationType; };
        std::string get_root_sugar_chain_id( ) { return chain_root_sugar_ID; };
        pybind11::dict get_root_info( ) { return rootSummary; };
        pybind11::dict get_protein_glycan_linkage_torsions( ) { return protein_glycan_linkage_torsion; };

        pybind11::dict get_glycan_summary( ) { return glycanSummary; };

        CarbohydrateStructure get_monosaccharide(const int glycanID);
        pybind11::list get_all_monosaccharides( ) { return sugars; };

        pybind11::dict query_glycomics_database( OfflineGlycomicsDatabase& importedDatabase, bool returnClosestMatches, bool returnAllPossiblePermutations, int nThreads );
        pybind11::list get_torsions_summary(OfflineTorsionsDatabase& importedDatabase);
        pybind11::dict get_SNFG_strings(bool includeClosestMatches);
        

        // pybind11::list return_permutations_of_glycan(bool returnAllPossiblePermutations, int nThreads) // Could be added under request. Right now don't see much use for it.

        void update_with_experimental_data(privateer::pyanalysis::XRayData& xray_data);
        void update_with_experimental_data(privateer::pyanalysis::CryoEMData& cryoem_data);
        bool check_if_updated_with_experimental_data() { return updatedWithExperimentalData; };
      private:
        privateer::pyanalysis::GlycosylationComposition parentGlycosylation;
        clipper::MGlycan glycan;
        std::vector<clipper::MSugar> sugars_in_glycan;
        
        int glycanID;
        int numberOfSugars;
        std::string glycanWURCS;
        pybind11::list uniqueMonosaccharides;
        int numberOfGlycosidicBonds;
        std::string glycosylationType;
        std::string chain_root_sugar_ID;
        pybind11::dict rootSummary;
        pybind11::dict protein_glycan_linkage_torsion;
        pybind11::dict glycanSummary;
        pybind11::list sugars;

        pybind11::dict glycoproteomicsDB;

        bool updatedWithExperimentalData;
        std::vector<std::pair<std::pair<clipper::MGlycan, std::vector<int>>,float>> outputGlycanPermutationContainer;
        
        void update_summary_of_glycan_after_dbquery()
        {
          auto tempGlycanSummary = glycanSummary;
          this->glycanSummary = pybind11::dict("GlycanID"_a=glycanID, "WURCS"_a=glycanWURCS, "GlycosylationType"_a=glycosylationType, "RootInfo"_a=rootSummary, "ProteinGlycanLinkageTorsion"_a=protein_glycan_linkage_torsion, "ExperimentalData"_a=updatedWithExperimentalData, "Glycoproteomics_DB_query"_a=glycoproteomicsDB);
        };
    };
  
    class CarbohydrateStructure 
    {
      public:
        CarbohydrateStructure() { };
        CarbohydrateStructure(clipper::MGlycan& mglycan, const int sugarID, const int glycanID, privateer::pyanalysis::GlycosylationComposition& parentGlycosylationComposition, privateer::pyanalysis::GlycanStructure& parentGlycanStructure){
          this->updatedWithExperimentalData = false;
          this->pyinit ( mglycan, sugarID, glycanID, parentGlycosylationComposition, parentGlycanStructure );
        };
        CarbohydrateStructure(clipper::MGlycan& mglycan, const int sugarID, const int glycanID, privateer::pyanalysis::GlycosylationComposition& parentGlycosylationComposition, privateer::pyanalysis::GlycanStructure& parentGlycanStructure, std::vector<clipper::MSugar>& list_of_sugars){
          this->updatedWithExperimentalData = true;
          this->pyinitWithExperimentalData ( mglycan, sugarID, glycanID, parentGlycosylationComposition, parentGlycanStructure, list_of_sugars );
        };
        CarbohydrateStructure(const int sugarID, std::vector<std::pair<clipper::String, clipper::MSugar>>& inputSugarList, privateer::pyanalysis::GlycosylationComposition& parentGlycosylationComposition, bool updatedWithExperimentalData){
          this->updatedWithExperimentalData = updatedWithExperimentalData;
          this->pyinitLigand( sugarID, inputSugarList, parentGlycosylationComposition );
        };
        ~CarbohydrateStructure() { };
        
        void pyinit (clipper::MGlycan& mglycan, const int sugarID, const int glycanID, privateer::pyanalysis::GlycosylationComposition& parentGlycosylationComposition, privateer::pyanalysis::GlycanStructure& parentGlycanStructure);
        void pyinitLigand ( const int sugarID, std::vector<std::pair<clipper::String, clipper::MSugar>>& inputSugarList, privateer::pyanalysis::GlycosylationComposition& parentGlycosylationComposition);
        void pyinitWithExperimentalData (clipper::MGlycan& mglycan, const int sugarID, const int glycanID, privateer::pyanalysis::GlycosylationComposition& parentGlycosylationComposition, privateer::pyanalysis::GlycanStructure& parentGlycanStructure, std::vector<clipper::MSugar>& list_of_sugars);

        void initialize_summary_of_sugar();

        bool operator==(const CarbohydrateStructure& inputSugar) const { return (sugar_pdb_id == inputSugar.get_sugar_pdb_id() && sugar_chain_id == inputSugar.get_chain_id()); }

        pybind11::dict get_sugar_summary( ) { return sugarSummary; };
        
        std::string get_chain_id() const { return sugar_chain_id; };
        int get_sugar_id( ) const { return sugarID; };
        int get_glycan_id( ) const { return glycanID; };
        int get_seqnum( ) const { return sugar_seqnum; };
        std::string get_sugar_pdb_id() const { return sugar_pdb_id; };

        std::string get_conformation_name() { return sugar_conformation_name; };
        std::string get_conformation_name_iupac() { return sugar_conformation_name_iupac; };
        float get_puckering_amplitude() { return sugar_puckering_amplitude; };
        std::string get_anomer() { return sugar_anomer; };
        std::string get_handedness() { return sugar_handedness; };
        std::string get_denomination() { return sugar_denomination; };
        int get_ring_cardinality() { return sugar_ring_cardinality; };
        pybind11::list get_cremer_pople_params() { return sugar_cremer_pople_params; };
        bool is_sane() { return sugar_sane; };
        std::string get_privateer_diagnostic() { return privateer_diagnostic; };
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
        float get_sugar_rscc() { return sugar_rscc; };
        float get_sugar_accum() { return sugar_accum; };
        bool get_sugar_occupancy_check() { return sugar_occupancy_check; }; 
        std::string get_glycosylation_context() { return sugar_context; };
        
        pybind11::list get_ring_atom_names() { return sugar_ring_atoms; };
        std::string get_wurcs_residue_code() { return clipper::data::convert_to_wurcs_residue_code(sugar.type().trim()); };
        int get_number_of_connections() { return parentGlycan.get_node(sugarID).number_of_connections(); };
        pybind11::list get_sugar_linkage_info() { return sugar_linkages; };

        void set_sugar_rscc(float input_sugar_rscc) { sugar_rscc = input_sugar_rscc; };
        void set_sugar_accum(float input_sugar_accum) { sugar_accum = input_sugar_accum; };
        void set_sugar_occupancy_check(bool input_sugar_occupancy_check) { sugar_occupancy_check = input_sugar_occupancy_check; };

        bool check_if_updated_with_experimental_data() { return updatedWithExperimentalData; };
      private:
        privateer::pyanalysis::GlycosylationComposition parentGlycosylation;
        privateer::pyanalysis::GlycanStructure parentGlycanStructure;
        clipper::MSugar sugar;
        clipper::MGlycan parentGlycan;
        clipper::MGlycan::Node sugarNode;

        pybind11::dict sugarSummary;

        int sugarID;
        int glycanID;
        int sugar_seqnum;
        std::string sugar_chain_id;
        std::string sugar_pdb_id;

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
        std::string privateer_diagnostic;
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
        float sugar_rscc = -1; // need to develop setter method as in privateer.cpp
        float sugar_accum = -1; // need to develop setter method as in privateer.cpp
        bool sugar_occupancy_check = false; // need to develop setter method as in privateer.cpp
        std::string sugar_context;

        pybind11::list sugar_ring_atoms;

        int sugar_connections;
        pybind11::list sugar_linkages;

        bool updatedWithExperimentalData;
    };

    class XRayData 
    {
      public:
        XRayData() { };
        XRayData(std::string& path_to_mtz_file, std::string& path_to_model_file, std::string& input_column_fobs_user, float ipradius, int nThreads, bool debug_output) {
          this->read_from_file ( path_to_mtz_file, path_to_model_file, input_column_fobs_user, ipradius, nThreads, debug_output);
        };
        ~XRayData() { };
        void read_from_file( std::string& path_to_mtz_file, std::string& path_to_model_file, std::string& input_column_fobs_user, float ipradius, int nThreads, bool debug_output);
        pybind11::list get_sugar_summary_with_experimental_data() { return sugar_summary_of_experimental_data; };
        pybind11::list get_ligand_summary_with_experimental_data() { return ligand_summary_of_experimental_data; };
        void print_cpp_console_output_summary() 
        { 
          int pos_slash = path_to_model_file.rfind("/");
          clipper::String path_to_model_file_clipper = path_to_model_file;
          privateer::util::print_monosaccharide_summary_python (false, false, pos_slash, false, finalLigandList, hklinfo, path_to_model_file_clipper); 
        };

        void write_output_maps();
        

        std::vector<std::pair< clipper::String , clipper::MSugar> > get_finalLigandList() { return finalLigandList; } // only c++
        std::vector<std::pair< clipper::String , clipper::MSugar> > get_finalLigandOnly() { return final_LigandsOnly; } // only c++

      private:
        bool debug_output;
        std::string path_to_model_file;
        clipper::MiniMol mmol;
        clipper::HKL_info hklinfo;
        clipper::String input_column_fobs;
        
        std::vector<std::pair< clipper::String , clipper::MSugar> > finalLigandList;
        std::vector<std::pair<clipper::String, clipper::MSugar>> final_LigandsOnly;
        pybind11::list sugar_summary_of_experimental_data;
        pybind11::list ligand_summary_of_experimental_data;

        // private methods
        pybind11::list generate_sugar_experimental_data_summary(std::vector<std::pair< clipper::String , clipper::MSugar>>& finalLigandList);
    };

    class CryoEMData 
    {
      public:
        CryoEMData() { };
        CryoEMData(std::string& path_to_mrc_file, std::string& path_to_model_file, float resolution, float ipradius, int nThreads, bool debug_output) {
          this->read_from_file ( path_to_mrc_file, path_to_model_file, resolution, ipradius, nThreads, debug_output);
        };
        ~CryoEMData() { };
        void read_from_file( std::string& path_to_mrc_file, std::string& path_to_model_file, float resolution, float ipradius, int nThreads, bool debug_output);
        pybind11::list get_sugar_summary_with_experimental_data() { return sugar_summary_of_experimental_data; };
        pybind11::list get_ligand_summary_with_experimental_data() { return ligand_summary_of_experimental_data; };
        void print_cpp_console_output_summary() 
        { 
          int pos_slash = path_to_model_file.rfind("/");
          clipper::String path_to_model_file_clipper = path_to_model_file;
          privateer::util::print_monosaccharide_summary_python (false, false, pos_slash, true, finalLigandList, hklinfo, path_to_model_file_clipper); 
        };
        

        std::vector<std::pair< clipper::String , clipper::MSugar> > get_finalLigandList() { return finalLigandList; } // only c++
        std::vector<std::pair< clipper::String , clipper::MSugar> > get_finalLigandOnly() { return final_LigandsOnly; } // only c++

      private:
        bool debug_output;
        std::string path_to_model_file;
        float resolution;
        clipper::MiniMol mmol;
        clipper::HKL_info hklinfo;
        
        std::vector<std::pair< clipper::String , clipper::MSugar> > finalLigandList;
        std::vector<std::pair<clipper::String, clipper::MSugar>> final_LigandsOnly;
        pybind11::list sugar_summary_of_experimental_data;
        pybind11::list ligand_summary_of_experimental_data;

        // private methods
        pybind11::list generate_sugar_experimental_data_summary(std::vector<std::pair< clipper::String , clipper::MSugar>>& finalLigandList);
    };

    class OfflineGlycomicsDatabase 
    {
      public:
        OfflineGlycomicsDatabase() { this->path_of_input_file = "nopath"; this->import_json_file(path_of_input_file); };
        OfflineGlycomicsDatabase(std::string& path_to_input_file) {
          this->import_json_file ( path_to_input_file);
        };
        void import_json_file( std::string& path_to_input_file );
        ~OfflineGlycomicsDatabase() { };

        std::vector<privateer::json::GlycomicsDatabase> return_imported_database() { return glytoucanglyconnectdatabase; };

      private:
        std::string path_of_input_file;
        std::vector<privateer::json::GlycomicsDatabase> glytoucanglyconnectdatabase;

    };

    class OfflineTorsionsDatabase 
    {
      public:
        OfflineTorsionsDatabase() { this->path_of_input_file = "nopath"; this->import_json_file(path_of_input_file); };
        OfflineTorsionsDatabase(std::string& path_to_input_file) {
          this->import_json_file ( path_to_input_file);
        };
        void import_json_file( std::string& path_to_input_file );
        ~OfflineTorsionsDatabase() { };

        std::vector<privateer::json::TorsionsDatabase> return_imported_database() { return torsionsdatabase; };

      private:
        std::string path_of_input_file;
        std::vector<privateer::json::TorsionsDatabase> torsionsdatabase;

    };

  }
}

#endif
