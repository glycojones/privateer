// Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
//    -> The name combines Mark Knopfler's 'Privateering' cracking album & Dr Cowtan's tradition of using maritime acronyms
//
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York



#include "privateer-json.h"
#include "clipper-glyco.h"
namespace privateer 
{
    namespace db
    {
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

        class OfflineTorsionsZScoreDatabase 
        {
        public:
            OfflineTorsionsZScoreDatabase() { this->path_of_input_file = "nopath"; this->import_json_file(path_of_input_file); };
            OfflineTorsionsZScoreDatabase(std::string& path_to_input_file) {
            this->import_json_file ( path_to_input_file);
            };
            void import_json_file( std::string& path_to_input_file );
            ~OfflineTorsionsZScoreDatabase() { };

            std::vector<privateer::json::TorsionsZScoreDatabase> return_imported_database() { return torsions_zscore_database; };
            void compute_z_score_for_glycan(std::vector<clipper::MGlycan::MGlycanTorsionSummary>& torsion_summary_of_glycan);

        private:
            std::string path_of_input_file;
            std::vector<privateer::json::TorsionsZScoreDatabase> torsions_zscore_database;
        };
    }
}