// Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
//    -> The name combines Mark Knopfler's 'Privateering' cracking album & Dr Cowtan's tradition of using maritime acronyms
//
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#include "privateer-db.h"

///////////////////////////////////////////////// Class OfflineGlycomicsDatabase ////////////////////////////////////////////////////////////////////
void privateer::db::OfflineGlycomicsDatabase::import_json_file(std::string &path_to_input_file)
{
    std::vector<privateer::json::GlycomicsDatabase> glycomics_database = privateer::json::read_json_file_for_glycomics_database(path_to_input_file);
    this->glytoucanglyconnectdatabase = glycomics_database;
}
///////////////////////////////////////////////// Class OfflineGlycomicsDatabase END ////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////// Class OfflineTorsionsDatabase ////////////////////////////////////////////////////////////////////
void privateer::db::OfflineTorsionsDatabase::import_json_file(std::string &path_to_input_file)
{
    std::vector<privateer::json::TorsionsDatabase> torsions_database = privateer::json::read_json_file_for_torsions_database(path_to_input_file);
    this->torsionsdatabase = torsions_database;
}

///////////////////////////////////////////////// Class OfflineTorsionsDatabase END ////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////// Class OfflineTorsionsZScoreDatabase ////////////////////////////////////////////////////////////////////
void privateer::db::OfflineTorsionsZScoreDatabase::import_json_file(std::string &path_to_input_file)
{
    privateer::json::GlobalTorsionZScore torsions_database = privateer::json::read_json_file_for_torsions_zscore_database(path_to_input_file);
    this->torsions_zscore_database = torsions_database;
}

void privateer::db::OfflineTorsionsZScoreDatabase::compute_z_score_for_glycan(std::vector<clipper::MGlycan::MGlycanTorsionSummary> &torsion_summary_of_glycan)
{
    std::cout << "Length of glycan torsions " << torsion_summary_of_glycan.size() << std::endl;
    std::cout << "Length of torsionsZScoreDatabase " << this->torsions_zscore_database.database_array.size() << std::endl;
}

///////////////////////////////////////////////// Class OfflineTorsionsZScoreDatabase END ///////////////////////////////////////////////////////////////////