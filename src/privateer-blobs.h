
// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2018 Haroldas Bagdonas & Kevin Cowtan & Jon Agirre
// York Structural Biology Laboratory
// The University of York
// mailto: hb1115@york.ac.uk
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk


#ifndef BLOBS_H_INCLUDED
#define BLOBS_H_INCLUDED

#include <iostream>
#include <fstream>
#include <iterator>
#include <regex>
#include <string>
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include <list>
#include "privateer-lib.h"
#include <clipper/clipper.h>
#include <clipper/clipper-cif.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <clipper/contrib/sfcalc_obs.h>
#include <clipper/minimol/minimol_utils.h>



extern "C"{
#include <stdlib.h>
}


struct GlycosylationMonomerMatch
{
	int PolymerID;
	int ResidueID;
};

struct PotentialGlycosylationSiteInfo
{
	int chainID;
	int monomerID;
	int typeOfGlycosylation;
};

struct GlycanToMiniMolIDs
{
	int proteinMiniMolID;
	int carbohydrateChainMiniMolID;
	int carbohydrateID;
};


std::vector<std::vector<GlycosylationMonomerMatch> > get_matching_monomer_positions(clipper::MiniMol& inputModel);
clipper::MiniMol get_model_without_waters(const clipper::String& ippdb);
bool check_glycosylation_presence(int chainID, int residueID, std::vector < clipper::MGlycan > glycanList);
clipper::Coord_orth getTargetPoint(clipper::Coord_orth& coord1, clipper::Coord_orth& coord2, int vectorShiftDistance);
void fillSearchArea(clipper::MiniMol& inputModel, clipper::Coord_orth& targetPos, clipper::Xmap<float>& sigmaa_dif_map, clipper::HKL_info& hklinfo, clipper::Map_stats& mapstats, int chainID, int monomerID);
void drawOriginPoint(clipper::MiniMol& inputModel, clipper::Coord_orth target, int chainID, int monomerID);
GlycanToMiniMolIDs getCarbohydrateRelationshipToMiniMol(clipper::MiniMol& inputModel, clipper::MSugar& carbohydrate, std::vector < clipper::MGlycan >& allSugars, int mglycanid, int sugaringlycanid);
double calculateMeanElectronDensityInTargetPosition(clipper::Coord_orth targetPos, clipper::Xmap<float>& sigmaa_dif_map, clipper::Map_stats& mapstats);
double calculateMeanElectronDensityForBiggerSphere(clipper::Coord_orth& targetPos, clipper::Xmap<float>& sigmaa_dif_map, clipper::Map_stats& mapstats, clipper::HKL_info& hklinfo);
std::vector<clipper::String> create_list_of_ignored_sugar_atoms(clipper::MSugar& carbohydrate);
std::vector<std::pair<PotentialGlycosylationSiteInfo, double> > get_electron_density_of_potential_glycosylation_sites(const std::vector<std::vector<GlycosylationMonomerMatch>>& informationVector, int vectorIndex, clipper::MiniMol& mmol, clipper::Xmap<float>& sigmaa_dif_map, clipper::HKL_info& hklinfo, std::vector < clipper::MGlycan >& glycanList, clipper::Map_stats& mapstats, float thresholdED, bool pdbexport = false);
std::vector<std::pair<GlycanToMiniMolIDs, double> > get_electron_density_of_potential_unmodelled_carbohydrate_monomers(std::vector < clipper::MSugar > glycanChain, clipper::MiniMol&inputModel, std::vector < clipper::MGlycan >& allSugars, int id, clipper::Xmap<float>& sigmaa_dif_map, clipper::HKL_info& hklinfo, clipper::Map_stats& mapstats, float thresholdED, bool pdbexport = false);



#endif