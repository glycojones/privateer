
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
	int FirstMMonomer;
	int LastMMonomer;
};

struct PotentialGlycosylationSiteInfo
{
	int chainID;
	int monomerID;
	int typeOfGlycosylation;
};

std::vector<std::vector<GlycosylationMonomerMatch> > get_matching_monomer_positions(const clipper::String& ippdb);
std::string get_HTML_output(const clipper::String& title, const clipper::String& ippdb);
clipper::MiniMol get_model_without_waters(const clipper::String& ippdb);
std::vector<std::pair<PotentialGlycosylationSiteInfo, double> > get_electron_density_of_potential_glycosylation_sites(const std::vector<std::vector<GlycosylationMonomerMatch>>& informationVector, int vectorIndex, clipper::MiniMol& mmol, clipper::Xmap<float>& sigmaa_dif_map, clipper::Grid_sampling& grid, clipper::HKL_info& hklinfo);

#endif