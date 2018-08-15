
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
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-minimol.h>
#include <iterator>
#include <regex>
#include <string>
#include <vector>
#include <utility>

extern "C"{
#include <stdlib.h>
}


struct GlycosylationMonomerMatch
{
	int PolymerID;
	int FirstMMonomer;
	int LastMMonomer;
};

std::vector<std::vector<GlycosylationMonomerMatch> > get_matching_mmonomer_positons(const clipper::String& ippdb);
std::string get_HTML_output(const clipper::String& title, const clipper::String& ippdb, clipper::String& HTMLopseq);
clipper::MiniMol get_model_without_waters(const clipper::String& ippdb);

#endif
