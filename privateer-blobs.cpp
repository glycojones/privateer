
// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2018 Haroldas Bagdonas & Kevin Cowtan & Jon Agirre
// York Structural Bioexpy Laboratory
// The University of York
// mailto: hb1115@york.ac.uk
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk



#include "privateer-blobs.h"

bool bestPointFinder(std::pair<clipper::Coord_orth, double> p1, std::pair<clipper::Coord_orth, double> p2) {
    return p1.second<p2.second;
}

std::vector<std::vector<GlycosylationMonomerMatch> > get_matching_monomer_positions(const clipper::String& ippdb)
{
	const int mmdbflags = ::mmdb::MMDBF_IgnoreBlankLines | ::mmdb::MMDBF_IgnoreDuplSeqNum | ::mmdb::MMDBF_IgnoreNonCoorPDBErrors | ::mmdb::MMDBF_IgnoreRemarks;
	clipper::MMDBfile mmdbwrk;
	clipper::MiniMol molwrk;
	mmdbwrk.SetFlag(mmdbflags);
	mmdbwrk.read_file(ippdb);
	mmdbwrk.import_minimol(molwrk);


	std::vector<GlycosylationMonomerMatch> NMMonomers;
	std::vector<GlycosylationMonomerMatch> CMMonomers;
	std::vector<GlycosylationMonomerMatch> OMMonomers;
	std::vector<GlycosylationMonomerMatch> SMMonomers;
	std::vector<GlycosylationMonomerMatch> NRemMMonomers;
	std::vector<std::vector<GlycosylationMonomerMatch> > MatchingMMonomers;


	std::regex isProteinNGlc("[N][A-Z^P][ST]|[N][A-Z][C]"); // N-Glycosylation = Asn-X-/Ser/Thr(X=anything except Proline) || Asn-X-Cys(X=Anything)
	std::regex isProteinCGlc("[W][A-Z][A-Z][W]|[W][ST][A-Z][C]"); // C-Glycosylation = Trp-X-X-Trp || Trp-Ser/Thr-X-Cys
	std::regex isProteinOGlc("[S]|[T]"); // temporary fix to loop through all Serines and Threonines. In Long term will be redesigned
	std::regex isProteinSGlc("[C]|[M]"); // temporary fix to loop through all Cysteines and Methionines. In Long term will be redesigned
	std::regex isProteinNGlcRemoved("[A][A-Z^P][ST]|[A][A-Z][C]|[Q][A-Z^P][ST]|[Q][A-Z][C]");
	std::smatch NMatch;
	std::smatch CMatch;
	std::smatch OMatch;
	std::smatch SMatch;
	std::smatch NRemMatch;

const char rtype1[21] =
  {  'A',  'R',  'N',  'D',  'C',  'Q',  'E',  'G',  'H',  'I',
	 'L',  'K',  'M',  'F',  'P',  'S',  'T',  'W',  'Y',  'V',
	 'M'};
const char rtype3[21][4] =
  {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
   "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
   "MSE"};

   const int ntype = sizeof( rtype1 ) / sizeof( rtype1[0] );
   clipper::String seqfile = "";


   for ( int c = 0; c < molwrk.size(); c++ ) {
	 clipper::String id = molwrk[c].id();
	 clipper::String seq;
	 for ( int r = 0; r < molwrk[c].size(); r++ ) {
	   char symbol = ' ';
	   for ( int t = 0; t < ntype; t++ )
	if ( molwrk[c][r].type() == rtype3[t] )
	{
	symbol = rtype1[t];
	}
	   if ( symbol != ' ' )
	   {
	   seq = seq + symbol;
	   }
	 }
	 if ( seq.length() > 0 )
	 {
		 std::string TEMPseq = seq;

		 auto NGlc_begin = std::sregex_iterator(seq.begin(), seq.end(), isProteinNGlc); // Establish sregex_iterator search range for isProteinNGlc
		 auto NGlc_end = std::sregex_iterator();
		 auto CGlc_begin = std::sregex_iterator(seq.begin(), seq.end(), isProteinCGlc); // same as above, just for C-glycosylation
		 auto CGlc_end = std::sregex_iterator();
		 auto OGlc_begin = std::sregex_iterator(seq.begin(), seq.end(), isProteinOGlc); // same as above, just for O-glycosylation
		 auto OGlc_end = std::sregex_iterator();
		 auto SGlc_begin = std::sregex_iterator(seq.begin(), seq.end(), isProteinSGlc); // same as above, just for S-glycosylation
		 auto SGlc_end = std::sregex_iterator();
		 auto NRemGlc_begin = std::sregex_iterator(seq.begin(), seq.end(), isProteinNGlcRemoved); // same as above, just for N-Removed-glycosylation
		 auto NRemGlc_end = std::sregex_iterator();

		if (std::regex_search(seq, isProteinNGlc))
		{
			int it = 0;
			for (std::sregex_iterator monomer = NGlc_begin; monomer != NGlc_end; ++monomer) // a loop to iterate through sregex_iterator to make sure that all matches of N-Glycosylation are found
			{
					NMatch = *monomer;
					int it2 = NMatch.length();
					for(; it < TEMPseq.length(); ++it)
					{
						if (it == NMatch.position())
						{
							NMMonomers.push_back({c, it, it + it2});
						}
						if (it == NMatch.position() + NMatch.length())
						{
							break;
						}
					}
				}
		}
		if (std::regex_search(seq, isProteinCGlc))
		{
			int it = 0;
			for (std::sregex_iterator monomer = CGlc_begin; monomer != CGlc_end; ++monomer) // a loop to iterate through sregex_iterator to make sure that all matches of N-Glycosylation are found
				{
				   CMatch = *monomer;
				   int it2 = CMatch.length();
				   for(; it < TEMPseq.length(); ++it)
					{
						if (it == CMatch.position())
						{
							CMMonomers.push_back({c, it, it + it2});
						}
						if (it == CMatch.position() + CMatch.length())
						{
							break;
						}
					}
				}
		}
		if (std::regex_search(seq, isProteinOGlc))
		{
			int it = 0;
			for (std::sregex_iterator monomer = OGlc_begin; monomer != OGlc_end; ++monomer) // a loop to iterate through sregex_iterator to make sure that all matches of N-Glycosylation are found
			{
				OMatch = *monomer;
				int it2 = OMatch.length();
				for(; it < TEMPseq.length(); ++it)
				{
					if (it == OMatch.position())
					{
						OMMonomers.push_back({c, it, it + it2});
					}
					if (it == OMatch.position() + OMatch.length())
					{
						break;
					}
				}
			}
		}
		if (std::regex_search(seq, isProteinSGlc))
		{
			int it = 0;
			for (std::sregex_iterator monomer = SGlc_begin; monomer != SGlc_end; ++monomer) // a loop to iterate through sregex_iterator to make sure that all matches of N-Glycosylation are found
			{
				SMatch = *monomer;
				int it2 = SMatch.length();
				for(; it < TEMPseq.length(); ++it)
				{
					if (it == SMatch.position())
					{
						SMMonomers.push_back({c, it, it + it2});
					}
					if (it == SMatch.position() + SMatch.length())
					{
						break;
					}
				}
			}
		}
		if (std::regex_search(seq, isProteinNGlcRemoved))
		{
			int it = 0;
			for (std::sregex_iterator monomer = NRemGlc_begin; monomer != NRemGlc_end; ++monomer) // a loop to iterate through sregex_iterator to make sure that all matches of N-Glycosylation are found
			{
				NRemMatch = *monomer;
				int it2 = NRemMatch.length();
				for(; it < TEMPseq.length(); ++it)
				{
					if (it == NRemMatch.position())
					{
						NRemMMonomers.push_back({c, it, it + it2});
					}
					if (it == NRemMatch.position() + NRemMatch.length())
					{
						break;
					}
				}
			}
		}
	}
  }
MatchingMMonomers.push_back({NMMonomers});
MatchingMMonomers.push_back({CMMonomers});
MatchingMMonomers.push_back({OMMonomers});
MatchingMMonomers.push_back({SMMonomers});
MatchingMMonomers.push_back({NRemMMonomers});


return MatchingMMonomers;
}



std::string get_HTML_output(const clipper::String& title, const clipper::String& ippdb)
{

		const int mmdbflags = ::mmdb::MMDBF_IgnoreBlankLines | ::mmdb::MMDBF_IgnoreDuplSeqNum | ::mmdb::MMDBF_IgnoreNonCoorPDBErrors | ::mmdb::MMDBF_IgnoreRemarks;
		clipper::MMDBfile mmdbwrk;
		clipper::MiniMol molwrk;
		mmdbwrk.SetFlag(mmdbflags);
		mmdbwrk.read_file(ippdb);
		mmdbwrk.import_minimol(molwrk);

		std::regex isProteinNGlc("[N][A-Z^P][ST]|[N][A-Z][C]"); // N-Glycosylation = Asn-X-/Ser/Thr(X=anything except Proline) || Asn-X-Cys(X=Anything)
		std::regex isProteinNGlcRemoved("[A][A-Z^P][ST]|[A][A-Z][C]|[Q][A-Z^P][ST]|[Q][A-Z][C]");
		std::regex isProteinCGlc("[W][A-Z][A-Z][W]|[W][ST][A-Z][C]"); // C-Glycosylation = Trp-X-X-Trp || Trp-Ser/Thr/X-Cys
		std::regex isProteinOGlc("[N][A-Z][T]|[N][A-Z][S]|[S][A-Z][A-Z][P]|[P][A-Z][T]|[T][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][P]|[T][A-Z][A-Z][P]");

		std::smatch NMatch;
		std::smatch CMatch;
		std::smatch OMatch;
		std::smatch NRemMatch;


		std::string HTMLhighlight = "<span style=\"background-color: #ffff00; color: #ff0000;\"><strong>";
		std::string HTMLhighlightEnd = "</strong></span>";


	const char rtype1[21] =
	  {  'A',  'R',  'N',  'D',  'C',  'Q',  'E',  'G',  'H',  'I',
		 'L',  'K',  'M',  'F',  'P',  'S',  'T',  'W',  'Y',  'V',
		 'M'};
	const char rtype3[21][4] =
	  {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
	   "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
	   "MSE"};

	   const int ntype = sizeof( rtype1 ) / sizeof( rtype1[0] );
	   clipper::String seqfile = "";
	   clipper::String HTMLfile;


	   for ( int c = 0; c < molwrk.size(); c++ ) {
		 clipper::String id = molwrk[c].id();
		 clipper::String seq;
		 for ( int r = 0; r < molwrk[c].size(); r++ ) {
		   char symbol = ' ';
		   for ( int t = 0; t < ntype; t++ )
		if ( molwrk[c][r].type() == rtype3[t] )
		{
		symbol = rtype1[t];
		}
		   if ( symbol != ' ' )
		   {
		   seq = seq + symbol;
		   }
		 }
		 if ( seq.length() > 0 )
		 {
			 std::string HTMLseq;
			 std::string TEMPseq = seq;
			 int NGlcounter = 0; //Reset the N-Glycosylation counter at the beginning of the loop
			 int CGlcounter = 0; //Reset the C-Glycosylation counter at the beginning of the loop
			 int OGlcounter = 0;
			 int NGlRemovedCounter = 0;
			 auto NGlc_begin = std::sregex_iterator(seq.begin(), seq.end(), isProteinNGlc); // Establish sregex_iterator search range for isProteinNGlc
			 auto NGlc_end = std::sregex_iterator();
			 auto CGlc_begin = std::sregex_iterator(seq.begin(), seq.end(), isProteinCGlc); // same as above, just for C-glycosylation
			 auto CGlc_end = std::sregex_iterator();
			 auto OGlc_begin = std::sregex_iterator(seq.begin(), seq.end(), isProteinOGlc); // same as above, just for C-glycosylation
			 auto OGlc_end = std::sregex_iterator();
			 auto NRemGlc_begin = std::sregex_iterator(seq.begin(), seq.end(), isProteinNGlcRemoved); // same as above, just for C-glycosylation
			 auto NRemGlc_end = std::sregex_iterator();


			if ((std::regex_search(seq, isProteinNGlc) == false) &&
				(std::regex_search(seq, isProteinCGlc) == false) &&
				(std::regex_search(seq, isProteinOGlc) == false) &&
				(std::regex_search(seq, isProteinNGlcRemoved) == false))
				{
				HTMLseq = HTMLseq + "<br" + "/>" + "No possible glycosylation sequence motiffs were detected." + "</p>";
				}
			else
			{
			if (std::regex_search(seq, isProteinNGlc))
			{
				int it = 0;
				for (std::sregex_iterator monomer = NGlc_begin; monomer != NGlc_end; ++monomer) // a loop to iterate through sregex_iterator to make sure that all matches of N-Glycosylation are found
				{
						NMatch = *monomer;
						int it2 = NMatch.length();
						++NGlcounter;
						for(; it < TEMPseq.length(); ++it)
						{
							if (it == NMatch.position())
							{
								HTMLseq.append(HTMLhighlight);
							}
							if (it == NMatch.position() + NMatch.length())
							{
								HTMLseq.append(HTMLhighlightEnd);
								break;
							}
						HTMLseq.push_back(seq[it]);
						}
					}
					for (; it < TEMPseq.length(); ++it)
					{
					HTMLseq.push_back(seq[it]);
					}
			HTMLseq = HTMLseq + "<br" + "/>" + "Detected N-linked glycosylation (s): " + std::to_string(NGlcounter) + "</p>";
			}
			if (std::regex_search(seq, isProteinCGlc))
			{
				int it = 0;
				for (std::sregex_iterator monomer = CGlc_begin; monomer != CGlc_end; ++monomer) // a loop to iterate through sregex_iterator to make sure that all matches of N-Glycosylation are found
					{
					   CMatch = *monomer;
					   int it2 = CMatch.length();
					   ++CGlcounter;
					   for(; it < TEMPseq.length(); ++it)
						{
							if (it == CMatch.position())
							{
								HTMLseq.append(HTMLhighlight);
							}
							if (it == CMatch.position() + CMatch.length())
							{
								HTMLseq.append(HTMLhighlightEnd);
								break;
							}
						HTMLseq.push_back(seq[it]);
						}
					}
					for (; it < TEMPseq.length(); ++it)
					{
					HTMLseq.push_back(seq[it]);
					}
			HTMLseq = HTMLseq + "<br" + "/>" + "Detected C-linked glycosylation (s): " + std::to_string(CGlcounter) + "</p>";
			}
			if (std::regex_search(seq, isProteinOGlc))
			{
				int it = 0;
				for (std::sregex_iterator monomer = OGlc_begin; monomer != OGlc_end; ++monomer) // a loop to iterate through sregex_iterator to make sure that all matches of N-Glycosylation are found
				{
					OMatch = *monomer;
					int it2 = OMatch.length();
					++OGlcounter;
					for(; it < TEMPseq.length(); ++it)
					{
						if (it == OMatch.position())
						{
							HTMLseq.append(HTMLhighlight);
						}
						if (it == OMatch.position() + OMatch.length())
						{
							HTMLseq.append(HTMLhighlightEnd);
							break;
						}
					HTMLseq.push_back(seq[it]);
					}
				}
				for (; it < TEMPseq.length(); ++it)
				{
				HTMLseq.push_back(seq[it]);
				}
			HTMLseq = HTMLseq + "<br" + "/>" + "Detected O-linked glycosylation (s): " + std::to_string(OGlcounter) + "</p>";
			}
			if (std::regex_search(seq, isProteinNGlcRemoved))
			{
				int it = 0;
				for (std::sregex_iterator monomer = NRemGlc_begin; monomer != NRemGlc_end; ++monomer) // a loop to iterate through sregex_iterator to make sure that all matches of N-Glycosylation are found
				{
					NRemMatch = *monomer;
					int it2 = NRemMatch.length();
					++NGlRemovedCounter;
					for(; it < TEMPseq.length(); ++it)
					{
						if (it == NRemMatch.position())
						{
							HTMLseq.append(HTMLhighlight);
						}
						if (it == NRemMatch.position() + NRemMatch.length())
						{
							HTMLseq.append(HTMLhighlightEnd);
							break;
						}
					HTMLseq.push_back(seq[it]);
					}
				}
				for (; it < TEMPseq.length(); ++it)
				{
				HTMLseq.push_back(seq[it]);
				}
			HTMLseq = HTMLseq + "<br" + "/>" + "Detected possible N-Glycosylation motiffs that could have been modified(mutation of ASN to ALA or GLN) to remove N-Glycosylation: "
			+ std::to_string(NGlRemovedCounter) + "</p>";
			}
		}
		HTMLfile = HTMLfile + "<p>&gt;" + title + "_" + id + "<br" + "/>" + HTMLseq + "</p>";
	   }


	}
	std::string HTMLoutput(HTMLfile.c_str());
	return HTMLoutput;
}


clipper::MiniMol get_model_without_waters(const clipper::String& ippdb)
{
	const int mmdbflags = ::mmdb::MMDBF_IgnoreBlankLines | ::mmdb::MMDBF_IgnoreDuplSeqNum | ::mmdb::MMDBF_IgnoreNonCoorPDBErrors | ::mmdb::MMDBF_IgnoreRemarks;
	clipper::MMDBfile mmdbwrk;
	clipper::MiniMol molwrk;
	mmdbwrk.SetFlag(mmdbflags);
	mmdbwrk.read_file(ippdb);
	mmdbwrk.import_minimol(molwrk);


	clipper::MiniMol molwrk_new( molwrk.spacegroup(), molwrk.cell() );
	for ( int c = 0; c < molwrk.size(); c++ )
	{
		clipper::MPolymer mp;
		for ( int r = 0; r < molwrk[c].size(); r++ )
		{
				if ( molwrk[c][r].type() != "HOH" )
				{
					mp.insert( molwrk[c][r] );
				}
				if ( r == (molwrk[c].size() - 1))
				{
					clipper::String chainID = molwrk[c].id();
					mp.set_id(chainID);
					molwrk_new.insert(mp);
					mp = clipper::MPolymer();
				}
		}
	 }

	return molwrk_new;
}

bool check_glycosylation_presence(clipper::String chainID, clipper::String residueID, std::vector < clipper::MGlycan > glycanList)
{
	for ( int i = 0 ; i < glycanList.size() ; i++ )
        {
			clipper::String glycanProteinResidueAttachedTo = glycanList[i].get_root().first.id().trim();
			clipper::String glycanProteinChainAttachedTo = glycanList[i].get_chain().substr(0,1);
		
			if((chainID == glycanProteinChainAttachedTo) && (residueID == glycanProteinResidueAttachedTo))
				return true;
		}
	return false;	
}

clipper::Coord_orth getTargetPoint(clipper::Coord_orth& coord1, clipper::Coord_orth& coord2, int vectorShiftDistance)
{
	clipper::Coord_orth coord; 

	clipper::Vec3<clipper::ftype> baseVector((coord1.x()-coord2.x()),(coord1.y()-coord2.y()), (coord1.z()-coord2.z()));
	// Create a 1A unit vector out of baseVector, to be used later in vector shifting
	clipper::Vec3<clipper::ftype> unitVector = baseVector.unit();

	// Obtain coordinates in the middle of suspected glycan density via 5A vector shift. This is the nearest glycan bonded via ND2 atom to ASN residue.
	coord = clipper::Coord_orth( (coord2.x()+(unitVector[0]*vectorShiftDistance)), (coord2.y()+(unitVector[1]*vectorShiftDistance)), (coord2.z()+(unitVector[2]*vectorShiftDistance)) );

	return coord;
}


void fillSearchArea(clipper::MiniMol& inputModel, clipper::Coord_orth& targetPos, clipper::Xmap<float>& sigmaa_dif_map, clipper::Grid_sampling& grid, clipper::HKL_info& hklinfo, int chainID, int monomerID)
{
	clipper::Coord_orth origin(targetPos.x()-1, targetPos.y()-1, targetPos.z()-1);
	clipper::Coord_orth destination(targetPos.x()+1, targetPos.y()+1, targetPos.z()+1);

	clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;

	i0 = clipper::Xmap_base::Map_reference_coord( sigmaa_dif_map, origin.coord_frac(hklinfo.cell()).coord_grid(grid) );

	for ( iu = i0; iu.coord().u() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).u(); iu.next_u() )
		for ( iv = iu; iv.coord().v() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).v(); iv.next_v() )
			for ( iw = iv; iw.coord().w() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).w(); iw.next_w() )
				{
						clipper::Coord_orth targetuvw = iw.coord_orth();
						clipper::Atom dummyAtom;
						dummyAtom.set_coord_orth(targetuvw);
						dummyAtom.set_element("H");
						clipper::MAtom dummyAtomExport(dummyAtom);
						inputModel[chainID][monomerID].insert(dummyAtomExport);
				}
}

//Need to improve this function to improve calculation of electron density values within the sphere. 
// refer to std::vector<clipper::Xmap_base::Map_reference_coord> in privateer.cpp
double calculateMeanElectronDensityInArea(clipper::Coord_orth& targetPos, clipper::Xmap<float>& sigmaa_dif_map, clipper::Grid_sampling& grid, clipper::HKL_info& hklinfo)
{
	double meanElectronDensity = 0.0;
	int n_points = 0;

	// Define origin and destination for drawing the sphere. Electron density data obtained from within the sphere later on.
	clipper::Coord_orth origin(targetPos.x()-1, targetPos.y()-1, targetPos.z()-1);
	clipper::Coord_orth destination(targetPos.x()+1, targetPos.y()+1, targetPos.z()+1);

	clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;

	if(!origin.is_null() && !destination.is_null())
	{
		i0 = clipper::Xmap_base::Map_reference_coord( sigmaa_dif_map, origin.coord_frac(hklinfo.cell()).coord_grid(grid) );

		for ( iu = i0; iu.coord().u() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).u(); iu.next_u() )
			for ( iv = iu; iv.coord().v() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).v(); iv.next_v() )
				for ( iw = iv; iw.coord().w() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).w(); iw.next_w() )
					{
						meanElectronDensity = meanElectronDensity + sigmaa_dif_map[iw];
						n_points++;
					}
	}
 
	meanElectronDensity = meanElectronDensity / n_points;

	return meanElectronDensity;
}


// Draw 2 vectors, one from CB to CA, another from CB to ND2. Return the max afterwards. 

std::vector<std::pair<PotentialGlycosylationSiteInfo, double> > get_electron_density_of_potential_glycosylation_sites(const std::vector<std::vector<GlycosylationMonomerMatch>>& informationVector, int vectorIndex, clipper::MiniMol& inputModel, clipper::Xmap<float>& sigmaa_dif_map, clipper::Grid_sampling& grid, clipper::HKL_info& hklinfo, std::vector < clipper::MGlycan >& glycanList, bool pdbexport) 
{

	std::vector<std::pair<PotentialGlycosylationSiteInfo, double> > finalVectorForBlobValues;
		if(vectorIndex == 0)
		{
			if (!informationVector[vectorIndex].empty())
			{
				std::list<clipper::String> ignoreAtomList = {""};
				int vectorShiftLimit = 5;
				for (int c = 0; c < informationVector[vectorIndex].size(); c++)
				{
					for(int r = informationVector[vectorIndex][c].FirstMMonomer; r < informationVector[vectorIndex][c].LastMMonomer; r++)
					{
						if (inputModel[informationVector[vectorIndex][c].PolymerID][r].type() == "ASN") // in N-Glycosylation a glycan is attached through ASN residue
                            {	
								bool siteAlreadyGlycosylated = check_glycosylation_presence(inputModel[informationVector[vectorIndex][c].PolymerID].id(), inputModel[informationVector[vectorIndex][c].PolymerID][r].id().trim(), glycanList);
								if(!siteAlreadyGlycosylated)
									{
										clipper::MAtom CBAtom;
										CBAtom = inputModel[informationVector[vectorIndex][c].PolymerID][r].find(" CB ", clipper::MM::ANY);

										clipper::Coord_orth CBCoordinate; // CB atom is used as a direction towards the glycan density
										CBCoordinate = CBAtom.coord_orth();

										std::vector<std::pair<clipper::Coord_orth, double>> pairs;

										for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][r].size(); natom++)
										{
											bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id()) != ignoreAtomList.end());
											if(!atomIgnored)
												{
													clipper::Coord_orth vectorOrigin;

													vectorOrigin = inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].coord_orth();
													
													for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
														{
															clipper::Coord_orth potentialTarget = getTargetPoint(vectorOrigin, CBCoordinate, vectorShift);

															double meanDensityValue = calculateMeanElectronDensityInArea(potentialTarget, sigmaa_dif_map, grid, hklinfo);
															std::pair<clipper::Coord_orth, double> tempDensityInfo(potentialTarget, meanDensityValue);
															pairs.push_back(tempDensityInfo);
														}
												}
										}
										const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
										clipper::Coord_orth bestTarget = bestPair->first;
										double bestDensityValue = bestPair->second;

										if(bestDensityValue > 0.08)
										{
										std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, r, vectorIndex}, bestDensityValue);
										finalVectorForBlobValues.push_back(densityInfo);
										fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, grid, hklinfo, informationVector[vectorIndex][c].PolymerID, r);
										}
									}
							}
					}
				}
			}
		}
		if(vectorIndex == 1)
		{
			if (!informationVector[vectorIndex].empty())
			{
				std::list<clipper::String> ignoreAtomList = {" O  ", " N  ", " C  "};
				int vectorShiftLimit = 10;
				for (int c = 0; c < informationVector[vectorIndex].size(); c++)
				{
					for(int r = informationVector[vectorIndex][c].FirstMMonomer; r < informationVector[vectorIndex][c].LastMMonomer; r++)
					{
						if (inputModel[informationVector[vectorIndex][c].PolymerID][r].type() == "TRP") // in C-Glycosylation a glycan is attached through TRP residue
                            {	
								bool siteAlreadyGlycosylated = check_glycosylation_presence(inputModel[informationVector[vectorIndex][c].PolymerID].id(), inputModel[informationVector[vectorIndex][c].PolymerID][r].id().trim(), glycanList);
								if(!siteAlreadyGlycosylated)
									{                                
										clipper::MAtom CD1Atom; // A glycan is attached to TRP residue via CD1 atom
										CD1Atom = inputModel[informationVector[vectorIndex][c].PolymerID][r].find(" CD1", clipper::MM::ANY);

										clipper::Coord_orth CD1Coordinate; // CD1 atom is used as a direction towards the glycan density
										CD1Coordinate = CD1Atom.coord_orth();

										std::vector<std::pair<clipper::Coord_orth, double>> pairs;

										for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][r].size(); natom++)
										{
											bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id()) != ignoreAtomList.end());
											if(!atomIgnored)
												{
													clipper::Coord_orth vectorOrigin;

													vectorOrigin = inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].coord_orth();
													
													for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
														{
															clipper::Coord_orth potentialTarget = getTargetPoint(CD1Coordinate, vectorOrigin, vectorShift);

															double meanDensityValue = calculateMeanElectronDensityInArea(potentialTarget, sigmaa_dif_map, grid, hklinfo);
															std::pair<clipper::Coord_orth, double> tempDensityInfo(potentialTarget, meanDensityValue);
															pairs.push_back(tempDensityInfo);
														}
												}
										}
										
										const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
										clipper::Coord_orth bestTarget = bestPair->first;
										double bestDensityValue = bestPair->second;

										if(bestDensityValue > 0.08)
										{
										std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, r, vectorIndex}, bestDensityValue);
										finalVectorForBlobValues.push_back(densityInfo);
										fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, grid, hklinfo, informationVector[vectorIndex][c].PolymerID, r);
										}
									}
							}
					}
				}
			}
		}
		if(vectorIndex == 2)
		{
			if (!informationVector[vectorIndex].empty())
			{
				for (int c = 0; c < informationVector[vectorIndex].size(); c++)
				{
					for(int r = informationVector[vectorIndex][c].FirstMMonomer; r < informationVector[vectorIndex][c].LastMMonomer; r++)
					{
						if (inputModel[informationVector[vectorIndex][c].PolymerID][r].type() == "SER") // in O-Glycosylation a glycan is attached through Thr or Ser residue
                            {	
								int vectorShiftLimit = 5; 
								bool siteAlreadyGlycosylated = check_glycosylation_presence(inputModel[informationVector[vectorIndex][c].PolymerID].id(), inputModel[informationVector[vectorIndex][c].PolymerID][r].id().trim(), glycanList);
								if(!siteAlreadyGlycosylated)
									{      								
										std::list<clipper::String> ignoreAtomList = {" O  "};
										clipper::MAtom OGAtom; // On Ser glycan attaches to OG, on Thr glycan attaches to OG1

										OGAtom = inputModel[informationVector[vectorIndex][c].PolymerID][r].find(" OG ", clipper::MM::ANY);

										clipper::Coord_orth OGCoordinate; // OG atom is used as a point of attachement of a glycan
										OGCoordinate = OGAtom.coord_orth();

										std::vector<std::pair<clipper::Coord_orth, double>> pairs;

										for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][r].size(); natom++)
										{
											bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id()) != ignoreAtomList.end());
											if(!atomIgnored)
												{
													clipper::Coord_orth vectorOrigin;

													vectorOrigin = inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].coord_orth();
													
													for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
														{
															clipper::Coord_orth potentialTarget = getTargetPoint(OGCoordinate, vectorOrigin, vectorShift);

															double meanDensityValue = calculateMeanElectronDensityInArea(potentialTarget, sigmaa_dif_map, grid, hklinfo);
															std::pair<clipper::Coord_orth, double> tempDensityInfo(potentialTarget, meanDensityValue);
															pairs.push_back(tempDensityInfo);
														}
												}
										}
										
										const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
										clipper::Coord_orth bestTarget = bestPair->first;
										double bestDensityValue = bestPair->second;

										if(bestDensityValue > 0.09)
										{
										std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, r, vectorIndex}, bestDensityValue);
										finalVectorForBlobValues.push_back(densityInfo);
										fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, grid, hklinfo, informationVector[vectorIndex][c].PolymerID, r);
										}
									}
							}
						if (inputModel[informationVector[vectorIndex][c].PolymerID][r].type() == "THR")
							{
								int vectorShiftLimit = 5;  
								bool siteAlreadyGlycosylated = check_glycosylation_presence(inputModel[informationVector[vectorIndex][c].PolymerID].id(), inputModel[informationVector[vectorIndex][c].PolymerID][r].id().trim(), glycanList);
								if(!siteAlreadyGlycosylated)
									{     
										               							
										std::list<clipper::String> ignoreAtomList = {" O  ", " N  ", " C  "};
										clipper::MAtom OG1Atom; // On Thr glycan attaches to OG1
										OG1Atom = inputModel[informationVector[vectorIndex][c].PolymerID][r].find(" OG1", clipper::MM::ANY);

										clipper::Coord_orth OG1Coordinate; // OG atom is used as a point of attachement of a glycan
										OG1Coordinate = OG1Atom.coord_orth();

										std::vector<std::pair<clipper::Coord_orth, double>> pairs;

										for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][r].size(); natom++)
										{
												bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id()) != ignoreAtomList.end());
												if(!atomIgnored)
												{
													clipper::Coord_orth vectorOrigin;

													vectorOrigin = inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].coord_orth();
													
													for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
														{
															clipper::Coord_orth potentialTarget = getTargetPoint(OG1Coordinate, vectorOrigin, vectorShift);

															double meanDensityValue = calculateMeanElectronDensityInArea(potentialTarget, sigmaa_dif_map, grid, hklinfo);
															std::pair<clipper::Coord_orth, double> tempDensityInfo(potentialTarget, meanDensityValue);
															pairs.push_back(tempDensityInfo);
														}
												}
										}
											
										const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
										clipper::Coord_orth bestTarget = bestPair->first;
										double bestDensityValue = bestPair->second;

										if(bestDensityValue > 0.09)
										{
										std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, r, vectorIndex}, bestDensityValue);
										finalVectorForBlobValues.push_back(densityInfo);
										fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, grid, hklinfo, informationVector[vectorIndex][c].PolymerID, r);
										}
									}
							}
					}
				}
			}
		}
		if(vectorIndex == 3)
		{
			if (!informationVector[vectorIndex].empty())
			{
				for (int c = 0; c < informationVector[vectorIndex].size(); c++)
				{
					for(int r = informationVector[vectorIndex][c].FirstMMonomer; r < informationVector[vectorIndex][c].LastMMonomer; r++)
					{
						if (inputModel[informationVector[vectorIndex][c].PolymerID][r].type() == "CYS") // in S-Glycosylation a glycan is attached through CYS or MET residue
                            {	
								int vectorShiftLimit = 3;  
								bool siteAlreadyGlycosylated = check_glycosylation_presence(inputModel[informationVector[vectorIndex][c].PolymerID].id(), inputModel[informationVector[vectorIndex][c].PolymerID][r].id().trim(), glycanList);
								if(!siteAlreadyGlycosylated)
									{                								
										std::list<clipper::String> ignoreAtomList = {" O  "};
										clipper::MAtom SGAtom; // On Ser glycan attaches to OG, on Thr glycan attaches to OG1

										SGAtom = inputModel[informationVector[vectorIndex][c].PolymerID][r].find(" SG ", clipper::MM::ANY);

										clipper::Coord_orth SGCoordinate; // SG atom is used as a point of attachement of a glycan
										SGCoordinate = SGAtom.coord_orth();

										std::vector<std::pair<clipper::Coord_orth, double>> pairs;

										for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][r].size(); natom++)
										{
											bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id()) != ignoreAtomList.end());
											if(!atomIgnored)
												{
													clipper::Coord_orth vectorOrigin;

													vectorOrigin = inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].coord_orth();
													
													for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
														{
															clipper::Coord_orth potentialTarget = getTargetPoint(SGCoordinate, vectorOrigin, vectorShift);

															double meanDensityValue = calculateMeanElectronDensityInArea(potentialTarget, sigmaa_dif_map, grid, hklinfo);
															std::pair<clipper::Coord_orth, double> tempDensityInfo(potentialTarget, meanDensityValue);
															pairs.push_back(tempDensityInfo);
														}
												}
										}
										
										const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
										clipper::Coord_orth bestTarget = bestPair->first;
										double bestDensityValue = bestPair->second;

										if(bestDensityValue > 0.09)
										{
										std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, r, vectorIndex}, bestDensityValue);
										finalVectorForBlobValues.push_back(densityInfo);
										fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, grid, hklinfo, informationVector[vectorIndex][c].PolymerID, r);
										}
									}
							}
						if (inputModel[informationVector[vectorIndex][c].PolymerID][r].type() == "MET") // in S-Glycosylation a glycan is attached through CYS or MET residue
							{
								int vectorShiftLimit = 3;
								bool siteAlreadyGlycosylated = check_glycosylation_presence(inputModel[informationVector[vectorIndex][c].PolymerID].id(), inputModel[informationVector[vectorIndex][c].PolymerID][r].id().trim(), glycanList);
								if(!siteAlreadyGlycosylated)
									{   
										               							
										std::list<clipper::String> ignoreAtomList = {""};
										clipper::MAtom SDAtom; // On Thr glycan attaches to OG1
										SDAtom = inputModel[informationVector[vectorIndex][c].PolymerID][r].find(" SD ", clipper::MM::ANY);

										clipper::Coord_orth SDCoordinate; // OG atom is used as a point of attachement of a glycan
										SDCoordinate = SDAtom.coord_orth();

										std::vector<std::pair<clipper::Coord_orth, double>> pairs;

										for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][r].size(); natom++)
										{
												bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id()) != ignoreAtomList.end());
												if(!atomIgnored)
												{
													clipper::Coord_orth vectorOrigin;

													vectorOrigin = inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].coord_orth();
													
													for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
														{
															clipper::Coord_orth potentialTarget = getTargetPoint(SDCoordinate, vectorOrigin, vectorShift);

															double meanDensityValue = calculateMeanElectronDensityInArea(potentialTarget, sigmaa_dif_map, grid, hklinfo);
															std::pair<clipper::Coord_orth, double> tempDensityInfo(potentialTarget, meanDensityValue);
															pairs.push_back(tempDensityInfo);
														}
												}
										}
											
										const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
										clipper::Coord_orth bestTarget = bestPair->first;
										double bestDensityValue = bestPair->second;

										if(bestDensityValue > 0.09)
										{
										std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, r, vectorIndex}, bestDensityValue);
										finalVectorForBlobValues.push_back(densityInfo);
										fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, grid, hklinfo, informationVector[vectorIndex][c].PolymerID, r);
										}
									}
							}
					}
				}
			}
		}
		if(vectorIndex == 4) // this one doesn't make sense, why would you find glycosylation if it was removed... especially in a crystal.
		{
			if (!informationVector[vectorIndex].empty())
			{
				for (int c = 0; c < informationVector[vectorIndex].size(); c++)
				{
					for(int r = informationVector[vectorIndex][c].FirstMMonomer; r < informationVector[vectorIndex][c].LastMMonomer; r++)
					{
						if (inputModel[informationVector[vectorIndex][c].PolymerID][r].type() == "ALA") 
                            {
								int vectorShiftLimit = 5;
								bool siteAlreadyGlycosylated = check_glycosylation_presence(inputModel[informationVector[vectorIndex][c].PolymerID].id(), inputModel[informationVector[vectorIndex][c].PolymerID][r].id().trim(), glycanList);
								if(!siteAlreadyGlycosylated)
									{                   								
										std::list<clipper::String> ignoreAtomList = {""};
										clipper::MAtom CBAtom; // CA to CB for Ala, CG to NE2 for GLN
										CBAtom = inputModel[informationVector[vectorIndex][c].PolymerID][r].find(" CB ", clipper::MM::ANY);

										clipper::Coord_orth CBCoordinate; // CB atom is used as a direction towards the glycan density
										CBCoordinate = CBAtom.coord_orth();

										std::vector<std::pair<clipper::Coord_orth, double>> pairs;

										for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][r].size(); natom++)
										{
											
											bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id()) != ignoreAtomList.end());
											if(!atomIgnored)
												{
													clipper::Coord_orth vectorOrigin;

													vectorOrigin = inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].coord_orth();
													
													for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
														{
															clipper::Coord_orth potentialTarget = getTargetPoint(CBCoordinate, vectorOrigin, vectorShift);

															double meanDensityValue = calculateMeanElectronDensityInArea(potentialTarget, sigmaa_dif_map, grid, hklinfo);
															std::pair<clipper::Coord_orth, double> tempDensityInfo(potentialTarget, meanDensityValue);
															pairs.push_back(tempDensityInfo);
														}
												}
										}
										const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
										clipper::Coord_orth bestTarget = bestPair->first;
										double bestDensityValue = bestPair->second;

										if(bestDensityValue > 0.09)
										{
										std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, r, vectorIndex}, bestDensityValue);
										finalVectorForBlobValues.push_back(densityInfo);
										fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, grid, hklinfo, informationVector[vectorIndex][c].PolymerID, r);
										}
									}	
							}
						if (inputModel[informationVector[vectorIndex][c].PolymerID][r].type() == "GLN") 
                            {	
								int vectorShiftLimit = 5;
								bool siteAlreadyGlycosylated = check_glycosylation_presence(inputModel[informationVector[vectorIndex][c].PolymerID].id(), inputModel[informationVector[vectorIndex][c].PolymerID][r].id().trim(), glycanList);
								if(!siteAlreadyGlycosylated)
									{                   								
										std::list<clipper::String> ignoreAtomList = {""};
										
										clipper::MAtom CGAtom; // CA to CB for Ala, CG to NE2 for GLN
										CGAtom = inputModel[informationVector[vectorIndex][c].PolymerID][r].find(" CG ", clipper::MM::ANY);

										clipper::Coord_orth CGCoordinate; // CB atom is used as a direction towards the glycan density
										CGCoordinate = CGAtom.coord_orth();

										std::vector<std::pair<clipper::Coord_orth, double>> pairs;

										for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][r].size(); natom++)
										{
											
											bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id()) != ignoreAtomList.end());
											if(!atomIgnored)
												{
													clipper::Coord_orth vectorOrigin;

													vectorOrigin = inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].coord_orth();
													
													for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
														{
															clipper::Coord_orth potentialTarget = getTargetPoint(vectorOrigin, CGCoordinate, vectorShift);

															double meanDensityValue = calculateMeanElectronDensityInArea(potentialTarget, sigmaa_dif_map, grid, hklinfo);
															std::pair<clipper::Coord_orth, double> tempDensityInfo(potentialTarget, meanDensityValue);
															pairs.push_back(tempDensityInfo);
														}
												}
										}
										const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
										clipper::Coord_orth bestTarget = bestPair->first;
										double bestDensityValue = bestPair->second;

										if(bestDensityValue > 0.09)
										{
										std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, r, vectorIndex}, bestDensityValue);
										finalVectorForBlobValues.push_back(densityInfo);
										fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, grid, hklinfo, informationVector[vectorIndex][c].PolymerID, r);
										}
									}
							}
					}
				}
			}
		}
	return finalVectorForBlobValues;
}
	


//Compiler settings: -I/y/people/hb1115/devtools/install/include -L/y/people/hb1115/devtools/install/lib -lclipper-minimol -lclipper-core -lclipper-mmdb -lmmdb2 -lclipper-ccp4

