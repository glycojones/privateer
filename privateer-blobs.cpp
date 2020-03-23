
// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2018 Haroldas Bagdonas & Kevin Cowtan & Jon Agirre
// York Structural Biology Laboratory
// The University of York
// mailto: hb1115@york.ac.uk
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk



#include "privateer-blobs.h"


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
	std::vector<GlycosylationMonomerMatch> NRemMMonomers;
	std::vector<std::vector<GlycosylationMonomerMatch> > MatchingMMonomers;


	std::regex isProteinNGlc("[N][A-Z^P][ST]|[N][A-Z][C]"); // N-Glycosylation = Asn-X-/Ser/Thr(X=anything except Proline) || Asn-X-Cys(X=Anything)
	std::regex isProteinNGlcRemoved("[A][A-Z^P][ST]|[A][A-Z][C]|[Q][A-Z^P][ST]|[Q][A-Z][C]");
	std::regex isProteinCGlc("[W][A-Z][A-Z][W]|[W][ST][A-Z][C]"); // C-Glycosylation = Trp-X-X-Trp || Trp-Ser/Thr-X-Cys
	std::regex isProteinOGlc("[N][A-Z][T]|[N][A-Z][S]|[S][A-Z][A-Z][P]|[P][A-Z][T]|[T][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][P]|[T][A-Z][A-Z][P]");
	// Need to patch O-glycosylation regex. Not enough sequence patterns, a lot of matching is missed on actual cases. 

	std::smatch NMatch;
	std::smatch CMatch;
	std::smatch OMatch;
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
		 auto OGlc_begin = std::sregex_iterator(seq.begin(), seq.end(), isProteinOGlc); // same as above, just for C-glycosylation
		 auto OGlc_end = std::sregex_iterator();
		 auto NRemGlc_begin = std::sregex_iterator(seq.begin(), seq.end(), isProteinNGlcRemoved); // same as above, just for C-glycosylation
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
					molwrk_new.insert(mp);
					mp = clipper::MPolymer();
				}
		}
	 }

	return molwrk_new;
}



std::vector<std::pair<PotentialGlycosylationSiteInfo, double> > get_electron_density_of_potential_glycosylation_sites(const std::vector<std::vector<GlycosylationMonomerMatch>>& informationVector, int vectorIndex, clipper::MiniMol& inputModel, clipper::Xmap<float>& sigmaa_dif_map, clipper::Grid_sampling& grid, clipper::HKL_info& hklinfo, bool pdbexport) 
{
	std::vector<std::pair<PotentialGlycosylationSiteInfo, double> > finalVectorForBlobValues;
		if(vectorIndex == 0)
		{
			if (!informationVector[vectorIndex].empty())
			{
				for (int c = 0; c < informationVector[vectorIndex].size(); c++)
				{
					for(int r = informationVector[vectorIndex][c].FirstMMonomer; r < informationVector[vectorIndex][c].LastMMonomer; r++)
					{
						if (inputModel[informationVector[vectorIndex][c].PolymerID][r].type() == "ASN") // in N-Glycosylation a glycan is attached through ASN residue
                            {	
                                clipper::Coord_orth ND2Coordinate; // A glycan is attached to ASN residue via ND2 atom
                                clipper::Coord_orth CBCoordinate; // CB atom is used as a direction towards the glycan density

                                // looping through atoms of ASN residue to find ND2 and CB atoms.
                                for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][r].size(); natom++)
                                    {
                                        if(inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id() == " ND2")
                                        ND2Coordinate = inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].coord_orth();

                                        if(inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id() == " CB ")
                                        CBCoordinate = inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].coord_orth();
                                    }



                                // Create a vector between ND2 and CB
                                clipper::Vec3<clipper::ftype> baseVector((ND2Coordinate.x()-CBCoordinate.x()),(ND2Coordinate.y()-CBCoordinate.y()), (ND2Coordinate.z()-CBCoordinate.z()));
                                // Create a 1A unit vector out of baseVector, to be used later in vector shifting
                                clipper::Vec3<clipper::ftype> unitVector = baseVector.unit();

                                // Obtain coordinates in the middle of suspected glycan density via 5A vector shift. This is the nearest glycan bonded via ND2 atom to ASN residue.
                                clipper::Coord_orth target( (CBCoordinate.x()+(unitVector[0]*5)), (CBCoordinate.y()+(unitVector[1]*5)), (CBCoordinate.z()+(unitVector[2]*5)) );

								if(pdbexport)
								{
									clipper::Atom dummyAtom;
									dummyAtom.set_coord_orth(target);
									dummyAtom.set_element("B");
									clipper::MAtom dummyAtomExport(dummyAtom);
									inputModel[informationVector[vectorIndex][c].PolymerID][r].insert(dummyAtomExport);
								}


                                double meanDensityExp = 0.0;
                                int n_points = 0;

                                // Define origin and destination for drawing the sphere. Electron density data obtained from within the sphere later on.
                                clipper::Coord_orth origin(target.x()-2, target.y()-2, target.z()-2);
                                clipper::Coord_orth destination(target.x()+2, target.y()+2, target.z()+2);

                                clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;


                                i0 = clipper::Xmap_base::Map_reference_coord( sigmaa_dif_map, origin.coord_frac(hklinfo.cell()).coord_grid(grid) );

								for ( iu = i0; iu.coord().u() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).u(); iu.next_u() )
                                    for ( iv = iu; iv.coord().v() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).v(); iv.next_v() )
                                        for ( iw = iv; iw.coord().w() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).w(); iw.next_w() )
                                            {
                                                meanDensityExp = meanDensityExp + sigmaa_dif_map[iw];
                                                n_points++;
                                            }

                                meanDensityExp = meanDensityExp / n_points;
								
								std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, r, vectorIndex}, meanDensityExp);

								finalVectorForBlobValues.push_back(densityInfo);
							}
					}
				}
			}
		}
		if(vectorIndex == 1)
		{
			if (!informationVector[vectorIndex].empty())
			{
				for (int c = 0; c < informationVector[vectorIndex].size(); c++)
				{
					for(int r = informationVector[vectorIndex][c].FirstMMonomer; r < informationVector[vectorIndex][c].LastMMonomer; r++)
					{
						if (inputModel[informationVector[vectorIndex][c].PolymerID][r].type() == "TRP") // // in C-Glycosylation a glycan is attached through TRP residue
                            {	
                                clipper::Coord_orth CD1Coordinate; // A glycan is attached to ASN residue via ND2 atom
                                clipper::Coord_orth CZ3Coordinate;

                                for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][r].size(); natom++)
                                    {
                                        if(inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id() == " CD1")
                                        CD1Coordinate = inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].coord_orth();

                                        if(inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id() == " CZ3")
                                        CZ3Coordinate = inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].coord_orth();
                                    }

                                    // Create a vector between ND2 and CB
                                clipper::Vec3<clipper::ftype> baseVector((CD1Coordinate.x()-CZ3Coordinate.x()),(CD1Coordinate.y()-CZ3Coordinate.y()), (CD1Coordinate.z()-CZ3Coordinate.z()));
                                    // Create a 1A unit vector out of baseVector, to be used later in vector shifting
                                clipper::Vec3<clipper::ftype> unitVector = baseVector.unit();

                                    // Create a target coordinate by applying 7A shift via unitVector, specific to Trp only.
                                clipper::Coord_orth target( (CZ3Coordinate.x()+(unitVector[0]*7)), (CZ3Coordinate.y()+(unitVector[1]*7)), (CZ3Coordinate.z()+(unitVector[2]*7)) );

								if(pdbexport)
								{
									clipper::Atom dummyAtom;
									dummyAtom.set_coord_orth(target);
									dummyAtom.set_element("B");
									clipper::MAtom dummyAtomExport(dummyAtom);
									inputModel[informationVector[vectorIndex][c].PolymerID][r].insert(dummyAtomExport);
								}

                                double meanDensityExp = 0.0;
                                int n_points = 0;

                                // Define origin and destination for drawing the sphere. Electron density data obtained from within the sphere later on.
                                clipper::Coord_orth origin(target.x()-2, target.y()-2, target.z()-2);
                                clipper::Coord_orth destination(target.x()+2, target.y()+2, target.z()+2);

                                clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;


                                i0 = clipper::Xmap_base::Map_reference_coord( sigmaa_dif_map, origin.coord_frac(hklinfo.cell()).coord_grid(grid) );

								for ( iu = i0; iu.coord().u() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).u(); iu.next_u() )
                                    for ( iv = iu; iv.coord().v() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).v(); iv.next_v() )
                                        for ( iw = iv; iw.coord().w() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).w(); iw.next_w() )
                                            {
                                                meanDensityExp = meanDensityExp + sigmaa_dif_map[iw];
                                                n_points++;
                                            }

                                meanDensityExp = meanDensityExp / n_points;

								std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, r, vectorIndex}, meanDensityExp);
								
								finalVectorForBlobValues.push_back(densityInfo);
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
						if (inputModel[informationVector[vectorIndex][c].PolymerID][r].type() == "SER" || inputModel[informationVector[vectorIndex][c].PolymerID][r].type() == "THR") // in O-Glycosylation a glycan is attached through Thr or Ser residue
                            {	
								clipper::Coord_orth CBorNE2Coordinate; // CA to CB for Ala, CG to NE2 for GLN
                                clipper::Coord_orth CAorCGCoordinate;

                                for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][r].size(); natom++)
                                    {
                                    if(inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id() == " CB " || inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id() == " NE2")
                                        CBorNE2Coordinate = inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].coord_orth();

                                    if(inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id() == " CA " || inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id() == " CG ")
                                        CAorCGCoordinate = inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].coord_orth();
                                    }

                                    // Create a vector between CA to CB for Ala and CG to NE2 for GLN
                                clipper::Vec3<clipper::ftype> baseVector((CBorNE2Coordinate.x()-CAorCGCoordinate.x()),(CBorNE2Coordinate.y()-CAorCGCoordinate.y()), (CBorNE2Coordinate.z()-CAorCGCoordinate.z()));
                                    // Create a 1A unit vector out of baseVector, to be used later in vector shifting
                                clipper::Vec3<clipper::ftype> unitVector = baseVector.unit();
                                    // Create a target coordinate by applying 5A shift via unitVector. Values may need adjusting as haven't found a proper pdb file to test with. Significant differences between Ala and Gln
                                clipper::Coord_orth target( (CAorCGCoordinate.x()+(unitVector[0]*5)), (CAorCGCoordinate.y()+(unitVector[1]*5)), (CAorCGCoordinate.z()+(unitVector[2]*5)) );

								if(pdbexport)
								{
									clipper::Atom dummyAtom;
									dummyAtom.set_coord_orth(target);
									dummyAtom.set_element("B");
									clipper::MAtom dummyAtomExport(dummyAtom);
									inputModel[informationVector[vectorIndex][c].PolymerID][r].insert(dummyAtomExport);
								}

                                double meanDensityExp = 0.0;
                                int n_points = 0;

                                clipper::Coord_orth origin(target.x()-2, target.y()-2, target.z()-2);
                                clipper::Coord_orth destination(target.x()+2, target.y()+2, target.z()+2);

                                clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;


                                i0 = clipper::Xmap_base::Map_reference_coord( sigmaa_dif_map, origin.coord_frac(hklinfo.cell()).coord_grid(grid) );

								for ( iu = i0; iu.coord().u() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).u(); iu.next_u() )
                                    for ( iv = iu; iv.coord().v() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).v(); iv.next_v() )
                                        for ( iw = iv; iw.coord().w() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).w(); iw.next_w() )
                                            {
                                                meanDensityExp = meanDensityExp + sigmaa_dif_map[iw];
                                                n_points++;
                                            }

								meanDensityExp = meanDensityExp / n_points;

								std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, r, vectorIndex}, meanDensityExp);
								
								finalVectorForBlobValues.push_back(densityInfo);
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
						if (inputModel[informationVector[vectorIndex][c].PolymerID][r].type() == "ALA" || inputModel[informationVector[vectorIndex][c].PolymerID][r].type() == "GLN") // if N-Glycosylation removed via pngase F or pngase A, then residue gets replaced to GLN or ALA. 
                            {	
                                clipper::Coord_orth OGCoordinate; // On Ser glycan attaches to OG, on Thr glycan attaches to OG1
                                clipper::Coord_orth CACoordinate; // CA coordinate used in both Ser and Thr to establish a proper direction for the vector.

                                for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][r].size(); natom++)
                                    {
                                        if(inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id() == " OG " || inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id() == " OG1")
                                        OGCoordinate = inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].coord_orth();

                                        if(inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].id() == " CA ")
                                        CACoordinate = inputModel[informationVector[vectorIndex][c].PolymerID][r][natom].coord_orth();
                                    }

                                    // Create a vector between OG/OG1 and CA for either Ser or Thr
                                clipper::Vec3<clipper::ftype> baseVector((OGCoordinate.x()-CACoordinate.x()),(OGCoordinate.y()-CACoordinate.y()), (OGCoordinate.z()-CACoordinate.z()));
                                    // Create a 1A unit vector out of baseVector, to be used later in vector shifting
                                clipper::Vec3<clipper::ftype> unitVector = baseVector.unit();

                                    // Create a target coordinate by applying 5A shift via unitVector
                                clipper::Coord_orth target( (CACoordinate.x()+(unitVector[0]*5)), (CACoordinate.y()+(unitVector[1]*5)), (CACoordinate.z()+(unitVector[2]*5)) );


								if(pdbexport)
								{
									clipper::Atom dummyAtom;
									dummyAtom.set_coord_orth(target);
									dummyAtom.set_element("B");
									clipper::MAtom dummyAtomExport(dummyAtom);
									inputModel[informationVector[vectorIndex][c].PolymerID][r].insert(dummyAtomExport);
								}

                                double meanDensityExp = 0.0;
                                int n_points = 0;

                                // Define origin and destination for drawing the sphere. Electron density data obtained from within the sphere later on.
                                clipper::Coord_orth origin(target.x()-2, target.y()-2, target.z()-2);
                                clipper::Coord_orth destination(target.x()+2, target.y()+2, target.z()+2);

                                clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;


                                i0 = clipper::Xmap_base::Map_reference_coord( sigmaa_dif_map, origin.coord_frac(hklinfo.cell()).coord_grid(grid) );

								for ( iu = i0; iu.coord().u() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).u(); iu.next_u() )
                                    for ( iv = iu; iv.coord().v() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).v(); iv.next_v() )
                                        for ( iw = iv; iw.coord().w() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).w(); iw.next_w() )
                                            {
                                                meanDensityExp = meanDensityExp + sigmaa_dif_map[iw];
                                                n_points++;
                                            }

								meanDensityExp = meanDensityExp / n_points;

								std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, r, vectorIndex}, meanDensityExp);

								finalVectorForBlobValues.push_back(densityInfo);
							}
					}
				}
			}
		}
	return finalVectorForBlobValues;
}
	


//Compiler settings: -I/y/people/hb1115/devtools/install/include -L/y/people/hb1115/devtools/install/lib -lclipper-minimol -lclipper-core -lclipper-mmdb -lmmdb2 -lclipper-ccp4
