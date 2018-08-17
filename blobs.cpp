
// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2018 Haroldas Bagdonas & Kevin Cowtan & Jon Agirre
// York Structural Biology Laboratory
// The University of York
// mailto: hb1115@york.ac.uk
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk



#include "blobs.h"


std::vector<std::vector<GlycosylationMonomerMatch> > get_matching_mmonomer_positons(const clipper::String& ippdb)
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
	std::regex isProteinCGlc("[W][A-Z][A-Z][W]|[W][ST][A-Z][C]"); // C-Glycosylation = Trp-X-X-Trp || Trp-Ser/Thr/X-Cys
	std::regex isProteinOGlc("[N][A-Z][T]|[N][A-Z][S]|[S][A-Z][A-Z][P]|[P][A-Z][T]|[T][A-Z][A-Z][A-Z][A-Z][A-Z][A-Z][P]|[T][A-Z][A-Z][P]");

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

//Compiler settings: -I/y/people/hb1115/devtools/install/include -L/y/people/hb1115/devtools/install/lib -lclipper-minimol -lclipper-core -lclipper-mmdb -lmmdb2 -lclipper-ccp4

int main(int argc, char** argv)
{

	CCP4Program prog( "GlycoSeqChecker", "0.3", "$Date: 2018/08/15" );

	// defaults
	clipper::String title 	 = "No Title";
	clipper::String ippdb    = "No File";


	// command input
	CCP4CommandInput args( argc, argv, true );
	int arg = 0;
	while ( ++arg < args.size() ) {
	  if ( args[arg] == "-title" ) {
		if ( ++arg < args.size() ) title = args[arg];
	  } else if ( args[arg] == "-pdbin" ) {
		if ( ++arg < args.size() ) ippdb = args[arg];
	  } else {
		std::cout << "Unrecognized:\t" << args[arg] << "\n";
		args.clear();
	  }
	}
	if ( args.size() <= 1 ) {
	  std::cout << "Usage: GlycoSeqChecker\n\t-pdbin <filename>\n\t\n-title <protein_title>\n";
	  exit(1);
	}

	std::cout << copy[0][1].LastMMonomer; // [0] = N-Glycosylation MMonomers std::vector, [1]=Second Match info of NGlc pattern in struct form.


}
