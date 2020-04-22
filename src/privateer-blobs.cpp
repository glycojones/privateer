
// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2018 Haroldas Bagdonas & Kevin Cowtan & Jon Agirre
// York Structural Bioexpy Laboratory
// The University of York
// mailto: hb1115@york.ac.uk
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk



#include "../include/privateer-blobs.h"

bool bestPointFinder(std::pair<clipper::Coord_orth, double> p1, std::pair<clipper::Coord_orth, double> p2) {
    return p1.second<p2.second;
}

std::vector<std::vector<GlycosylationMonomerMatch> > get_matching_monomer_positions(clipper::MiniMol& inputModel)
{

	std::vector<GlycosylationMonomerMatch> NMMonomers;
	std::vector<GlycosylationMonomerMatch> CMMonomers;
	std::vector<GlycosylationMonomerMatch> OMMonomers;
	std::vector<GlycosylationMonomerMatch> SMMonomers;
	std::vector<GlycosylationMonomerMatch> NRemMMonomers;
	std::vector<std::vector<GlycosylationMonomerMatch> > MatchingMMonomers;



   for ( int pol = 0; pol < inputModel.size(); pol++ ) 
   {
	 for ( int mon = 0; mon < inputModel[pol].size(); mon++ ) 
	 {
		if (inputModel[pol][mon].type() == "ASN" || inputModel[pol][mon].type() == "ARG")
		{
			NMMonomers.push_back({pol, mon});
		}
		
		if (inputModel[pol][mon].type() == "TRP")
		{
			CMMonomers.push_back({pol, mon});
		}

		if (inputModel[pol][mon].type() == "SER" || inputModel[pol][mon].type() == "THR")
		{
			OMMonomers.push_back({pol, mon});
		}

		if (inputModel[pol][mon].type() == "CYS" || inputModel[pol][mon].type() == "MET")
		{
			SMMonomers.push_back({pol, mon});
		}

		if (inputModel[pol][mon].type() == "ALA" || inputModel[pol][mon].type() == "GLN")
		{
			NRemMMonomers.push_back({pol, mon});
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


void fillSearchArea(clipper::MiniMol& inputModel, clipper::Coord_orth& targetPos, clipper::Xmap<float>& sigmaa_dif_map, clipper::HKL_info& hklinfo, clipper::Map_stats& mapstats, int chainID, int monomerID)
{
	// Define origin and destination for drawing the sphere. Electron density data obtained from within the sphere later on.
	clipper::Coord_orth origin(targetPos.x()-2, targetPos.y()-2, targetPos.z()-2);
	clipper::Coord_orth destination(targetPos.x()+2, targetPos.y()+2, targetPos.z()+2);

	clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;

	clipper::Grid_sampling grid( hklinfo.spacegroup(), hklinfo.cell(), hklinfo.resolution() );  // define grid

	if(!origin.is_null() && !destination.is_null())
	{
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
 } 


void drawOriginPoint(clipper::MiniMol& inputModel, clipper::Coord_orth target, int chainID, int monomerID)
{
		clipper::U_aniso_orth null(0);

		clipper::Atom dummyAtom;
		dummyAtom.set_coord_orth(target);
		dummyAtom.set_element("H");
		dummyAtom.set_occupancy(0);
		dummyAtom.set_u_aniso_orth(null);
		dummyAtom.set_u_iso(0);
		clipper::MAtom dummyAtomExport(dummyAtom);
		dummyAtomExport.set_id(" DUM");
		inputModel[chainID][monomerID].insert(dummyAtomExport);
}


GlycanToMiniMolIDs getCarbohydrateRelationshipToMiniMol(clipper::MiniMol& inputModel, clipper::MSugar& carbohydrate, std::vector < clipper::MGlycan >& allSugars, int mglycanid, int sugaringlycanid)
{
	clipper::String glycanAttachedTo = allSugars[mglycanid].get_chain();
	int backboneID;
	int carbohydrateID;
	GlycanToMiniMolIDs output;

	for(int chain = 0; chain < inputModel.size(); chain++)
	{
		if(inputModel[chain].id() == glycanAttachedTo) 
		{
			backboneID = chain;
			break;
		}
	}
	
	for(int residue = 0; residue < inputModel[backboneID].size(); residue++)
	{
		if(inputModel[backboneID][residue].id() == carbohydrate.id())
		{
			carbohydrateID = residue;
			break;
		}
	}

	output = {backboneID, carbohydrateID, sugaringlycanid};
	return output;
}

// FIX ME: I mistook arguments in coot code. box_radius = radius, not contour_level. contour_level is tIsoLevel. Need to rewrite this bit.
double calculateMeanElectronDensityInTargetPosition(clipper::Coord_orth targetPos, clipper::Xmap<float>& sigmaa_dif_map, clipper::Map_stats& mapstats)
{
	double meanElectronDensity = 0.0;
	int n_points = 0;

	float map_sigma = mapstats.std_dev();
	float box_radius = 5.00 * map_sigma;

	int isample_step = 1;

		// Define origin and destination for drawing the sphere. Electron density data obtained from within the sphere later on.
	clipper::Coord_orth origin(targetPos.x()-0.8, targetPos.y()-0.8, targetPos.z()-0.8);
	clipper::Coord_orth destination(targetPos.x()+0.8, targetPos.y()+0.8, targetPos.z()+0.8);

	clipper::Coord_frac originref = origin.coord_frac(sigmaa_dif_map.cell());
	clipper::Coord_frac destinationref = destination.coord_frac(sigmaa_dif_map.cell());
	
	clipper::Coord_frac origin0(
			    originref.u() - box_radius/sigmaa_dif_map.cell().descr().a(),
			    originref.v() - box_radius/sigmaa_dif_map.cell().descr().b(),
			    originref.w() - box_radius/sigmaa_dif_map.cell().descr().c() );
   	clipper::Coord_frac origin1(
			    originref.u() + box_radius/sigmaa_dif_map.cell().descr().a(),
			    originref.v() + box_radius/sigmaa_dif_map.cell().descr().b(),
			    originref.w() + box_radius/sigmaa_dif_map.cell().descr().c() );

	clipper::Coord_frac destination0(
			    destinationref.u() - box_radius/sigmaa_dif_map.cell().descr().a(),
			    destinationref.v() - box_radius/sigmaa_dif_map.cell().descr().b(),
			    destinationref.w() - box_radius/sigmaa_dif_map.cell().descr().c() );
   	clipper::Coord_frac destination1(
			    destinationref.u() + box_radius/sigmaa_dif_map.cell().descr().a(),
			    destinationref.v() + box_radius/sigmaa_dif_map.cell().descr().b(),
			    destinationref.w() + box_radius/sigmaa_dif_map.cell().descr().c() );

	clipper::Grid_range gridorigin(origin0.coord_grid(sigmaa_dif_map.grid_sampling()),
			origin1.coord_grid(sigmaa_dif_map.grid_sampling()));

	clipper::Grid_range griddestination(destination0.coord_grid(sigmaa_dif_map.grid_sampling()),
			destination1.coord_grid(sigmaa_dif_map.grid_sampling()));
	

    //  std::cout << "    box_radius " << box_radius << std::endl;
    //  std::cout << "    centre_point: " << targetPos.format() << std::endl;
    //  std::cout << "    origin0: " << origin0.format() << std::endl;
    //  std::cout << "    origin1: " << origin1.format() << std::endl;
    //  std::cout << "    destination0: " << destination0.format() << std::endl;
    //  std::cout << "    destination1: " << destination1.format() << std::endl;
    //  std::cout << "    origin0.coord_orth(): " << origin0.coord_orth(sigmaa_dif_map.cell()).format() << std::endl;
    //  std::cout << "    origin1.coord_orth(): " << origin1.coord_orth(sigmaa_dif_map.cell()).format() << std::endl;
    //  std::cout << "    destination0.coord_orth(): " << destination0.coord_orth(sigmaa_dif_map.cell()).format() << std::endl;
    //  std::cout << "    destination1.coord_orth(): " << destination1.coord_orth(sigmaa_dif_map.cell()).format() << std::endl;
    //  std::cout << "    grid: " << grid.format() << std::endl;
	//  std::cout << "    gridorigin: " << gridorigin.format() << std::endl;
	//  std::cout << "    griddestination: " << griddestination.format() << std::endl;
	if(!originref.is_null() && !destinationref.is_null())
	{
		clipper::Xmap_base::Map_reference_coord ix( sigmaa_dif_map );
		int u, v, w, ii;
		for (w = gridorigin.min().w(); w <= griddestination.max().w(); w+=isample_step ) {
			for (v = gridorigin.min().v(); v <= griddestination.max().v(); v+=isample_step) {
				ix.set_coord(clipper::Coord_grid( gridorigin.min().u(), v, w ));
				for (u = gridorigin.min().u(); u <= griddestination.max().u(); u+=isample_step) 
					{
						// std::cout << "    ix: " << ix.coord_orth().format() << std::endl;
						meanElectronDensity = meanElectronDensity + sigmaa_dif_map[ix];
						n_points++;
						for(ii=0; ii<isample_step; ii++) ix.next_u(); 
					} 
			}
		}
	}


	meanElectronDensity = meanElectronDensity / n_points;

	return meanElectronDensity;
}

// FIX ME: I mistook arguments in coot code. box_radius = radius, not contour_level. contour_level is tIsoLevel. Need to rewrite this bit. 
double calculateMeanElectronDensityForBiggerSphere(clipper::Coord_orth& targetPos, clipper::Xmap<float>& sigmaa_dif_map, clipper::Map_stats& mapstats, clipper::HKL_info& hklinfo)
{
	double meanElectronDensity = 0.0;
	int n_points = 0;

	float map_sigma = mapstats.std_dev();
	float rmsdLimit = 3.00 * map_sigma;

	clipper::Grid_sampling grid( hklinfo.spacegroup(), hklinfo.cell(), hklinfo.resolution() );  // define grid

	// Define origin and destination for drawing the sphere. Electron density data obtained from within the sphere later on.
	clipper::Coord_orth origin(targetPos.x()-2, targetPos.y()-2, targetPos.z()-2);
	clipper::Coord_orth destination(targetPos.x()+2, targetPos.y()+2, targetPos.z()+2);

	clipper::Xmap_base::Map_reference_coord i0, iu, iv, iw;

	if(!origin.is_null() && !destination.is_null())
	{
		i0 = clipper::Xmap_base::Map_reference_coord( sigmaa_dif_map, origin.coord_frac(hklinfo.cell()).coord_grid(grid) );

		for ( iu = i0; iu.coord().u() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).u(); iu.next_u() )
			for ( iv = iu; iv.coord().v() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).v(); iv.next_v() )
				for ( iw = iv; iw.coord().w() <= destination.coord_frac(hklinfo.cell()).coord_grid(grid).w(); iw.next_w() )
					{
						if(sigmaa_dif_map[iw] > rmsdLimit)
						meanElectronDensity = meanElectronDensity + sigmaa_dif_map[iw];
						n_points++;
					}
	}
 
	meanElectronDensity = meanElectronDensity / n_points;

	return meanElectronDensity;
}
 




std::vector<clipper::String> create_list_of_ignored_sugar_atoms(clipper::MSugar& carbohydrate)
{
	std::vector<clipper::String> ignoreAtomList;
	std::vector<clipper::MAtom> ringMembers = carbohydrate.ring_members();
	

	for(int i = 0; i < ringMembers.size(); i++)
	{
		clipper::String atomID = ringMembers[i].id();
		ignoreAtomList.push_back(atomID);
	}

	for(int atom = 0; atom < carbohydrate.size(); atom++)
	{
		bool atomAlreadyIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), carbohydrate[atom].id()) != ignoreAtomList.end());
		if(!atomAlreadyIgnored)
			if(carbohydrate[atom].element() != " O") ignoreAtomList.push_back(carbohydrate[atom].id());
	}
	return ignoreAtomList;
}


// TO DO after: possible improvements, after determining best point, expand the cube at that point to get all electron density and see whether there would be discernible difference between false positives and true positives. 
std::vector<std::pair<PotentialGlycosylationSiteInfo, double> > get_electron_density_of_potential_glycosylation_sites(const std::vector<std::vector<GlycosylationMonomerMatch>>& informationVector, int vectorIndex, clipper::MiniMol& inputModel, clipper::Xmap<float>& sigmaa_dif_map, clipper::HKL_info& hklinfo, std::vector < clipper::MGlycan >& glycanList, clipper::Map_stats& mapstats, float thresholdED, bool pdbexport) 
{
	float thresholdEDBestBlob = 0.070;
	std::vector<std::pair<PotentialGlycosylationSiteInfo, double> > finalVectorForBlobValues;
		if(vectorIndex == 0)
		{
			if (!informationVector[vectorIndex].empty())
			{
				std::list<clipper::String> ignoreAtomList = {" ND2"};
				int vectorShiftLimit = 5;
				for (int c = 0; c < informationVector[vectorIndex].size(); c++)
				{
						if (inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() == "ASN") // in N-Glycosylation a glycan is attached through ASN residue
                            {	
								bool siteAlreadyGlycosylated = check_glycosylation_presence(inputModel[informationVector[vectorIndex][c].PolymerID].id(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id().trim(), glycanList);
								if(!siteAlreadyGlycosylated)
									{
										clipper::MAtom ND2Atom;
										
										try {
										ND2Atom = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].find(" ND2", clipper::MM::ANY);
										} catch (const clipper::Message_fatal& error) {
										std::cerr << "Unable to find necessary ND2 atom for residue" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id() << "-" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() << " in Chain " << inputModel[informationVector[vectorIndex][c].PolymerID].id() << "\n" << "\n";
										}
										
										clipper::Coord_orth ND2Coordinate; // ND2 atom is used as a direction towards the glycan density
										ND2Coordinate = ND2Atom.coord_orth();

										std::vector<std::pair<clipper::Coord_orth, double>> pairs;

										for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].size(); natom++)
										{
											bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].id()) != ignoreAtomList.end());
											if(!atomIgnored)
												{
													clipper::Coord_orth vectorOrigin;

													vectorOrigin = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].coord_orth();
													
													for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
														{
															clipper::Coord_orth potentialTarget = getTargetPoint(ND2Coordinate, vectorOrigin, vectorShift);

															double meanDensityValue = calculateMeanElectronDensityInTargetPosition(potentialTarget, sigmaa_dif_map, mapstats);
															std::pair<clipper::Coord_orth, double> tempDensityInfo(potentialTarget, meanDensityValue);
															pairs.push_back(tempDensityInfo);
														}
												}
										}
										const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
										clipper::Coord_orth bestTarget = bestPair->first;
										double bestDensityValue = bestPair->second;

										if(bestDensityValue > thresholdEDBestBlob)
										{
											double meanDensityValueBiggerArea = calculateMeanElectronDensityForBiggerSphere(bestTarget, sigmaa_dif_map, mapstats, hklinfo);
											if(meanDensityValueBiggerArea > thresholdED)
											{
												std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID, vectorIndex}, meanDensityValueBiggerArea);
												finalVectorForBlobValues.push_back(densityInfo);

													if(pdbexport)
													{
														drawOriginPoint(inputModel, bestTarget, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
														// fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, hklinfo, mapstats, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
													}
											}
										}
									}
							}
				}
			}
			if (!informationVector[vectorIndex].empty())
			{
				std::list<clipper::String> ignoreAtomList = {" NH2", " NH1"};
				int vectorShiftLimit = 5;
				for (int c = 0; c < informationVector[vectorIndex].size(); c++)
				{
						if (inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() == "ARG") // in N-Glycosylation a glycan is attached through ARG residue
                            {	
								bool siteAlreadyGlycosylated = check_glycosylation_presence(inputModel[informationVector[vectorIndex][c].PolymerID].id(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id().trim(), glycanList);
								if(!siteAlreadyGlycosylated)
									{
										clipper::MAtom NH1Atom;
										clipper::MAtom NH2Atom;
										
										try {
										NH1Atom = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].find(" NH1", clipper::MM::ANY);
										} catch (const clipper::Message_fatal& error) {
										std::cerr << "Unable to find necessary NH1 atom for residue" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id() << "-" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() << " in Chain " << inputModel[informationVector[vectorIndex][c].PolymerID].id() << "\n" << "\n";
										}

										try {
										NH2Atom = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].find(" NH2", clipper::MM::ANY);
										} catch (const clipper::Message_fatal& error) {
										std::cerr << "Unable to find necessary NH2 atom for residue" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id() << "-" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() << " in Chain " << inputModel[informationVector[vectorIndex][c].PolymerID].id() << "\n" << "\n";
										}
										
										clipper::Coord_orth NH1Coordinate; // NH1 atom is used as a direction towards the glycan density
										NH1Coordinate = NH1Atom.coord_orth();

										clipper::Coord_orth NH2Coordinate; // NH1 atom is used as a direction towards the glycan density
										NH2Coordinate = NH2Atom.coord_orth();

										std::vector<std::pair<clipper::Coord_orth, double>> pairs;

										for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].size(); natom++)
										{
											bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].id()) != ignoreAtomList.end());
											if(!atomIgnored)
												{
													clipper::Coord_orth vectorOrigin;

													vectorOrigin = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].coord_orth();
													
													for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
														{
															clipper::Coord_orth potentialTargetNH1 = getTargetPoint(NH1Coordinate, vectorOrigin, vectorShift);
															clipper::Coord_orth potentialTargetNH2 = getTargetPoint(NH1Coordinate, vectorOrigin, vectorShift);

															double meanDensityValueNH1 = calculateMeanElectronDensityInTargetPosition(potentialTargetNH1, sigmaa_dif_map, mapstats);
															double meanDensityValueNH2 = calculateMeanElectronDensityInTargetPosition(potentialTargetNH1, sigmaa_dif_map, mapstats);
															std::pair<clipper::Coord_orth, double> tempDensityInfoNH1(potentialTargetNH1, meanDensityValueNH1);
															std::pair<clipper::Coord_orth, double> tempDensityInfoNH2(potentialTargetNH2, meanDensityValueNH2);
															pairs.push_back(tempDensityInfoNH1);
															pairs.push_back(tempDensityInfoNH2);
														}
												}
										}
										const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
										clipper::Coord_orth bestTarget = bestPair->first;
										double bestDensityValue = bestPair->second;

										if(bestDensityValue > thresholdEDBestBlob)
										{
											double meanDensityValueBiggerArea = calculateMeanElectronDensityForBiggerSphere(bestTarget, sigmaa_dif_map, mapstats, hklinfo);
											if(meanDensityValueBiggerArea > thresholdED)
											{
												std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID, vectorIndex}, meanDensityValueBiggerArea);
												finalVectorForBlobValues.push_back(densityInfo);

													if(pdbexport)
													{
														drawOriginPoint(inputModel, bestTarget, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
														// fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, hklinfo, mapstats, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
													}
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
				std::list<clipper::String> ignoreAtomList = {" O  ", " N  ", " C  ", " CD1"};
				int vectorShiftLimit = 10;
				for (int c = 0; c < informationVector[vectorIndex].size(); c++)
				{
						if (inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() == "TRP") // in C-Glycosylation a glycan is attached through TRP residue
                            {	
								bool siteAlreadyGlycosylated = check_glycosylation_presence(inputModel[informationVector[vectorIndex][c].PolymerID].id(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id().trim(), glycanList);
								if(!siteAlreadyGlycosylated)
									{                                
										clipper::MAtom CD1Atom; // A glycan is attached to TRP residue via CD1 atom
										
										try {
										CD1Atom = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].find(" CD1", clipper::MM::ANY);
										} catch (const clipper::Message_fatal& error) {
										std::cerr << "Unable to find necessary CD1 atom for residue" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id() << "-" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() << " in Chain " << inputModel[informationVector[vectorIndex][c].PolymerID].id() << "\n" << "\n";
										}
										
										clipper::Coord_orth CD1Coordinate; // CD1 atom is used as a direction towards the glycan density
										CD1Coordinate = CD1Atom.coord_orth();

										std::vector<std::pair<clipper::Coord_orth, double>> pairs;

										for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].size(); natom++)
										{
											bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].id()) != ignoreAtomList.end());
											if(!atomIgnored)
												{
													clipper::Coord_orth vectorOrigin;

													vectorOrigin = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].coord_orth();
													
													for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
														{
															clipper::Coord_orth potentialTarget = getTargetPoint(CD1Coordinate, vectorOrigin, vectorShift);

															double meanDensityValue = calculateMeanElectronDensityInTargetPosition(potentialTarget, sigmaa_dif_map, mapstats);
															std::pair<clipper::Coord_orth, double> tempDensityInfo(potentialTarget, meanDensityValue);
															pairs.push_back(tempDensityInfo);
														}
												}
										}
										
										const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
										clipper::Coord_orth bestTarget = bestPair->first;
										double bestDensityValue = bestPair->second;

										if(bestDensityValue > thresholdEDBestBlob)
										{
											double meanDensityValueBiggerArea = calculateMeanElectronDensityForBiggerSphere(bestTarget, sigmaa_dif_map, mapstats, hklinfo);
											if(meanDensityValueBiggerArea > thresholdED)
											{
												std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID, vectorIndex}, meanDensityValueBiggerArea);
												finalVectorForBlobValues.push_back(densityInfo);
													if(pdbexport)
													{
														drawOriginPoint(inputModel, bestTarget, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
														// fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, hklinfo, mapstats, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
													}
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
						if (inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() == "SER") // in O-Glycosylation a glycan is attached through Thr or Ser residue
                            {	
								int vectorShiftLimit = 5; 
								bool siteAlreadyGlycosylated = check_glycosylation_presence(inputModel[informationVector[vectorIndex][c].PolymerID].id(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id().trim(), glycanList);
								if(!siteAlreadyGlycosylated)
									{
										std::list<clipper::String> ignoreAtomList = {" O  ", " OG "};
										clipper::MAtom OGAtom; // On Ser glycan attaches to OG, on Thr glycan attaches to OG1
										
										try {
										OGAtom = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].find(" OG ", clipper::MM::ANY);
										} catch (const clipper::Message_fatal& error) {
										std::cerr << "Unable to find necessary OG atom for residue" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id() << "-" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() << " in Chain " << inputModel[informationVector[vectorIndex][c].PolymerID].id() << "\n" << "\n";
										}
										

										clipper::Coord_orth OGCoordinate; // OG atom is used as a point of attachement of a glycan
										OGCoordinate = OGAtom.coord_orth();

										std::vector<std::pair<clipper::Coord_orth, double>> pairs;

										for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].size(); natom++)
										{
											bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].id()) != ignoreAtomList.end());
											if(!atomIgnored)
												{
													clipper::Coord_orth vectorOrigin;

													vectorOrigin = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].coord_orth();
													
													for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
														{
															clipper::Coord_orth potentialTarget = getTargetPoint(OGCoordinate, vectorOrigin, vectorShift);

															double meanDensityValue = calculateMeanElectronDensityInTargetPosition(potentialTarget, sigmaa_dif_map, mapstats);
															std::pair<clipper::Coord_orth, double> tempDensityInfo(potentialTarget, meanDensityValue);
															pairs.push_back(tempDensityInfo);
														}
												}
										}
										
										const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
										clipper::Coord_orth bestTarget = bestPair->first;
										double bestDensityValue = bestPair->second;

										if(bestDensityValue > thresholdEDBestBlob)
										{
											double meanDensityValueBiggerArea = calculateMeanElectronDensityForBiggerSphere(bestTarget, sigmaa_dif_map, mapstats, hklinfo);
											if(meanDensityValueBiggerArea > thresholdED)
											{
												std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID, vectorIndex}, meanDensityValueBiggerArea);
												finalVectorForBlobValues.push_back(densityInfo);
													if(pdbexport)
													{
														drawOriginPoint(inputModel, bestTarget, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
														// fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, hklinfo, mapstats, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
													}
											}
										}
									}
							}
						if (inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() == "THR")
							{
								int vectorShiftLimit = 5;  
								bool siteAlreadyGlycosylated = check_glycosylation_presence(inputModel[informationVector[vectorIndex][c].PolymerID].id(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id().trim(), glycanList);
								if(!siteAlreadyGlycosylated)
									{   
										std::list<clipper::String> ignoreAtomList = {" O  ", " N  ", " C  ", " OG1"};
										clipper::MAtom OG1Atom; // On Thr glycan attaches to OG1

										try {
										OG1Atom = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].find(" OG1", clipper::MM::ANY);
										} catch (const clipper::Message_fatal& error) {
										std::cerr << "Unable to find necessary OG1 atom for residue" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id() << "-" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() << " in Chain " << inputModel[informationVector[vectorIndex][c].PolymerID].id() << "\n" << "\n";
										}

										clipper::Coord_orth OG1Coordinate; // OG atom is used as a point of attachement of a glycan
										OG1Coordinate = OG1Atom.coord_orth();

										std::vector<std::pair<clipper::Coord_orth, double>> pairs;

										for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].size(); natom++)
										{
												bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].id()) != ignoreAtomList.end());
												if(!atomIgnored)
												{
													clipper::Coord_orth vectorOrigin;

													vectorOrigin = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].coord_orth();
													
													for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
														{
															clipper::Coord_orth potentialTarget = getTargetPoint(OG1Coordinate, vectorOrigin, vectorShift);

															double meanDensityValue = calculateMeanElectronDensityInTargetPosition(potentialTarget, sigmaa_dif_map, mapstats);
															std::pair<clipper::Coord_orth, double> tempDensityInfo(potentialTarget, meanDensityValue);
															pairs.push_back(tempDensityInfo);
														}
												}
										}
											
										const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
										clipper::Coord_orth bestTarget = bestPair->first;
										double bestDensityValue = bestPair->second;

										if(bestDensityValue > thresholdEDBestBlob)
										{
											double meanDensityValueBiggerArea = calculateMeanElectronDensityForBiggerSphere(bestTarget, sigmaa_dif_map, mapstats, hklinfo);
											if(meanDensityValueBiggerArea > thresholdED)
											{
												std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID, vectorIndex}, meanDensityValueBiggerArea);
												finalVectorForBlobValues.push_back(densityInfo);
													if(pdbexport)
													{
														drawOriginPoint(inputModel, bestTarget, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
														// fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, hklinfo, mapstats, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
													}
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
						if (inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() == "CYS") // in S-Glycosylation a glycan is attached through CYS or MET residue
                            {	
								int vectorShiftLimit = 3;  
								bool siteAlreadyGlycosylated = check_glycosylation_presence(inputModel[informationVector[vectorIndex][c].PolymerID].id(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id().trim(), glycanList);
								if(!siteAlreadyGlycosylated)
									{                								
										std::list<clipper::String> ignoreAtomList = {" O  ", " SG "};
										clipper::MAtom SGAtom; // On Ser glycan attaches to OG, on Thr glycan attaches to OG1

										try {
										SGAtom = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].find(" SG ", clipper::MM::ANY);
										} catch (const clipper::Message_fatal& error) {
										std::cerr << "Unable to find necessary SG atom for residue" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id() << "-" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() << " in Chain " << inputModel[informationVector[vectorIndex][c].PolymerID].id() << "\n" << "\n";
										}

										clipper::Coord_orth SGCoordinate; // SG atom is used as a point of attachement of a glycan
										SGCoordinate = SGAtom.coord_orth();

										std::vector<std::pair<clipper::Coord_orth, double>> pairs;

										for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].size(); natom++)
										{
											bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].id()) != ignoreAtomList.end());
											if(!atomIgnored)
												{
													clipper::Coord_orth vectorOrigin;

													vectorOrigin = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].coord_orth();
													
													for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
														{
															clipper::Coord_orth potentialTarget = getTargetPoint(SGCoordinate, vectorOrigin, vectorShift);

															double meanDensityValue = calculateMeanElectronDensityInTargetPosition(potentialTarget, sigmaa_dif_map, mapstats);
															std::pair<clipper::Coord_orth, double> tempDensityInfo(potentialTarget, meanDensityValue);
															pairs.push_back(tempDensityInfo);
														}
												}
										}
										
										const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
										clipper::Coord_orth bestTarget = bestPair->first;
										double bestDensityValue = bestPair->second;

										if(bestDensityValue > thresholdEDBestBlob)
										{
											double meanDensityValueBiggerArea = calculateMeanElectronDensityForBiggerSphere(bestTarget, sigmaa_dif_map, mapstats, hklinfo);
											if(meanDensityValueBiggerArea > thresholdED)
											{
												std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID, vectorIndex}, meanDensityValueBiggerArea);
												finalVectorForBlobValues.push_back(densityInfo);
													if(pdbexport)
													{
														drawOriginPoint(inputModel, bestTarget, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
														// fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, hklinfo, mapstats, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
													}
											}
										}
									}
							}
						if (inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() == "MET") // in S-Glycosylation a glycan is attached through CYS or MET residue
							{
								int vectorShiftLimit = 3;
								bool siteAlreadyGlycosylated = check_glycosylation_presence(inputModel[informationVector[vectorIndex][c].PolymerID].id(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id().trim(), glycanList);
								if(!siteAlreadyGlycosylated)
									{   
										               							
										std::list<clipper::String> ignoreAtomList = {" SD "};
										clipper::MAtom SDAtom; // On Thr glycan attaches to OG1
										
										try {
										SDAtom = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].find(" SD ", clipper::MM::ANY);
										} catch (const clipper::Message_fatal& error) {
										std::cerr << "Unable to find necessary SD atom for residue" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id() << "-" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() << " in Chain " << inputModel[informationVector[vectorIndex][c].PolymerID].id() << "\n" << "\n";
										}

										clipper::Coord_orth SDCoordinate; // OG atom is used as a point of attachement of a glycan
										SDCoordinate = SDAtom.coord_orth();

										std::vector<std::pair<clipper::Coord_orth, double>> pairs;

										for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].size(); natom++)
										{
												bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].id()) != ignoreAtomList.end());
												if(!atomIgnored)
												{
													clipper::Coord_orth vectorOrigin;

													vectorOrigin = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].coord_orth();
													
													for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
														{
															clipper::Coord_orth potentialTarget = getTargetPoint(SDCoordinate, vectorOrigin, vectorShift);

															double meanDensityValue = calculateMeanElectronDensityInTargetPosition(potentialTarget, sigmaa_dif_map, mapstats);
															std::pair<clipper::Coord_orth, double> tempDensityInfo(potentialTarget, meanDensityValue);
															pairs.push_back(tempDensityInfo);
														}
												}
										}
											
										const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
										clipper::Coord_orth bestTarget = bestPair->first;
										double bestDensityValue = bestPair->second;

										if(bestDensityValue > thresholdEDBestBlob)
										{
											double meanDensityValueBiggerArea = calculateMeanElectronDensityForBiggerSphere(bestTarget, sigmaa_dif_map, mapstats, hklinfo);
											if(meanDensityValueBiggerArea > thresholdED)
											{
												std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID, vectorIndex}, meanDensityValueBiggerArea);
												finalVectorForBlobValues.push_back(densityInfo);
													if(pdbexport)
													{
														drawOriginPoint(inputModel, bestTarget, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
														// fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, hklinfo, mapstats, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
													}
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
						if (inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() == "ALA") 
                            {
								int vectorShiftLimit = 5;
								bool siteAlreadyGlycosylated = check_glycosylation_presence(inputModel[informationVector[vectorIndex][c].PolymerID].id(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id().trim(), glycanList);
								if(!siteAlreadyGlycosylated)
									{                   								
										std::list<clipper::String> ignoreAtomList = {" CB "};
										clipper::MAtom CBAtom; // CA to CB for Ala, CG to NE2 for GLN

										try {
										CBAtom = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].find(" CB ", clipper::MM::ANY);
										} catch (const clipper::Message_fatal& error) {
										std::cerr << "Unable to find necessary CB atom for residue" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id() << "-" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() << " in Chain " << inputModel[informationVector[vectorIndex][c].PolymerID].id() << "\n" << "\n";
										}

										clipper::Coord_orth CBCoordinate; // CB atom is used as a direction towards the glycan density
										CBCoordinate = CBAtom.coord_orth();

										std::vector<std::pair<clipper::Coord_orth, double>> pairs;

										for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].size(); natom++)
										{
											
											bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].id()) != ignoreAtomList.end());
											if(!atomIgnored)
												{
													clipper::Coord_orth vectorOrigin;

													vectorOrigin = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].coord_orth();
													
													for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
														{
															clipper::Coord_orth potentialTarget = getTargetPoint(CBCoordinate, vectorOrigin, vectorShift);

															double meanDensityValue = calculateMeanElectronDensityInTargetPosition(potentialTarget, sigmaa_dif_map, mapstats);
															std::pair<clipper::Coord_orth, double> tempDensityInfo(potentialTarget, meanDensityValue);
															pairs.push_back(tempDensityInfo);
														}
												}
										}
										const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
										clipper::Coord_orth bestTarget = bestPair->first;
										double bestDensityValue = bestPair->second;

										if(bestDensityValue > thresholdEDBestBlob)
										{
											double meanDensityValueBiggerArea = calculateMeanElectronDensityForBiggerSphere(bestTarget, sigmaa_dif_map, mapstats, hklinfo);
											if(meanDensityValueBiggerArea > thresholdED) 
											{
												std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID, vectorIndex}, meanDensityValueBiggerArea);
												finalVectorForBlobValues.push_back(densityInfo);
													if(pdbexport)
													{
														drawOriginPoint(inputModel, bestTarget, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
														// fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, hklinfo, mapstats, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
													}
											}
										}
									}	
							}
						if (inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() == "GLN") 
                            {	
								int vectorShiftLimit = 5;
								bool siteAlreadyGlycosylated = check_glycosylation_presence(inputModel[informationVector[vectorIndex][c].PolymerID].id(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id().trim(), glycanList);
								if(!siteAlreadyGlycosylated)
									{                   								
										std::list<clipper::String> ignoreAtomList = {" CG "};
										
										clipper::MAtom CGAtom; // CA to CB for Ala, CG to NE2 for GLN

										try {
										CGAtom = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].find(" CG ", clipper::MM::ANY);
										} catch (const clipper::Message_fatal& error) {
										std::cerr << "Unable to find necessary CG atom for residue" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].id() << "-" << inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].type() << " in Chain " << inputModel[informationVector[vectorIndex][c].PolymerID].id() << "\n" << "\n";
										}

										
										clipper::Coord_orth CGCoordinate; // CB atom is used as a direction towards the glycan density
										CGCoordinate = CGAtom.coord_orth();

										std::vector<std::pair<clipper::Coord_orth, double>> pairs;

										for (int natom = 0; natom < inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID].size(); natom++)
										{
											
											bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].id()) != ignoreAtomList.end());
											if(!atomIgnored)
												{
													clipper::Coord_orth vectorOrigin;

													vectorOrigin = inputModel[informationVector[vectorIndex][c].PolymerID][informationVector[vectorIndex][c].ResidueID][natom].coord_orth();
													
													for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
														{
															clipper::Coord_orth potentialTarget = getTargetPoint(vectorOrigin, CGCoordinate, vectorShift);

															double meanDensityValue = calculateMeanElectronDensityInTargetPosition(potentialTarget, sigmaa_dif_map, mapstats);
															std::pair<clipper::Coord_orth, double> tempDensityInfo(potentialTarget, meanDensityValue);
															pairs.push_back(tempDensityInfo);
														}
												}
										}
										const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
										clipper::Coord_orth bestTarget = bestPair->first;
										double bestDensityValue = bestPair->second;

										if(bestDensityValue > thresholdEDBestBlob)
										{
											double meanDensityValueBiggerArea = calculateMeanElectronDensityForBiggerSphere(bestTarget, sigmaa_dif_map, mapstats, hklinfo);
											if(meanDensityValueBiggerArea > thresholdED)
											{
												std::pair<PotentialGlycosylationSiteInfo, double> densityInfo(PotentialGlycosylationSiteInfo{informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID, vectorIndex}, meanDensityValueBiggerArea);
												finalVectorForBlobValues.push_back(densityInfo);
													if(pdbexport)
													{
														drawOriginPoint(inputModel, bestTarget, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
														// fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, hklinfo, mapstats, informationVector[vectorIndex][c].PolymerID, informationVector[vectorIndex][c].ResidueID);
													}
											}
										}
									}
							}
				}
			}
		}
	return finalVectorForBlobValues;
}

std::vector<std::pair<GlycanToMiniMolIDs, double> > get_electron_density_of_potential_unmodelled_carbohydrate_monomers(std::vector < clipper::MSugar > glycanChain, clipper::MiniMol&inputModel, std::vector < clipper::MGlycan >& allSugars, int id, clipper::Xmap<float>& sigmaa_dif_map, clipper::HKL_info& hklinfo, clipper::Map_stats& mapstats, float thresholdED, bool pdbexport)
{
	float thresholdEDBestBlob = 0.070;
	std::vector<std::pair<GlycanToMiniMolIDs, double> > finalVectorForBlobValues;
	int vectorShiftLimit = 5;
	for (int monomer = 0; monomer < glycanChain.size(); monomer++)
	{
		std::vector<clipper::String> ignoreAtomList = create_list_of_ignored_sugar_atoms(glycanChain[monomer]);

		for (int atom = 0; atom < glycanChain[monomer].size(); atom++)
		{
			std::vector<std::pair<clipper::Coord_orth, double>> pairs;
			clipper::Coord_orth sugarCentre = glycanChain[monomer].ring_centre();
			bool atomIgnored = (std::find(ignoreAtomList.begin(), ignoreAtomList.end(), glycanChain[monomer][atom].id()) != ignoreAtomList.end());
			if (!atomIgnored) {
				clipper::Coord_orth linkageAtomLocation = glycanChain[monomer][atom].coord_orth();


				for(int vectorShift = 1; vectorShift <= vectorShiftLimit; vectorShift++)
					{
						clipper::Coord_orth potentialTarget = getTargetPoint(linkageAtomLocation, sugarCentre, vectorShift);

						double meanDensityValue = calculateMeanElectronDensityInTargetPosition(potentialTarget, sigmaa_dif_map, mapstats);
						std::pair<clipper::Coord_orth, double> tempDensityInfo(potentialTarget, meanDensityValue);
						pairs.push_back(tempDensityInfo);
					}
			
			const auto bestPair = max_element(pairs.begin(), pairs.end(), bestPointFinder);
			clipper::Coord_orth bestTarget = bestPair->first;
			double bestDensityValue = bestPair->second;

			if(bestDensityValue > thresholdEDBestBlob)
			{
				GlycanToMiniMolIDs identification = getCarbohydrateRelationshipToMiniMol(inputModel, glycanChain[monomer], allSugars, id, monomer);
				double meanDensityValueBiggerArea = calculateMeanElectronDensityForBiggerSphere(bestTarget, sigmaa_dif_map, mapstats, hklinfo);
				if(meanDensityValueBiggerArea > thresholdED)
				{
					std::pair<GlycanToMiniMolIDs, double> densityInfo(GlycanToMiniMolIDs{identification.proteinMiniMolID, identification.carbohydrateChainMiniMolID, identification.carbohydrateID}, meanDensityValueBiggerArea);
					finalVectorForBlobValues.push_back(densityInfo);
						if(pdbexport) 
						{
							drawOriginPoint(inputModel, bestTarget, identification.proteinMiniMolID, identification.carbohydrateChainMiniMolID);
							// fillSearchArea(inputModel, bestTarget, sigmaa_dif_map, mapstats, identification.proteinMiniMolID, identification.carbohydrateChainMiniMolID);
						}
				}

			}
			}
		}
	}
	return finalVectorForBlobValues;
}
	




