/*! \file clipper-glyco.cpp
  Implementation file for sugar data */

// clipper-glyco.cc: a set of tools for handling sugars
// version  0.6.0
// 2013 Jon Agirre & Kevin Cowtan, The University of York
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk
//
// Compile with 'make DEBUG=-DDUMP' to enable debugging messages
//
//L  This library is free software and is distributed under the terms
//L  and conditions of version 2.1 of the GNU Lesser General Public
//L  Licence (LGPL) with the following additional clause:
//L
//L     `You may also combine or link a "work that uses the Library" to
//L     produce a work containing portions of the Library, and distribute
//L     that work under terms of your choice, provided that you give
//L     prominent notice with each copy of the work that the specified
//L     version of the Library is used in it, and that you include or
//L     provide public access to the complete corresponding
//L     machine-readable source code for the Library including whatever
//L     changes were used in the work. (i.e. If you make changes to the
//L     Library you must distribute those, but you do not need to
//L     distribute source or object code to those portions of the work
//L     not covered by this licence.)'
//L
//L  Note that this clause grants an additional right and does not impose
//L  any additional restriction, and so does not affect compatibility
//L  with the GNU General Public Licence (GPL). If you wish to negotiate
//L  other terms, please contact the maintainer.
//L
//L  You can redistribute it and/or modify the library under the terms of
//L  the GNU Lesser General Public License as published by the Free Software
//L  Foundation; either version 2.1 of the License, or (at your option) any
//L  later version.
//L
//L  This library is distributed in the hope that it will be useful, but
//L  WITHOUT ANY WARRANTY; without even the implied warranty of
//L  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//L  Lesser General Public License for more details.
//L
//L  You should have received a copy of the CCP4 licence and/or GNU
//L  Lesser General Public License along with this library; if not, write
//L  to the CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
//L  The GNU Lesser General Public can also be obtained by writing to the
//L  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//L  MA 02111-1307 USA


#include <fstream>
#include <algorithm>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <clipper/minimol/minimol_utils.h>
#include "clipper-glyco.h"

#define DBG std::cout << "[" << __FUNCTION__ << "] - "


using namespace clipper;


/*! Constructor: Empty constructor for later initialisation. */

MSugar::MSugar( )
{

}

/*! Constructor: create a new sugar object from a standard MMonomer.
	  If reference data for the sugar cannot be found in the database, the members of the ring will be determined using a recursive version of Fleury's algorithm for finding eulerian cycles in undirected graphs
	\param ml A MiniMol object containing this and neighbouring sugars, which may affect the stereochemistry
	\param mm A MMonomer object that will be extended into a sugar
	\return The MSugar object, which will contain cremer-pople parameters, conformation code, anomer, handedness and linkage information */

MSugar::MSugar(const clipper::MiniMol& ml, const clipper::MMonomer& mm) 
{

	// we calculate the non-bond object first, then continue with normal creation
	const clipper::MAtomNonBond& nb = MAtomNonBond (ml, 5.0);
	MSugar(ml, mm, nb);	
}


/*! Constructor: create a new sugar object from a standard MMonomer.
	If reference data for the sugar cannot be found in the database, the members of the ring will be determined using a recursive version of Fleury's algorithm for finding eulerian cycles in undirected graphs
	\param ml A MiniMol object containing this and neighbouring sugars
	\param mm An MMonomer object that will be extended into a sugar
	\param nb An MAtomNonBond object to be used for the determination of the stereochemistry
	\return The MSugar object, which will contain cremer-pople parameters, conformation code, anomer, handedness and linkage information */

MSugar::MSugar(const clipper::MiniMol& ml, const clipper::MMonomer& mm, const clipper::MAtomNonBond& nb)
{	

    copy(mm,clipper::MM::COPY_MPC);	// import_data from MMonomer

    this->sugar_supported = true; 
    this->sugar_parent_molecule = &ml;
    this->sugar_parent_molecule_nonbond = &nb; // store pointers
    this->sugar_index = db_not_checked;
    this->sugar_index = 9999; // default value for "not found in database".
    this->sugar_alternate_confcode = " "; // initially, we would like to suppose this	
    this->sugar_context = "";
    
    #ifdef DUMP
        DBG << std::endl << std::endl << "looking for " << this->id() << " " << this->type().trim() << " on the database..." << std::endl;
    #endif

    this->sugar_found_db = lookup_database(this->type().trim());

    sugar_bfactor = 0.0;

    for (int i=0; i < this->size(); i++)
    {
        MSugar mstmp= *this;
        sugar_bfactor += mstmp[i].u_iso();
    }

    sugar_bfactor /= this->size();
    sugar_bfactor = clipper::Util::u2b(sugar_bfactor);

    if ( this->sugar_found_db )
    {
        #ifdef DUMP
            DBG << "found it! " << std::endl;
        #endif

	std::vector<clipper::String> buffer = clipper::data::sugar_database[sugar_index].ring_atoms.trim().split(" ");
	for (int i=0 ; i < buffer.size() ; i++) 
	{
	    int index_atom = 0;
	    index_atom = this->lookup(buffer[i],clipper::MM::ANY);
	    if (index_atom == -1) 
	    {
		this->sugar_supported = false;
		this->sugar_sane = false; // we don't support cyclic sugars with less or more than 5-6 ring atoms
		this->sugar_denomination = "    unsupported    ";
		this->sugar_anomer = "X";
		this->sugar_handedness = "X";
		return;
	    }
	    else if ((*this)[index_atom].occupancy() < 1.0) // we will need to grab conformation A
	    {
		if ( get_altconf((*this)[index_atom]) != ' ' )
		{ 
		    index_atom = this->lookup(buffer[i].trim()+" :A",clipper::MM::UNIQUE);
				
		    if (index_atom == -1) // same case as above
		        index_atom = this->lookup(buffer[i].trim()+" :B",clipper::MM::UNIQUE); // try grabbing conformation B instead 
		    else 
                        sugar_alternate_confcode = " :A";
	
		    if (index_atom == -1) // we've tried A and B and it still fails... so we're going to give up for now 
		    {
			this->sugar_supported = false;
			this->sugar_sane = false;
			this->sugar_denomination = "    unsupported    ";
			this->sugar_anomer = "X";
		        this->sugar_handedness = "X";
			return;
		    }
		    else sugar_alternate_confcode = " :B";
		}
 	    }
	    sugar_ring_elements.push_back((*this)[index_atom]);
	}
    }
    else
    {
	this->sugar_ring_elements = this->ringMembers();
    }

	///#ifdef DUMP
	//	DBG << "Ring members successfully determined." << std::endl;
	//#endif

	if (this->sugar_ring_elements.size() == 5)
	{
		this->cremerPople_furanose(*this->sugar_parent_molecule, mm);	
		this->sugar_conformation = conformationFuranose(this->sugar_cremer_pople_params[1]);
		
		#ifdef DUMP
			DBG << "After checking the conformation..." << std::endl;
		#endif

	}
	else if (this->sugar_ring_elements.size() == 6)
	{
		this->cremerPople_pyranose(*this->sugar_parent_molecule, mm);
		this->sugar_conformation = conformationPyranose(this->sugar_cremer_pople_params[1], this->sugar_cremer_pople_params[2]);
	}
	else
	{
		this->sugar_supported = false;
		this->sugar_sane = false; // we don't support cyclic sugars with less or more than 5-6 ring atoms
		this->sugar_denomination = "    unsupported    ";
		this->sugar_anomer = "X";
		this->sugar_handedness = "X";
	}

	if (sugar_ring_elements[1].name().trim().find("C1") != std::string::npos) this->sugar_denomination = "aldo";
	else this->sugar_denomination = "keto";

	if (sugar_ring_elements.size() == 5) this->sugar_denomination = this->sugar_denomination + "furanose";
	else this->sugar_denomination = this->sugar_denomination + "pyranose";

	this->sugar_denomination = clipper::String( this->sugar_anomer + "-" + this->sugar_handedness + "-" + this->sugar_denomination );

	// sanity check:

	this->sugar_sane = false;

	#ifdef DUMP
		DBG << "Just before examining the ring..." << std::endl;
	#endif

	
	if ( examine_ring() ) sugar_diag_ring = true; else sugar_diag_ring = false;

        clipper::String ref_conformation;
        clipper::ftype  ref_puckering;   
        clipper::ftype  ref_bonds_rmsd;
        clipper::ftype  ref_angles_rmsd;

	if (this->sugar_found_db)
	{

            ref_conformation = clipper::data::sugar_database[sugar_index].ref_conformation;
            ref_puckering    = clipper::data::sugar_database[sugar_index].ref_puckering;
            ref_bonds_rmsd   = clipper::data::sugar_database[sugar_index].ref_bonds_rmsd;
            ref_angles_rmsd  = clipper::data::sugar_database[sugar_index].ref_angles_rmsd;

	    if ( ( ( sugar_handedness != "D" ) && ( clipper::data::sugar_database[sugar_index].handedness != "D" ) ) 
		     || ( ( sugar_handedness != "L" ) && (clipper::data::sugar_database[sugar_index].handedness != "L" ) ) )
	        sugar_diag_chirality = true; 
	    else sugar_diag_chirality = false;


	    if ( ( ( sugar_anomer == "alpha") && ( clipper::data::sugar_database[sugar_index].anomer != "B" ) )
		     || ( ( sugar_anomer == "beta") && ( clipper::data::sugar_database[sugar_index].anomer != "A" ) ) )
	        sugar_diag_anomer = true;
	    else sugar_diag_anomer = false;

            if ( ref_conformation == conformation_name() )
                sugar_diag_conformation = true;
            else 
            {
                if ( ( conformation_name() == "4c1" ) && ( sugar_handedness != "L" ))
                    sugar_diag_conformation = true;
                else if ( ( conformation_name() == "1c4" ) && ( sugar_handedness != "D" ))
                    sugar_diag_conformation = true;
                else if ( ring_cardinality() < 6) 
                   sugar_diag_conformation = true;
                else  sugar_diag_conformation = false;
            }

            if ( sugar_diag_conformation )
            {    
                if (( puckering_amplitude() > ref_puckering - 0.15 ) && (puckering_amplitude() < ref_puckering + 0.15))
                    sugar_diag_puckering = true;
                else sugar_diag_puckering = false;
            
                if ( sugar_ring_bond_rmsd < ( ref_bonds_rmsd + 0.039 ) ) 
                    sugar_diag_bonds_rmsd = true; 
	        else sugar_diag_bonds_rmsd = false; 

	        if ( sugar_ring_angle_rmsd < ( ref_angles_rmsd + 3.0 ) )
	            sugar_diag_angles_rmsd = true; 
	        else sugar_diag_angles_rmsd = false; 
            }
            else 
            {
                if (( puckering_amplitude() > 0.9 ) || ( puckering_amplitude() < 0.4 ))
                    sugar_diag_puckering = false;
                else sugar_diag_puckering = true;

                sugar_diag_bonds_rmsd = sugar_diag_angles_rmsd = true;
            }

	    if ( sugar_diag_puckering && sugar_diag_anomer && sugar_diag_chirality && sugar_diag_ring )
            sugar_sane = true;
	}
	else
	{
	    sugar_diag_anomer=true;
	    sugar_diag_chirality=true;  // perform a generic test based on rough ideal values

	    if (sugar_ring_elements.size() == 5)
	    {
                if (sugar_ring_bond_rmsd < 0.040 ) 
                    sugar_diag_bonds_rmsd = true; 
		else sugar_diag_bonds_rmsd = false; 
	    }
	    else 
	    {		
	        if (sugar_ring_bond_rmsd < 0.035 ) 
		    sugar_diag_bonds_rmsd = true; 
		else sugar_diag_bonds_rmsd = false; 
	    }

	    if (sugar_ring_elements.size() ==  5)
	    {
	        if ((sugar_ring_angle_rmsd > 4.0 ) && (sugar_ring_angle_rmsd < 8.0)) 
		    sugar_diag_angles_rmsd = true; 
		else sugar_diag_angles_rmsd = false;
	    }
	    else 
	    {
	        if (sugar_ring_angle_rmsd < 4.0 ) 
		    sugar_diag_angles_rmsd = true; 
	        else sugar_diag_angles_rmsd = false;
	    }

            if ( ( conformation_name() == "4c1" ) && ( sugar_handedness != "L" ))
                sugar_diag_conformation = true;
            else if ( ( conformation_name() == "1c4" ) && ( sugar_handedness != "D" ))
                sugar_diag_conformation = true;
            else sugar_diag_conformation = false;

            if (( puckering_amplitude() > 0.9 ) || ( puckering_amplitude() < 0.4 ))
                sugar_diag_puckering = false;
            else sugar_diag_puckering = true;
            
	    if ( sugar_diag_puckering && sugar_diag_anomer && sugar_diag_chirality && sugar_diag_ring )
            sugar_sane = true;

	}

	#ifdef DUMP
	    DBG << "Just after examining the ring, exiting the constructor, good job!" << std::endl;
	#endif
}




/*! Constructor: create a new sugar object from a standard clipper::MMonomer, providing external validation data.
	If reference data for the sugar cannot be found in the database, the members of the ring will be determined using a recursive version of Fleury's algorithm for finding eulerian cycles in undirected graphs
	\param ml A MiniMol object containing this and neighbouring sugars
	\param mm An MMonomer object that will be extended into a sugar
	\param nb An MAtomNonBond object to be used for the determination of the stereochemistry
	\return The MSugar object, which will contain cremer-pople parameters, conformation code, anomer, handedness and linkage information */

MSugar::MSugar(const clipper::MiniMol& ml, const clipper::MMonomer& mm, const clipper::MAtomNonBond& nb, clipper::data::sugar_database_entry& validation_data)
{	
    copy(mm,clipper::MM::COPY_MPC);	// import_data from MMonomer

    this->sugar_supported = true; 
    this->sugar_parent_molecule = &ml;
    this->sugar_parent_molecule_nonbond = &nb; // store pointers
    this->sugar_index = db_not_checked;
    this->sugar_index = 9999; // default value for "not found in database".
    this->sugar_alternate_confcode = " "; // initially, we would like to suppose this	
    this->sugar_context = "";
    
    #ifdef DUMP
        DBG << std::endl << std::endl << "looking for " << this->id() << " " << this->type().trim() << " on the database..." << std::endl;
    #endif

    this->sugar_found_db = true;

    sugar_bfactor = 0.0;

    for (int i=0; i < this->size(); i++)
    {
        MSugar mstmp= *this;
        sugar_bfactor += mstmp[i].u_iso();
    }

    sugar_bfactor /= this->size();
    sugar_bfactor = clipper::Util::u2b(sugar_bfactor);

    #ifdef DUMP
	DBG << "found it! " << std::endl;
    #endif

    std::vector <clipper::String> buffer = validation_data.ring_atoms.trim().split(" ");
    for (int i=0 ; i < buffer.size() ; i++) 
    {
	int index_atom = 0;
	index_atom = this->lookup(buffer[i],clipper::MM::ANY);
	if (index_atom == -1) 
	{
	    this->sugar_supported = false;
	    this->sugar_sane = false; // we don't support cyclic sugars with less or more than 5-6 ring atoms
	    this->sugar_denomination = "    unsupported    ";
	    this->sugar_anomer = "X";
	    this->sugar_handedness = "X";
	    return;
	}
	else if ((*this)[index_atom].occupancy() < 1.0) // we will need to grab conformation A
	{
	    if ( get_altconf((*this)[index_atom]) != ' ' )
	    { 
		index_atom = this->lookup(buffer[i].trim()+" :A",clipper::MM::UNIQUE);
				
		if (index_atom == -1) // same case as above
		    index_atom = this->lookup(buffer[i].trim()+" :B",clipper::MM::UNIQUE); // try grabbing conformation B instead 
		else 
                    sugar_alternate_confcode = " :A";
	
		if (index_atom == -1) // we've tried A and B and it still fails... so we're going to give up for now 
		{
		    this->sugar_supported = false;
		    this->sugar_sane = false;
		    this->sugar_denomination = "    unsupported    ";
		    this->sugar_anomer = "X";
		    this->sugar_handedness = "X";
		    return;
		}
		else sugar_alternate_confcode = " :B";
	    }
 	}
	sugar_ring_elements.push_back((*this)[index_atom]);
    }

    if (this->sugar_ring_elements.size() == 5)
    {
	this->cremerPople_furanose(*this->sugar_parent_molecule, mm);	
	this->sugar_conformation = conformationFuranose(this->sugar_cremer_pople_params[1]);
		
	#ifdef DUMP
	    DBG << "After checking the conformation..." << std::endl;
	#endif

    }
    else if (this->sugar_ring_elements.size() == 6)
    {
	this->cremerPople_pyranose(*this->sugar_parent_molecule, mm);
	this->sugar_conformation = conformationPyranose(this->sugar_cremer_pople_params[1], this->sugar_cremer_pople_params[2]);
    }
    else
    {
	this->sugar_supported = false;
	this->sugar_sane = false; // we don't support cyclic sugars with less or more than 5-6 ring atoms
	this->sugar_denomination = "    unsupported    ";
	this->sugar_anomer = "X";
	this->sugar_handedness = "X";
    }

    if (sugar_ring_elements[1].name().trim().find("C1") != std::string::npos) this->sugar_denomination = "aldo";
    else this->sugar_denomination = "keto";

    if (sugar_ring_elements.size() == 5) this->sugar_denomination = this->sugar_denomination + "furanose";
    else this->sugar_denomination = this->sugar_denomination + "pyranose";

    this->sugar_denomination = clipper::String( this->sugar_anomer + "-" + this->sugar_handedness + "-" + this->sugar_denomination );

    // sanity check:

    this->sugar_sane = false;

    #ifdef DUMP
	DBG << "Just before examining the ring..." << std::endl;
    #endif

	
    if ( examine_ring() ) sugar_diag_ring = true; else sugar_diag_ring = false;

    clipper::String ref_conformation;
    clipper::ftype  ref_puckering;   
    clipper::ftype  ref_bonds_rmsd;
    clipper::ftype  ref_angles_rmsd;


    ref_conformation = validation_data.ref_conformation;
    ref_puckering    = validation_data.ref_puckering;
    ref_bonds_rmsd   = validation_data.ref_bonds_rmsd;
    ref_angles_rmsd  = validation_data.ref_angles_rmsd;

    if ( ( ( sugar_handedness != "D" ) && ( validation_data.handedness != "D" ) ) 
	    || ( ( sugar_handedness != "L" ) && ( validation_data.handedness != "L" ) ) )
	sugar_diag_chirality = true; 
    else 
        sugar_diag_chirality = false;


    if ( ( ( sugar_anomer == "alpha") && ( validation_data.anomer != "B" ) )
	    || ( ( sugar_anomer == "beta") && ( validation_data.anomer != "A" ) ) )
	sugar_diag_anomer = true;
    else 
        sugar_diag_anomer = false;

    if ( ref_conformation == conformation_name() )
        sugar_diag_conformation = true;
    else 
    {
        if ( ( conformation_name() == "4c1" ) && ( sugar_handedness != "L" ))
            sugar_diag_conformation = true;
        else if ( ( conformation_name() == "1c4" ) && ( sugar_handedness != "D" ))
            sugar_diag_conformation = true;
        else if ( ring_cardinality() < 6) 
            sugar_diag_conformation = true;
        else  sugar_diag_conformation = false;
    }

    if ( sugar_diag_conformation )
    {    
        if (( puckering_amplitude() > ref_puckering - 0.15 ) && (puckering_amplitude() < ref_puckering + 0.15))
            sugar_diag_puckering = true;
        else sugar_diag_puckering = false;
            
        if ( sugar_ring_bond_rmsd < ( ref_bonds_rmsd + 0.039 ) ) 
            sugar_diag_bonds_rmsd = true; 
	else sugar_diag_bonds_rmsd = false; 

	if ( sugar_ring_angle_rmsd < ( ref_angles_rmsd + 3.0 ) )
	    sugar_diag_angles_rmsd = true; 
	else sugar_diag_angles_rmsd = false; 
    }
    else 
    {
        if (( puckering_amplitude() > 0.9 ) || ( puckering_amplitude() < 0.4 ))
            sugar_diag_puckering = false;
        else sugar_diag_puckering = true;

        sugar_diag_bonds_rmsd = sugar_diag_angles_rmsd = true;
    }

    if ( sugar_diag_puckering && sugar_diag_anomer && sugar_diag_chirality && sugar_diag_ring )
        sugar_sane = true;

    #ifdef DUMP
	DBG << "Just after examining the ring, exiting the constructor, good job!" << std::endl;
    #endif

}









/*! Checks if the sugar is in the database of sugars. If found, it stores
		its index and returns true.
	\param name Three-letter code for the sugar
	\return True if found, false otherwise */

bool MSugar::lookup_database(clipper::String name)
{
	for (int i = 0; i < clipper::data::sugar_database_size ; i++)
	{
		if (name.trim() == clipper::data::sugar_database[i].name_short.trim())
		{
			this->sugar_index = i;
			this->sugar_found_db = true;
			return true;
		}
	}

	this->sugar_index = db_not_found;
	this->sugar_found_db = false;
	return false;

}



/*! Internal function, not for public use
		\param mmol The parent MiniMol object
		\return A vector of real values containing the cremer-pople parameters
*/

std::vector<clipper::ftype> MSugar::cremerPople_pyranose(const clipper::MiniMol& mmol, clipper::MMonomer mm) // modifies object, uses copy of mm
{
    clipper::ftype centre_x, centre_y, centre_z;
    centre_x = centre_y = centre_z = 0.0;

	bool lurd_reverse = false;	// to be used when the configurational carbon is lower-ranked than the carbon making the link,
								// which is involved in a configurational switch that effectively reverses the LURD mnemonic

	std::vector<clipper::MAtom> ring_atoms = ring_members(); // this has already been determined in the constructor

	clipper::String nz1 = ring_atoms[0].name().trim(); // In-ring oxygen
    clipper::String nz2 = ring_atoms[1].name().trim(); // Anomeric carbon
    clipper::String nz3 = ring_atoms[2].name().trim(); // The rest of the in-ring carbons...
    clipper::String nz4 = ring_atoms[3].name().trim();
    clipper::String nz5 = ring_atoms[4].name().trim();
    clipper::String nz6 = ring_atoms[5].name().trim(); // End of ring

	#ifdef DUMP
		DBG << "getting the stereochemistry..." << std::endl;
	#endif

	stereochemistry_pairs stereo = get_stereochemistry(mmol);

	#ifdef DUMP
		DBG << "done." << std::endl;
	#endif

	this->sugar_anomeric_carbon = stereo.first.first;
	this->sugar_anomeric_substituent = stereo.first.second;
	this->sugar_configurational_carbon = stereo.second.first;
	this->sugar_configurational_substituent = stereo.second.second;

	clipper::String anomeric_carbon = stereo.first.first.name().trim();
	clipper::String anomeric_substituent = stereo.first.second.name().trim();
	clipper::String configurational_carbon = stereo.second.first.name().trim();
	clipper::String configurational_substituent = stereo.second.second.name().trim();

	if ( configurational_carbon != "XXX" )
		if ( ( ring_atoms[5].name().trim() == configurational_carbon) || !is_part_of_ring( sugar_configurational_carbon, ring_atoms ) ) // we have to take into account the rotation performed prior to closing the ring
			lurd_reverse = true;

    // find geometrical centre

    centre_x += (*this)[this->lookup(nz1,clipper::MM::ANY)].coord_orth().x();
    centre_y += (*this)[this->lookup(nz1,clipper::MM::ANY)].coord_orth().y();
    centre_z += (*this)[this->lookup(nz1,clipper::MM::ANY)].coord_orth().z();

    centre_x += (*this)[this->lookup(nz2,clipper::MM::ANY)].coord_orth().x();
    centre_y += (*this)[this->lookup(nz2,clipper::MM::ANY)].coord_orth().y();
    centre_z += (*this)[this->lookup(nz2,clipper::MM::ANY)].coord_orth().z();

    centre_x += (*this)[this->lookup(nz3,clipper::MM::ANY)].coord_orth().x();
    centre_y += (*this)[this->lookup(nz3,clipper::MM::ANY)].coord_orth().y();
    centre_z += (*this)[this->lookup(nz3,clipper::MM::ANY)].coord_orth().z();

    centre_x += (*this)[this->lookup(nz4,clipper::MM::ANY)].coord_orth().x();
    centre_y += (*this)[this->lookup(nz4,clipper::MM::ANY)].coord_orth().y();
    centre_z += (*this)[this->lookup(nz4,clipper::MM::ANY)].coord_orth().z();

    centre_x += (*this)[this->lookup(nz5,clipper::MM::ANY)].coord_orth().x();
    centre_y += (*this)[this->lookup(nz5,clipper::MM::ANY)].coord_orth().y();
    centre_z += (*this)[this->lookup(nz5,clipper::MM::ANY)].coord_orth().z();

    centre_x += (*this)[this->lookup(nz6,clipper::MM::ANY)].coord_orth().x();
    centre_y += (*this)[this->lookup(nz6,clipper::MM::ANY)].coord_orth().y();
    centre_z += (*this)[this->lookup(nz6,clipper::MM::ANY)].coord_orth().z();

    centre_x = centre_x / 6;
    centre_y = centre_y / 6;
    centre_z = centre_z / 6;

    clipper::Coord_orth centre(centre_x,centre_y,centre_z);
    this->sugar_centre = centre;

    //#ifdef DUMP
	//		DBG << "Ring centre: " << centre.format() << std::endl;
	//#endif

    clipper::RTop_orth shift(clipper::Mat33<>::identity(), (- centre)); //clipper::Vec3<>::null()

    mm.transform(shift); // recentre the sugar coordinates

	stereo.first.first.transform(shift);
	stereo.first.second.transform(shift);
	stereo.second.first.transform(shift);
	stereo.second.second.transform(shift);

    //#ifdef DUMP
	//		DBG << shift.format() << std::endl;
	//#endif

    clipper::Vec3<clipper::ftype> rPrime(0.0, 0.0, 0.0);
    clipper::Vec3<clipper::ftype> r2Prime(0.0, 0.0, 0.0);
    clipper::Vec3<clipper::ftype> n(0.0, 0.0, 0.0);

    double argument = 0.0; // (j-1) = 0 for the first case
    rPrime +=  mm[mm.lookup(nz1,clipper::MM::ANY)].coord_orth() * sin(argument);
    r2Prime += mm[mm.lookup(nz1,clipper::MM::ANY)].coord_orth() * cos(argument);

    argument = (2.0 * clipper::Util::pi())/6.0;
    rPrime +=  mm[mm.lookup(nz2,clipper::MM::ANY)].coord_orth() * sin(argument);
    r2Prime += mm[mm.lookup(nz2,clipper::MM::ANY)].coord_orth() * cos(argument);

    argument = (2.0 * clipper::Util::pi() * 2)/6.0;
    rPrime +=  mm[mm.lookup(nz3,clipper::MM::ANY)].coord_orth() * sin(argument);
    r2Prime += mm[mm.lookup(nz3,clipper::MM::ANY)].coord_orth() * cos(argument);

    argument = (2.0 * clipper::Util::pi() * 3)/6.0;
    rPrime +=  mm[mm.lookup(nz4,clipper::MM::ANY)].coord_orth() * sin(argument);
    r2Prime += mm[mm.lookup(nz4,clipper::MM::ANY)].coord_orth() * cos(argument);

    argument = (2.0 * clipper::Util::pi() * 4)/6.0;
    rPrime +=  mm[mm.lookup(nz5,clipper::MM::ANY)].coord_orth() * sin(argument);
    r2Prime += mm[mm.lookup(nz5,clipper::MM::ANY)].coord_orth() * cos(argument);

    argument = (2.0 * clipper::Util::pi() * 5)/6.0;
    rPrime +=  mm[mm.lookup(nz6,clipper::MM::ANY)].coord_orth() * sin(argument);
    r2Prime += mm[mm.lookup(nz6,clipper::MM::ANY)].coord_orth() * cos(argument);

    n = (clipper::Vec3<clipper::ftype>::cross(rPrime, r2Prime)).unit();

    clipper::ftype z1, z2, z3, z4, z5, z6, z_anomeric_carbon, z_anomeric_substituent, z_configurational_carbon, z_configurational_substituent;

    z1 = clipper::Vec3<clipper::ftype>::dot(mm[mm.lookup(nz1,clipper::MM::ANY)].coord_orth(), n);
    z2 = clipper::Vec3<clipper::ftype>::dot(mm[mm.lookup(nz2,clipper::MM::ANY)].coord_orth(), n);
    z3 = clipper::Vec3<clipper::ftype>::dot(mm[mm.lookup(nz3,clipper::MM::ANY)].coord_orth(), n);
    z4 = clipper::Vec3<clipper::ftype>::dot(mm[mm.lookup(nz4,clipper::MM::ANY)].coord_orth(), n);
    z5 = clipper::Vec3<clipper::ftype>::dot(mm[mm.lookup(nz5,clipper::MM::ANY)].coord_orth(), n);
    z6 = clipper::Vec3<clipper::ftype>::dot(mm[mm.lookup(nz6,clipper::MM::ANY)].coord_orth(), n);

		///////// handedness //////////////

		clipper::MAtom nz6_substituent;
		clipper::ftype z6_substituent;

		nz6_substituent.set_id("XXX");

		const std::vector<clipper::MAtomIndexSymmetry> neighbourhood = this->sugar_parent_molecule_nonbond->atoms_near(ring_atoms[5].coord_orth(), 1.2);

		for ( int i = 0 ; i < neighbourhood.size() ; i++ )
		{
			if  (( mmol.atom( neighbourhood[i] ).element().trim() != "H" ) && (mmol.atom(neighbourhood[i]).name().trim() != ring_atoms[5].name().trim())) 
			// the target substituent could be anything apart from H, in-ring oxygen or in-ring carbon
			{
				if ( bonded ( mmol.atom(neighbourhood[i]), ring_atoms[5]) ) // check link
					if ( !is_part_of_ring(mmol.atom(neighbourhood[i]), ring_atoms) && (ring_atoms[5].occupancy() == mmol.atom(neighbourhood[i]).occupancy() ) ) // eliminate ring neighbours
					{
						if (mmol.atom(neighbourhood[i]).element().trim() != "C") nz6_substituent = mmol.atom(neighbourhood[i]); 
							// don't grab a carbon as substituent unless there's really no other option
						else if (nz6_substituent.name().trim() == "XXX") nz6_substituent = mmol.atom(neighbourhood[i]);
					}
			}
		}

		nz6_substituent.transform(shift); // we still need to recentre the atom
		z6_substituent = clipper::Vec3<clipper::ftype>::dot(nz6_substituent.coord_orth(), n);

		#ifdef DUMP
			DBG << "last in-ring carbon has occupancy " << ring_atoms[5].occupancy() << " and it's substituent is " 
			<< nz6_substituent.name().trim() << " with occupancy " << nz6_substituent.occupancy() << std::endl;
		#endif

	///////// stereochemistry ///////////

    z_anomeric_carbon = clipper::Vec3<clipper::ftype>::dot( stereo.first.first.coord_orth(), n );
    z_anomeric_substituent = clipper::Vec3<clipper::ftype>::dot( stereo.first.second.coord_orth(), n );

    z_configurational_carbon = clipper::Vec3<clipper::ftype>::dot( stereo.second.first.coord_orth(), n );
    z_configurational_substituent = clipper::Vec3<clipper::ftype>::dot( stereo.second.second.coord_orth(), n );

    clipper::ftype totalPuckering = sqrt(pow(z1,2) + pow(z2,2) + pow(z3,2) + pow(z4,2) + pow(z5,2) + pow(z6,2));

    clipper::ftype phi2, theta, q2, q3, angleCos, argASin;

    q3 = sqrt(1.0/6.0) * ( z3 + z5 + z1 - z2 - z4 - z6);

    theta = acos(q3 / totalPuckering);
    q2 = totalPuckering * sin (theta);

    theta *= (180.0/clipper::Util::pi()); // convert theta to degrees for sharing data

    angleCos = acos ( sqrt(1.0/3.0) * (z1 + z2*cos(4.0*clipper::Util::pi()/6.0) 
                + z3*cos(8.0*clipper::Util::pi()/6.0) + z4*cos(12.0*clipper::Util::pi()/6.0) 
                + z5*cos(16.0*clipper::Util::pi()/6.0) + z6*cos(20.0*clipper::Util::pi()/6.0)) / q2 );
    
    argASin = -sqrt(1.0/3.0) * (z2*sin(4.0*clipper::Util::pi()/6.0) + z3*sin(8.0*clipper::Util::pi()/6.0) 
                + z4*sin(12.0*clipper::Util::pi()/6.0) + z5*sin(16.0*clipper::Util::pi()/6.0) 
                + z6*sin(20.0*clipper::Util::pi()/6.0));

    // two possible values for phi, check other equivalence with sin (eqn 13 in the cremer-pople paper)

    if (float(q2*sin(angleCos)) == float(argASin)) phi2 = angleCos; // we compute internally in ftype (double precision) and then lower the precision for a reasonable comparison
    else phi2 = 2*clipper::Util::pi() - angleCos;

    phi2 *= (180.0/clipper::Util::pi()); // convert phi2 to degrees as well

    std::vector<clipper::ftype> cpParams;

	this->sugar_cremer_pople_params.push_back(totalPuckering);
	this->sugar_cremer_pople_params.push_back(phi2);
	this->sugar_cremer_pople_params.push_back(theta);

	cpParams.push_back(totalPuckering);
    cpParams.push_back(phi2);
    cpParams.push_back(theta);
    cpParams.push_back(q2);
    cpParams.push_back(q3);
    cpParams.push_back(this->sugar_conformation);

    if (( (z_anomeric_substituent > z_anomeric_carbon) && (z_configurational_substituent > z_configurational_carbon) ) || ( (z_anomeric_substituent < z_anomeric_carbon) && (z_configurational_substituent < z_configurational_carbon) ) )
		{
			if (lurd_reverse)
			{
				cpParams.push_back(anomer_beta);
				this->sugar_anomer="beta";
			}
			else
			{
				cpParams.push_back(anomer_alpha);
				this->sugar_anomer="alpha";
			}
		}
		else
		{
			if (lurd_reverse)
			{
				cpParams.push_back(anomer_alpha);
				this->sugar_anomer="alpha";
			}
			else
			{
				cpParams.push_back(anomer_beta);
				this->sugar_anomer="beta";
			}
		}

		#ifdef DUMP
			DBG << "an_c= " << anomeric_carbon << "/" << z_anomeric_carbon << " - an_subs= " << anomeric_substituent << "/" << z_anomeric_substituent << std::endl;
			DBG << "conf_c= " << configurational_carbon << "/" << z_configurational_carbon << " - conf_subs= " << configurational_substituent << "/" << z_configurational_substituent << std::endl;
			DBG << "z6= " << z6 << " z6_subs= " << z6_substituent << std::endl;
		#endif

    cpParams.push_back( z6 - z6_substituent );

    if ((is_part_of_ring(nz6_substituent, this->sugar_ring_elements)) || (nz6_substituent.name().trim() == "XXX")) this->sugar_handedness = "N";
    else if ((z6 - z6_substituent) < 0) this->sugar_handedness = "D"; else this->sugar_handedness = "L";

	#ifdef DUMP
		DBG << "Finished Cremer-Pople analysis, returning to caller..." << std::endl;
	#endif

    return cpParams;

}


    /*! Internal function, not for public use
        \param mmol The parent MiniMol object
	\return A vector of real values containing the cremer-pople parameters
    */

    std::vector<clipper::ftype> MSugar::cremerPople_furanose(const clipper::MiniMol& mmol, clipper::MMonomer mm) // modifies object, copies mm
    {
        clipper::ftype centre_x, centre_y, centre_z;
        centre_x = centre_y = centre_z = 0.0;

	MSugar sugar = *this;

	std::vector<clipper::MAtom> ring_atoms = ring_members();

	bool lurd_reverse = false; // to be used when the configurational carbon is ranked lower than the carbon making the link,
				   // which is involved in a configurational switch that effectively reverses the LURD mnemonic

	clipper::String nz1 = ring_atoms[0].name().trim(); // determine the relevant atoms
        clipper::String nz2 = ring_atoms[1].name().trim();
        clipper::String nz3 = ring_atoms[2].name().trim();
        clipper::String nz4 = ring_atoms[3].name().trim();
        clipper::String nz5 = ring_atoms[4].name().trim();

	stereochemistry_pairs stereo = get_stereochemistry( mmol );
	clipper::String anomeric_carbon = stereo.first.first.name().trim();
	clipper::String anomeric_substituent = stereo.first.second.name().trim();
	clipper::String configurational_carbon = stereo.second.first.name().trim();
	clipper::String configurational_substituent = stereo.second.second.name().trim();

	this->sugar_anomeric_carbon = stereo.first.first;
	this->sugar_anomeric_substituent = stereo.first.second;
	this->sugar_configurational_carbon = stereo.second.first;
	this->sugar_configurational_substituent = stereo.second.second;

	#ifdef DUMP
		DBG << "After getting the stereochemistry" << std::endl;
	#endif 

	
	if ( configurational_carbon != "XXX" )
	    if ( ( ring_atoms[4].name().trim() == configurational_carbon) || !is_part_of_ring( sugar_configurational_carbon, ring_atoms ) ) 
                // we have to take into account the rotation performed prior to closing the ring
	        lurd_reverse = true;

        // find geometrical centre
        centre_x += sugar[sugar.lookup(nz1,clipper::MM::ANY)].coord_orth().x();
        centre_y += sugar[sugar.lookup(nz1,clipper::MM::ANY)].coord_orth().y();
        centre_z += sugar[sugar.lookup(nz1,clipper::MM::ANY)].coord_orth().z();

        centre_x += sugar[sugar.lookup(nz2,clipper::MM::ANY)].coord_orth().x();
        centre_y += sugar[sugar.lookup(nz2,clipper::MM::ANY)].coord_orth().y();
        centre_z += sugar[sugar.lookup(nz2,clipper::MM::ANY)].coord_orth().z();

        centre_x += sugar[sugar.lookup(nz3,clipper::MM::ANY)].coord_orth().x();
        centre_y += sugar[sugar.lookup(nz3,clipper::MM::ANY)].coord_orth().y();
        centre_z += sugar[sugar.lookup(nz3,clipper::MM::ANY)].coord_orth().z();

        centre_x += sugar[sugar.lookup(nz4,clipper::MM::ANY)].coord_orth().x();
        centre_y += sugar[sugar.lookup(nz4,clipper::MM::ANY)].coord_orth().y();
        centre_z += sugar[sugar.lookup(nz4,clipper::MM::ANY)].coord_orth().z();

        centre_x += sugar[sugar.lookup(nz5,clipper::MM::ANY)].coord_orth().x();
        centre_y += sugar[sugar.lookup(nz5,clipper::MM::ANY)].coord_orth().y();
        centre_z += sugar[sugar.lookup(nz5,clipper::MM::ANY)].coord_orth().z();

        centre_x = centre_x / 5;
        centre_y = centre_y / 5;
        centre_z = centre_z / 5;

        clipper::Coord_orth centre(centre_x,centre_y,centre_z);
	this->sugar_centre = centre;

        #ifdef DUMP
	    DBG << "Ring centre: " << centre.format() << std::endl;
	#endif

        clipper::RTop_orth shift(clipper::Mat33<>::identity(), (- centre)); //clipper::Vec3<>::null()
        sugar.transform(shift); // recentre the sugar coordinates

	stereo.first.first.transform(shift);
	stereo.first.second.transform(shift);

	stereo.second.first.transform(shift);
	stereo.second.second.transform(shift);

        clipper::Vec3<clipper::ftype> rPrime(0.0, 0.0, 0.0);
        clipper::Vec3<clipper::ftype> r2Prime(0.0, 0.0, 0.0);
        clipper::Vec3<clipper::ftype> n(0.0, 0.0, 0.0);

        double argument = 0.0; // (j-1) = 0 for the first case
        rPrime +=  sugar[sugar.lookup(nz1,clipper::MM::ANY)].coord_orth() * sin(argument);
        r2Prime += sugar[sugar.lookup(nz1,clipper::MM::ANY)].coord_orth() * cos(argument);

        argument = (2.0 * clipper::Util::pi())/5.0;
        rPrime +=  sugar[sugar.lookup(nz2,clipper::MM::ANY)].coord_orth() * sin(argument);
        r2Prime += sugar[sugar.lookup(nz2,clipper::MM::ANY)].coord_orth() * cos(argument);

        argument = (2.0 * clipper::Util::pi() * 2)/5.0;
        rPrime +=  sugar[sugar.lookup(nz3,clipper::MM::ANY)].coord_orth() * sin(argument);
        r2Prime += sugar[sugar.lookup(nz3,clipper::MM::ANY)].coord_orth() * cos(argument);

        argument = (2.0 * clipper::Util::pi() * 3)/5.0;
        rPrime +=  sugar[sugar.lookup(nz4,clipper::MM::ANY)].coord_orth() * sin(argument);
        r2Prime += sugar[sugar.lookup(nz4,clipper::MM::ANY)].coord_orth() * cos(argument);

        argument = (2.0 * clipper::Util::pi() * 4)/5.0;
        rPrime +=  sugar[sugar.lookup(nz5,clipper::MM::ANY)].coord_orth() * sin(argument);
        r2Prime += sugar[sugar.lookup(nz5,clipper::MM::ANY)].coord_orth() * cos(argument);

        n = (clipper::Vec3<clipper::ftype>::cross(rPrime, r2Prime)).unit();

        clipper::ftype z1, z2, z3, z4, z5, z_anomeric_carbon, z_anomeric_substituent, z_configurational_carbon, z_configurational_substituent;

        z1 = clipper::Vec3<clipper::ftype>::dot(sugar[sugar.lookup(nz1,clipper::MM::ANY)].coord_orth(), n);
        z2 = clipper::Vec3<clipper::ftype>::dot(sugar[sugar.lookup(nz2,clipper::MM::ANY)].coord_orth(), n);
        z3 = clipper::Vec3<clipper::ftype>::dot(sugar[sugar.lookup(nz3,clipper::MM::ANY)].coord_orth(), n);
        z4 = clipper::Vec3<clipper::ftype>::dot(sugar[sugar.lookup(nz4,clipper::MM::ANY)].coord_orth(), n);
        z5 = clipper::Vec3<clipper::ftype>::dot(sugar[sugar.lookup(nz5,clipper::MM::ANY)].coord_orth(), n);

	///////// handedness //////////////

	clipper::MAtom nz5_substituent;
	clipper::ftype z5_substituent;

	nz5_substituent.set_id("XXX");

	const std::vector<clipper::MAtomIndexSymmetry> neighbourhood = this->sugar_parent_molecule_nonbond->atoms_near(ring_atoms[4].coord_orth(), 1.5);

	for ( int i = 0 ; i < neighbourhood.size() ; i++ )
	{
	    #ifdef DUMP
	        DBG << "Neighbour found for " << ring_atoms[4].name() << ": " << (mmol.atom(neighbourhood[i]).name()) 
                    << " at distance " << clipper::Coord_orth::length ( mmol.atom(neighbourhood[i]).coord_orth(), ring_atoms[4].coord_orth() ) << std::endl;
	    #endif
		
	    if (( mmol.atom(neighbourhood[i]).element().trim() != "H" ) && 
                (mmol.atom(neighbourhood[i]).name().trim() != ring_atoms[4].name().trim())) // the target substituent could be anything apart from H, in-ring oxygen or in-ring carbon
	    {
	        if ( bonded(mmol.atom(neighbourhood[i]), ring_atoms[4]) ) // check link
		    if ( !is_part_of_ring ( mmol.atom(neighbourhood[i]), ring_atoms ) ) // eliminate ring neighbours
		    {
		        #ifdef DUMP
                            DBG << "Analysing " <<  mmol.atom(neighbourhood[i]).name() << " with ID " <<  mmol.atom(neighbourhood[i]).id() << " as substituent for C4" << std::endl;
			#endif  
			if ( mmol.atom(neighbourhood[i]).element().trim() != "C" ) 
                            nz5_substituent = mmol.atom ( neighbourhood[i] ); // don't grab a carbon as substituent unless there's really no other option
			else if ( nz5_substituent.name().trim() == "XXX" ) 
                            nz5_substituent = mmol.atom ( neighbourhood[i] );
		    }
	    }
	}

	#ifdef DUMP
		DBG << "substituent at last in-ring carbon: " << nz5_substituent.name().trim() << std::endl;
	#endif

	nz5_substituent.transform(shift); // we still need to recentre the atom
	z5_substituent = clipper::Vec3<clipper::ftype>::dot(nz5_substituent.coord_orth(), n);

	///////// stereochemistry ///////////

        z_anomeric_carbon = clipper::Vec3<clipper::ftype>::dot( stereo.first.first.coord_orth(), n );
        z_anomeric_substituent = clipper::Vec3<clipper::ftype>::dot( stereo.first.second.coord_orth(), n );

        z_configurational_carbon = clipper::Vec3<clipper::ftype>::dot( stereo.second.first.coord_orth(), n );
        z_configurational_substituent = clipper::Vec3<clipper::ftype>::dot( stereo.second.second.coord_orth(), n );

	#ifdef DUMP
            DBG << "Projections done." << std::endl;
	#endif

        clipper::ftype totalPuckering = sqrt(pow(z1,2) + pow(z2,2) + pow(z3,2) + pow(z4,2) + pow(z5,2));

        clipper::ftype phi2, q2, argACos, argASin;

        argACos =  sqrt(2.0/5.0) * (z1 + z2*cos(4.0*clipper::Util::pi()/5.0) 
                    + z3*cos(8.0*clipper::Util::pi()/5.0) + z4*cos(12.0*clipper::Util::pi()/5.0) 
                    + z5*cos(16.0*clipper::Util::pi()/5.0));
    
        argASin = -sqrt(2.0/5.0) * (z2*sin(4.0*clipper::Util::pi()/5.0) + z3*sin(8.0*clipper::Util::pi()/5.0) 
                    + z4*sin(12.0*clipper::Util::pi()/5.0) + z5*sin(16.0*clipper::Util::pi()/5.0));
    
        phi2 = (clipper::Util::pi()/2) - atan( argACos / argASin); // we want arccotan but there's no such thing in C math I'm afraid
        clipper::ftype phi3 = atan (argASin / argACos );
        //phi2 += clipper::Util::pi()/2; // atan(w) result is given in the [-pi/2,+pi/2] range

        q2 = (argACos / cos(phi2));
        
//        if (float((q2*sin(phi2)) != float(argASin) )) phi2 = 2*clipper::Util::pi() - phi2;
        
//        std::cout << " phi2: " << phi2 << " phi3: " << phi3 << " q2: " << q2 << std::endl;
//        std::cout << " argASin: " << argASin << " q2*sin(phi2): " << float(q2*sin(phi2)) << " q2*sin(phi3): " << float(q2*sin(phi3)) << std::endl; 
//        std::cout << " argACos: " << argACos << " q2*cos(phi2): " << float(q2*cos(phi2)) << " q2*cos(phi3): " << float(q2*cos(phi3)) << std::endl;

        q2 < 0 ? phi2 = 180.0 + phi2 * ( 180.0/clipper::Util::pi() ) :  phi2 *= (180.0/clipper::Util::pi()); // convert phi2 to degrees as well
        
        if ( q2 < 0 ) q2 *= -1;

        std::vector<clipper::ftype> cpParams;

	this->sugar_cremer_pople_params.push_back(totalPuckering);
	this->sugar_cremer_pople_params.push_back(phi2);    
	this->sugar_cremer_pople_params.push_back(-1); // there's no theta for furanoses

	cpParams.push_back(totalPuckering);
        cpParams.push_back(phi2);
        cpParams.push_back(-1); // theta not defined for furanoses
        cpParams.push_back(q2);
        cpParams.push_back(-1); // m=2, so there's no q3
    
        if (( (z_anomeric_substituent > z_anomeric_carbon) && (z_configurational_substituent > z_configurational_carbon) ) || 
            ( (z_anomeric_substituent < z_anomeric_carbon) && (z_configurational_substituent < z_configurational_carbon) ) )
	{
	    if (lurd_reverse)
	    {
		cpParams.push_back(anomer_beta);
		this->sugar_anomer="beta";
	    }
	    else
	    {
                cpParams.push_back(anomer_alpha);
                this->sugar_anomer="alpha";
	    }
	}
	else
	{
            if (lurd_reverse)
            {
                cpParams.push_back(anomer_alpha);
                this->sugar_anomer="alpha";
            }
            else
            {
                cpParams.push_back(anomer_beta);
                this->sugar_anomer="beta";
            }
        }

	#ifdef DUMP
            DBG << "an_c= " << anomeric_carbon << "/" << z_anomeric_carbon << " - an_subs= " << anomeric_substituent << "/" << z_anomeric_substituent << std::endl;
            DBG << "conf_c= " << configurational_carbon << "/" << z_configurational_carbon << " - conf_subs= " 
                << configurational_substituent << "/" << z_configurational_substituent << std::endl;
	    DBG << "z5= " << z5 << " z5_subs= " << z5_substituent << std::endl;
	#endif

        cpParams.push_back( z5 - z5_substituent );

        if ((is_part_of_ring(nz5_substituent, this->sugar_ring_elements)) || (nz5_substituent.name().trim() == "XXX")) this->sugar_handedness = "N";
        else if ((z5 - z5_substituent) < 0) this->sugar_handedness = "D"; else this->sugar_handedness = "L";

        return cpParams;
}


/*! Internal function, not for public use
	\param phi Angle1 of Cremer-Pople calculations
	\param theta Angle2 of Cremer-Pople calculations
	\return An integer value describing the conformation
*/

int MSugar::conformationPyranose(const clipper::ftype& phi, const clipper::ftype& theta) const
{
    int confCode = 0;

    if (theta <= 22.5) confCode = conf_pyranose_4C1; // canonical chair
    else if ((theta > 22.5) && (theta <= 67.5)) // envelopes and half-chairs
             {
                 if ((phi > 15.0) && (phi <= 45.0)) confCode = conf_pyranose_OH1;
                 else if ((phi > 45.0) && (phi <= 75.0)) confCode = conf_pyranose_E1;
                 else if ((phi > 75.0) && (phi <= 105.0)) confCode = conf_pyranose_2H1;
                 else if ((phi > 105.0) && (phi <= 135.0)) confCode = conf_pyranose_2E;
                 else if ((phi > 135.0) && (phi <= 165.0)) confCode = conf_pyranose_2H3;
                 else if ((phi > 165.0) && (phi <= 195.0)) confCode = conf_pyranose_E3;
                 else if ((phi > 195.0) && (phi <= 225.0)) confCode = conf_pyranose_4H3;
                 else if ((phi > 225.0) && (phi <= 255.0)) confCode = conf_pyranose_4E;
                 else if ((phi > 255.0) && (phi <= 285.0)) confCode = conf_pyranose_4H5;
                 else if ((phi > 285.0) && (phi <= 315.0)) confCode = conf_pyranose_E5;
                 else if ((phi > 315.0) && (phi <= 345.0)) confCode = conf_pyranose_OH5;
                 else if ((phi > 345.0) || (phi <= 15.0)) confCode = conf_pyranose_OE;
             }
    else if ((theta > 67.5) && (theta <= 112.5)) // boats and skew boats
             {
                 if ((phi > 15.0) && (phi <= 45.0)) confCode = conf_pyranose_3S1;
                 else if ((phi > 45.0) && (phi <= 75.0)) confCode = conf_pyranose_B14;
                 else if ((phi > 75.0) && (phi <= 105.0)) confCode = conf_pyranose_5S1;
                 else if ((phi > 105.0) && (phi <= 135.0)) confCode = conf_pyranose_25B;
                 else if ((phi > 135.0) && (phi <= 165.0)) confCode = conf_pyranose_2SO;
                 else if ((phi > 165.0) && (phi <= 195.0)) confCode = conf_pyranose_B3O;
                 else if ((phi > 195.0) && (phi <= 225.0)) confCode = conf_pyranose_1S3;
                 else if ((phi > 225.0) && (phi <= 255.0)) confCode = conf_pyranose_14B;
                 else if ((phi > 255.0) && (phi <= 285.0)) confCode = conf_pyranose_1S5;
                 else if ((phi > 285.0) && (phi <= 315.0)) confCode = conf_pyranose_B25;
                 else if ((phi > 315.0) && (phi <= 345.0)) confCode = conf_pyranose_OS2;
                 else if ((phi > 345.0) || (phi <= 15.0)) confCode = conf_pyranose_3OB;

             }
    else if ((theta > 112.5) && (theta <= 157.5)) // envelopes and half-chairs
             {
                 if ((phi > 15.0) && (phi <= 45.0)) confCode = conf_pyranose_3H4;
                 else if ((phi > 45.0) && (phi <= 75.0)) confCode = conf_pyranose_E4;
                 else if ((phi > 75.0) && (phi <= 105.0)) confCode = conf_pyranose_5H4;
                 else if ((phi > 105.0) && (phi <= 135.0)) confCode = conf_pyranose_5E;
                 else if ((phi > 135.0) && (phi <= 165.0)) confCode = conf_pyranose_5HO;
                 else if ((phi > 165.0) && (phi <= 195.0)) confCode = conf_pyranose_EO;
                 else if ((phi > 195.0) && (phi <= 225.0)) confCode = conf_pyranose_1HO;
                 else if ((phi > 225.0) && (phi <= 255.0)) confCode = conf_pyranose_1E;
                 else if ((phi > 255.0) && (phi <= 285.0)) confCode = conf_pyranose_1H2;
                 else if ((phi > 285.0) && (phi <= 315.0)) confCode = conf_pyranose_E2;
                 else if ((phi > 315.0) && (phi <= 345.0)) confCode = conf_pyranose_3H2;
                 else if ((phi > 345.0) || (phi <= 15.0)) confCode = conf_pyranose_3E;

             }
    else if (theta >= 157.5) confCode = conf_pyranose_1C4; // canonical chair

    return confCode;
}


/*! Internal function, not for public use
	\param phi Angle1 of Cremer-Pople calculations
	\return An integer value describing the conformation
*/

int MSugar::conformationFuranose(const clipper::ftype& phi) const
{
    int confCode = 0;
    
    if ((phi > 9.0) && (phi <= 27.0)) confCode = conf_furanose_OT1;
    else if ((phi > 27.0)  && (phi <=  45.0)) confCode = conf_furanose_EV1;
    else if ((phi > 45.0)  && (phi <=  63.0)) confCode = conf_furanose_2T1;
    else if ((phi > 63.0)  && (phi <=  81.0)) confCode = conf_furanose_2EV;
    else if ((phi > 81.0)  && (phi <=  99.0)) confCode = conf_furanose_2T3;
    else if ((phi > 99.0)  && (phi <= 117.0)) confCode = conf_furanose_EV3;
    else if ((phi >117.0)  && (phi <= 135.0)) confCode = conf_furanose_4T3;
    else if ((phi >135.0)  && (phi <= 153.0)) confCode = conf_furanose_4EV;
    else if ((phi >153.0)  && (phi <= 171.0)) confCode = conf_furanose_4TO;
    else if ((phi > 171.0) && (phi <= 189.0)) confCode = conf_furanose_EVO;
    else if ((phi > 189.0) && (phi <= 207.0)) confCode = conf_furanose_1TO;
    else if ((phi > 207.0) && (phi <= 225.0)) confCode = conf_furanose_1EV;
    else if ((phi > 225.0) && (phi <= 243.0)) confCode = conf_furanose_1T2;
    else if ((phi > 243.0) && (phi <= 261.0)) confCode = conf_furanose_EV2;
    else if ((phi > 261.0) && (phi <= 279.0)) confCode = conf_furanose_3T2;
    else if ((phi > 279.0) && (phi <= 297.0)) confCode = conf_furanose_3EV;
    else if ((phi > 297.0) && (phi <= 315.0)) confCode = conf_furanose_3T4;
    else if ((phi > 315.0) && (phi <= 333.0)) confCode = conf_furanose_EV4;
    else if ((phi > 333.0) && (phi <= 351.0)) confCode = conf_furanose_OT4;
    else if ((phi > 351.0) || (phi <=   9.0)) confCode = conf_furanose_OEV;

    return confCode;
}


/*! Internal function, not for public use. Uses a recursive version of Fleury's algorithm for finding eulerian cycles in undirected graphs
	\return A vector of MAtoms containing the members of the sugar ring
*/

std::vector<clipper::MAtom> MSugar::ringMembers() const
{

	MSugar mm = *this;

	MSugar::visited_arcs background;
	background = std::vector<std::pair<clipper::MAtom, clipper::MAtom> >();

	std::vector<clipper::MAtom> buffer = findPath(mm, 0, background);


	// now, we must filter and reorder the result:

	int index = 1;
	if (buffer.size() > 2)
		while ( (buffer[0].name().trim() != buffer[index].name().trim()) && (index < buffer.size()) ) index++;

	buffer.erase( buffer.begin() );
	buffer.resize( index, clipper::MAtom() );

	std::vector<clipper::MAtom> result;

	#ifdef DUMP
		DBG << "Dumping ring contents... " << std::endl;
		for (int runner = 0 ; runner < buffer.size() ; runner++ ) DBG << buffer[runner].name().trim() << " with occupancy = " << buffer[runner].occupancy() << std::endl;
	#endif

	for ( int i = 0 ; i < buffer.size() ; i++ ) if ( buffer[i].element().trim() == "O" ) result.push_back( buffer[i] ); // find the oxygen and assign it to the first position

	for ( int i = 0 ; i < buffer.size() ; i++ )
	{
		if ( buffer[i].element().trim() == "C" )
		{
			if ( result.size() > 1 )
			{
				// grab carbon ranks and compare them

				std::stringstream stream_1(buffer[i].name().trim().split("C")[0]);
				int origin;
				stream_1 >> origin;

				std::stringstream stream_2(result[1].name().trim().split("C")[0]);
				int destination;
				stream_2 >> destination;

				//#ifdef DUMP
				//	DBG << buffer[i].name().trim() << " with index " << origin << " VS " << result[1].name().trim() << " with index " << destination << std::endl;
				//#endif

				int subindex = 1;

				while ( ( origin > destination ) && (subindex < result.size() ) )
				{

				//#ifdef DUMP
				//	DBG << "Subindex is " << subindex << " and result[subindex] is " << result[subindex].name() << std::endl;
				//#endif

					if ( ++subindex < result.size() )
					{
						std::stringstream stream_3( result[subindex].name().trim().split("C")[0] );
						stream_3 >> destination;
					}
				}

				//#ifdef DUMP
				//	DBG << "Inserting " << origin << " which is smaller than " << destination << " at position " << subindex << " and result.size() is " << result.size() << std::endl;
				//#endif

				result.insert( result.begin() + subindex , buffer[i] );
			}
			else
			{
				result.push_back( buffer[i] );
				//#ifdef DUMP
				//	DBG << "Inserting " << buffer[i].name().trim() << " as first item after the oxygen..." << std::endl;
				//#endif
			}
		}
	}

	#ifdef DUMP
		DBG << "Successfully determined ring members!" << std::endl;
	#endif

	return result;
}


/*! Internal function, not for public use
	\param mm A monomer containing a sugar ring
	\param current_atom The index we're checking now
	\param background A visited_arcs data structure containing those pairs previously visited
	\return A vector of atoms which are within reach of the input atom
*/

const std::vector<clipper::MAtom> MSugar::findPath(const clipper::MMonomer& mm, int current_atom, MSugar::visited_arcs& background) const
{
	std::vector<clipper::MAtom> result;
	std::vector<clipper::MAtom> available_paths;

	result = std::vector<clipper::MAtom>();
	available_paths = std::vector<clipper::MAtom>();

	available_paths = findBonded(mm[current_atom], background);

	if (available_paths.size() == 0)
	{
		return result; // return empty vector in case we have explored all the available bonds
	}

	for ( int i = 0 ; i < available_paths.size() ; i++ )
	{

		background.push_back( std::pair<clipper::MAtom, clipper::MAtom> (mm[current_atom], available_paths[i] ) );

		if ( closes_ring( available_paths[i], background ) && result.empty())
		{
			result.push_back( available_paths[i] ); // we'll check for this atom at a later stage
			result.push_back( mm[current_atom] );
			return result;
		}
		else if (!result.empty())
		{
			result.push_back( mm[current_atom] );
			return result;
		}
		else
		{
			std::vector<clipper::MAtom> tmpResult = findPath(mm, mm.lookup(available_paths[i].name().trim(),clipper::MM::ANY), background);
			result.insert(result.end(), tmpResult.begin(), tmpResult.end());
		}
	}

	if (!result.empty()) result.push_back(mm[current_atom]);

	return result;

}

/*! Internal function, not for public use
	\param ma An atom which may eventually close the ring
	\param background A visited_arcs data structure containing those pairs previously visited
	\return True if the atom closes the ring, false otherwise
*/

bool MSugar::closes_ring(const clipper::MAtom& ma, MSugar::visited_arcs& background) const
{
	for ( int i = 0 ; i < background.size() ; i++ )
	{
		if ( background[i].first.name().trim() == ma.name().trim() )
		{
			return true;
		}
	}
	return false;
}

/*! Internal function for getting bonded atoms from the same MMonomer. Doesn't check other MMonomers
	\param ma The atom whose bonded neighbours we intent to retrieve
	\param background A struct visited_arcs to check for already visited nodes
	\return An std::vector of MAtoms containing the bonded atoms
*/

const std::vector<clipper::MAtom> MSugar::findBonded(const clipper::MAtom& ma, MSugar::visited_arcs& background ) const
{
	std::vector<clipper::MAtom> result;

	MSugar mm = *this;

	for (int i = 0; i < mm.atom_list().size() ; i++)
	{
		clipper::ftype bond_length = clipper::Coord_orth::length(mm[i].coord_orth(),ma.coord_orth());

		if ((bond_length > 0.5) && (bond_length < 1.61)) // any type of bond
			if (!lookup_visited(background, std::pair<clipper::MAtom, clipper::MAtom>(ma, mm[i])))
				if (get_altconf(ma) == get_altconf(mm[i])) // we need to stay on the same conformation
					result.push_back(mm[i]);
	}

	return result;
}

/*! Internal function for checking if an arc (bond) of a graph (molecule) has already been visited
	\param background A visited_arcs data structure containing the graph-crawling history
	\param pairs An std::pair containing the arc to be visited, composed of two clipper::MAtom
	\return True if the arc has been visited, false otherwise
*/

bool MSugar::lookup_visited(const MSugar::visited_arcs& va, const std::pair<clipper::MAtom, clipper::MAtom> pairs) const
{
	for (int i = 0 ; i < va.size() ; i++) if (((pairs.first.name().trim() == va[i].first.name().trim()) && ((pairs.second.name().trim() == va[i].second.name().trim()))) || ((pairs.first.name().trim() == va[i].second.name().trim()) && ((pairs.second.name().trim() == va[i].first.name().trim())))) return true;
	return false;
}

/*! Internal function for getting the two carbon atoms required for the determination of the anomer type (alpha/beta)
	\param mmol A clipper::MiniMol object containing the parent molecule
	\return A stereochemistry_pairs data structure containing the anomeric carbon and the highest ranked stereocentre with their respective substituents.
					In case of not finding an atom: id='XXX'.
*/

MSugar::stereochemistry_pairs MSugar::get_stereochemistry(const clipper::MiniMol& mmol) // returns a pair of (C,substituent)-containing pairs of atoms for anomer type determination
{

	MSugar mm = *this;

	std::pair < std::pair<clipper::MAtom, clipper::MAtom>, std::pair<clipper::MAtom, clipper::MAtom > > result;
	std::vector<clipper::MAtom> ring_atoms = ring_members();

	clipper::MAtom anomeric_carbon, anomeric_substituent, configurational_carbon, configurational_substituent;

	anomeric_carbon = clipper::MAtom();
	anomeric_substituent = clipper::MAtom();
	configurational_carbon = clipper::MAtom();
	configurational_substituent = clipper::MAtom();

	anomeric_carbon.set_id("XXX");
	anomeric_substituent.set_id("XXX");

	configurational_carbon.set_id("XXX");
	configurational_substituent.set_id("XXX");

	if ( !ring_atoms.empty() )
	{
            if ( ring_atoms[1].element().trim() == "C" )        // the atom in position 1 is the anomeric carbon, let's identify it's substituent:
            {
                anomeric_carbon = ring_atoms[1]; // we're always getting the anomeric carbon here, given the way we detect/order the ring

                const std::vector<clipper::MAtomIndexSymmetry> neighbourhood = this->sugar_parent_molecule_nonbond->atoms_near(ring_atoms[1].coord_orth(), 1.2);

	        for ( int i = 0 ; i < neighbourhood.size() ; i++ )
                {
                    if (( mmol.atom(neighbourhood[i]).element().trim() != "H" ) 
                    && (mmol.atom(neighbourhood[i]).name().trim() != anomeric_carbon.name().trim()) )
		    //&& (( get_altconf(mmol.atom(neighbourhood[i])) == ' ' ) 
		    //|| ( get_altconf(mmol.atom(neighbourhood[i])) == 'A' )))
		    // the target substituent could be anything apart from H, in-ring oxygen or in-ring carbon
		    {
		        if ( bonded(neighbourhood[i], ring_atoms[1] ) )
                        {
                            if ( !is_part_of_ring(mmol.atom(neighbourhood[i]), ring_atoms))
                            {
                                if ( get_altconf(mmol.atom(neighbourhood[i])) == get_altconf(ring_atoms[1]) )
                                {
                                    if (mmol.atom(neighbourhood[i]).element().trim() != "C") anomeric_substituent = mmol.atom(neighbourhood[i]); 
                                    // don't grab a carbon as substituent unless there's really no other option
                                    else if (anomeric_substituent.name().trim() == "XXX") anomeric_substituent = mmol.atom(neighbourhood[i]);

                                }
                                //}
                            else
                            {
                                if (mmol.atom(neighbourhood[i]).element().trim() != "C") anomeric_substituent = mmol.atom(neighbourhood[i]); 
                                // don't grab a carbon as substituent unless there's really no other option
                                else if (anomeric_substituent.name().trim() == "XXX") anomeric_substituent = mmol.atom(neighbourhood[i]);
                            }
                        }
                    }
                }
            }
        }
    }
	
	#ifdef DUMP
		DBG << "Anomeric carbon: " << anomeric_carbon.id() << "  Substituent: " << anomeric_substituent.id() << std::endl;
	#endif
 
	result.first.first = anomeric_carbon;
	result.first.second = anomeric_substituent;

	for ( int i = 2 ; i < ring_atoms.size() ; i++ ) // we start checking for the highest ranked stereocentre at the in-ring carbon next (clockwise) to the anomeric carbon
		if (ring_atoms[i].element().trim() == "C")
			if ( is_stereocentre(ring_atoms[i], mmol) )
			{
				configurational_carbon = ring_atoms[i];
				
				// get the configurational carbon's target substituent
				const std::vector<clipper::MAtomIndexSymmetry> neighbourhood = this->sugar_parent_molecule_nonbond->atoms_near(configurational_carbon.coord_orth(), 1.2); // was 1.4

				for ( int i2 = 0 ; i2 < neighbourhood.size() ; i2++ )
				{
					if ( !is_part_of_ring(mmol.atom(neighbourhood[i2]),ring_atoms) 
						&& (mmol.atom(neighbourhood[i2]).element().trim() != "H" )
						&& (get_altconf(mmol.atom(neighbourhood[i2])) == get_altconf(anomeric_carbon)))
			//			&& (( get_altconf(mmol.atom(neighbourhood[i2])) == ' ' )  
			//			|| ( get_altconf(mmol.atom(neighbourhood[i2])) == 'A' ) ) ) // the target substituent could be anything apart from H, in-ring oxygen or in-ring carbon
					{
						if ( bonded (mmol.atom(neighbourhood[i2]), configurational_carbon) ) // check link
							if (( mmol.atom(neighbourhood[i2]).name().trim() != ring_atoms[i-1].name().trim()) 
								&& ( mmol.atom(neighbourhood[i2]).name().trim() != ring_atoms[0].name().trim() ) 
								&& ( mmol.atom(neighbourhood[i2]).name().trim() != configurational_carbon.name().trim() )) 
									// eliminate ring neighbours
									configurational_substituent = mmol.atom(neighbourhood[i2]);
					}
				}
			}


	#ifdef DUMP
		DBG << "(in-ring) configurational carbon: " << configurational_carbon.id() << "  substituent: " << configurational_substituent.id() << std::endl;
	#endif

	// we've recorded the highest ranked in-ring carbon atom & substituent in configurational_*

	clipper::MAtom next_carbon = configurational_substituent; 

	while ( is_stereocentre( next_carbon, mmol ) && (next_carbon.id().trim() != configurational_carbon.name().trim() ) ) ////////////////////// CONTINUE HERE 
	{
		configurational_carbon = next_carbon;
		const std::vector<clipper::MAtomIndexSymmetry> neighbourhood = this->sugar_parent_molecule_nonbond->atoms_near(configurational_carbon.coord_orth(), 1.2);
	
		for ( int i2 = 0 ; i2 < neighbourhood.size() ; i2++ )
		{
			if ( !is_part_of_ring(mmol.atom(neighbourhood[i2]),ring_atoms) 
				&& (mmol.atom(neighbourhood[i2]).element().trim() != "H" ) 
				&& (get_altconf(mmol.atom(neighbourhood[i2])) == get_altconf(configurational_carbon)))
//				&& (( get_altconf(mmol.atom(neighbourhood[i2])) == ' ' ) 
//				|| ( get_altconf(mmol.atom(neighbourhood[i2])) == 'A' ) ) ) // the target substituent could be anything apart from H, in-ring oxygen or in-ring carbon
			{   // set next carbon and configurational (subs, carbon)
				if ( bonded (mmol.atom(neighbourhood[i2]), configurational_carbon) ) // check link
				{	
						if ( (mmol.atom(neighbourhood[i2])).element().trim() == "C")
						{	if ( clipper::Coord_orth::length( mmol.atom( neighbourhood[i2] ).coord_orth(), ring_atoms[ring_atoms.size()-1].coord_orth() ) 
							   > clipper::Coord_orth::length( configurational_carbon.coord_orth(), ring_atoms[ring_atoms.size()-1].coord_orth() ))  // A carbon, linked to this carbon and further away from ring 
									next_carbon = mmol.atom(neighbourhood[i2]);
						}
						else configurational_substituent = mmol.atom(neighbourhood[i2]);
						
						/*if (( (mmol.atom(neighbourhood[i2])).element().trim() == "C") && 
							 ( clipper::Coord_orth::length( mmol.atom( neighbourhood[i2] ).coord_orth(), ring_atoms[ring_atoms.size()-1].coord_orth() )) 
							   > clipper::Coord_orth::length( configurational_carbon.coord_orth(), ring_atoms[ring_atoms.size()-1].coord_orth() ) ) // A carbon, linked to this carbon and further away from ring 
							next_carbon = mmol.atom(neighbourhood[i2]);
						else configurational_substituent = mmol.atom(neighbourhood[i2]); */
				}							

			} //check distance to analyse carbons further away from last carbon
		}
	}

	result.second.first = configurational_carbon;
	result.second.second = configurational_substituent;
	
	#ifdef DUMP
		DBG << "Configurational carbon: " << configurational_carbon.id() << "  Substituent: " << configurational_substituent.id() << std::endl;
	#endif

	return result;

}

/*! Internal function for checking whether a carbon atom qualifies as a stereocentre
	\param ma A clipper::MAtom object containing the atom to check
	\param mmol A clipper::MiniMol object for getting the neighbours
	\return A logic value
*/

bool MSugar::is_stereocentre(const clipper::MAtom& ma, const clipper::MiniMol& mmol)
{
	std::vector< clipper::MAtom > substituent_list; 
	std::vector< clipper::MAtom > ring_atoms = ring_members();
		
	if ( ma.element().trim() != "C" ) return false;

	const std::vector<clipper::MAtomIndexSymmetry> neighbourhood = this->sugar_parent_molecule_nonbond->atoms_near(ma.coord_orth(), 1.2); // tried and tested with 1.4
		
	for ( int k = 0 ; k < neighbourhood.size() ; k++ )
	{

		clipper::ftype distance = 0.0;

		if ( neighbourhood[k].symmetry() == 0 )
		{
			distance = clipper::Coord_orth::length(mmol.atom(neighbourhood[k]).coord_orth(), ma.coord_orth());
		}
		else // this neighbour is actually a symmetry mate, so let's find the symmetry operator that places a copy of the molecule close to the anomeric carbon
		{
			clipper::Spacegroup spgr = mmol.spacegroup();
			clipper::Coord_frac f1 = mmol.atom( neighbourhood[k] ).coord_orth().coord_frac( mmol.cell() );
			clipper::Coord_frac f2 = ma.coord_orth().coord_frac( mmol.cell() );
			f1 = spgr.symop(neighbourhood[k].symmetry()) * f1;
			f1 = f1.lattice_copy_near( f2 );
			distance = sqrt(( f2 - f1 ).lengthsq( mmol.cell() ));
		}

		if ( distance < 2.0 ) // check link
			if (( mmol.atom(neighbourhood[k]).element().trim() != "H" ) && (mmol.atom(neighbourhood[k]).name().trim() != ma.name().trim() ) && (get_altconf(mmol.atom(neighbourhood[k])) == get_altconf(ma)) ) 
			{
				//#ifdef DUMP
				//	DBG << "Counting " << mmol.atom(neighbourhood[k]).id() << " as substituent from a total of " << neighbourhood.size() << " atoms with symop " << neighbourhood[k].symmetry() << std::endl;
				//#endif

				bool found = false;
				
				for (int subs_no = 0 ; subs_no < substituent_list.size() ; subs_no++ ) 
					if (( substituent_list[subs_no].element().trim() != "C") && ( substituent_list[subs_no].element().trim() == mmol.atom(neighbourhood[k]).element().trim() ) )
						found = true;
						
				if ( !found || (mmol.atom(neighbourhood[k]).name().trim() == ring_atoms[0].name().trim()))
				{	
					substituent_list.push_back(mmol.atom(neighbourhood[k]));
				}

			}
	}

	//#ifdef DUMP
	//	DBG << "Number of substituents: " << substituent_list.size() << std::endl;
	//#endif

	if ( substituent_list.size() > 2 ) return true;
	else return false;

}

/*! Internal function for checking whether an atom is part of the ring
	\param ma A clipper::MAtom object containing the atom to check.
	\param ring_atoms An std::vector of clipper::MAtom's containing the sugar ring
	\return A boolean value with the obvious answer
*/

bool MSugar::is_part_of_ring(const clipper::MAtom& ma, const std::vector<clipper::MAtom> ring_atoms) const
{
	for (int i = 0; i < ring_atoms.size(); i++)
		if ( ring_atoms[i].name().trim() == ma.name().trim() ) return true;

	return false;
}


/*! Internal function for checking whether two atoms are bonded or not
	\param ma_one A clipper::MAtom object
	\param ma_two A clipper::MAtom object
	\return A boolean value with the obvious answer
*/

bool MSugar::bonded(const clipper::MAtomIndexSymmetry& ma_one, const clipper::MAtom& ma_two) const
{
    ////////////// CRITICAL: check symm mates ///////////////

    const int poly = ma_one.polymer();
    const int mono = ma_one.monomer();
    const int atom = ma_one.atom();

    double distance = 0.0;

    if ( ma_one.symmetry() == 0 )
    {
    	distance = clipper::Coord_orth::length(sugar_parent_molecule->atom(ma_one).coord_orth(), ma_two.coord_orth());
    }
    else // this neighbour is actually a symmetry mate, so let's find the symmetry operator that places a copy of the molecule close to the anomeric carbon
    {
    	clipper::Spacegroup spgr = sugar_parent_molecule->spacegroup();
    	//clipper::Coord_frac f1 = sugar_parent_molecule[poly][mono][atom].coord_orth().coord_frac(sugar_parent_molecule->cell());
    	clipper::Coord_frac f1 = sugar_parent_molecule->atom(ma_one).coord_orth().coord_frac(sugar_parent_molecule->cell());
    	clipper::Coord_frac f2 = ma_two.coord_orth().coord_frac(sugar_parent_molecule->cell());
    	f1 = spgr.symop(ma_one.symmetry()) * f1;
    	f1 = f1.lattice_copy_near( f2 );
		
	distance = sqrt(( f2 - f1 ).lengthsq( sugar_parent_molecule->cell() ));

    }
	
    if ( sugar_parent_molecule->atom(ma_one).element().trim() == "C" )
    {
        if ( ma_two.element().trim() == "C" )
	    if ((distance >  1.18 ) && ( distance < 1.60 )) return true; // C-C or C=C
	    else return false;
	else if ( ma_two.element().trim() == "N" )
	    if ((distance >  1.24 ) && ( distance < 1.52 )) return true; // C-N or C=N
	    else return false;
	else if ( ma_two.element().trim() == "O" )
	    if ((distance >  1.16 ) && ( distance < 1.50 )) return true; // C-O or C=O
	    else return false;
        else if ( ma_two.element().trim() == "S" )
            if ((distance >  1.50 ) && ( distance < 2.00 )) return true; // C-S or C=S
            else return false;
	else if ( ma_two.element().trim() == "H" )
	{
	    if ((distance >  0.96 ) && ( distance < 1.14 )) return true; // C-H
	    else return false;
	}
    }
    else if ( sugar_parent_molecule->atom(ma_one).element().trim() == "N" )
    {
        if ( ma_two.element().trim() == "C" )
	    if ((distance >  1.24 ) && ( distance < 1.52 )) return true; // N-C or N=C
	    else return false;
        else if ( ma_two.element().trim() == "H" )
	{
	    if ((distance >  0.90 ) && ( distance < 1.10 )) return true; // N-H
	    else return false;
	}
    }
    else if ( sugar_parent_molecule->atom(ma_one).element().trim() == "O" )
    {
	if ( ma_two.element().trim() == "C" )
	    if ((distance >  1.16 ) && ( distance < 1.50 )) return true; // O-C or O=C
	    else return false;
        else if ( ma_two.element().trim() == "H" )
	{
	    if ((distance >  0.88 ) && ( distance < 1.04 )) return true; // O-H
	    else return false;
	}
    }
    else if ( sugar_parent_molecule->atom(ma_one).element().trim() == "S" )
	if ( ma_two.element().trim() == "C" )
	    if ((distance >  1.26 ) && ( distance < 2.00 )) return true; // O-C or O=C
	    else return false;
        else if ( ma_two.element().trim() == "H" )
	{
	    if ((distance >  0.78 ) && ( distance < 1.44 )) return true; // O-H
	    else return false;
	}
    else if (( distance > 1.2) && (distance < 1.8)) return true; // unknown bond
	else return false;

    return false; // in case we haven't found any match
}


/*! Internal function for getting the alternate conformation code
	\param ma A clipper::MAtom object
	\return A character containing the code
*/

const char MSugar::get_altconf(const clipper::MAtom& ma) const
{
	clipper::String identifier = ma.id();
	if (identifier.size() > 5) return identifier[5];
	else return ' ';                                    // The alternate conformation code is the fifth character in the complete identificator.
}                                                           // We will return a blank character if there is no code present or if it is, but is blank


/*! Internal function for getting the alternate conformation code
*/

bool MSugar::examine_ring()
{
    int i;

    { // first element: angle ( [n-1] - O - anomeric C ), bond ( [n-1] - O ), torsion ( [n-1] - O - anomeric C - next C )
	clipper::Vec3<clipper::ftype> vec_a(sugar_ring_elements[sugar_ring_elements.size()-1].coord_orth().x() - sugar_ring_elements[0].coord_orth().x(),
					    sugar_ring_elements[sugar_ring_elements.size()-1].coord_orth().y() - sugar_ring_elements[0].coord_orth().y(),
					    sugar_ring_elements[sugar_ring_elements.size()-1].coord_orth().z() - sugar_ring_elements[0].coord_orth().z() );

	clipper::Vec3<clipper::ftype> vec_b(sugar_ring_elements[1].coord_orth().x() - sugar_ring_elements[0].coord_orth().x(),
					    sugar_ring_elements[1].coord_orth().y() - sugar_ring_elements[0].coord_orth().y(),
					    sugar_ring_elements[1].coord_orth().z() - sugar_ring_elements[0].coord_orth().z() );

	clipper::ftype norm_a = sqrt(pow(vec_a[0],2) + pow(vec_a[1],2) + pow(vec_a[2],2));
	clipper::ftype norm_b = sqrt(pow(vec_b[0],2) + pow(vec_b[1],2) + pow(vec_b[2],2));
	
        sugar_ring_angles.push_back( ( clipper::Util::rad2d(acos(clipper::Vec3<clipper::ftype>::dot(vec_a, vec_b)/ (norm_a *norm_b) ))) );
	sugar_ring_bonds.push_back( clipper::Coord_orth::length(sugar_ring_elements[sugar_ring_elements.size()-1].coord_orth(),sugar_ring_elements[0].coord_orth()) ); 
		
	sugar_ring_torsion.push_back ( clipper::Util::rad2d ( clipper::Coord_orth::torsion ( sugar_ring_elements[sugar_ring_elements.size()-1].coord_orth(),
													      sugar_ring_elements[0].coord_orth(),
													      sugar_ring_elements[1].coord_orth(),
													      sugar_ring_elements[2].coord_orth() ) ) );
    }

    for (i = 1 ; i < sugar_ring_elements.size() -1 ; i++ )  // calculate bond angles (all), lengths (all) and torsions (minus last one)
    {
        clipper::Vec3<clipper::ftype> vec_a(sugar_ring_elements[i-1].coord_orth().x() - sugar_ring_elements[i].coord_orth().x(),
                                                                                        sugar_ring_elements[i-1].coord_orth().y() - sugar_ring_elements[i].coord_orth().y(),
                                                                                        sugar_ring_elements[i-1].coord_orth().z() - sugar_ring_elements[i].coord_orth().z() );

        clipper::Vec3<clipper::ftype> vec_b(sugar_ring_elements[i+1].coord_orth().x() - sugar_ring_elements[i].coord_orth().x(),
                                                                                        sugar_ring_elements[i+1].coord_orth().y() - sugar_ring_elements[i].coord_orth().y(),
                                                                                        sugar_ring_elements[i+1].coord_orth().z() - sugar_ring_elements[i].coord_orth().z() );

        clipper::ftype norm_a = sqrt(pow(vec_a[0],2) + pow(vec_a[1],2) + pow(vec_a[2],2));
        clipper::ftype norm_b = sqrt(pow(vec_b[0],2) + pow(vec_b[1],2) + pow(vec_b[2],2));

        sugar_ring_angles.push_back( clipper::Util::rad2d((acos(clipper::Vec3<clipper::ftype>::dot(vec_a, vec_b)/ (norm_a * norm_b) ))) );
        sugar_ring_bonds.push_back( clipper::Coord_orth::length(sugar_ring_elements[i+1].coord_orth(),sugar_ring_elements[i].coord_orth()) ); 

        if ( i != (sugar_ring_elements.size()-2 ) )
	    sugar_ring_torsion.push_back (clipper::Util::rad2d(clipper::Coord_orth::torsion(sugar_ring_elements[i-1].coord_orth(),
                                                                                            sugar_ring_elements[i].coord_orth(),
                                                                                            sugar_ring_elements[i+1].coord_orth(),
                                                                                            sugar_ring_elements[i+2].coord_orth() ) ));
        else 
	    sugar_ring_torsion.push_back ((clipper::Util::rad2d( clipper::Coord_orth::torsion(sugar_ring_elements[i-1].coord_orth(),
                                                                                              sugar_ring_elements[i].coord_orth(),
                                                                                              sugar_ring_elements[i+1].coord_orth(),
                                                                                              sugar_ring_elements[0].coord_orth() )) ))  ;
    }
    
    // calculate last torsion angle
    
    sugar_ring_torsion.push_back (clipper::Util::rad2d(clipper::Coord_orth::torsion(sugar_ring_elements[sugar_ring_elements.size()-2].coord_orth(),
                                                                                    sugar_ring_elements[sugar_ring_elements.size()-1].coord_orth(),
                                                                                    sugar_ring_elements[0].coord_orth(),
                                                                                    sugar_ring_elements[1].coord_orth() ) ));


    clipper::Vec3<clipper::ftype> vec_a(sugar_ring_elements[i-1].coord_orth().x() - sugar_ring_elements[i].coord_orth().x(),
                                                                                    sugar_ring_elements[i-1].coord_orth().y() - sugar_ring_elements[i].coord_orth().y(),
                                                                                    sugar_ring_elements[i-1].coord_orth().z() - sugar_ring_elements[i].coord_orth().z() );

    clipper::Vec3<clipper::ftype> vec_b(sugar_ring_elements[0].coord_orth().x() - sugar_ring_elements[i].coord_orth().x(),
                                                                                  sugar_ring_elements[0].coord_orth().y() - sugar_ring_elements[i].coord_orth().y(),
                                                                                  sugar_ring_elements[0].coord_orth().z() - sugar_ring_elements[i].coord_orth().z() );

    clipper::ftype norm_a = sqrt(pow(vec_a[0],2) + pow(vec_a[1],2) + pow(vec_a[2],2));
    clipper::ftype norm_b = sqrt(pow(vec_b[0],2) + pow(vec_b[1],2) + pow(vec_b[2],2));

    sugar_ring_angles.push_back( clipper::Util::rad2d(acos(clipper::Vec3<clipper::ftype>::dot(vec_a, vec_b)/ (norm_a *norm_b) )));
    sugar_ring_bonds.push_back( clipper::Coord_orth::length(sugar_ring_elements[0].coord_orth(),sugar_ring_elements[i].coord_orth()) ); 

    clipper::ftype rmsd_bonds, rmsd_angles;
    rmsd_bonds = 0.0;
    rmsd_angles = 0.0;

    for (int j = 0 ; j < sugar_ring_bonds.size() ; j++ )
    {
        if (( j == 0) || (j == sugar_ring_bonds.size()-1)) rmsd_bonds = rmsd_bonds + pow((sugar_ring_bonds[j] - 1.430), 2);
	else rmsd_bonds = rmsd_bonds + pow((sugar_ring_bonds[j] - 1.530), 2);  
		
	if (( j == 0) || (j == sugar_ring_angles.size()-1)) rmsd_angles = rmsd_angles + pow((sugar_ring_angles[j] - 112.0), 2);
	else rmsd_angles = rmsd_angles + pow((sugar_ring_angles[j] - 109.0), 2);  
    }

    rmsd_bonds = rmsd_bonds / sugar_ring_bonds.size();
    rmsd_bonds = sqrt(rmsd_bonds);

    rmsd_angles = rmsd_angles / sugar_ring_angles.size();
    rmsd_angles = sqrt(rmsd_angles);
	
    sugar_ring_angle_rmsd = rmsd_angles;
    sugar_ring_bond_rmsd = rmsd_bonds;

    for (i = 0 ; i < sugar_ring_elements.size() -1 ; i++)
        if (!bonded(sugar_ring_elements[i], sugar_ring_elements[i+1])) 
	    return false; 

        if (!bonded(sugar_ring_elements[i], sugar_ring_elements[0])) 
	    return false;
        else 
	    return true;
}

/*! Internal function for checking whether two atoms are bonded or not
 * 	\param ma_one A clipper::MAtom object
 * 		\param ma_two A clipper::MAtom object
 * 			\return A boolean value with the obvious answer
 * 			*/

bool MSugar::bonded(const clipper::MAtom& ma_one, const clipper::MAtom& ma_two) const
{
    clipper::ftype distance = clipper::Coord_orth::length( ma_one.coord_orth(), ma_two.coord_orth() );
	
    if ( ma_one.element().trim() == "C" )
    {
        if ( ma_two.element().trim() == "C" )
	    if ((distance >  1.18 ) && ( distance < 1.60 )) return true; // C-C or C=C
	    else return false;
	else if ( ma_two.element().trim() == "N" )
	    if ((distance >  1.24 ) && ( distance < 1.52 )) return true; // C-N or C=N
	    else return false;
        else if ( ma_two.element().trim() == "O" )
            if ((distance >  1.16 ) && ( distance < 1.50 )) return true; // C-O or C=O
            else return false;
        else if ( ma_two.element().trim() == "S" )
            if (( distance > 1.2 ) && ( distance < 2.0 )) return true; // C-S or C=S
            else return false;
        else if ( ma_two.element().trim() == "H" )
            if ((distance >  0.96 ) && ( distance < 1.14 )) return true; // C-H				
            else return false;
    }
    else if ( ma_one.element().trim() == "N" )
    {
        if ( ma_two.element().trim() == "C" )
            if ((distance >  1.24 ) && ( distance < 1.52 )) return true; // N-C or N=C
            else return false;
        else if ( ma_two.element().trim() == "H" )
            if ((distance >  0.90 ) && ( distance < 1.10 )) return true; // N-H
            else return false;
    }
    else if ( ma_one.element().trim() == "O" )
    {
        if ( ma_two.element().trim() == "C" )
            if ((distance >  1.16 ) && ( distance < 1.50 )) return true; // O-C or O=C
            else return false;
        else if ( ma_two.element().trim() == "H" )
            if ((distance >  0.88 ) && ( distance < 1.04 )) return true; // O-H
            else return false;
    }
    else if ( ma_one.element().trim() == "S" )
    {
        if ( ma_two.element().trim() == "C" )
            if ((distance >  1.26 ) && ( distance < 2.00 )) return true; // O-C or O=C
            else return false;
        else if ( ma_two.element().trim() == "H" )
            if ((distance >  0.78 ) && ( distance < 1.24 )) return true; // O-H
            else return false;
    }
    else if (( distance > 1.2) && (distance < 2.0)) return true; // unknown bond
    else return false;

    return false; // in case we haven't found any match
}


MDisaccharide::MDisaccharide ( clipper::MiniMol& mmol, const clipper::MAtomNonBond& manb, clipper::MMonomer& mm )
{
    int index = search_disaccharides ( mm.type().c_str() ); // we know beforehand that this is a known disaccharide, no need to re-check
        
    clipper::data::sugar_database_entry val_string_one = clipper::data::disaccharide_database[index].sugar_one;
    clipper::data::sugar_database_entry val_string_two = clipper::data::disaccharide_database[index].sugar_two;

    sugar_one = clipper::MSugar ( mmol, mm, manb, val_string_one );
    sugar_two = clipper::MSugar ( mmol, mm, manb, val_string_two ); 
    
    sugar_one.set_type ( clipper::String( sugar_one.type().trim() + "[" + clipper::data::disaccharide_database[index].sugar_one.name_short + "]" ));
    sugar_two.set_type ( clipper::String( sugar_two.type().trim() + "[" + clipper::data::disaccharide_database[index].sugar_two.name_short + "]" ));
}


MGlycan::MGlycan ( clipper::String chain, const clipper::MMonomer& root_aa, clipper::MSugar& root_sugar )
{
    root.first = root_aa;
    t_node first_node;
    sugars.push_back ( root_sugar ); 
    first_node.sugar = sugars.back();
    node_list.push_back ( first_node );
    root.second = node_list[0];
    this->chain = chain;
}


clipper::String MGlycan::carbname_of(clipper::String name) const
{
    clipper::String new_name = name;
    
    if ( name == "NAG" ) new_name = "GlcNac";
    else if ( name == "NGA" ) new_name = "GalNac";
    else if ( name == "MAN" ) new_name = "Man"   ;
    else if ( name == "BMA" ) new_name = "Man"   ;
    else if ( name == "SIA" ) new_name = "SiaAc" ;
    else if ( name == "GAL" ) new_name = "Gal"   ;
    else if ( name == "GLA" ) new_name = "Gal"   ;
    else if ( name == "FUC" ) new_name = "Fuc"   ;
    else if ( name == "FCB" ) new_name = "Fuc"   ; 
    else if ( name == "GLC" ) new_name = "Glc"   ;
    else if ( name == "BGC" ) new_name = "Glc"   ;
    else if ( name == "DAN" ) new_name = "NeuAc"   ; 

    return new_name;
}



clipper::String MGlycan::print_linear ( const bool print_info, const bool html_format, const bool translate )
{
    clipper::String buffer;
    
    if ( html_format ) buffer.insert ( 0, "</sub>" ); 
    buffer.insert ( 0, root.first.id().trim() );
    if ( html_format ) buffer.insert ( 0, "<sub>" );
    
    buffer.insert( 0, this->get_root().type().c_str() );
    
    html_format ? buffer.insert ( 0, "&#8722;" ) : buffer.insert( 0, "-" );
    clipper::String anomer;
    if ( html_format ) root.second.sugar.anomer() == "alpha" ? anomer = "&#945;" : anomer = "&#946;";
    else root.second.sugar.anomer() == "alpha" ? anomer = "a" : anomer = "b"; 
    
    html_format ? buffer.insert ( 0, anomer ) : buffer.insert ( 0, anomer );
    html_format ? buffer.insert ( 0, " &#8592; " ) : buffer.insert( 0, "-" );
        
    clipper::MSugar msug = node_list[0].sugar;    

    if ( print_info )
    {
        if ( html_format ) buffer.insert ( 0, "</sub>" ); 
        buffer.insert ( 0, msug.id().trim() );
        if ( html_format ) buffer.insert ( 0, "<sub>" ); 
    }
    if ( html_format ) buffer.insert ( 0, "</span>" );     
    translate ? buffer.insert ( 0, carbname_of ( msug.type().trim() ).c_str() ) : buffer.insert( 0, msug.type().c_str() ); 
    if ( html_format ) 
    {
        buffer.insert ( 0, "\">" );
//        std::ostringstream rscc;
//        rscc <<  msug.get_rscc();
//        buffer.insert ( 0, rscc.str() );
        buffer.insert ( 0, msug.type() + " " + msug.id().trim() + " in " + msug.conformation_name() );
        buffer.insert ( 0, "<span title=\"" );      
    }

    if ( node_list.size() < 2 ) return buffer;
    else
    {
        if ( html_format ) 
            node_list[0].connections[0].type == "alpha" ? anomer = "&#945;" : anomer = "&#946;";
        else 
            node_list[0].connections[0].type == "alpha" ? anomer = "a" : anomer = "b";
        
        std::ostringstream s;
        s << node_list[0].connections[0].index;
        
        html_format ? buffer.insert ( 0, "&#8722;" ) : buffer.insert( 0, "-" );
        buffer.insert( 0, s.str() );
        buffer.insert( 0, anomer );
        html_format ? buffer.insert ( 0, "&#8592;" ) : buffer.insert( 0, "-" );
    }   

    bool branching = false;
    t_node branched_from;

    for ( int i = 1 ; i < node_list.size() ; i++ )
    {
        msug = node_list[i].sugar;    
        
        if ( print_info )
        {
            if ( html_format ) 
                buffer.insert ( 0, "</sub>" ); 
            
            buffer.insert ( 0, msug.id().trim() );
            
            if ( html_format ) 
                buffer.insert ( 0, "<sub>" ); 
        }
    
        if ( html_format ) buffer.insert ( 0, "</span>" );      
        translate ? buffer.insert ( 0, carbname_of ( msug.type().trim() ).c_str() ) : buffer.insert( 0, msug.type().c_str() ); 
        if ( html_format ) 
        {
            buffer.insert ( 0, "\">" );
            buffer.insert ( 0, msug.type() + " " + msug.id().trim() + " in " + msug.conformation_name() );
            buffer.insert ( 0, "<span title=\"" );      
        }

        if (( node_list[i].connections.size() == 0 ) && branching )
        {
            if ( html_format ) 
                branched_from.connections[1].type == "alpha" ? anomer = "&#945;" : anomer = "&#946;";
            else 
                branched_from.connections[1].type == "alpha" ? anomer = "a" : anomer = "b";
            
            std::ostringstream s;
            s << branched_from.connections[1].index;
            
            branching = false;
            buffer.insert ( 0, "(" );
            html_format ? buffer.insert ( 0, "&#8722;" ) : buffer.insert ( 0, "-" );
            buffer.insert( 0, s.str() );
            buffer.insert( 0, anomer );
            html_format ? buffer.insert ( 0, "&#8592;" ) : buffer.insert( 0, "-" );        
        }    
        else if ( node_list[i].connections.size() > 1 )
        {
            branching = true;
            buffer.insert ( 0, ")" );
            branched_from = node_list[i];
        }
        
        if ( node_list[i].connections.size() > 0 )
        {
            if ( html_format ) 
                node_list[i].connections[0].type == "alpha" ? anomer = "&#945;" : anomer = "&#946;";
            else 
                node_list[i].connections[0].type == "alpha" ? anomer = "a" : anomer = "b";
            
            std::ostringstream s;
            s << node_list[i].connections[0].index;
            
            html_format ? buffer.insert ( 0, "&#8722;" ) : buffer.insert( 0, "-" );
            buffer.insert( 0, s.str() );
            buffer.insert( 0, anomer );
            html_format ? buffer.insert ( 0, "&#8592;" ) : buffer.insert( 0, "-" );
        }
        
    }
    return buffer;
}   


bool MGlycan::link_sugars ( int link, clipper::MSugar& first_sugar, clipper::MSugar& next_sugar )
{
    int index = 0;
    bool found;

    for ( int i = 0 ; i < node_list.size() ; i++ )
        if ( node_list[i].sugar.id().c_str() == first_sugar.id().c_str() )
        {
            found = true;
            index = i;
            break;
        }
    
    if (!found) return true;

    t_node new_node;
    new_node.sugar = next_sugar; // create a new node with the next sugar

    t_connection new_connection;  // create a new connection pointing to the next sugar, to be added to the previous one
    new_connection.index = link; // temporary fix to enable early testing, substitute with an elaborate thing
    new_connection.type = next_sugar.anomer();
    new_connection.node = &new_node; 
    sugars.push_back ( next_sugar ); 
    node_list[index].connections.push_back ( new_connection ); // add the new connection to the previous node
    node_list.push_back ( new_node ); // add the new sugar to the node list
    
    return false;
}


/*! Internal function for getting the alternate conformation code
	\param ma A clipper::MAtom object
	\return A character containing the code
*/

const char MGlycology::get_altconf(const clipper::MAtom& ma) const
{
	clipper::String identifier = ma.id();
	if (identifier.size() > 5) return identifier[5];
	else return ' ';                                    // The alternate conformation code is the fifth character in the complete identificator.
}                                                           // We will return a blank character if there is no code present or if it is, but is blank


void MGlycology::extend_tree ( clipper::MGlycan& mg, clipper::MSugar& msug )
{
    std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > contacts = get_contacts ( msug );

    const clipper::MiniMol& tmpmol = *this->mmol;

    for (int i = 0 ; i < contacts.size() ; i++ )
    {
        clipper::MSugar tmpsug = clipper::MSugar ( *this->mmol, tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()], *this->manb );
        mg.link_sugars ( parse_order ( contacts[i].first.id() ), msug, tmpsug );
        //std::cout << "id: " << contacts[i].first.id() << " parsed: " << parse_order(contacts[i].first.id()) << std::endl;
        extend_tree ( mg, tmpsug );
    }
}

MGlycology::MGlycology ( const clipper::MiniMol& mmol )
{
    const clipper::MAtomNonBond nb = MAtomNonBond ( mmol, 1.0 );
    MGlycology( mmol, nb );
}


MGlycology::MGlycology ( const clipper::MiniMol& mmol, const clipper::MAtomNonBond& manb)
{

    this->manb = &manb;
    this->mmol = &mmol;

    std::vector < clipper::MMonomer > potential_n_roots;
    std::vector < clipper::MMonomer > potential_o_roots;

    for ( int pol = 0; pol < mmol.size() ; pol++ )
        for ( int mon = 0 ; mon < mmol[pol].size() ; mon++ )
        {
            // Will need to keep this list up to date with the latest discoveries
            // To do: include check on other ligands, such as lipids (ceramide O-glycosylation, for instance)
            
            if ( mmol[pol][mon].type() == "ASN" ) potential_n_roots.push_back ( mmol[pol][mon] ); // n-linked GlcNAc ?
            else if ( mmol[pol][mon].type() == "ARG" ) potential_n_roots.push_back ( mmol[pol][mon] ); // Arginine rhamnosylation?
            else if ( mmol[pol][mon].type() == "THR" ) potential_o_roots.push_back ( mmol[pol][mon] ); // o-linked stuff ?
            else if ( mmol[pol][mon].type() == "SER" ) potential_o_roots.push_back ( mmol[pol][mon] );
            else if ( mmol[pol][mon].type() == "LYS" ) potential_o_roots.push_back ( mmol[pol][mon] );
            else if ( mmol[pol][mon].type() == "TYR" ) potential_o_roots.push_back ( mmol[pol][mon] );
            else if ( mmol[pol][mon].type() == "PRO" ) potential_o_roots.push_back ( mmol[pol][mon] );
        }

    for ( int i = 0 ; i < potential_n_roots.size() ; i++ )  // create n-glycan roots with first sugar
    {        
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_n_roots[i] ) ;    
        
        for ( int j = 0 ; j < linked.size() ; j++ )
        {
            const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()]; 
            
            if ( tmpmon.type().trim() == "NAG" ) 
            { 
                if ( clipper::MSugar::search_database( tmpmon.type().c_str() ) )
                {
                    clipper::MSugar sugar( mmol, tmpmon, manb );
                    list_of_sugars.push_back ( sugar );
                    clipper::MGlycan mg = clipper::MGlycan ( clipper::MPolymer::id_tidy(mmol[linked[j].second.polymer()].id()), 
                                                             potential_n_roots[i], list_of_sugars[list_of_sugars.size()-1] );
                    mg.set_kind_of_glycan ( "n-glycan" );
                    list_of_glycans.push_back ( mg );
                    break;
                }
            }
        }
    }


    for ( int i = 0 ; i < potential_o_roots.size() ; i++ )  // create o-glycan roots with first sugar
    {        
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_o_roots[i] ) ;    
        
        for ( int j = 0 ; j < linked.size() ; j++ )
        {
            const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()]; 
            
            if (( tmpmon.type().trim() == "NGA" ) || (tmpmon.type().trim() == "FUC" ) || (tmpmon.type().trim() == "RAM" ) || (tmpmon.type().trim() == "BGC" ) || (tmpmon.type().trim() == "BMA" ) || (tmpmon.type().trim() == "NAG" ))
            { 
                if ( clipper::MSugar::search_database( tmpmon.type().c_str() ) )
                {
                    clipper::MSugar sugar( mmol, tmpmon, manb );
                    list_of_sugars.push_back ( sugar );
                    clipper::MGlycan mg = clipper::MGlycan ( clipper::MPolymer::id_tidy(mmol[linked[j].second.polymer()].id()), 
                                                             potential_o_roots[i], list_of_sugars[list_of_sugars.size()-1] );
                    mg.set_kind_of_glycan ( "o-glycan" );
                    list_of_glycans.push_back ( mg );
                    break;
                }
            }
        }
    }

    
    for ( int i = 0 ; i < list_of_glycans.size() ; i++ )
    {
        std::vector < clipper::MSugar > sugar_list = list_of_glycans[i].get_sugars();
        extend_tree ( list_of_glycans[i] , sugar_list[0] );
    }
    

}

const std::vector < std::pair< clipper::MAtom, clipper::MAtomIndexSymmetry > > MGlycology::get_contacts ( const clipper::MMonomer& mm )
{
    // mm can be aminoacid (typically ASN) or sugar: ND2, O2, O3, O4, O6  

    std::vector < clipper::MAtom > candidates; // look for the possibilities
    std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > tmpresults;

    if ( mm.type().trim() == "ASN" ) 
    {
        int id = mm.lookup ( "ND2", clipper::MM::ANY );
        
        if ( id != -1 )
            candidates.push_back ( mm[id] );
        else return tmpresults; // empty result

    }
    else if ( mm.type().trim() == "SER" )
    {
        int id = mm.lookup ( "OG", clipper::MM::ANY );
        
        if ( id != -1 )
            candidates.push_back ( mm[id] );
        else return tmpresults; // empty result

    }
    else if ( mm.type().trim() == "THR" )
    {
        int id = mm.lookup ( "OG1", clipper::MM::ANY );
        
        if ( id != -1 )
            candidates.push_back ( mm[id] );
        else return tmpresults; // empty result

    }
    else 
    {

        int id = mm.lookup ( "O2", clipper::MM::ANY );
        
        if ( id != -1 )
            candidates.push_back ( mm[id] );
        
        id = mm.lookup ( "O3", clipper::MM::ANY );
        
        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O4", clipper::MM::ANY );
        
        if ( id != -1 )
            candidates.push_back ( mm[id] );
        
        id = mm.lookup ( "O6", clipper::MM::ANY );
        
        if ( id != -1 )
            candidates.push_back ( mm[id] );
    }
 
    if ( candidates.size() == 0 ) 
        return tmpresults;  // empty result
    else 
    {
        const clipper::MiniMol& tmpmol = *mmol;

        for ( int i = 0 ; i < candidates.size() ; i++ )
        {
            std::vector < clipper::MAtomIndexSymmetry > contacts = this->manb->atoms_near ( candidates[i].coord_orth(), 2.0 );

            for (int j = 0 ; j < contacts.size() ; j++ ) 
            {
                if ((tmpmol[contacts[j].polymer()][contacts[j].monomer()].id().trim() != mm.id().trim())
                && ( clipper::Coord_orth::length ( tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()].coord_orth(), candidates[i].coord_orth() ) < 1.61 ))
                {   
                    if (( get_altconf(tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()]) == 'A')
                    || (get_altconf(tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()]) == ' '))
                    {    
                        std::pair < clipper::MAtom , clipper::MAtomIndexSymmetry > link_tmp;
                        link_tmp.first = candidates[i];                   
                        link_tmp.second = contacts[j];
                        tmpresults.push_back ( link_tmp );
                    }
                }  
            }
        }
        return tmpresults;
    }
}







