/*! \file clipper-glyco.cpp
  Implementation file for sugar data */

// clipper-glyco.cpp: a set of tools for handling sugars
// version  0.9.1
// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York




#include "clipper-glyco.h"

#define DBG std::cout << "[" << __FUNCTION__ << "] - "



using namespace clipper;


///////////////////////// MSugar ///////////////////////////////

/*! Constructor: Empty constructor for later initialisation. */
MSugar::MSugar( )
{

}

/*! Constructor: create a new sugar object from a standard MMonomer.
	  If reference data for the sugar cannot be found in the database, the members of the ring will be determined using a recursive version of Fleury's algorithm for finding eulerian cycles in undirected graphs
	\param ml A MiniMol object containing this and neighbouring sugars, which may affect the stereochemistry
	\param mm A MMonomer object that will be extended into a sugar
	\return The MSugar object, which will contain cremer-pople parameters, conformation code, anomer, handedness and linkage information */

MSugar::MSugar(const clipper::MiniMol& ml, const clipper::String chainID, const clipper::MMonomer& mm, bool& debug_output, char alt_conf)
{

	// we calculate the non-bond object first, then continue with normal creation
    this->debug_output = debug_output;
	const clipper::MAtomNonBond& nb = MAtomNonBond (ml, 5.0);
	MSugar(ml, chainID, mm, nb, debug_output, alt_conf);
}
/*
ok_with_ring() = True
ok_with_bonds_rmsd() = True
ok_with_angles_rmsd() = True
ok_with_anomer() = True
ok_with_chirality() = True
ok_with_conformation() = True

In MSugar steps where if statement that returns incomplete MSugar, set these to false - the way that I did with cremer_pople_params and sugar_conformation_code
*/

/*! Constructor: create a new sugar object from a standard MMonomer.
	If reference data for the sugar cannot be found in the database, the members of the ring will be determined using a recursive version of Fleury's algorithm for finding eulerian cycles in undirected graphs
	\param ml A MiniMol object containing this and neighbouring sugars
	\param mm An MMonomer object that will be extended into a sugar
	\param nb An MAtomNonBond object to be used for the determination of the stereochemistry
	\return The MSugar object, which will contain cremer-pople parameters, conformation code, anomer, handedness and linkage information */

MSugar::MSugar(const clipper::MiniMol& ml, const clipper::String chainID, const clipper::MMonomer& mm, const clipper::MAtomNonBond& nb, bool& debug_output, char alt_conf )
{

    copy(mm,clipper::MM::COPY_MPC);	// import_data from MMonomer

    this->sugar_chain_id = chainID;
    this->debug_output = debug_output;
    this->sugar_supported = true;
    this->sugar_parent_molecule = &ml;
    this->sugar_parent_molecule_nonbond = &nb; // store pointers
    this->sugar_index = db_not_checked;
    this->sugar_index = 9999; // default value for "not found in database".
    this->sugar_alternate_confcode = " "; // initially, we would like to suppose this
    this->sugar_context = "";
    this->sugar_pdb_id = mm.id();
    this->sugar_seqnum = mm.seqnum();
    this->sugar_name_short = mm.type().trim();

    if(debug_output)
    {
        std::cout << std::endl ;
        DBG << "looking for " << this->id() << " " << this->type().trim() << " on the database..." << std::endl;
        // DBG << "with the coords of = " << this->find("C1", MM::ANY).coord_orth().format() << " on the database..." << std::endl;
        // alt_conf != ' ' ? DBG << "Alternate locator supplied: " << alt_conf << std::endl : true;
    }

    if(debug_output)
    {
        std::cout << "Size of (*this).size() " << (*this).size() << std::endl;
        for(int i = 0; i < (*this).size(); i++)
        {
            std::cout << "Atom ID: (*this)[" << i << "].id " << (*this)[i].id() << std::endl;
        }
    }

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
        if(debug_output)
        {
            DBG << "found it! " << std::endl;
        }

        std::vector<clipper::String> buffer = clipper::data::sugar_database[sugar_index].ring_atoms.trim().split(" ");

        this->sugar_name_full = clipper::data::sugar_database[sugar_index].name_long;

        for (int i=0 ; i < buffer.size() ; i++)
        {
            int index_atom = 0;

            if ( alt_conf != ' ' ) // alternate conformations are present
            {
                alt_conf == 'A' ? sugar_alternate_confcode = " :A" : sugar_alternate_confcode = " :B";

                index_atom = this->lookup(buffer[i].trim()+sugar_alternate_confcode,clipper::MM::UNIQUE);

                if(debug_output)
                {
                    DBG << "index_atom in line 137" << index_atom << std::endl;
                }

                if (index_atom == -1) // we've tried A and B and it still fails... so we're going to give up for now
                {
                    this->sugar_supported = false;
                    this->sugar_sane = false;
                    this->sugar_denomination = "    unsupported    ";
                    this->sugar_anomer = "X";
                    this->sugar_handedness = "X";
                    
                    for(int j = 0; j < 3; j++)
                        this->sugar_cremer_pople_params.push_back(-1);

                    this->sugar_conformation = 0;

                    return;
                }
            }
            else
            {

                index_atom = this->lookup(buffer[i],clipper::MM::ANY);


                if(debug_output)
                {
                    DBG << "index_atom in line 164 = " << index_atom << " value of buffer[" << i << "] =" << buffer[i] << std::endl;
                }

                if (index_atom == -1)
                {
                    this->sugar_supported = false;
                    this->sugar_sane = false; // we don't support cyclic sugars with less or more than 5-6 ring atoms
                    this->sugar_denomination = "    unsupported    ";
                    this->sugar_anomer = "X";
                    this->sugar_handedness = "X";

                    for(int j = 0; j < 3; j++)
                       this->sugar_cremer_pople_params.push_back(-1);

                    this->sugar_conformation = 0;

                    return;
                }
            }
            if(debug_output)
            {
                DBG << "trying to push (*this)[index_atom] in line 182" << (*this)[index_atom].id() << std::endl;
            }
            sugar_ring_elements.push_back((*this)[index_atom]);
        }
    }
    else
    {
        if(debug_output)
        {
            std::cout << std::endl ;
            DBG << "Unable to find the " << this->id() << " in the database. Initiating ringMembers() function" << std::endl;
        }
        this->sugar_ring_elements = this->ringMembers();
    }

	if (this->sugar_ring_elements.size() == 5)
	{
		this->cremerPople_furanose(*this->sugar_parent_molecule, mm);
		this->sugar_conformation = conformationFuranose(this->sugar_cremer_pople_params[1]);

		if(debug_output)
        {
			DBG << "After checking the conformation..." << std::endl;
		}
	}
	else if (this->sugar_ring_elements.size() == 6)
	{
		this->cremerPople_pyranose(*this->sugar_parent_molecule, mm);
		this->sugar_conformation = conformationPyranose(this->sugar_cremer_pople_params[1], this->sugar_cremer_pople_params[2]);

        if(debug_output)
        {
			DBG << "After checking the conformation..." << std::endl;
		}
	}
	else
	{
		this->sugar_supported = false;
		this->sugar_sane = false; // we don't support cyclic sugars with less or more than 5-6 ring atoms
		this->sugar_denomination = "    unsupported    ";
		this->sugar_anomer = "X";
		this->sugar_handedness = "X";

        for(int j = 0; j < 3; j++)
            this->sugar_cremer_pople_params.push_back(-1);
        
        this->sugar_conformation = 0;
        
	}

    if (sugar_ring_elements[1].name().trim().find("C1") != std::string::npos)
    {
        this->sugar_type = "aldose";
        this->sugar_denomination = "aldo";
    }
    else
    {
        this->sugar_type = "ketose";
        this->sugar_denomination = "keto";
    }

	if (sugar_ring_elements.size() == 5)
        this->sugar_denomination = this->sugar_denomination + "furanose";
	else
        this->sugar_denomination = this->sugar_denomination + "pyranose";

	this->sugar_denomination = clipper::String( this->sugar_anomer + "-" + this->sugar_handedness + "-" + this->sugar_denomination );

	// sanity check:

	this->sugar_sane = false;

	if(debug_output)
    {
		DBG << "Just before examining the ring..." << std::endl;
	}


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
	    else
            sugar_diag_chirality = false;


	    if ( ( ( sugar_anomer == "alpha") && ( clipper::data::sugar_database[sugar_index].anomer != "B" ) )
		     || ( ( sugar_anomer == "beta") && ( clipper::data::sugar_database[sugar_index].anomer != "A" ) ) )
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
            if (( puckering_amplitude() > ref_puckering - 0.18 ) && (puckering_amplitude() < ref_puckering + 0.15 ))
                sugar_diag_puckering = true;
            else
                sugar_diag_puckering = false;

            if ( sugar_ring_bond_rmsd < ( ref_bonds_rmsd + 0.039 ) )
                sugar_diag_bonds_rmsd = true;
            else
                sugar_diag_bonds_rmsd = false;

            if ( sugar_ring_angle_rmsd < ( ref_angles_rmsd + 3.0 ) )
                sugar_diag_angles_rmsd = true;
            else
                sugar_diag_angles_rmsd = false;
        }
        else
        {
            if (( puckering_amplitude() > 0.9 ) || ( puckering_amplitude() < 0.42 ))
                sugar_diag_puckering = false;
            else
                sugar_diag_puckering = true;

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
            else
                sugar_diag_bonds_rmsd = false;
	    }
	    else
	    {
	        if (sugar_ring_bond_rmsd < 0.035 )
                sugar_diag_bonds_rmsd = true;
            else
                sugar_diag_bonds_rmsd = false;
	    }

	    if (sugar_ring_elements.size() ==  5)
	    {
	        if ((sugar_ring_angle_rmsd > 4.0 ) && (sugar_ring_angle_rmsd < 8.0))
                sugar_diag_angles_rmsd = true;
            else
                sugar_diag_angles_rmsd = false;
	    }
	    else
	    {
	        if (sugar_ring_angle_rmsd < 4.0 )
                sugar_diag_angles_rmsd = true;
	        else
                sugar_diag_angles_rmsd = false;
	    }

        if ( ( conformation_name() == "4c1" ) && ( sugar_handedness != "L" ))
            sugar_diag_conformation = true;
        else if ( ( conformation_name() == "1c4" ) && ( sugar_handedness != "D" ))
            sugar_diag_conformation = true;
        else
            sugar_diag_conformation = false;

        if (( puckering_amplitude() > 0.9 ) || ( puckering_amplitude() < 0.42 ))
            sugar_diag_puckering = false;
        else
            sugar_diag_puckering = true;

	    if ( sugar_diag_puckering && sugar_diag_anomer && sugar_diag_chirality && sugar_diag_ring )
            sugar_sane = true;

	}

	if(debug_output)
    {
	    DBG << "Just after examining the ring, exiting the constructor, good job!" << std::endl;
	}
}




/*! Constructor: create a new sugar object from a standard clipper::MMonomer, providing external validation data.
	If reference data for the sugar cannot be found in the database, the members of the ring will be determined using a recursive version of Fleury's algorithm for finding eulerian cycles in undirected graphs
	\param ml A MiniMol object containing this and neighbouring sugars
	\param mm An MMonomer object that will be extended into a sugar
	\param nb An MAtomNonBond object to be used for the determination of the stereochemistry
    \param validation_data A data structure containing validation data that overrides whatever there is in the database
	\return The MSugar object, which will contain cremer-pople parameters, conformation code, anomer, handedness and linkage information */

MSugar::MSugar(const clipper::MiniMol& ml, const clipper::String chainID, const clipper::MMonomer& mm, const clipper::MAtomNonBond& nb, clipper::data::sugar_database_entry& validation_data, bool& debug_output, char alt_conf )
{
    copy(mm,clipper::MM::COPY_MPC);	// import_data from MMonomer

    this->sugar_chain_id = chainID;
    this->debug_output = debug_output;
    this->sugar_supported = true;
    this->sugar_parent_molecule = &ml;
    this->sugar_parent_molecule_nonbond = &nb; // store pointers
    this->sugar_index = db_not_checked;
    this->sugar_index = 9999; // default value for "not found in database".
    this->sugar_alternate_confcode = " "; // initially, we would like to suppose this
    this->sugar_context = "";
    this->sugar_pdb_id = mm.id();
    this->sugar_seqnum = mm.seqnum();
    this->sugar_name_short = mm.type().trim();
    this->sugar_name_full = validation_data.name_long;

    if(debug_output)
    {
        std::cout << std::endl ;
        DBG << "looking for " << this->id() << " " << this->type().trim() << " on the database..." << std::endl;
        // DBG << "with the coords of = " << this->find("C1", MM::ANY).coord_orth().format() << " on the database..." << std::endl;
        // alt_conf != ' ' ? DBG << "Alternate locator supplied: " << alt_conf << std::endl : true;
    }

    this->sugar_found_db = true;

    sugar_bfactor = 0.0;

    for (int i=0; i < this->size(); i++)
    {
        MSugar mstmp= *this;
        sugar_bfactor += mstmp[i].u_iso();
    }

    sugar_bfactor /= this->size();
    sugar_bfactor = clipper::Util::u2b(sugar_bfactor);

    if(debug_output)
    {
        DBG << "found it! " << std::endl;
    }

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

            for(int j = 0; j < 3; j++)
                this->sugar_cremer_pople_params.push_back(-1);

            this->sugar_conformation = 0;

            return;
        }
        else if ( alt_conf != ' ' ) // alternate conformations are present
        {
            alt_conf == 'A' ? sugar_alternate_confcode = " :A" : sugar_alternate_confcode = " :B";

            index_atom = this->lookup(buffer[i].trim()+sugar_alternate_confcode,clipper::MM::UNIQUE);

            if (index_atom == -1) // we've tried A and B and it still fails... so we're going to give up for now
            {
                this->sugar_supported = false;
                this->sugar_sane = false;
                this->sugar_denomination = "    unsupported    ";
                this->sugar_anomer = "X";
                this->sugar_handedness = "X";
                return;
            }
        }
        sugar_ring_elements.push_back((*this)[index_atom]);
    }

    if (this->sugar_ring_elements.size() == 5)
    {
        this->cremerPople_furanose(*this->sugar_parent_molecule, mm);
        this->sugar_conformation = conformationFuranose(this->sugar_cremer_pople_params[1]);

        if(debug_output)
        {
            DBG << "After checking the conformation..." << std::endl;
        }
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

        for(int j = 0; j < 3; j++)
            this->sugar_cremer_pople_params.push_back(-1);

        this->sugar_conformation = 0;

    }

    if (sugar_ring_elements[1].name().trim().find("C1") != std::string::npos)
    {
        this->sugar_type = "aldose";
        this->sugar_denomination = "aldo";
    }
    else
    {
        this->sugar_type = "ketose";
        this->sugar_denomination = "keto";
    }

    if (sugar_ring_elements.size() == 5)
        this->sugar_denomination = this->sugar_denomination + "furanose";
    else
        this->sugar_denomination = this->sugar_denomination + "pyranose";

    this->sugar_denomination = clipper::String( this->sugar_anomer + "-" + this->sugar_handedness + "-" + this->sugar_denomination );

    // sanity check:
    this->sugar_sane = false;

    if(debug_output)
    {
        DBG << "Just before examining the ring..." << std::endl;
    }

    if ( examine_ring() )
        sugar_diag_ring = true;
    else
        sugar_diag_ring = false;

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
        else
            sugar_diag_conformation = false;
    }

    if ( sugar_diag_conformation )
    {
        if (( puckering_amplitude() > ref_puckering - 0.18 ) && (puckering_amplitude() < ref_puckering + 0.15))
            sugar_diag_puckering = true;
        else
            sugar_diag_puckering = false;

        if ( sugar_ring_bond_rmsd < ( ref_bonds_rmsd + 0.039 ) )
            sugar_diag_bonds_rmsd = true;
        else
            sugar_diag_bonds_rmsd = false;

        if ( sugar_ring_angle_rmsd < ( ref_angles_rmsd + 3.0 ) )
            sugar_diag_angles_rmsd = true;
        else
            sugar_diag_angles_rmsd = false;
    }
    else
    {
        if (( puckering_amplitude() > 0.9 ) || ( puckering_amplitude() < 0.42 ))
            sugar_diag_puckering = false;
        else sugar_diag_puckering = true;

        sugar_diag_bonds_rmsd = sugar_diag_angles_rmsd = true;
    }

    if ( sugar_diag_puckering && sugar_diag_anomer && sugar_diag_chirality && sugar_diag_ring )
        sugar_sane = true;

    if(debug_output)
    {
	    DBG << "Just after examining the ring, exiting the constructor, good job!" << std::endl;
    }

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

	if(debug_output)
    {
		DBG << "getting the stereochemistry..." << std::endl;
	}

	stereochemistry_pairs stereo = get_stereochemistry(mmol);

	if(debug_output)
    {
		DBG << "done." << std::endl;
	}

	this->sugar_anomeric_carbon = stereo.first.first;
	this->sugar_anomeric_substituent = stereo.first.second;
	this->sugar_configurational_carbon = stereo.second.first;
	this->sugar_configurational_substituent = stereo.second.second;

	clipper::String anomeric_carbon = stereo.first.first.name().trim();
	clipper::String anomeric_substituent = stereo.first.second.name().trim();
	clipper::String configurational_carbon = stereo.second.first.name().trim();
	clipper::String configurational_substituent = stereo.second.second.name().trim();

	if ( configurational_carbon != "XXX" )
		if ( ( ring_atoms[5].name().trim() == configurational_carbon) ||
               !is_part_of_ring( sugar_configurational_carbon, ring_atoms ) )
			lurd_reverse = true; // we have to take into account the rotation performed prior to closing the ring

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

    clipper::RTop_orth shift(clipper::Mat33<>::identity(), (- centre)); //clipper::Vec3<>::null()

    mm.transform(shift); // recentre the sugar coordinates

/*	stereo.first.first.transform(shift);
	stereo.first.second.transform(shift);
	stereo.second.first.transform(shift);
	stereo.second.second.transform(shift); */

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
    sugar_mean_plane = n;

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
            {
                if ( !is_part_of_ring(mmol.atom(neighbourhood[i]), ring_atoms) && (ring_atoms[5].occupancy() == mmol.atom(neighbourhood[i]).occupancy() ) ) // eliminate ring neighbours
                {
                    if (mmol.atom(neighbourhood[i]).element().trim() != "C") nz6_substituent = mmol.atom(neighbourhood[i]);
                        // don't grab a carbon as substituent unless there's really no other option
                    else if (nz6_substituent.name().trim() == "XXX") nz6_substituent = mmol.atom(neighbourhood[i]);
                }
            }
        }
    }

    nz6_substituent.transform(shift); // we still need to recentre the atom
    z6_substituent = clipper::Vec3<clipper::ftype>::dot(nz6_substituent.coord_orth(), n);

    if(debug_output)
    {
        DBG << "last in-ring carbon has occupancy " << ring_atoms[5].occupancy() << " and its substituent is "
            << nz6_substituent.name().trim() << " with occupancy " << nz6_substituent.occupancy() << std::endl;
    }

	///////// stereochemistry ///////////

    clipper::Vec3<clipper::ftype> anomeric_plane;


    clipper::Vec3<clipper::ftype> vec1_anomer ( sugar_ring_elements[0].coord_orth().x() - sugar_ring_elements[1].coord_orth().x(),
                                                sugar_ring_elements[0].coord_orth().y() - sugar_ring_elements[1].coord_orth().y(),
                                                sugar_ring_elements[0].coord_orth().z() - sugar_ring_elements[1].coord_orth().z() );

    clipper::Vec3<clipper::ftype> vec2_anomer ( sugar_ring_elements[2].coord_orth().x() - sugar_ring_elements[1].coord_orth().x(),
                                                sugar_ring_elements[2].coord_orth().y() - sugar_ring_elements[1].coord_orth().y(),
                                                sugar_ring_elements[2].coord_orth().z() - sugar_ring_elements[1].coord_orth().z() );

    anomeric_plane = clipper::Vec3<clipper::ftype>::cross ( vec1_anomer, vec2_anomer ); // orthogonal, pointing up

    z_anomeric_carbon = clipper::Vec3<clipper::ftype>::dot( stereo.first.first.coord_orth(), anomeric_plane );
    z_anomeric_substituent = clipper::Vec3<clipper::ftype>::dot( stereo.first.second.coord_orth(), anomeric_plane );



    if ( stereo.second.first.name().trim() == sugar_ring_elements[5].name().trim() ) // most general case
    {
        clipper::Vec3<clipper::ftype> config_plane;


        clipper::Vec3<clipper::ftype> vec1_config ( sugar_ring_elements[4].coord_orth().x() - sugar_ring_elements[5].coord_orth().x(),
                                                    sugar_ring_elements[4].coord_orth().y() - sugar_ring_elements[5].coord_orth().y(),
                                                    sugar_ring_elements[4].coord_orth().z() - sugar_ring_elements[5].coord_orth().z() );

        clipper::Vec3<clipper::ftype> vec2_config ( sugar_ring_elements[0].coord_orth().x() - sugar_ring_elements[5].coord_orth().x(),
                                                    sugar_ring_elements[0].coord_orth().y() - sugar_ring_elements[5].coord_orth().y(),
                                                    sugar_ring_elements[0].coord_orth().z() - sugar_ring_elements[5].coord_orth().z() );

        config_plane = clipper::Vec3<clipper::ftype>::cross ( vec1_config, vec2_config ); // orthogonal, pointing up

        z_configurational_carbon = clipper::Vec3<clipper::ftype>::dot( stereo.second.first.coord_orth(), config_plane );
        z_configurational_substituent = clipper::Vec3<clipper::ftype>::dot( stereo.second.second.coord_orth(), config_plane );

        if ( z_configurational_substituent > z_configurational_carbon )
            this->sugar_handedness = "D";
        else
            this->sugar_handedness = "L";

    }
    else
    {
        //stereo.second.first.transform(shift);
        //stereo.second.second.transform(shift);

        z_configurational_carbon = clipper::Vec3<clipper::ftype>::dot( stereo.second.first.coord_orth(), n );
        z_configurational_substituent = clipper::Vec3<clipper::ftype>::dot( stereo.second.second.coord_orth(), n );

        if ((is_part_of_ring(nz6_substituent, this->sugar_ring_elements)) || (nz6_substituent.name().trim() == "XXX"))
            this->sugar_handedness = "N";
        else if ((z6 - z6_substituent) < 0)
            this->sugar_handedness = "D";
        else
            this->sugar_handedness = "L";
    }



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

    if ( stereo.second.first.name().trim() == "C8") { // Simplified bit for sialic acids
      if (z_anomeric_substituent > z_anomeric_carbon) {
        cpParams.push_back(anomer_beta);
				this->sugar_anomer="beta";
      }
      else {
        cpParams.push_back(anomer_alpha);
				this->sugar_anomer="alpha";
      }
    }
    else if (( (z_anomeric_substituent > z_anomeric_carbon) && (z_configurational_substituent > z_configurational_carbon) ) || ( (z_anomeric_substituent < z_anomeric_carbon) && (z_configurational_substituent < z_configurational_carbon) ) )
		{
			if (lurd_reverse)
			{
				cpParams.push_back(anomer_beta);
				this->sugar_anomer="beta";
      }
      else {
        cpParams.push_back(anomer_alpha);
				this->sugar_anomer="alpha";
      }
    }
    else if (( (z_anomeric_substituent > z_anomeric_carbon) && (z_configurational_substituent > z_configurational_carbon) ) || ( (z_anomeric_substituent < z_anomeric_carbon) && (z_configurational_substituent < z_configurational_carbon) ) )
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

    if(debug_output)
    {
        DBG << "an_c= " << anomeric_carbon << "/" << z_anomeric_carbon << " - an_subs= " << anomeric_substituent << "/" << z_anomeric_substituent << std::endl;
        DBG << "conf_c= " << configurational_carbon << "/" << z_configurational_carbon << " - conf_subs= " << configurational_substituent << "/" << z_configurational_substituent << std::endl;
        DBG << "z6= " << z6 << " z6_subs= " << z6_substituent << std::endl;
    }

    cpParams.push_back( z6 - z6_substituent );

	if(debug_output)
    {
		DBG << "Finished Cremer-Pople analysis, returning to caller..." << std::endl;
	}

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

	if(debug_output)
    {
		DBG << "After getting the stereochemistry" << std::endl;
	}


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

        if(debug_output)
        {
            DBG << "Ring centre: " << centre.format() << std::endl;
        }

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
        sugar_mean_plane = n;

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
            if(debug_output)
            {
                DBG << "Neighbour found for " << ring_atoms[4].name()
                    << ": " << (mmol.atom(neighbourhood[i]).name())
                    << " at distance "
                    << clipper::Coord_orth::length ( mmol.atom(neighbourhood[i]).coord_orth(), ring_atoms[4].coord_orth() )
                    << std::endl;
            }

            if (( mmol.atom(neighbourhood[i]).element().trim() != "H" ) &&
                (mmol.atom(neighbourhood[i]).name().trim() != ring_atoms[4].name().trim())) // the target substituent could be anything apart from H, in-ring oxygen or in-ring carbon
            {
                if ( bonded(mmol.atom(neighbourhood[i]), ring_atoms[4]) ) // check link
                    if ( !is_part_of_ring ( mmol.atom(neighbourhood[i]), ring_atoms ) ) // eliminate ring neighbours
                    {

                        if ( mmol.atom(neighbourhood[i]).element().trim() != "C" )
                            nz5_substituent = mmol.atom ( neighbourhood[i] ); // don't grab a carbon as substituent unless there's really no other option
                        else if ( nz5_substituent.name().trim() == "XXX" )
                            nz5_substituent = mmol.atom ( neighbourhood[i] );
                    }
            }
        }

	if(debug_output)
    {
		DBG << "substituent at last in-ring carbon: " << nz5_substituent.name().trim() << std::endl;
	}

	nz5_substituent.transform(shift); // we still need to recentre the atom
	z5_substituent = clipper::Vec3<clipper::ftype>::dot(nz5_substituent.coord_orth(), n);



	///////// stereochemistry ///////////

    clipper::Vec3<clipper::ftype> anomeric_plane;


    clipper::Vec3<clipper::ftype> vec1_anomer ( sugar_ring_elements[0].coord_orth().x() - sugar_ring_elements[1].coord_orth().x(),
                                                sugar_ring_elements[0].coord_orth().y() - sugar_ring_elements[1].coord_orth().y(),
                                                sugar_ring_elements[0].coord_orth().z() - sugar_ring_elements[1].coord_orth().z() );

    clipper::Vec3<clipper::ftype> vec2_anomer ( sugar_ring_elements[2].coord_orth().x() - sugar_ring_elements[1].coord_orth().x(),
                                                sugar_ring_elements[2].coord_orth().y() - sugar_ring_elements[1].coord_orth().y(),
                                                sugar_ring_elements[2].coord_orth().z() - sugar_ring_elements[1].coord_orth().z() );

    anomeric_plane = clipper::Vec3<clipper::ftype>::cross ( vec1_anomer, vec2_anomer ); // orthogonal, pointing up

    //float angle_anomer = clipper::Util::rad2d ( get_angle ( vector_anomer, vector_config ) );

    z_anomeric_carbon = clipper::Vec3<clipper::ftype>::dot( stereo.first.first.coord_orth(), anomeric_plane );
    z_anomeric_substituent = clipper::Vec3<clipper::ftype>::dot( stereo.first.second.coord_orth(), anomeric_plane );

    if ( stereo.second.first.coord_orth() == sugar_ring_elements[4].coord_orth() )
    {
        clipper::Vec3<clipper::ftype> config_plane;


        clipper::Vec3<clipper::ftype> vec1_config ( sugar_ring_elements[3].coord_orth().x() - sugar_ring_elements[4].coord_orth().x(),
                                                    sugar_ring_elements[3].coord_orth().y() - sugar_ring_elements[4].coord_orth().y(),
                                                    sugar_ring_elements[3].coord_orth().z() - sugar_ring_elements[4].coord_orth().z() );

        clipper::Vec3<clipper::ftype> vec2_config ( sugar_ring_elements[0].coord_orth().x() - sugar_ring_elements[4].coord_orth().x(),
                                                    sugar_ring_elements[0].coord_orth().y() - sugar_ring_elements[4].coord_orth().y(),
                                                    sugar_ring_elements[0].coord_orth().z() - sugar_ring_elements[4].coord_orth().z() );

        config_plane = clipper::Vec3<clipper::ftype>::cross ( vec1_config, vec2_config ); // orthogonal, pointing up

        z_configurational_carbon = clipper::Vec3<clipper::ftype>::dot( stereo.second.first.coord_orth(), config_plane );
        z_configurational_substituent = clipper::Vec3<clipper::ftype>::dot( stereo.second.second.coord_orth(), config_plane );

    }
    else
    {
        z_configurational_carbon = clipper::Vec3<clipper::ftype>::dot( stereo.second.first.coord_orth(), n );
        z_configurational_substituent = clipper::Vec3<clipper::ftype>::dot( stereo.second.second.coord_orth(), n );
    }

    clipper::ftype totalPuckering = sqrt(pow(z1,2) + pow(z2,2) + pow(z3,2) + pow(z4,2) + pow(z5,2));

    clipper::ftype phi2, q2, argACos, argASin;

    argACos =  sqrt(2.0/5.0) * (z1 + z2*cos(4.0*clipper::Util::pi()/5.0)
                + z3*cos(8.0*clipper::Util::pi()/5.0) + z4*cos(12.0*clipper::Util::pi()/5.0)
                + z5*cos(16.0*clipper::Util::pi()/5.0));

    argASin = -sqrt(2.0/5.0) * (z2*sin(4.0*clipper::Util::pi()/5.0) + z3*sin(8.0*clipper::Util::pi()/5.0)
                + z4*sin(12.0*clipper::Util::pi()/5.0) + z5*sin(16.0*clipper::Util::pi()/5.0));

    phi2 = (clipper::Util::pi()/2) - atan( argACos / argASin); // we want arccotan but there's no such thing in C math I'm afraid
    clipper::ftype phi3 = atan (argASin / argACos );

    q2 = (argACos / cos(phi2));

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

	if(debug_output)
    {
        DBG << "an_c= " << anomeric_carbon << "/" << z_anomeric_carbon << " - an_subs= " << anomeric_substituent << "/" << z_anomeric_substituent << std::endl;
        DBG << "conf_c= " << configurational_carbon << "/" << z_configurational_carbon << " - conf_subs= "
            << configurational_substituent << "/" << z_configurational_substituent << std::endl;
	    DBG << "z5= " << z5 << " z5_subs= " << z5_substituent << std::endl;
	}

    cpParams.push_back( z5 - z5_substituent );

    if ((is_part_of_ring(nz5_substituent, this->sugar_ring_elements)) || (nz5_substituent.name().trim() == "XXX"))
        this->sugar_handedness = "N";
    else if ((z5 - z5_substituent) < 0)
        this->sugar_handedness = "D";
    else
        this->sugar_handedness = "L";

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
    else if (phi == -1 && theta == -1) confCode = 0;

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
    else if (phi == -1) confCode = 0;

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


    if(buffer.size() > 0)
	    buffer.erase( buffer.begin() );

	buffer.resize( index, clipper::MAtom() );

	std::vector<clipper::MAtom> result;

	if(debug_output)
    {
		DBG << "Dumping ring contents... " << std::endl;
		for (int runner = 0 ; runner < buffer.size() ; runner++ ) DBG << buffer[runner].name().trim() << " with occupancy = " << buffer[runner].occupancy() << std::endl;
	}

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

				//if(debug_output)
                //{
				//	DBG << buffer[i].name().trim() << " with index " << origin << " VS " << result[1].name().trim() << " with index " << destination << std::endl;
				//}

				int subindex = 1;

				while ( ( origin > destination ) && (subindex < result.size() ) )
				{

				//if(debug_output)
                //{
				//	DBG << "Subindex is " << subindex << " and result[subindex] is " << result[subindex].name() << std::endl;
				//}

					if ( ++subindex < result.size() )
					{
						std::stringstream stream_3( result[subindex].name().trim().split("C")[0] );
						stream_3 >> destination;
					}
				}

				//if(debug_output)
                //{
				//	DBG << "Inserting " << origin << " which is smaller than " << destination << " at position " << subindex << " and result.size() is " << result.size() << std::endl;
				//}

				result.insert( result.begin() + subindex , buffer[i] );
			}
			else
			{
				result.push_back( buffer[i] );
				//if(debug_output)
                //{
				//	DBG << "Inserting " << buffer[i].name().trim() << " as first item after the oxygen..." << std::endl;
				//}
			}
		}
	}

	if(debug_output)
    {
		DBG << "Successfully determined ring members!" << std::endl;
	}

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
				if (altconf_compatible(get_altconf(ma), get_altconf(mm[i]))) // we need to stay on the same conformation
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
        if ( ring_atoms[1].element().trim() == "C" )        // the atom in position 1 is the anomeric carbon, let's identify its substituent:
        {
            anomeric_carbon = ring_atoms[1]; // we're always getting the anomeric carbon here, given the way we detect/order the ring
            const std::vector<clipper::MAtomIndexSymmetry> neighbourhood = this->sugar_parent_molecule_nonbond->atoms_near(ring_atoms[1].coord_orth(), 1.2); // 1.2

            for ( int i = 0 ; i < neighbourhood.size() ; i++ )
            {
                if (( mmol.atom(neighbourhood[i]).element().trim() != "H" )
                 && ( mmol.atom(neighbourhood[i]).name().trim() != anomeric_carbon.name().trim()) )
                {
                    if ( bonded(neighbourhood[i], ring_atoms[1] ) )
                    {
                        if ( !is_part_of_ring(mmol.atom(neighbourhood[i]), ring_atoms))
                        {
                            if ( altconf_compatible(get_altconf ( mmol.atom(neighbourhood[i]) ), get_altconf ( ring_atoms[1] ) ))
                            {
                                int polymerID = neighbourhood[i].polymer(); // needed for ugly fix down below.
                                int monomerID = neighbourhood[i].monomer();
                                if (mmol.atom(neighbourhood[i]).element().trim() != "C" ) // O6 is part of the ring if only checking by name!!!
                                    anomeric_substituent = mmol.atom(neighbourhood[i]);
                                else if (mmol[polymerID][monomerID].type().trim() == "TRP") // really ugly fix for C/TRP-mannosylation...
                                    if(mmol.atom(neighbourhood[i]).element().trim() == "C")
                                        anomeric_substituent = mmol.atom(neighbourhood[i]);
                            }
                            /*else if ( get_altconf ( mmol.atom(neighbourhood[i]) ) == 'A' && get_altconf ( ring_atoms[1] ) != 'B')
                            {
                                if (mmol.atom(neighbourhood[i]).element().trim() != "C" ) // O6 is part of the ring if only checking by name!!!
                                    anomeric_substituent = mmol.atom(neighbourhood[i]);
                            }
                            else if ( get_altconf ( ring_atoms[1] ) == 'A' && get_altconf ( mmol.atom(neighbourhood[i]) ) !='B' )
                            {
                                if (mmol.atom(neighbourhood[i]).element().trim() != "C" ) // O6 is part of the ring if only checking by name!!!
                                    anomeric_substituent = mmol.atom(neighbourhood[i]);
                            }*/
                        }
                    }
                }
            }
        }
    }

	if(debug_output)
    {
		DBG << "Anomeric carbon: " << anomeric_carbon.id() << "  Substituent: " << anomeric_substituent.id() << std::endl;
	}

	result.first.first = anomeric_carbon;
	result.first.second = anomeric_substituent;

	for ( int i = 2 ; i < ring_atoms.size() ; i++ ) // we start checking for the highest ranked stereocentre at the in-ring carbon next (clockwise) to the anomeric carbon
		if (ring_atoms[i].element().trim() == "C")
			if ( is_stereocentre(ring_atoms[i], mmol) )
			{
				configurational_carbon = ring_atoms[i];

                if(debug_output)
                {
                    DBG << "i: " << i << "\tring_atoms.size() " << ring_atoms.size() << std::endl;
                }

				// get the configurational carbon's target substituent
				const std::vector<clipper::MAtomIndexSymmetry> neighbourhood = this->sugar_parent_molecule_nonbond->atoms_near(configurational_carbon.coord_orth(), 1.2); // was 1.4

				for ( int i2 = 0 ; i2 < neighbourhood.size() ; i2++ )
				{
                    if(debug_output && i2 == 1 || debug_output && i2 == neighbourhood.size()-2)
                    {
                        DBG << "i2: " << i2 << "\tneighbourhood.size() " << neighbourhood.size() << std::endl;
                    }
					if ( !is_part_of_ring(mmol.atom(neighbourhood[i2]),ring_atoms)
						&& (mmol.atom(neighbourhood[i2]).element().trim() != "H" )
						&& altconf_compatible(get_altconf(mmol.atom(neighbourhood[i2])), get_altconf(anomeric_carbon)))
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


	if(debug_output)
    {
		DBG << "(in-ring) configurational carbon: " << configurational_carbon.id() << "  substituent: " << configurational_substituent.id() << std::endl;
	}

	// we've recorded the highest ranked in-ring carbon atom & substituent in configurational_*

	clipper::MAtom next_carbon = configurational_substituent;

	while ( is_stereocentre( next_carbon, mmol ) && (next_carbon.id().trim() != configurational_carbon.id().trim() ) ) ////////////////////// CONTINUE HERE
	{
		configurational_carbon = next_carbon;
		const std::vector<clipper::MAtomIndexSymmetry> neighbourhood = this->sugar_parent_molecule_nonbond->atoms_near(configurational_carbon.coord_orth(), 1.2);

		for ( int i2 = 0 ; i2 < neighbourhood.size() ; i2++ )
		{
            if(debug_output && i2 == 1 || debug_output && i2 == neighbourhood.size()-2)
            {
                DBG << "i2 in while loop: " << i2 << "\tneighbourhood.size() " << neighbourhood.size() << std::endl;
            }
			if ( !is_part_of_ring(mmol.atom(neighbourhood[i2]),ring_atoms)
				&& (mmol.atom(neighbourhood[i2]).element().trim() != "H" )
				&& altconf_compatible(get_altconf(mmol.atom(neighbourhood[i2])), get_altconf(configurational_carbon)))
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

	if(debug_output)
    {
		DBG << "Configurational carbon: " << configurational_carbon.id() << "  Substituent: " << configurational_substituent.id() << std::endl;
	}

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
		else // neighbour is a symmetry mate, so let's find the symmetry operator that places a copy close to the anomeric carbon
		{
			clipper::Spacegroup spgr = mmol.spacegroup();
			clipper::Coord_frac f1 = mmol.atom( neighbourhood[k] ).coord_orth().coord_frac( mmol.cell() );
			clipper::Coord_frac f2 = ma.coord_orth().coord_frac( mmol.cell() );
			f1 = spgr.symop(neighbourhood[k].symmetry()) * f1;
			f1 = f1.lattice_copy_near( f2 );
			distance = sqrt(( f2 - f1 ).lengthsq( mmol.cell() ));
		}

		if ( distance < 2.0 ) // check link
			if (( mmol.atom(neighbourhood[k]).element().trim() != "H" ) && (mmol.atom(neighbourhood[k]).name().trim() != ma.name().trim() ) && altconf_compatible(get_altconf(mmol.atom(neighbourhood[k])), get_altconf(ma)) )
			{
				//if(debug_output)
                //{
				//	DBG << "Counting " << mmol.atom(neighbourhood[k]).id() << " as substituent from a total of " << neighbourhood.size() << " atoms with symop " << neighbourhood[k].symmetry() << std::endl;
				//}

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

	//if(debug_output)
    //{
	//	DBG << "Number of substituents: " << substituent_list.size() << std::endl;
	//}

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
	for (int i = 0; i < ring_atoms.size(); i++) // check by coords seems to be a much better option than checking by name!
		if ( ring_atoms[i].coord_orth() == ma.coord_orth() ) return true;

	return false;
}


/*! Internal function for checking whether two atoms are bonded or not
	\param ma_one A clipper::MAtom object
	\param ma_two A clipper::MAtom object
	\return A boolean value with the answer
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
        {
            if ((distance >  1.18 ) && ( distance < 1.62 ))
                return true; // C-C or C=C
            else
                return false;
        }
        else if ( ma_two.element().trim() == "N" )
        {
            if ((distance >  1.24 ) && ( distance < 1.62 ))
                return true; // C-N or C=N
            else
                return false;
        }
        else if ( ma_two.element().trim() == "O" )
        {
            if ((distance >  1.16 ) && ( distance < 1.60 ))
                return true; // C-O or C=O
            else
                return false;
        }
        else if ( ma_two.element().trim() == "S" )
        {
            if ((distance >  1.50 ) && ( distance < 2.00 ))
                return true; // C-S or C=S
            else
                return false;
        }
        else if ( ma_two.element().trim() == "H" )
        {
            if ((distance >  0.96 ) && ( distance < 1.14 ))
                return true; // C-H
            else
                return false;
        }
    }
    else if ( sugar_parent_molecule->atom(ma_one).element().trim() == "N" )
    {
        if ( ma_two.element().trim() == "C" )
        {
            if ((distance >  1.24 ) && ( distance < 1.62 ))
                return true; // N-C or N=C
            else
                return false;
        }
        else if ( ma_two.element().trim() == "H" )
        {
            if ((distance >  0.90 ) && ( distance < 1.10 ))
                return true; // N-H
            else
                return false;
        }
    }
    else if ( sugar_parent_molecule->atom(ma_one).element().trim() == "O" )
    {
        if ( ma_two.element().trim() == "C" )
        {
            if ((distance >  1.16 ) && ( distance < 1.60 ))
                return true; // O-C or O=C
            else
                return false;
        }
        else if ( ma_two.element().trim() == "H" )
        {
            if ((distance >  0.88 ) && ( distance < 1.04 ))
                return true; // O-H
            else
                return false;
        }
    }
    else if ( sugar_parent_molecule->atom(ma_one).element().trim() == "S" )
    {
        if ( ma_two.element().trim() == "C" )
        {
            if ((distance >  1.26 ) && ( distance < 2.00 ))
                return true; // S-C or S=C
            else return false;
        }
        else if ( ma_two.element().trim() == "H" )
        {
            if ((distance >  0.78 ) && ( distance < 1.44 ))
                return true; // S-H
            else return false;
        }
    }
    else     /// unknown bond
        if (( distance > 1.2) && (distance < 1.8))
            return true;
        else
            return false;

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
        {
            if(debug_output)
            {
                DBG << std::endl << "Returning false - sugar_ring_elements[" << i << "].id().trim() "  << sugar_ring_elements[i].id().trim() << "\t\t\tsugar_ring_elements[" << i+1 << "].id().trim() " << sugar_ring_elements[i+1].id().trim() << std::endl;
            }
            return false;
        }
        if (!bonded(sugar_ring_elements[i], sugar_ring_elements[0]))
        {
            if(debug_output)
            {
                DBG << std::endl << "Returning false - sugar_ring_elements[" << i << "].id().trim() "  << sugar_ring_elements[i].id().trim() << "\t\t\tsugar_ring_elements[" << 0 << "].id().trim() " << sugar_ring_elements[0].id().trim() << std::endl;
            }
            return false;
        }
        else
        {
            return true;
        }
	    
}

/*! Internal function for checking whether two atoms are bonded or not
 * 	\param ma_one A clipper::MAtom object
 * 	\param ma_two A clipper::MAtom object
 * 	\return A boolean value with the obvious answer
 */

bool MSugar::bonded(const clipper::MAtom& ma_one, const clipper::MAtom& ma_two) const
{
    clipper::ftype distance = clipper::Coord_orth::length( ma_one.coord_orth(), ma_two.coord_orth() );
    if(debug_output)
    {
        DBG << std::endl << "Received - ma_one.id().trim() " << ma_one.id().trim() << "\t\t\tma_two.id().trim() " << ma_two.id().trim() << "\tdistance = " << distance << std::endl;
    }
    if ( ma_one.element().trim() == "C" )
    {
        if ( ma_two.element().trim() == "C" )
        {
            if ((distance >  1.18 ) && ( distance < 1.62 ))
                return true; // C-C or C=C
            else
                return false;
        }
        else if ( ma_two.element().trim() == "N" )
        {
            if ((distance >  1.24 ) && ( distance < 1.62 ))
                return true; // C-N or C=N
            else
                return false;
        }
        else if ( ma_two.element().trim() == "O" )
        {
            if ((distance >  1.16 ) && ( distance < 1.60 ))
                return true; // C-O or C=O
            else
                return false;
        }
        else if ( ma_two.element().trim() == "S" )
        {
            if (( distance > 1.2 ) && ( distance < 2.0 ))
                return true; // C-S or C=S
            else
                return false;
        }
        else if ( ma_two.element().trim() == "H" )
        {
            if ((distance >  0.96 ) && ( distance < 1.14 ))
                return true; // C-H
            else
                return false;
        }
    }
    else if ( ma_one.element().trim() == "N" )
    {
        if ( ma_two.element().trim() == "C" )
        {
            if ((distance >  1.24 ) && ( distance < 1.62 ))
                return true; // N-C or N=C
            else
                return false;
        }
        else if ( ma_two.element().trim() == "H" )
        {
            if ((distance >  0.90 ) && ( distance < 1.10 ))
                return true; // N-H
            else
                return false;
        }
    }
    else if ( ma_one.element().trim() == "O" )
    {
        if ( ma_two.element().trim() == "C" )
        {
            if ((distance >  1.16 ) && ( distance < 1.60 ))
                return true; // O-C or O=C
            else
                return false;
        }
        else if ( ma_two.element().trim() == "H" )
        {
            if ((distance >  0.88 ) && ( distance < 1.04 ))
                return true; // O-H
            else
                return false;
        }
    }
    else if ( ma_one.element().trim() == "S" )
    {
        if ( ma_two.element().trim() == "C" )
        {
            if ((distance >  1.26 ) && ( distance < 2.00 ))
                return true; // O-C or O=C
            else
                return false;
        }
        else if ( ma_two.element().trim() == "H" )
        {
            if ((distance >  0.78 ) && ( distance < 1.24 ))
                return true; // O-H
            else
                return false;
        }
    }
    else if (( distance > 1.2) && (distance < 2.0))
        return true; // unknown bond
    else
        return false;

    return false; // in case we haven't found any match
}


std::vector < std::pair< clipper::MAtomIndexSymmetry, clipper::ftype > > MSugar::get_stacked_residues ( ) const
{
    clipper::MAtom ma;
    clipper::Coord_orth centre_apolar;

    std::vector < std::pair< clipper::MAtomIndexSymmetry, clipper::ftype > > results;

    if(ring_members().size() < 6 || type_of_sugar().trim() == "unsupported")
        return results;

    if ( this->handedness() == "D" )
        centre_apolar = clipper::Coord_orth((ring_members()[1].coord_orth().x() +
                                             ring_members()[3].coord_orth().x() +
                                             ring_members()[5].coord_orth().x() ) / 3.0,
                                            (ring_members()[1].coord_orth().y() +
                                             ring_members()[3].coord_orth().y() +
                                             ring_members()[5].coord_orth().y() ) / 3.0,
                                            (ring_members()[1].coord_orth().z() +
                                             ring_members()[3].coord_orth().z() +
                                             ring_members()[5].coord_orth().z() ) / 3.0 );
    else
        centre_apolar = clipper::Coord_orth((ring_members()[1].coord_orth().x() +
                                             ring_members()[5].coord_orth().x() ) / 2.0,
                                            (ring_members()[1].coord_orth().y() +
                                             ring_members()[5].coord_orth().y() ) / 2.0,
                                            (ring_members()[1].coord_orth().z() +
                                             ring_members()[5].coord_orth().z() ) / 2.0 );

    const std::vector<clipper::MAtomIndexSymmetry> neighbourhood = this->sugar_parent_molecule_nonbond->atoms_near(centre_apolar, 5.0);

    const clipper::MiniMol& mmol  = *sugar_parent_molecule;

	for ( int k = 0 ; k < neighbourhood.size() ; k++ )
	{
        const clipper::MMonomer& mmon = mmol[neighbourhood[k].polymer()][neighbourhood[k].monomer()];

        if (( mmon.type() != "TRP" ) &&
            ( mmon.type() != "TYR" ) &&
            ( mmon.type() != "PHE" ) &&
            ( mmon.type() != "HIS" ))
            continue;

		clipper::ftype distance = 0.0;

		if ( neighbourhood[k].symmetry() == 0 )
		{
			distance = clipper::Coord_orth::length(mmol.atom(neighbourhood[k]).coord_orth(), centre_apolar);
		}
		else // this neighbour is actually a symmetry mate
		{
			clipper::Spacegroup spgr = mmol.spacegroup();
			clipper::Coord_frac f1 = mmol.atom( neighbourhood[k] ).coord_orth().coord_frac( mmol.cell() );
			clipper::Coord_frac f2 = centre_apolar.coord_frac( mmol.cell() );
			f1 = spgr.symop(neighbourhood[k].symmetry()) * f1;
			f1 = f1.lattice_copy_near( f2 );
			distance = sqrt(( f2 - f1 ).lengthsq( mmol.cell() ));
		}

		if ( distance < 4.0 ) // check distance, compliant with Hudson et al., JACS 2015
        {
            clipper::Vec3<ftype> aromatic_plane = find_aromatic_plane ( mmon );
            clipper::ftype angle = get_angle ( aromatic_plane, ring_mean_plane () );

            if ( angle > 2.75 ) // 0 < angle < Pi, angle must be 30deg at most, also compliant with Hudson et al., JACS 2015
            {
                std::pair < clipper::MAtomIndexSymmetry, ftype > individual_result;
                individual_result.first = neighbourhood[k];
                individual_result.second = angle;

                int residue;

                for ( residue = 0 ; residue < results.size(); residue++ )
                {
                    if ( mmol[results[residue].first.polymer()][results[residue].first.monomer()].id() == mmol[neighbourhood[k].polymer()][neighbourhood[k].monomer()].id() )
                        break;
                }

                if ( residue == results.size() ) // not found in the results vector
                    results.push_back ( individual_result );
            }
         }
	}
    return results;
}



///////////////////////// MDisaccharide ///////////////////////////////


MDisaccharide::MDisaccharide ( clipper::MiniMol& mmol, const clipper::MAtomNonBond& manb, const clipper::String chainID, clipper::MMonomer& mm, bool& debug_output )
{
    this->debug_output = debug_output;
    int index = search_disaccharides ( mm.type().c_str() ); // we know beforehand that this is a known disaccharide, no need to re-check

    clipper::data::sugar_database_entry val_string_one = clipper::data::disaccharide_database[index].sugar_one;
    clipper::data::sugar_database_entry val_string_two = clipper::data::disaccharide_database[index].sugar_two;

    sugar_one = clipper::MSugar ( mmol, chainID, mm, manb, val_string_one, debug_output );
    sugar_two = clipper::MSugar ( mmol, chainID, mm, manb, val_string_two, debug_output );

    sugar_one.set_type ( clipper::String( sugar_one.type().trim() + "[" + clipper::data::disaccharide_database[index].sugar_one.name_short + "]" ));
    sugar_two.set_type ( clipper::String( sugar_two.type().trim() + "[" + clipper::data::disaccharide_database[index].sugar_two.name_short + "]" ));
}



///////////////////////// MGlycan ///////////////////////////////




MGlycan::MGlycan ( clipper::String chain, clipper::MMonomer& root_aa, clipper::MSugar& root_sugar, clipper::String& root_sugar_chain_id, bool& debug_output, std::string expression_system )
{
    this->debug_output = debug_output;
    root.second = clipper::MSugar(root_sugar);
    sugars.push_back ( root.second );
    Node first_node( root.second );

    node_list.push_back ( first_node );

    root.first = root_aa;

    this->chain = chain;
    this->chain_root_sugar = root_sugar_chain_id;

    if(debug_output)
    {
        DBG << "root.first: " << root.first.type() << "; root.second: " << root.second.type() << std::endl;
    }

    /*if ( expression_system != "undefined" )
        set_annotations ( expression_system );*/
}

MGlycan::MGlycan ( clipper::String chain, clipper::MSugar& root_sugar, clipper::String& root_sugar_chain_id, bool& debug_output, std::string expression_system )
{
    this->debug_output = debug_output;
    root.second = clipper::MSugar(root_sugar);
    sugars.push_back ( root.second );
    Node first_node( root.second );

    node_list.push_back ( first_node );

    root.first = clipper::MMonomer();

    this->chain = chain;
    this->chain_root_sugar = root_sugar_chain_id;

    if(debug_output)
    {
        DBG << "root.first: " << root.first.type() << "; root.second: " << root.second.type() << std::endl;
    }

    /*if ( expression_system != "undefined" )
        set_annotations ( expression_system );*/
}


clipper::String MGlycan::print_linear ( const bool print_info, const bool html_format, const bool translate )
{
    clipper::String buffer = "";

    if(debug_output)
    {
        DBG << "Glycan length: " << sugars.size() << std::endl;
    }

    if ( html_format ) buffer.insert ( 0, "</sub>" );
    buffer.insert ( 0, root.first.id().trim() );
    if ( html_format ) buffer.insert ( 0, "<sub>" );

    if(debug_output)
    {
        DBG << "Accessed the root, which contains this sugar: " << root.second.type() <<  std::endl;
    }

    buffer.insert( 0, root.first.type().c_str() );

    html_format ? buffer.insert ( 0, "&#8722;" ) : buffer.insert( 0, "-" );

    clipper::String anomer;

    if ( html_format ) root.second.anomer() == "alpha" ? anomer = "&#945;" : anomer = "&#946;";
    else root.second.anomer() == "alpha" ? anomer = "a" : anomer = "b";

    html_format ? buffer.insert ( 0, anomer ) : buffer.insert ( 0, anomer );
    html_format ? buffer.insert ( 0, " &#8592; " ) : buffer.insert( 0, "-" );

    clipper::MSugar msug = node_list.front().get_sugar();

        if(debug_output)
        {
            DBG << "Node list size: " << node_list.size() << std::endl;
            DBG << "Accessed the first sugar!" << std::endl;
        }

    if ( print_info )
    {
        if ( html_format ) buffer.insert ( 0, "</sub>" );
        buffer.insert ( 0, msug.id().trim() );
        if ( html_format ) buffer.insert ( 0, "<sub>" );
    }
    if ( html_format ) buffer.insert ( 0, "</span>" );
    translate ? buffer.insert ( 0, clipper::data::carbname_of ( msug.type().trim() ).c_str() ) : buffer.insert( 0, msug.type().c_str() );
    if ( html_format )
    {
        buffer.insert ( 0, "\">" );
        buffer.insert ( 0, msug.type() + " " + msug.id().trim() + " in " + msug.conformation_name() );
        buffer.insert ( 0, "<span title=\"" );
    }

        if(debug_output)
        {
            DBG << "Node list size: " << node_list.size() << std::endl;
        }

    if ( node_list.size() < 2 ) return buffer;
    else
    {
        if ( html_format )
            node_list[0].get_connection(0).get_anomericity() == "alpha" ? anomer = "&#945;" : anomer = "&#946;";
        else
            node_list[0].get_connection(0).get_anomericity() == "alpha" ? anomer = "a" : anomer = "b";

        std::ostringstream s;
        s << node_list[0].get_connection(0).get_order();

        html_format ? buffer.insert ( 0, "&#8722;" ) : buffer.insert( 0, "-" );
        buffer.insert( 0, s.str() );
        buffer.insert( 0, anomer );
        html_format ? buffer.insert ( 0, "&#8592;" ) : buffer.insert( 0, "-" );
    }

    bool branching = false;
    Node branched_from;

    for ( int i = 1 ; i < node_list.size() ; i++ ) // loop over all sugars except the first one (already processed)
    {
        msug = node_list[i].get_sugar();

        if ( print_info )
        {
            if ( html_format )
                buffer.insert ( 0, "</sub>" );

            buffer.insert ( 0, msug.id().trim() );

            if ( html_format )
                buffer.insert ( 0, "<sub>" );
        }

        if ( html_format ) buffer.insert ( 0, "</span>" );
        translate ? buffer.insert ( 0, clipper::data::carbname_of ( msug.type().trim() ).c_str() ) : buffer.insert( 0, msug.type().c_str() );
        if ( html_format )
        {
            buffer.insert ( 0, "\">" );
            buffer.insert ( 0, msug.type() + " " + msug.id().trim() + " in " + msug.conformation_name() );
            buffer.insert ( 0, "<span title=\"" );
        }

        if (( node_list[i].number_of_connections() == 0 ) && branching ) // close branch
        {
            if ( html_format )
                branched_from.get_connection(1).get_anomericity() == "alpha" ? anomer = "&#945;" : anomer = "&#946;";
            else
                branched_from.get_connection(1).get_anomericity() == "alpha" ? anomer = "a" : anomer = "b";

            std::ostringstream s;
            s << branched_from.get_connection(1).get_order();

            branching = false;
            buffer.insert ( 0, "(" );
            html_format ? buffer.insert ( 0, "&#8722;" ) : buffer.insert ( 0, "-" );
            buffer.insert( 0, s.str() );
            buffer.insert( 0, anomer );
            html_format ? buffer.insert ( 0, "&#8592;" ) : buffer.insert( 0, "-" );
        }
        else if ( node_list[i].number_of_connections() > 1 ) // open branch
        {
            branching = true;
            buffer.insert ( 0, ")" );
            branched_from = node_list[i];
        }

        if ( node_list[i].number_of_connections() > 0 ) // draw the actual connection
        {
            if ( html_format )
                node_list[i].get_connection(0).get_anomericity() == "alpha" ? anomer = "&#945;" : anomer = "&#946;";
            else
                node_list[i].get_connection(0).get_anomericity() == "alpha" ? anomer = "a" : anomer = "b";

            std::ostringstream s;
            s << node_list[i].get_connection(0).get_order();

            html_format ? buffer.insert ( 0, "&#8722;" ) : buffer.insert( 0, "-" );
            buffer.insert( 0, s.str() );
            buffer.insert( 0, anomer );
            html_format ? buffer.insert ( 0, "&#8592;" ) : buffer.insert( 0, "-" );
        }

    }
    return buffer;
}

// Fix me: need to add support for 2-8 links (Sialic Acids)
// TO DO: review atom assignments regarding furanoses, as the changes here were made when I misunderstood some chemistry
bool MGlycan::link_sugars ( int link, clipper::MSugar& first_sugar, clipper::MSugar& next_sugar, clipper::MAtom& donorAtom, clipper::MAtom& acceptorAtom )
{
    int index = 0;
    bool found = false;

    if(debug_output)
    {
        DBG << "Linking " << first_sugar.type() << " with " << next_sugar.type() << std::endl;
    }

    for ( int i = 0 ; i < node_list.size() ; i++ )
        if ( strcmp( node_list[i].get_sugar().id().c_str(), first_sugar.id().c_str()) == 0)
        {
            found = true;
            index = i;
            break;
        }
    

    if (!found)
    {
        if(debug_output)
        {
            DBG << "We haven't found a match for the first sugar. This is bad news." << std::endl;
        }

        return true;
    }

    Node new_node( next_sugar );      // create a new node with the next sugar
    node_list.push_back ( new_node ); // add the new sugar to the node list

    if(debug_output)
    {
        DBG << "Creating new connection.\nlink: " << link << "\tnext_sugar.anomer(): " << next_sugar.anomer() << "\tnode_list.size()-1: " << node_list.size()-1 << std::endl;
        DBG << "Anomeric substituent: " << next_sugar.anomeric_substituent().id().trim() << " located at: " << next_sugar.anomeric_substituent().coord_orth().format() << std::endl;
    }

    Linkage new_connection ( link, next_sugar.anomer(), node_list.size()-1 );

    clipper::ftype  omega_nine  = 0.0, 
                    omega_eight = 0.0, 
                    omega_seven = 0.0, 
                    omega_six   = 0.0, 
                    omega       = 0.0, 
                    psi         = 0.0, 
                    phi         = 0.0,
                    phi_cone_ctwo_oeight_ceight = 0.0;

    clipper::MAtom actual_c1, first_sugar_c1, c1, o1, c2, o2, c3, o3, c4, o4, c5, o5_next_sugar, o5, next_sugar_ring_oxygen, c6, o6, first_sugar_ring_oxygen, c7, o7, c8, o8, c9, o9;
    if ( link == 8 )
    {
        // next_sugar 5 SIA - MAN
        // first_sugar 4 SIA - BMA
        if(next_sugar.ring_members().size() == 6)
        {
            o6 = next_sugar.ring_members()[0];              // O6
            c2 = next_sugar.ring_members()[1];              // C2
            c1 = next_sugar[next_sugar.lookup("C1",clipper::MM::ANY)];
            o8 = next_sugar.anomeric_substituent();         // O8 usually
        }
        else if(next_sugar.ring_members().size() == 5)
        {
            o6 = next_sugar.ring_members()[0];              // O5
            c1 = next_sugar[next_sugar.lookup("C1",clipper::MM::ANY)]; 
            c2 = next_sugar.ring_members()[1];              // C2
            o8 = next_sugar.anomeric_substituent();         // O8 usually
        }
        else
        {
            std::cout << "ERROR: Unsupported ring size, expecting either a 5 membered or 6 membered ring." << std::endl;
            std::cout << "\tRing size received: " << next_sugar.ring_members().size() << std::endl;
            std::cout << "\tnext_sugar info: " << next_sugar.type().trim() << "-" << next_sugar.id().trim() << std::endl;
            throw std::runtime_error("Unsupported ring size");
        }

        if(first_sugar.ring_members().size() == 6)
        {
            c6 = first_sugar.ring_members()[5];
            first_sugar_ring_oxygen = first_sugar.ring_members()[0]; // O6 for omega_seven
            c7 = first_sugar[first_sugar.lookup("C7",clipper::MM::ANY)];
            o7 = first_sugar[first_sugar.lookup("O7",clipper::MM::ANY)];
            c8 = first_sugar[first_sugar.lookup("C8",clipper::MM::ANY)];
            
            // Temporary workaround until more general solution can be generated...
            if(first_sugar.type().trim() == "SIA" || first_sugar.type().trim() == "SLB")
            {
                c9 = first_sugar[first_sugar.lookup("C9",clipper::MM::ANY)];
                o9 = first_sugar[first_sugar.lookup("O9",clipper::MM::ANY)];
            }

        }
        else if(first_sugar.ring_members().size() == 5)
        {
            c6 = first_sugar.ring_members()[4];
            first_sugar_ring_oxygen = first_sugar.ring_members()[0]; // O6 for omega_seven
            c7 = first_sugar[first_sugar.lookup("C7",clipper::MM::ANY)];
            o7 = first_sugar[first_sugar.lookup("O7",clipper::MM::ANY)];
            c8 = first_sugar[first_sugar.lookup("C8",clipper::MM::ANY)];


            if(first_sugar.type().trim() == "SIA" || first_sugar.type().trim() == "SLB")
            {
                c9 = first_sugar[first_sugar.lookup("C9",clipper::MM::ANY)];
                o9 = first_sugar[first_sugar.lookup("O9",clipper::MM::ANY)];
            }
        }
        else
        {
            std::cout << "ERROR: Unsupported ring size, expecting either a 5 membered or 6 membered ring." << std::endl;
            std::cout << "\tRing size received: " << first_sugar.ring_members().size() << std::endl;
            std::cout << "\tfirst_sugar info: " << first_sugar.type().trim() << "-" << first_sugar.id().trim() << std::endl;
            throw std::runtime_error("Unsupported ring size");
        }

        phi   = clipper::Coord_orth::torsion ( o6.coord_orth(),
                                               c2.coord_orth(),
                                               o8.coord_orth(),
                                               c8.coord_orth() );
        
        phi_cone_ctwo_oeight_ceight   = clipper::Coord_orth::torsion (  c1.coord_orth(),
                                                                        c2.coord_orth(),
                                                                        o8.coord_orth(),
                                                                        c8.coord_orth() );

        psi   = clipper::Coord_orth::torsion ( c2.coord_orth(),
                                               o8.coord_orth(),
                                               c8.coord_orth(),
                                               c7.coord_orth() );

        omega_seven = clipper::Coord_orth::torsion (    o7.coord_orth(),
                                                        c7.coord_orth(),
                                                        c6.coord_orth(),
                                                        first_sugar_ring_oxygen.coord_orth() );

        omega_eight = clipper::Coord_orth::torsion (    o8.coord_orth(),
                                                        c8.coord_orth(),
                                                        c7.coord_orth(),
                                                        o7.coord_orth() );


        if(first_sugar.type().trim() == "SIA")
        {
            omega_nine = clipper::Coord_orth::torsion (     o9.coord_orth(),
                                                            c9.coord_orth(),
                                                            c8.coord_orth(),
                                                            o8.coord_orth() );
        }
        
        new_connection.set_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), clipper::Util::rad2d(omega), clipper::Util::rad2d(omega_six), clipper::Util::rad2d(omega_seven), clipper::Util::rad2d(omega_eight), clipper::Util::rad2d(omega_nine),  clipper::Util::rad2d(phi_cone_ctwo_oeight_ceight) );
        if(debug_output)
        {
            std::vector<clipper::ftype32> torsions = new_connection.get_torsions();
            DBG << "Torsions for link = " << link << ", phi = " << torsions[0] << "\t\tphi_c1c2o8c8 = " << torsions[5] << "\t\tpsi = " << torsions[1] << "\t\tomega_seven = " << torsions[2] << "\t\tomega_eight = " << torsions[3] << "\t\tomega_nine = " << torsions[4] << std::endl;
        }
        new_connection.set_linkage_atoms(donorAtom, acceptorAtom);
        add_torsions_for_plots(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom);
    }
    else if ( link == 7 )
    {
        // next_sugar 6 GMH - MAN
        // first_sugar 3 GM0 - BMA
        if(next_sugar.ring_members().size() == 6)
        {
            o5 = next_sugar.ring_members()[0];              // O5
            c1 = next_sugar.ring_members()[1];              // C1
            o7 = next_sugar.anomeric_substituent();         // O7 usually
        }
        else if(next_sugar.ring_members().size() == 5)
        {
            o5 = next_sugar.ring_members()[0];              // O5
            // c1 = next_sugar[next_sugar.lookup("C1",clipper::MM::ANY)]; // not going to be in the ring, but going to form the glycosidic bond regardless
            c1 = next_sugar.ring_members()[1];              // C2
            o7 = next_sugar.anomeric_substituent();         // O7 usually
        }
        else
        {
            std::cout << "ERROR: Unsupported ring size, expecting either a 5 membered or 6 membered ring." << std::endl;
            std::cout << "\tRing size received: " << next_sugar.ring_members().size() << std::endl;
            std::cout << "\tnext_sugar info: " << next_sugar.type().trim() << "-" << next_sugar.id().trim() << std::endl;
            throw std::runtime_error("Unsupported ring size");
        }

        if(first_sugar.ring_members().size() == 6)
        {
            c5 = first_sugar.ring_members()[5];
            c6 = first_sugar[first_sugar.lookup("C6",clipper::MM::ANY)];
            o6 = first_sugar[first_sugar.lookup("O6",clipper::MM::ANY)];
            c7 = first_sugar[first_sugar.lookup("C7",clipper::MM::ANY)];
            first_sugar_ring_oxygen = first_sugar.ring_members()[0]; // O5 for omega
        }
        else if(first_sugar.ring_members().size() == 5)
        {
            c5 = first_sugar.ring_members()[4];
            c6 = first_sugar[first_sugar.lookup("C6",clipper::MM::ANY)];
            o6 = first_sugar[first_sugar.lookup("O6",clipper::MM::ANY)];
            c7 = first_sugar[first_sugar.lookup("C7",clipper::MM::ANY)];
            first_sugar_ring_oxygen = first_sugar.ring_members()[0]; // O5 for omega
        }
        else
        {
            std::cout << "ERROR: Unsupported ring size, expecting either a 5 membered or 6 membered ring." << std::endl;
            std::cout << "\tRing size received: " << first_sugar.ring_members().size() << std::endl;
            std::cout << "\tfirst_sugar info: " << first_sugar.type().trim() << "-" << first_sugar.id().trim() << std::endl;
            throw std::runtime_error("Unsupported ring size");
        }

        phi   = clipper::Coord_orth::torsion ( o5.coord_orth(),
                                               c1.coord_orth(),
                                               o7.coord_orth(),
                                               c7.coord_orth() );

        psi   = clipper::Coord_orth::torsion ( c1.coord_orth(),
                                               o7.coord_orth(),
                                               c7.coord_orth(),
                                               c6.coord_orth() );

        // 6-7 = 7
        // 7-8 = 8
        // 8-9 = 9

        omega_six = clipper::Coord_orth::torsion (  o6.coord_orth(),
                                                    c6.coord_orth(),
                                                    c5.coord_orth(),
                                                    first_sugar_ring_oxygen.coord_orth() );

        omega_seven = clipper::Coord_orth::torsion (    o7.coord_orth(),
                                                        c7.coord_orth(),
                                                        c6.coord_orth(),
                                                        o6.coord_orth() );


        new_connection.set_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), clipper::Util::rad2d(omega), clipper::Util::rad2d(omega_six), clipper::Util::rad2d(omega_seven), clipper::Util::rad2d(omega_eight), clipper::Util::rad2d(omega_nine),  clipper::Util::rad2d(phi_cone_ctwo_oeight_ceight) );
        if(debug_output)
        {
            std::vector<clipper::ftype32> torsions = new_connection.get_torsions();
            DBG << "Torsions for link = " << link << ", phi = " << torsions[0] << "\t\tpsi = " << torsions[1] << "\t\tomega_six = " << torsions[2] << "\t\tomega_seven = " << torsions[3] << std::endl; 
        }
        new_connection.set_linkage_atoms(donorAtom, acceptorAtom);
        add_torsions_for_plots(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom);
    }
    else if ( link == 6 )
    {
        if(next_sugar.ring_members().size() == 6)
        {
            o5 = next_sugar.ring_members()[0];              // O5
            c1 = next_sugar.ring_members()[1];              // C1
            o6 = next_sugar.anomeric_substituent();         // O6 usually
        }
        else if(next_sugar.ring_members().size() == 5)
        {
            o5 = next_sugar.ring_members()[0];              // O5
            // c1 = next_sugar[next_sugar.lookup("C1",clipper::MM::ANY)]; // not going to be in the ring, but going to form the glycosidic bond regardless
            c1 = next_sugar.ring_members()[1];              // C2
            o6 = next_sugar.anomeric_substituent();         // O6 usually
        }
        else
        {
            std::cout << "ERROR: Unsupported ring size, expecting either a 5 membered or 6 membered ring." << std::endl;
            std::cout << "\tRing size received: " << next_sugar.ring_members().size() << std::endl;
            std::cout << "\tnext_sugar info: " << next_sugar.type().trim() << "-" << next_sugar.id().trim() << std::endl;
            throw std::runtime_error("Unsupported ring size");
        }

        if(first_sugar.ring_members().size() == 6)
        {
            c5 = first_sugar.ring_members()[5];             // C5 - for pyranose
            c6 = first_sugar.configurational_substituent(); // C6
            first_sugar_ring_oxygen = first_sugar.ring_members()[0]; // O5 for omega
        }
        else if(first_sugar.ring_members().size() == 5)
        {
            c5 = first_sugar.ring_members()[4];             // C5 - for furanose
            c6 = first_sugar.configurational_substituent(); // C6
            first_sugar_ring_oxygen = first_sugar.ring_members()[0]; // O5 for omega
        }
        else
        {
            std::cout << "ERROR: Unsupported ring size, expecting either a 5 membered or 6 membered ring." << std::endl;
            std::cout << "\tRing size received: " << first_sugar.ring_members().size() << std::endl;
            std::cout << "\tfirst_sugar info: " << first_sugar.type().trim() << "-" << first_sugar.id().trim() << std::endl;
            throw std::runtime_error("Unsupported ring size");
        }

        phi   = clipper::Coord_orth::torsion ( o5.coord_orth(),
                                               c1.coord_orth(),
                                               o6.coord_orth(),
                                               c6.coord_orth() );

        psi   = clipper::Coord_orth::torsion ( c1.coord_orth(),
                                               o6.coord_orth(),
                                               c6.coord_orth(),
                                               c5.coord_orth() );

        omega = clipper::Coord_orth::torsion (  first_sugar_ring_oxygen.coord_orth(),
                                                c5.coord_orth(),
                                                c6.coord_orth(),
                                                o6.coord_orth() );

        new_connection.set_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), clipper::Util::rad2d(omega), clipper::Util::rad2d(omega_six), clipper::Util::rad2d(omega_seven), clipper::Util::rad2d(omega_eight), clipper::Util::rad2d(omega_nine),  clipper::Util::rad2d(phi_cone_ctwo_oeight_ceight) );
        if(debug_output)
        {
            std::vector<clipper::ftype32> torsions = new_connection.get_torsions();
            DBG << "Torsions for link = " << link << ", phi = " << torsions[0] << "\t\tpsi = " << torsions[1] << "\t\tomega = " << torsions[2] << std::endl; 
        }
        new_connection.set_linkage_atoms(donorAtom, acceptorAtom);
        add_torsions_for_plots(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom);
    }
    else if ( link == 5 )
    {
        // next_sugar GMH - BMA
        // first_sugar KDO - NAG
        if(next_sugar.ring_members().size() == 6)
        {
            next_sugar_ring_oxygen = next_sugar.ring_members()[0];  // O5
            c1 = next_sugar.ring_members()[1];              // C1
            o5 = next_sugar.anomeric_substituent();         // O5 usually
        }
        else if(next_sugar.ring_members().size() == 5)
        {
            next_sugar_ring_oxygen = next_sugar.ring_members()[0];  // O5
            // c1 = next_sugar[next_sugar.lookup("C1",clipper::MM::ANY)]; // not going to be in the ring, but going to form the glycosidic bond regardless
            c1 = next_sugar.ring_members()[1];              // C2
            o5 = next_sugar.anomeric_substituent();         // O4 usually
        }
        else
        {
            std::cout << "ERROR: Unsupported ring size, expecting either a 5 membered or 6 membered ring." << std::endl;
            std::cout << "\tRing size received: " << next_sugar.ring_members().size() << std::endl;
            std::cout << "\tnext_sugar info: " << next_sugar.type().trim() << "-" << next_sugar.id().trim() << std::endl;
            throw std::runtime_error("Unsupported ring size");
        }

        if(first_sugar.ring_members().size() == 6)
        {
            c5 = first_sugar.ring_members()[4];             // C5
            c6 = first_sugar.ring_members()[5];             // C6
        }
        else if(first_sugar.ring_members().size() == 5)
        {
            c5 = first_sugar.ring_members()[3];             // C5
            c6 = first_sugar.ring_members()[4];             // C6
        }
        else
        {
            std::cout << "ERROR: Unsupported ring size, expecting either a 5 membered or 6 membered ring." << std::endl;
            std::cout << "\tRing size received: " << first_sugar.ring_members().size() << std::endl;
            std::cout << "\tfirst_sugar info: " << first_sugar.type().trim() << "-" << first_sugar.id().trim() << std::endl;
            throw std::runtime_error("Unsupported ring size");
        }

        phi   = clipper::Coord_orth::torsion ( next_sugar_ring_oxygen.coord_orth(),
                                               c1.coord_orth(),
                                               o5.coord_orth(),
                                               c5.coord_orth() );

        psi   = clipper::Coord_orth::torsion ( c1.coord_orth(),
                                               o5.coord_orth(),
                                               c5.coord_orth(),
                                               c6.coord_orth() );

        new_connection.set_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), clipper::Util::rad2d(omega), clipper::Util::rad2d(omega_six), clipper::Util::rad2d(omega_seven), clipper::Util::rad2d(omega_eight), clipper::Util::rad2d(omega_nine),  clipper::Util::rad2d(phi_cone_ctwo_oeight_ceight) );
        if(debug_output)
        {
            std::vector<clipper::ftype32> torsions = new_connection.get_torsions();
            DBG << "Torsions for link = " << link << ", phi = " << torsions[0] << "\t\tpsi = " << torsions[1] << std::endl;  
        }
        new_connection.set_linkage_atoms(donorAtom, acceptorAtom);
        add_torsions_for_plots(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom);
    }
    else if ( link == 4 )
    {
        if(next_sugar.ring_members().size() == 6)
        {
            o5 = next_sugar.ring_members()[0];              // O5
            c1 = next_sugar.ring_members()[1];              // C1
            o4 = next_sugar.anomeric_substituent();         // O4 usually
        }
        else if(next_sugar.ring_members().size() == 5)
        {
            o5 = next_sugar.ring_members()[0];              // O5
            // c1 = next_sugar[next_sugar.lookup("C1",clipper::MM::ANY)]; // not going to be in the ring, but going to form the glycosidic bond regardless
            c1 = next_sugar.ring_members()[1];              // C2
            o4 = next_sugar.anomeric_substituent();         // O4 usually
        }
        else
        {
            std::cout << "ERROR: Unsupported ring size, expecting either a 5 membered or 6 membered ring." << std::endl;
            std::cout << "\tRing size received: " << next_sugar.ring_members().size() << std::endl;
            std::cout << "\tnext_sugar info: " << next_sugar.type().trim() << "-" << next_sugar.id().trim() << std::endl;
            throw std::runtime_error("Unsupported ring size");
        }

        if(first_sugar.ring_members().size() == 6)
        {
            c4 = first_sugar.ring_members()[4];             // C4
            c5 = first_sugar.ring_members()[5];             // C5
        }
        else if(first_sugar.ring_members().size() == 5)
        {
            c4 = first_sugar.ring_members()[3];             // C4
            c5 = first_sugar.ring_members()[4];             // C5
        }
        else
        {
            std::cout << "ERROR: Unsupported ring size, expecting either a 5 membered or 6 membered ring." << std::endl;
            std::cout << "\tRing size received: " << first_sugar.ring_members().size() << std::endl;
            std::cout << "\tfirst_sugar info: " << first_sugar.type().trim() << "-" << first_sugar.id().trim() << std::endl;
            throw std::runtime_error("Unsupported ring size");
        }

        phi   = clipper::Coord_orth::torsion ( o5.coord_orth(),
                                               c1.coord_orth(),
                                               o4.coord_orth(),
                                               c4.coord_orth() );

        psi   = clipper::Coord_orth::torsion ( c1.coord_orth(),
                                               o4.coord_orth(),
                                               c4.coord_orth(),
                                               c5.coord_orth() );

        new_connection.set_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), clipper::Util::rad2d(omega), clipper::Util::rad2d(omega_six), clipper::Util::rad2d(omega_seven), clipper::Util::rad2d(omega_eight), clipper::Util::rad2d(omega_nine),  clipper::Util::rad2d(phi_cone_ctwo_oeight_ceight) );
        if(debug_output)
        {
            std::vector<clipper::ftype32> torsions = new_connection.get_torsions();
            DBG << "Torsions for link = " << link << ", phi = " << torsions[0] << "\t\tpsi = " << torsions[1] << std::endl;  
        }
        new_connection.set_linkage_atoms(donorAtom, acceptorAtom);
        add_torsions_for_plots(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom);
    }
    else if ( link == 3 )
    {
        if(next_sugar.ring_members().size() == 6)
        {
            o5 = next_sugar.ring_members()[0];              // O5
            c1 = next_sugar.ring_members()[1];              // C1
            o3 = next_sugar.anomeric_substituent();         // O3 usually
        }
        else if(next_sugar.ring_members().size() == 5)
        {
            o5 = next_sugar.ring_members()[0];              // O5
            // c1 = next_sugar[next_sugar.lookup("C1",clipper::MM::ANY)]; // not going to be in the ring, but going to form the glycosidic bond regardless
            c1 = next_sugar.ring_members()[1];              // C2
            o3 = next_sugar.anomeric_substituent();         // O3 usually
        }
        else
        {
            std::cout << "ERROR: Unsupported ring size, expecting either a 5 membered or 6 membered ring." << std::endl;
            std::cout << "\tRing size received: " << next_sugar.ring_members().size() << std::endl;
            std::cout << "\tnext_sugar info: " << next_sugar.type().trim() << "-" << next_sugar.id().trim() << std::endl;
            throw std::runtime_error("Unsupported ring size");
        }

        if(first_sugar.ring_members().size() == 6)
        {
            c3 = first_sugar.ring_members()[3];             // C3
            c4 = first_sugar.ring_members()[4];             // C4
        }
        else if(first_sugar.ring_members().size() == 5)
        {
            c3 = first_sugar.ring_members()[2];             // C3
            c4 = first_sugar.ring_members()[3];             // C4
        }
        else
        {
            std::cout << "ERROR: Unsupported ring size, expecting either a 5 membered or 6 membered ring." << std::endl;
            std::cout << "\tRing size received: " << first_sugar.ring_members().size() << std::endl;
            std::cout << "\tfirst_sugar info: " << first_sugar.type().trim() << "-" << first_sugar.id().trim() << std::endl;
            throw std::runtime_error("Unsupported ring size");
        }


        phi   = clipper::Coord_orth::torsion ( o5.coord_orth(),
                                               c1.coord_orth(),
                                               o3.coord_orth(),
                                               c3.coord_orth() );

        psi   = clipper::Coord_orth::torsion ( c1.coord_orth(),
                                               o3.coord_orth(),
                                               c3.coord_orth(),
                                               c4.coord_orth() );

        new_connection.set_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), clipper::Util::rad2d(omega), clipper::Util::rad2d(omega_six), clipper::Util::rad2d(omega_seven), clipper::Util::rad2d(omega_eight), clipper::Util::rad2d(omega_nine),  clipper::Util::rad2d(phi_cone_ctwo_oeight_ceight) );
        if(debug_output)
        {
            std::vector<clipper::ftype32> torsions = new_connection.get_torsions();
            DBG << "Torsions for link = " << link << ", phi = " << torsions[0] << "\t\tpsi = " << torsions[1] << std::endl;
        }
        new_connection.set_linkage_atoms(donorAtom, acceptorAtom);
        add_torsions_for_plots(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom);
    }
    else if ( link == 2 )
    {
        if(next_sugar.ring_members().size() == 6)
        {
            o5 = next_sugar.ring_members()[0];              // O5
            c1 = next_sugar.ring_members()[1];              // C1
            o2 = next_sugar.anomeric_substituent();         // O2 usually
        }
        else if(next_sugar.ring_members().size() == 5)
        {
            o5 = next_sugar.ring_members()[0];              // O5
            // c1 = next_sugar[next_sugar.lookup("C1",clipper::MM::ANY)]; // not going to be in the ring, but going to form the glycosidic bond regardless
            c1 = next_sugar.ring_members()[1];              // C2
            o2 = next_sugar.anomeric_substituent();         // O2 usually
        }
        else
        {
            std::cout << "ERROR: Unsupported ring size, expecting either a 5 membered or 6 membered ring." << std::endl;
            std::cout << "\tRing size received: " << next_sugar.ring_members().size() << std::endl;
            std::cout << "\tnext_sugar info: " << next_sugar.type().trim() << "-" << next_sugar.id().trim() << std::endl;
            throw std::runtime_error("Unsupported ring size");
        }

        if(first_sugar.ring_members().size() == 6)
        {
            c2 = first_sugar.ring_members()[2];             // C2
            c3 = first_sugar.ring_members()[3];             // C3
        }
        else if(first_sugar.ring_members().size() == 5)
        {
            c2 = first_sugar.ring_members()[1];             // C2
            c3 = first_sugar.ring_members()[2];             // C3
        }
        else
        {
            std::cout << "ERROR: Unsupported ring size, expecting either a 5 membered or 6 membered ring." << std::endl;
            std::cout << "\tRing size received: " << first_sugar.ring_members().size() << std::endl;
            std::cout << "\tfirst_sugar info: " << first_sugar.type().trim() << "-" << first_sugar.id().trim() << std::endl;
            throw std::runtime_error("Unsupported ring size");
        }

        phi   = clipper::Coord_orth::torsion ( o5.coord_orth(),
                                               c1.coord_orth(),
                                               o2.coord_orth(),
                                               c2.coord_orth() );

        psi   = clipper::Coord_orth::torsion ( c1.coord_orth(),
                                               o2.coord_orth(),
                                               c2.coord_orth(),
                                               c3.coord_orth() );
        new_connection.set_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), clipper::Util::rad2d(omega), clipper::Util::rad2d(omega_six), clipper::Util::rad2d(omega_seven), clipper::Util::rad2d(omega_eight), clipper::Util::rad2d(omega_nine),  clipper::Util::rad2d(phi_cone_ctwo_oeight_ceight) );
        if(debug_output)
        {
            std::vector<clipper::ftype32> torsions = new_connection.get_torsions();
            DBG << "Torsions for link = " << link << ", phi = " << torsions[0] << "\t\tpsi = " << torsions[1] << std::endl;
        }
        new_connection.set_linkage_atoms(donorAtom, acceptorAtom);
        add_torsions_for_plots(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom);
    }
    else if ( link == 1 )
    {
        // next_sugar 5 GLC - 4 MAN
        // first_sugar 4 FRU - 3 BMA
        if(next_sugar.ring_members().size() == 6)
        {
            o5 = next_sugar.ring_members()[0];              // O5
            c1 = next_sugar.ring_members()[1];              // C1
            o1 = next_sugar.anomeric_substituent();         // O1 usually
        }
        else if(next_sugar.ring_members().size() == 5)
        {
            o5 = first_sugar.ring_members()[0];              // O5
            o5_next_sugar = next_sugar.ring_members()[0];
            // c1 = next_sugar[next_sugar.lookup("C1",clipper::MM::ANY)]; // not going to be in the ring, but going to form the glycosidic bond regardless
            c1 = next_sugar.ring_members()[1];              // C2
            actual_c1 = next_sugar[next_sugar.lookup("C1",clipper::MM::ANY)];
            o1 = next_sugar.anomeric_substituent();         // O1 usually
            c3 = next_sugar.ring_members()[2];             // C3
        }
        else
        {
            std::cout << "ERROR: Unsupported ring size, expecting either a 5 membered or 6 membered ring." << std::endl;
            std::cout << "\tRing size received: " << next_sugar.ring_members().size() << std::endl;
            std::cout << "\tnext_sugar info: " << next_sugar.type().trim() << "-" << next_sugar.id().trim() << std::endl;
            throw std::runtime_error("Unsupported ring size");
        }

        if(first_sugar.ring_members().size() == 6)
        {
            first_sugar_c1 = first_sugar.ring_members()[1];  // C1
            c2 = first_sugar.ring_members()[2];             // C2
        }
        else if(first_sugar.ring_members().size() == 5)
        {
            first_sugar_c1 = first_sugar.ring_members()[1];  // C2
            c2 = first_sugar.ring_members()[2];             // C3
        }
        else
        {
            std::cout << "ERROR: Unsupported ring size, expecting either a 5 membered or 6 membered ring." << std::endl;
            std::cout << "\tRing size received: " << first_sugar.ring_members().size() << std::endl;
            std::cout << "\tfirst_sugar info: " << first_sugar.type().trim() << "-" << first_sugar.id().trim() << std::endl;
            throw std::runtime_error("Unsupported ring size");
        }

        if(next_sugar.ring_members().size() == 6)
        {

            phi   = clipper::Coord_orth::torsion (  o5.coord_orth(),
                                                    c1.coord_orth(),
                                                    o1.coord_orth(),
                                                    first_sugar_c1.coord_orth() );

            psi   = clipper::Coord_orth::torsion (  c1.coord_orth(),
                                                    o1.coord_orth(),
                                                    first_sugar_c1.coord_orth(),
                                                    c2.coord_orth() );

            new_connection.set_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), clipper::Util::rad2d(omega), clipper::Util::rad2d(omega_six), clipper::Util::rad2d(omega_seven), clipper::Util::rad2d(omega_eight), clipper::Util::rad2d(omega_nine),  clipper::Util::rad2d(phi_cone_ctwo_oeight_ceight) );
            if(debug_output)
            {
                std::vector<clipper::ftype32> torsions = new_connection.get_torsions();
                DBG << "Torsions for when next_sugar has 6 members, link = " << link << ", phi = " << torsions[0] << "\t\tpsi = " << torsions[1] << std::endl;
            }
            new_connection.set_linkage_atoms(donorAtom, acceptorAtom);
            add_torsions_for_plots(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom);
        }
        else if(next_sugar.ring_members().size() == 5)
        {
            phi   = clipper::Coord_orth::torsion (  o5.coord_orth(),
                                                    first_sugar_c1.coord_orth(),
                                                    o1.coord_orth(),
                                                    actual_c1.coord_orth() );

            // Can't figure this shit out, to come back later. 
            psi   = clipper::Coord_orth::torsion (  c2.coord_orth(),
                                                    first_sugar_c1.coord_orth(),
                                                    o1.coord_orth(),
                                                    actual_c1.coord_orth() );
            
            omega = clipper::Coord_orth::torsion (  o1.coord_orth(),
                                                    actual_c1.coord_orth(),
                                                    c1.coord_orth(),
                                                    c3.coord_orth() );

            new_connection.set_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), clipper::Util::rad2d(omega), clipper::Util::rad2d(omega_six), clipper::Util::rad2d(omega_seven), clipper::Util::rad2d(omega_eight), clipper::Util::rad2d(omega_nine),  clipper::Util::rad2d(phi_cone_ctwo_oeight_ceight) );
            if(debug_output)
            {
                std::vector<clipper::ftype32> torsions = new_connection.get_torsions();
                DBG << "Torsions for when next_sugar has 5 members, link = " << link << ", phi = " << torsions[0] << "\t\tpsi = " << torsions[1] << "\t\tomega = " << torsions[2] << std::endl;
            }
            new_connection.set_linkage_atoms(donorAtom, acceptorAtom);
            add_torsions_for_plots(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom);
        }
    }

    sugars.push_back ( next_sugar );
    node_list[index].add_connection ( new_connection ); // add the new connection to the previous node

    return false;
}

void MGlycan::add_torsions_for_plots(float Phi, float Psi, clipper::String first_residue_name, clipper::MAtom first_atom, clipper::String second_residue_name, clipper::MAtom second_atom)
{
    if(!all_torsions_within_mglycan.empty())
    {
        auto search_result = std::find_if(all_torsions_within_mglycan.begin(), all_torsions_within_mglycan.end(), [first_residue_name, second_residue_name](MGlycanTorsionSummary& element)
        {
            return first_residue_name == element.first_residue_name && second_residue_name == element.second_residue_name;
        });

        if(search_result != std::end(all_torsions_within_mglycan))
        {
            MGlycanTorsionSummary& found_torsion_description = *search_result;
            found_torsion_description.torsions.push_back(std::make_pair(Phi, Psi));
            found_torsion_description.atoms.push_back(std::make_pair(first_atom, second_atom));
        }
        else
        {
            MGlycanTorsionSummary new_torsion;
            std::string type;
            if(clipper::data::is_amino_acid(first_residue_name) || clipper::data::is_amino_acid(second_residue_name))
                type = "protein-sugar";
            else
                type = "sugar-sugar";
            
            new_torsion.type = type;
            new_torsion.first_residue_name = first_residue_name;
            new_torsion.second_residue_name = second_residue_name;
            new_torsion.atoms.push_back(std::make_pair(first_atom, second_atom));
            new_torsion.torsions.push_back(std::make_pair(Phi, Psi));

            all_torsions_within_mglycan.push_back(new_torsion);
        }
    }
    else
    {
        MGlycanTorsionSummary first_torsion;
        std::string type;
        if(clipper::data::is_amino_acid(first_residue_name) || clipper::data::is_amino_acid(second_residue_name))
            type = "protein-sugar";
        else
            type = "sugar-sugar";
        
        first_torsion.type = type;
        first_torsion.first_residue_name = first_residue_name;
        first_torsion.second_residue_name = second_residue_name;
        first_torsion.atoms.push_back(std::make_pair(first_atom, second_atom));
        first_torsion.torsions.push_back(std::make_pair(Phi, Psi));

        all_torsions_within_mglycan.push_back(first_torsion);
    }
}


void MGlycan::set_annotations ( std::string expression_system )
{
    // supported expression systems: fungal, yeast, plant, insect, mammalian and human

    // ROOT CHECKS

    if ( get_type() == "n-glycan" )
    {
        if ( clipper::data::carbname_of(node_list[0].get_sugar().type().trim()) == "GlcNAc" )
        {
            if ( node_list[0].get_sugar().anomer() != "beta" )
                add_link_annotation ( " Warning: n-GlcNAc linkage should be beta! " );
        }
    }
    else if ( get_type() == "o-glycan" )
    {
        if ( clipper::data::carbname_of(node_list[0].get_sugar().type().trim()) == "GalNAc" )
        {
            if ( node_list[0].get_sugar().anomer() != "alpha" )
                add_link_annotation ( " Warning: o-GalNAc linkage should be alpha! " );
        }
        else if ( clipper::data::carbname_of(node_list[0].get_sugar().type().trim()) == "GlcNAc" )
        {
            if ( node_list[0].get_sugar().anomer() != "beta" )
                add_link_annotation ( " Warning: o-GlcNAc linkage should be beta! " );
        }
        else if ( get_type() == "c-glycan" )
        {
          if ( clipper::data::carbname_of(node_list[0].get_sugar().type().trim()) == "Man" )
          {
              if ( node_list[0].get_sugar().anomer() != "alpha" )
                  add_link_annotation ( " Warning: Man-Trp linkage should be alpha! " );
          }
        }
        return; // no more o-glycosylation validation for now; will come back to this later
    }


    // LEVEL 1 CHECKS

    int n_conn = node_list[0].number_of_connections();
    Node second_glcnac;


    for ( int i = 0 ; i < n_conn; i++ )
    {
        Node actual_node = node_list[node_list[0].get_connection(i).get_linked_node_id()];

        if ( clipper::data::carbname_of(actual_node.get_sugar().type().trim()) == "GlcNAc" )
        {
            second_glcnac = actual_node;
            if ( node_list[0].get_connection(i).get_order() != 4 || node_list[0].get_connection(i).get_anomericity() != "beta" )
                    node_list[0].get_connection(i).add_annotation ( " Warning: this GlcNAc-GlcNAc linkage should be beta 1-4 " );
        }
        else if ( clipper::data::carbname_of(actual_node.get_sugar().type().trim()) == "Fuc" )
        {
            if ( expression_system == "mammalian" || expression_system == "human" )
            {
                if ( node_list[0].get_connection(i).get_order() != 6 || node_list[0].get_connection(i).get_anomericity() != "alpha" )
                {
                    if(debug_output)
                    {
                        DBG << std::endl << "Wrong core linkage" << std::endl;
                    }
                    node_list[0].get_connection(i).add_annotation ( " Warning: this GlcNAc-Fuc linkage should be alpha 1-6 " );
                }
            }
            else if ( expression_system == "plant" || expression_system == "insect" )
            {
                if ( node_list[0].get_connection(i).get_order() != 3 || node_list[0].get_connection(i).get_anomericity() != "alpha" )
                {
                    node_list[0].get_connection(i).add_annotation ( " Warning: this GlcNAc-Fuc linkage should be alpha 1-3 " );
                }
            }
        }
        else // unexpected sugar linked to the core
        {
            node_list[0].get_connection(i).add_annotation ( " Warning: unusual core-linked sugar. Use mass spectrometry to confirm!" );
            actual_node.add_annotation ( " Warning: unusual core-linked sugar. Use mass spectrometry to confirm!" );
        }
    }

    Node branching_mannose;

    if ( second_glcnac.is_initialised() ) // move to the second GlcNAc
    {
        int n_conn = second_glcnac.number_of_connections();

        for ( int i = 0 ; i < n_conn; i++ )
        {
            Node actual_node = node_list[second_glcnac.get_connection(i).get_linked_node_id()];

            if ( clipper::data::carbname_of(actual_node.get_sugar().type().trim()) == "Man" )
            {
                branching_mannose = actual_node;
                if ( second_glcnac.get_connection(i).get_order() != 4 || second_glcnac.get_connection(i).get_anomericity() != "beta" )
                    second_glcnac.get_connection(i).add_annotation ( " Warning: this GlcNAc-Man linkage should be beta 1-4 " );
            }
            else // unexpected sugar linked to the core
            {
                second_glcnac.get_connection(i).add_annotation ( " Warning: unusual core-linked sugar. Use mass spectrometry to confirm!" );
                actual_node.add_annotation ( " Warning: unusual core-linked sugar. Use mass spectrometry to confirm!" );
            }
        }
    }
}

/*
Added by Haroldas Bagdonas (hb1115@york.ac.uk) on 03/01/2020
Function used to convert an int number to a char letter. Needed for generation of WURCS strings, specifically to conver residue IDs to letters and leave numbers
for linkage information.

Last modified on: 03/01/2020
*/
char MGlycan::convertNumberToLetter(int number)
{
    std::string alphabet("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
    return alphabet.at(number % alphabet.size());
}

/*
Added by Haroldas Bagdonas (hb1115@york.ac.uk) on 03/01/2020
Function used to obtain a std::vector populated with WURCS residue descriptions of unique carbohydrate monomers in the glycan chain

Performance considerations: uses a non ideal implementation of making sure that only unique strings are added to std::vector
                            via std::find function. This gives O(N*N) complexity. However, this method is used regardless for the
                            following reasons:
                            1.) This std::vector is unlikely to ever be significantly huge, thus unlikely to cause a huge performance penalty.
                            2.) WURCS requires that the sequence of items pushed into vector are maintained - we don't want
                            to rearrange the order of individual items within the std::vector.

Last modified on: 03/01/2020
*/
std::vector < std::string > MGlycan::obtain_unique_WURCS_residues()
{
    std::vector < std::string > uniqueResidues;

    for(int i = 0; i < node_list.size(); i++)
    {
        clipper::MSugar msug;
        std::string msug_wurcs_string;

        msug = node_list[i].get_sugar();
        msug_wurcs_string = clipper::data::convert_to_wurcs_residue_code ( msug.type().trim() );

        if (std::find(uniqueResidues.begin(), uniqueResidues.end(), msug_wurcs_string) == uniqueResidues.end()) {
            uniqueResidues.push_back(msug_wurcs_string);
        }
    }
    return uniqueResidues;
}

std::vector < std::string > MGlycan::obtain_unique_residue_codes()
{
    std::vector < std::string > uniqueResidues;
    
    for(int i = 0; i < node_list.size(); i++)
    {   
        clipper::MSugar msug;
        std::string msug_code_string;

        msug = node_list[i].get_sugar();
        msug_code_string =  msug.type().trim(); 

        if (std::find(uniqueResidues.begin(), uniqueResidues.end(), msug_code_string) == uniqueResidues.end()) {
            uniqueResidues.push_back(msug_code_string);
        }
    }
    return uniqueResidues;
}

/*
Added by Haroldas Bagdonas (hb1115@york.ac.uk) on 03/01/2020
Function used to obtain a total number of glycosidic bonds within the chain. Needed for generation of WURCS strings.

Last modified on: 03/01/2020
*/
const int MGlycan::obtain_total_number_of_glycosidic_bonds()
{
    int totalConnections = 0;

    for(int i = 0; i < node_list.size(); i++)
    {
        int numOfConnectionsPerResidue = node_list[i].number_of_connections();
        totalConnections += numOfConnectionsPerResidue;
    }
    return totalConnections;
}

/*

Added by Haroldas Bagdonas (hb1115@york.ac.uk) on 02/01/2020
Function that produces a WURCS string describing a Glycan sequence in the Glycoprotein Chain.

WURCS=Version/Unique Residue Count, Chain Length Count, Number of linkages between monomers/[Unique Monomer Description 1][Unique Monomer Description 2]/Unique Monomer 1 at position A-Unique Monomer 1 at position B/Glycosidic bond between Monomer A from position 4 connecting to Monomer B at position 1.
WURCS=2.0/5,9,8/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O]/1-1-2-3-1-3-1-4-5/a4-b1_b4-c1_c3-d1_c6-f1_d4-e1_f4-g1_g4-h1_h6-i2


TEST CASE:
PDB ID:                                 3v8x
GlyTouCan ID for glycan on Chain A:     G98736SM
GlyTouCan ID for glycan on Chain B:     G74608QW

WURCS2GTC DEMO:
Glycan on Chain A:
https://api.glycosmos.org/glytoucan/sparql/wurcs2gtcids?wurcs=WURCS=2.0/5,10,9/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O]/1-1-2-3-1-4-5-3-1-4/a4-b1_b4-c1_c3-d1_c6-h1_d4-e1_e4-f1_f6-g2_h4-i1_i4-j1
Chain on Chain B:
https://api.glycosmos.org/glytoucan/sparql/wurcs2gtcids?wurcs=WURCS=2.0/5,9,8/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O]/1-1-2-3-1-3-1-4-5/a4-b1_b4-c1_c3-d1_c6-f1_d4-e1_f4-g1_g4-h1_h6-i2

12/02/2020: All previously noted issues have been addressed. The current version of the code should be much more robust and not have any significant bugs.

Last modified on: 12/02/2020

*/
clipper::String MGlycan::generate_wurcs()
{
    int glycanLength = sugars.size();
    const int numberOfConnections = obtain_total_number_of_glycosidic_bonds();
    std::vector<std::string> uniqueResidueList = obtain_unique_WURCS_residues();

    clipper::String wurcs_string = "WURCS=";

    // WURCS version
    wurcs_string += "2.0";
    wurcs_string += "/";

    // Unit Count(Unique Residue Count, Chain Length Count, Number of linkages between monomers)

    wurcs_string += std::to_string(uniqueResidueList.size()) + "," + std::to_string(glycanLength) + "," + std::to_string(numberOfConnections) + "/";

    // Unique Monomer description loop.
    for (int i = 0; i < uniqueResidueList.size(); i++)
    {
        wurcs_string += "[" + uniqueResidueList[i] + "]";
    }
    wurcs_string += "/";

    // Add first monomer ID.

    wurcs_string += "1";

    if (node_list.size() < 2)
    {
        wurcs_string += "/";
        return wurcs_string; // if only a single monomer, linkage information does not exist.
    } // If the glycan chain has more than 1 monomer, add other monomers and include linkage information between monomers.
    else
    {
        clipper::MSugar msug;
        int connectedToNodeID;
        std::ostringstream linkagePositionInitial;

        wurcs_string += "-";

        // Add sequence information, by relating to monomer descriptions contained within std::vector < std::string > uniqueResidueList.
        for (int i = 1; i < node_list.size(); i++)
        {
            std::string msug_wurcs_string;

            msug = node_list[i].get_sugar();
            msug_wurcs_string = clipper::data::convert_to_wurcs_residue_code(msug.type().trim());

            std::vector<std::string>::iterator residueAssigner = std::find(uniqueResidueList.begin(), uniqueResidueList.end(), msug_wurcs_string);

            int residueID = std::distance(uniqueResidueList.begin(), residueAssigner);

            wurcs_string += std::to_string(residueID + 1);
            if (i < (node_list.size() - 1))
                wurcs_string += "-";
        }

        if(debug_output)
        {
            DBG << "Type of sugar via ::MSugar.full_type() = " << msug.full_type() << std::endl;
            DBG << "Number of connections for msug/node_list[0]: " << node_list[0].number_of_connections() << std::endl
                << std::endl;
        }

        wurcs_string += "/";

        // Add linkage information

        if (node_list[0].number_of_connections() > 0)
        {
            for (int j = 0; j < node_list[0].number_of_connections(); j++ )
            {
            std::ostringstream linkagePosition;
            connectedToNodeID = node_list[0].get_connection(j).get_linked_node_id();
            msug = node_list[connectedToNodeID].get_sugar();
            linkagePosition << node_list[0].get_connection(j).get_order();
            std::pair<clipper::MAtom, clipper::MAtom> linkage_atoms = node_list[0].get_connection(j).get_linkage_atoms();

            wurcs_string += convertNumberToLetter(0);
            wurcs_string += linkagePosition.str();

            wurcs_string += "-";

            wurcs_string += convertNumberToLetter(connectedToNodeID);
            if (msug.full_type() == "ketose")
            {
                wurcs_string += "2";
                if(linkage_atoms.first.element().trim() == "S" || linkage_atoms.first.element().trim() == "F" || linkage_atoms.first.element().trim() == "N")
                    wurcs_string += "*" + linkage_atoms.first.element().trim() + "*";
            }
            else
            {
                wurcs_string += "1";
                if(linkage_atoms.first.element().trim() == "S" || linkage_atoms.first.element().trim() == "F" || linkage_atoms.first.element().trim() == "N")
                    wurcs_string += "*" + linkage_atoms.first.element().trim() + "*";
            }
            

            wurcs_string += "_";
            }
        }

        // Describe the rest of the linkages.
        for (int i = 1; i < node_list.size(); i++)
        {
            msug = node_list[i].get_sugar();

            if(debug_output)
            {
                DBG << "Type of sugar via ::MSugar.full_type() = " << msug.full_type() << std::endl;
                DBG << "Number of connections for msug/node_list[" << i << "]: " << node_list[i].number_of_connections() << std::endl
                    << std::endl;
                DBG << "Residue code via ::MSugar.type().trim() = " << msug.type().trim() << std::endl;
            }

            if (node_list[i].number_of_connections() > 0)
            {
                for (int j = 0; j < node_list[i].number_of_connections(); j++)
                {
                    if(debug_output)
                    {
                        DBG << "Connection: " << j+1 << " out of " << node_list[i].number_of_connections() << std::endl;
                    }
                    
                    std::ostringstream linkagePosition;
                    connectedToNodeID = node_list[i].get_connection(j).get_linked_node_id();
                    
                    if(debug_output)
                    {
                        DBG << "connectedToNodeID: " << connectedToNodeID << std::endl;
                    }
                    msug = node_list[connectedToNodeID].get_sugar();
                    linkagePosition << node_list[i].get_connection(j).get_order();
                    std::pair<clipper::MAtom, clipper::MAtom> linkage_atoms = node_list[i].get_connection(j).get_linkage_atoms();
                    
                    if(debug_output)
                    {
                        DBG << "linkagePosition: " << linkagePosition.str() << std::endl;
                    }

                    wurcs_string += convertNumberToLetter(i);
                    wurcs_string += linkagePosition.str();

                    wurcs_string += "-";

                    wurcs_string += convertNumberToLetter(connectedToNodeID);
                    if (msug.full_type() == "ketose")
                    {
                        wurcs_string += "2";
                        if(linkage_atoms.first.element().trim() == "S" || linkage_atoms.first.element().trim() == "F" || linkage_atoms.first.element().trim() == "N")
                            wurcs_string += "*" + linkage_atoms.first.element().trim() + "*";
                    }
                    else
                    {
                        wurcs_string += "1";
                        if(linkage_atoms.first.element().trim() == "S" || linkage_atoms.first.element().trim() == "F" || linkage_atoms.first.element().trim() == "N")
                            wurcs_string += "*" + linkage_atoms.first.element().trim() + "*";
                    }
            

                    if (i < (node_list.size() - 2))
                        wurcs_string += "_";
                }
            }
        }
    }

    if (!wurcs_string.empty() && wurcs_string.back() == '_')
    {
        wurcs_string.pop_back();
    }

    // this->wurcs = wurcs_string;
    return wurcs_string;
}

std::string MGlycan::write_ring_ext_restraints ( float weight ) {

  std::string buffer = "";

  std::vector<clipper::MSugar> sugar_list = this->get_sugars();
  for ( int i = 0; i < sugar_list.size(); i++ ) {
    buffer += "# " + sugar_list[i].type() + " " + sugar_list[i].id() + "\n";
    std::string residue = sugar_list[i].id();
    std::string chain = this->get_chain();
    if ( this->kind_of_glycan == "c-glycan" ) { // needs 1C4 restraints
      buffer += "external torsion first chain " + chain + " residue " + residue + " atom O5 next" +
                                      " chain " + chain + " residue " + residue + " atom C1 next" +
                                      " chain " + chain + " residue " + residue + " atom C2 next" +
                                      " chain " + chain + " residue " + residue + " atom C3 value -55.71 sigma 0.1 period 1\n";

      buffer += "external torsion first chain " + chain + " residue " + residue + " atom C1 next" +
                                      " chain " + chain + " residue " + residue + " atom C2 next" +
                                      " chain " + chain + " residue " + residue + " atom C3 next" +
                                      " chain " + chain + " residue " + residue + " atom C4 value  51.72 sigma 0.1 period 1\n";

      buffer += "external torsion first chain " + chain + " residue " + residue + " atom C2 next" +
                                      " chain " + chain + " residue " + residue + " atom C3 next" +
                                      " chain " + chain + " residue " + residue + " atom C4 next" +
                                      " chain " + chain + " residue " + residue + " atom C5 value -47.55 sigma 0.1 period 1\n";

      buffer += "external torsion first chain " + chain + " residue " + residue + " atom C3 next" +
                                      " chain " + chain + " residue " + residue + " atom C4 next" +
                                      " chain " + chain + " residue " + residue + " atom C5 next" +
                                      " chain " + chain + " residue " + residue + " atom O5 value  45.67 sigma 0.1 period 1\n";

      buffer += "external torsion first chain " + chain + " residue " + residue + " atom C4 next" +
                                      " chain " + chain + " residue " + residue + " atom C5 next" +
                                      " chain " + chain + " residue " + residue + " atom O5 next" +
                                      " chain " + chain + " residue " + residue + " atom C1 value -51.06 sigma 0.1 period 1\n";

      buffer += "external torsion first chain " + chain + " residue " + residue + " atom C5 next" +
                                      " chain " + chain + " residue " + residue + " atom O5 next" +
                                      " chain " + chain + " residue " + residue + " atom C1 next" +
                                      " chain " + chain + " residue " + residue + " atom C2 value 56.33 sigma 0.1 period 1\n\n";
    }
    else {
      buffer += "external torsion first chain " + chain + " residue " + residue + " atom O5 next" +
                                      " chain " + chain + " residue " + residue + " atom C1 next" +
                                      " chain " + chain + " residue " + residue + " atom C2 next" +
                                      " chain " + chain + " residue " + residue + " atom C3 value  55.71 sigma 0.1 period 1\n";

      buffer += "external torsion first chain " + chain + " residue " + residue + " atom C1 next" +
                                      " chain " + chain + " residue " + residue + " atom C2 next" +
                                      " chain " + chain + " residue " + residue + " atom C3 next" +
                                      " chain " + chain + " residue " + residue + " atom C4 value -51.72 sigma 0.1 period 1\n";

      buffer += "external torsion first chain " + chain + " residue " + residue + " atom C2 next" +
                                      " chain " + chain + " residue " + residue + " atom C3 next" +
                                      " chain " + chain + " residue " + residue + " atom C4 next" +
                                      " chain " + chain + " residue " + residue + " atom C5 value  47.55 sigma 0.1 period 1\n";

      buffer += "external torsion first chain " + chain + " residue " + residue + " atom C3 next" +
                                      " chain " + chain + " residue " + residue + " atom C4 next" +
                                      " chain " + chain + " residue " + residue + " atom C5 next" +
                                      " chain " + chain + " residue " + residue + " atom O5 value -45.67 sigma 0.1 period 1\n";

      buffer += "external torsion first chain " + chain + " residue " + residue + " atom C4 next" +
                                      " chain " + chain + " residue " + residue + " atom C5 next" +
                                      " chain " + chain + " residue " + residue + " atom O5 next" +
                                      " chain " + chain + " residue " + residue + " atom C1 value  51.06 sigma 0.1 period 1\n";

      buffer += "external torsion first chain " + chain + " residue " + residue + " atom C5 next" +
                                      " chain " + chain + " residue " + residue + " atom O5 next" +
                                      " chain " + chain + " residue " + residue + " atom C1 next" +
                                      " chain " + chain + " residue " + residue + " atom C2 value -56.33 sigma 0.1 period 1\n\n";
    }
  }
  return buffer;
}

std::string MGlycan::write_link_ext_restraints ( float weight ) {

  std::string buffer = "";
  return buffer;
}


void MGlycan::remove_node_at_index ( int index )
{
    if (index>node_list.size()-1)
    {
        int lastElementID = node_list.size() - 1;


        for (int i = 0; i < node_list.size(); i++)
        {
            if (node_list[i].number_of_connections() > 0)
            {
                int lastNodeToDelete = -1;
                for (int j = 0; j < node_list[i].number_of_connections(); j++ )
                {
                    int connectedToNodeID = node_list[i].get_connection(j).get_linked_node_id();
                    if(connectedToNodeID > index)
                        {
                            node_list[i].get_connection(j).modify_linked_node_id(connectedToNodeID - 1);
                            continue;
                        }
                    if(connectedToNodeID == index) lastNodeToDelete = connectedToNodeID;
                }
                if (lastNodeToDelete != -1) node_list[i].remove_connection(lastNodeToDelete);
            }
        }
        node_list.pop_back();
        sugars.pop_back();
    }
    else
    {
        for (int i = 0; i < node_list.size(); i++)
        {
            if (node_list[i].number_of_connections() > 0)
            {
                int lastNodeToDelete = -1;
                for (int j = 0; j < node_list[i].number_of_connections(); j++ )
                {
                    int connectedToNodeID = node_list[i].get_connection(j).get_linked_node_id();
                    if(connectedToNodeID > index)
                        {
                            node_list[i].get_connection(j).modify_linked_node_id(connectedToNodeID - 1);
                            continue;
                        }
                    if(connectedToNodeID == index) lastNodeToDelete = connectedToNodeID;
                }
                if (lastNodeToDelete != -1) node_list[i].remove_connection(lastNodeToDelete);
            }
        }
        node_list.erase(node_list.begin() + index);
        sugars.erase(sugars.begin() + index);
    }
}

void MGlycan::replace_sugar_at_index ( int index, clipper::MSugar& donor )
{
    node_list[index].set_sugar(donor);
}

void MGlycan::update_msugar_in_root ( clipper::MSugar& newmsug )
{
    clipper::MMonomer cmonomer = root.first;

    root = std::make_pair(cmonomer, newmsug);
}

///////////////////////// MGlycology ///////////////////////////////



MGlycology::MGlycology ( const clipper::MiniMol& mmol, bool debug_output, std::string expression_system )
{
    const clipper::MAtomNonBond nb = MAtomNonBond ( mmol, 1.0 );

    this->debug_output = debug_output;
    pyinit( mmol, nb, debug_output, expression_system );
}


MGlycology::MGlycology ( const clipper::MiniMol& mmol, const clipper::MAtomNonBond& manb, bool debug_output, std::string expression_system )
{
    this->debug_output = debug_output;
    this->manb = &manb;
    this->mmol = &mmol;

    this->expression_system = expression_system;

    std::vector<std::pair<clipper::MMonomer, clipper::String>> potential_n_roots;
    std::vector<std::pair<clipper::MMonomer, clipper::String>> potential_o_roots;
    std::vector<std::pair<clipper::MMonomer, clipper::String>> potential_s_roots;
    std::vector<std::pair<clipper::MMonomer, clipper::String>> potential_c_roots;
    std::vector<std::pair<clipper::MMonomer, clipper::String>> potential_rootless_polysaccharides;
    std::vector<std::pair<clipper::MMonomer, clipper::String>> special_rootless_polysaccharides;
    

    for ( int pol = 0; pol < mmol.size() ; pol++ )
        for ( int mon = 0 ; mon < mmol[pol].size() ; mon++ )
        {
            std::vector<clipper::MAtom> KDOlikerecipe;
            for (int atom = 0; atom < mmol[pol][mon].size(); atom++)
            {
                clipper::MAtom tmpAtom = mmol[pol][mon][atom];
                if(get_altconf(tmpAtom) != ' ')
                {
                    const char altconf = get_altconf(tmpAtom);
                    
                    std::string altConfSymbol(1, altconf);
                    std::string O1_with_altconf = "O1 :" + altConfSymbol;
                    std::string C1_with_altconf = "C1 :" + altConfSymbol;
                    std::string OA_with_altconf = "O1A :" + altConfSymbol;
                    std::string OB_with_altconf = "O1B :" + altConfSymbol;

                    if (tmpAtom.id().trim() == O1_with_altconf)
                        potential_rootless_polysaccharides.push_back( std::make_pair(mmol[pol][mon], mmol[pol].id()) );
                    
                    if (tmpAtom.id().trim() == C1_with_altconf || tmpAtom.id().trim() == OA_with_altconf || tmpAtom.id().trim() == OB_with_altconf)
                        KDOlikerecipe.push_back(tmpAtom);
                }
                else
                    if (tmpAtom.id().trim() == "O1") 
                        potential_rootless_polysaccharides.push_back( std::make_pair(mmol[pol][mon], mmol[pol].id()) );
                    
                    if (tmpAtom.id().trim() == "C1" || tmpAtom.id().trim() == "O1A" || tmpAtom.id().trim() == "O1B")
                        KDOlikerecipe.push_back(tmpAtom);
            }
            if(KDOlikerecipe.size() == 3)
                special_rootless_polysaccharides.push_back( std::make_pair(mmol[pol][mon], mmol[pol].id()) );

            // Will need to keep this list up to date with the latest discoveries
            // To do: include check on other ligands, such as lipids (e.g. ceramide O-glycosylation)
            if ( mmol[pol][mon].type() == "ASN" ) potential_n_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // n-linked GlcNAc ?
            else if ( mmol[pol][mon].type() == "ARG" ) potential_n_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // Arginine rhamnosylation?
            else if ( mmol[pol][mon].type() == "THR" ) potential_o_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // o-linked stuff ?
            else if ( mmol[pol][mon].type() == "SER" ) potential_o_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) );
            else if ( mmol[pol][mon].type() == "LYS" ) potential_o_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) );
            else if ( mmol[pol][mon].type() == "TYR" ) potential_o_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) );
            else if ( mmol[pol][mon].type() == "CYS" ) potential_s_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // s-linked stuff ?
            else if ( mmol[pol][mon].type() == "TRP" ) potential_c_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // C-linked stuff for C/TRP-mannosylation
        }
    for(int i = 0; i < special_rootless_polysaccharides.size(); i++)
        potential_rootless_polysaccharides.push_back(special_rootless_polysaccharides[i]);
    for ( int i = 0 ; i < potential_n_roots.size() ; i++ )  // create n-glycan roots with first sugar
    {
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_n_roots[i].first ) ;

        for ( int j = 0 ; j < linked.size() ; j++ )
        {
            const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()];

            if ( tmpmon.type().trim() == "NAG" )
            {
                if ( clipper::MSugar::search_database( tmpmon.type().c_str() ) )
                {
                    clipper::MSugar sugar( mmol, mmol[linked[j].second.polymer()].id().trim(), tmpmon, manb, debug_output );
                    if(sugar.ring_members().size() == 5 || sugar.ring_members().size() == 6)
                    {
                        list_of_sugars.push_back ( sugar );

                        if(debug_output)
                        {
                            DBG << "Created the MSugar object" << std::endl;
                        }

                        if(debug_output)
                        {
                            DBG << "potential n roots is " << potential_n_roots[i].first.type() << std::endl;
                            DBG << "sugar is " << sugar.type() << std::endl;
                            DBG << "id is " << mmol[linked[j].second.polymer()].id().trim() << std::endl;
                        }
                        clipper::String root_sugar_chain_id = mmol[linked[j].second.polymer()].id().trim().substr(0,1);
                        clipper::MGlycan mg (   potential_n_roots[i].second,
                                                potential_n_roots[i].first,
                                                list_of_sugars.back(),
                                                root_sugar_chain_id,
                                                debug_output,
                                                this->expression_system );

                        if(debug_output)
                        {
                            DBG << "Exited the glycan constructor!" << std::endl;
                        }

                        mg.set_kind_of_glycan ( "n-glycan" );

                        if ( linked[j].second.monomer()+2 < mmol[linked[j].second.polymer()].size() )
                        {
                            if ( mmol[linked[j].second.polymer()][linked[j].second.monomer()+2].type().trim() != "THR" &&
                                mmol[linked[j].second.polymer()][linked[j].second.monomer()+2].type().trim() != "SER" )
                                // this is not a consensus ASN-GlcNAc glycosylation point
                                mg.add_root_annotation ( " Warning: this glycosylation point does not follow the Asn-X-Thr/Ser consensus sequence. ");
                        }
                        
                        clipper::MAtom o5 = sugar.ring_members()[0];              // O5
                        clipper::MAtom c1 = sugar.ring_members()[1];              // C1
                        clipper::MAtom nd2= sugar.anomeric_substituent();         // ND2 usually
                        clipper::MAtom cg = potential_n_roots[i].first.find("CG", clipper::MM::ANY);      // CG
                        clipper::MAtom cb = potential_n_roots[i].first.find("CB", clipper::MM::ANY);      // CB
                        clipper::ftype phi, psi;

                        phi   = clipper::Coord_orth::torsion (  o5.coord_orth(),
                                                                c1.coord_orth(),
                                                                nd2.coord_orth(),
                                                                cg.coord_orth() );

                        psi   = clipper::Coord_orth::torsion (  c1.coord_orth(),
                                                                nd2.coord_orth(),
                                                                cg.coord_orth(),
                                                                cb.coord_orth() );

                        if ( psi < 0 )
                            psi = clipper::Util::twopi() + psi;

                        mg.set_glycosylation_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi) );
                        mg.add_torsions_for_plots(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), potential_n_roots[i].first.type().trim(), nd2, sugar.type().trim(), c1);

                        list_of_glycans.push_back ( mg );
                        break;
                    }
                }
            }
            else if ( tmpmon.type().trim() == "NDG" )
            {
                if ( clipper::MSugar::search_database( tmpmon.type().c_str() ) )
                {
                    clipper::MSugar sugar( mmol, mmol[linked[j].second.polymer()].id().trim(), tmpmon, manb, debug_output );
                    if(sugar.ring_members().size() == 5 || sugar.ring_members().size() == 6)
                    {
                        list_of_sugars.push_back ( sugar );

                        if(debug_output)
                        {
                            DBG << "Created the MSugar object" << std::endl;
                        }

                        if(debug_output)
                        {
                            DBG << "potential n roots is " << potential_n_roots[i].first.type() << std::endl;
                            DBG << "sugar is " << sugar.type() << std::endl;
                            DBG << "id is " << mmol[linked[j].second.polymer()].id().trim() << std::endl;
                        }

                        clipper::String root_sugar_chain_id = mmol[linked[j].second.polymer()].id().trim().substr(0,1);
                        clipper::MGlycan mg (   potential_n_roots[i].second,
                                                potential_n_roots[i].first,
                                                list_of_sugars.back(),
                                                root_sugar_chain_id,
                                                debug_output,
                                                this->expression_system );
                        mg.set_kind_of_glycan ( "n-glycan" );

                        clipper::MAtom o5 = sugar.ring_members()[0];              // O5
                        clipper::MAtom c1 = sugar.ring_members()[1];              // C1
                        clipper::MAtom nd2= sugar.anomeric_substituent();         // ND2 usually
                        clipper::MAtom cg = potential_n_roots[i].first.find("CG", clipper::MM::ANY);      // CG
                        clipper::MAtom cb = potential_n_roots[i].first.find("CB", clipper::MM::ANY);      // CB
                        clipper::ftype phi, psi;

                        phi   = clipper::Coord_orth::torsion (  o5.coord_orth(),
                                                                c1.coord_orth(),
                                                                nd2.coord_orth(),
                                                                cg.coord_orth() );

                        psi   = clipper::Coord_orth::torsion (  c1.coord_orth(),
                                                                nd2.coord_orth(),
                                                                cg.coord_orth(),
                                                                cb.coord_orth() );

                        mg.set_glycosylation_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi) );
                        mg.add_torsions_for_plots(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), potential_n_roots[i].first.type().trim(), nd2, sugar.type().trim(), c1);


                        if ( linked[j].second.monomer()+2 < mmol[linked[j].second.polymer()].size() )
                        {
                            if ( mmol[linked[j].second.polymer()][linked[j].second.monomer()+2].type().trim() != "THR" &&
                                mmol[linked[j].second.polymer()][linked[j].second.monomer()+2].type().trim() != "SER" )
                                // this is not a consensus ASN-GlcNAc glycosylation point
                                mg.add_root_annotation ( " Warning: this glycosylation point does not follow the Asn-X-Thr/Ser consensus sequence. ");
                        }

                        list_of_glycans.push_back ( mg );

                        break;
                    }
                }
            }
        }
    }


    for ( int i = 0 ; i < potential_o_roots.size() ; i++ )  // create o-glycan roots with first sugar
    {
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_o_roots[i].first ) ;

        for ( int j = 0 ; j < linked.size() ; j++ )
        {
            const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()];

            if (( tmpmon.type().trim() == "NGA" ) || ( tmpmon.type().trim() == "A2G" ) || (tmpmon.type().trim() == "FUC" ) || (tmpmon.type().trim() == "RAM" ) || (tmpmon.type().trim() == "BGC" ) || (tmpmon.type().trim() == "BMA" ) || (tmpmon.type().trim() == "NAG" )
                || (tmpmon.type().trim() == "NDG" ))
            {
                if ( clipper::MSugar::search_database( tmpmon.type().c_str() ) )
                {
                    if(debug_output)
                    {
                        DBG << "Creating the MSugar object from potential_o_roots" << std::endl;
                    }
                    
                    clipper::MSugar sugar( mmol, mmol[linked[j].second.polymer()].id().trim(), tmpmon, manb, debug_output );
                    if(sugar.ring_members().size() == 5 || sugar.ring_members().size() == 6)
                    {
                        list_of_sugars.push_back ( sugar );

                        if(debug_output)
                        {
                            DBG << "Created the MSugar object" << std::endl;
                        }

                        if(debug_output)
                        {
                            DBG << "potential o roots is " << potential_o_roots[i].first.type() << std::endl;
                            DBG << "sugar is " << sugar.type() << std::endl;
                            DBG << "id is " << mmol[linked[j].second.polymer()].id().trim() << std::endl;
                        }

                        clipper::String root_sugar_chain_id = mmol[linked[j].second.polymer()].id().trim().substr(0,1);
                        clipper::MGlycan mg (   potential_o_roots[i].second,
                                                potential_o_roots[i].first,
                                                list_of_sugars.back(),
                                                root_sugar_chain_id,
                                                debug_output,
                                                this->expression_system );

                        mg.set_kind_of_glycan ( "o-glycan" );

                        clipper::MAtom o5 = sugar.ring_members()[0];              // O5
                        clipper::MAtom c1 = sugar.ring_members()[1];              // C1
                        clipper::MAtom og1= sugar.anomeric_substituent();         // OG/OG1 SER/THR
                        clipper::MAtom cg = potential_o_roots[i].first.find("CB", clipper::MM::ANY);      // CB
                        clipper::MAtom cb = potential_o_roots[i].first.find("CA", clipper::MM::ANY);      // CA
                        clipper::ftype phi, psi;

                        phi   = clipper::Coord_orth::torsion (  o5.coord_orth(),
                                                                c1.coord_orth(),
                                                                og1.coord_orth(),
                                                                cg.coord_orth() );

                        psi   = clipper::Coord_orth::torsion (  c1.coord_orth(),
                                                                og1.coord_orth(),
                                                                cg.coord_orth(),
                                                                cb.coord_orth() );

                        mg.set_glycosylation_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi) );
                        mg.add_torsions_for_plots(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), potential_o_roots[i].first.type().trim(), og1, sugar.type().trim(), c1);


                        list_of_glycans.push_back ( mg );
                        break;
                    }
                }
            }
        }
    }


    for ( int i = 0 ; i < potential_s_roots.size() ; i++ )  // create o-glycan roots with first sugar
    {
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_s_roots[i].first ) ;

        for ( int j = 0 ; j < linked.size() ; j++ )
        {
            const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()];

            if (( tmpmon.type().trim() == "NGA" ) || ( tmpmon.type().trim() == "A2G" ) || (tmpmon.type().trim() == "FUC" ) || (tmpmon.type().trim() == "RAM" ) || (tmpmon.type().trim() == "BGC" ) || (tmpmon.type().trim() == "BMA" ) || (tmpmon.type().trim() == "NAG" )
                || (tmpmon.type().trim() == "NDG" ))
            {
                if ( clipper::MSugar::search_database( tmpmon.type().c_str() ) )
                {
                    clipper::MSugar sugar( mmol, mmol[linked[j].second.polymer()].id().trim(), tmpmon, manb, debug_output );
                    if(sugar.ring_members().size() == 5 || sugar.ring_members().size() == 6)
                    {
                        list_of_sugars.push_back ( sugar );

                        if(debug_output)
                        {  
                            DBG << "Created the MSugar object" << std::endl;
                        }

                        if(debug_output)
                        {
                            DBG << "potential s roots is " << potential_s_roots[i].first.type() << std::endl;
                            DBG << "sugar is " << sugar.type() << std::endl;
                            DBG << "id is " << mmol[linked[j].second.polymer()].id().trim() << std::endl;
                        }
                        
                        clipper::String root_sugar_chain_id = mmol[linked[j].second.polymer()].id().trim().substr(0,1);
                        clipper::MGlycan mg (   potential_s_roots[i].second,
                                                potential_s_roots[i].first,
                                                list_of_sugars.back(),
                                                root_sugar_chain_id,
                                                debug_output,
                                                this->expression_system );
                        mg.set_kind_of_glycan ( "s-glycan" );
                        list_of_glycans.push_back ( mg );
                        break;
                    }
                }
            }
        }
    }

    for ( int i = 0 ; i < potential_c_roots.size() ; i++ )  // create c-glycan roots with first sugar
    {
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_c_roots[i].first ) ;

        for ( int j = 0 ; j < linked.size() ; j++ )
        {
            const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()];

            if (( tmpmon.type().trim() == "MAN") || (tmpmon.type().trim() == "BMA") || (tmpmon.type().trim() == "7D1") || (tmpmon.type().trim() == "M1P") || (tmpmon.type().trim() == "M6P") || (tmpmon.type().trim() == "MBF") || (tmpmon.type().trim() == "MMA") || (tmpmon.type().trim() == "OPM") || (tmpmon.type().trim() == "M6D"))
            {
                if ( clipper::MSugar::search_database( tmpmon.type().c_str() ) )
                {
                    clipper::MSugar sugar( mmol, mmol[linked[j].second.polymer()].id().trim(), tmpmon, manb, debug_output );
                    if(sugar.ring_members().size() == 5 || sugar.ring_members().size() == 6)
                    {
                        list_of_sugars.push_back ( sugar );

                        if(debug_output)
                        {
                            DBG << "Created the MSugar object" << std::endl;
                        }

                        if(debug_output)
                        {
                            DBG << "potential c roots is " << potential_s_roots[i].first.type() << std::endl;
                            DBG << "sugar is " << sugar.type() << std::endl;
                            DBG << "id is " << mmol[linked[j].second.polymer()].id().trim() << std::endl;
                        }

                        clipper::String root_sugar_chain_id = mmol[linked[j].second.polymer()].id().trim().substr(0,1);
                        clipper::MGlycan mg (   potential_c_roots[i].second,
                                                potential_c_roots[i].first,
                                                list_of_sugars.back(),
                                                root_sugar_chain_id,
                                                debug_output,
                                                this->expression_system );
                        mg.set_kind_of_glycan ( "c-glycan" );

                        // Check 1C4 conformation on ring - see John et al, Nature Chemical Biology 2021 (17):428–437
                        if ( sugar.conformation_name() == "1c4" ){
                        sugar.override_conformation_diag ( true );
                        }

                        clipper::MAtom o5 = sugar.ring_members()[0];              // O5
                        clipper::MAtom c1 = sugar.ring_members()[1];              // C1
                        clipper::MAtom cd1= sugar.anomeric_substituent();         // CD1 TRP
                        clipper::MAtom cg = potential_c_roots[i].first.find("CG", clipper::MM::ANY);      // CG
                        clipper::MAtom cb = potential_c_roots[i].first.find("CB", clipper::MM::ANY);      // CB
                        clipper::ftype phi, psi;

                        phi   = clipper::Coord_orth::torsion (  o5.coord_orth(),
                                                                c1.coord_orth(),
                                                                cd1.coord_orth(),
                                                                cg.coord_orth() );

                        psi   = clipper::Coord_orth::torsion (  c1.coord_orth(),
                                                                cd1.coord_orth(),
                                                                cg.coord_orth(),
                                                                cb.coord_orth() );

                        mg.set_glycosylation_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi) );
                        mg.add_torsions_for_plots(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), potential_c_roots[i].first.type().trim(), cd1, sugar.type().trim(), c1);

                        if ( linked[j].second.monomer()+3 < mmol[linked[j].second.polymer()].size() )
                        // Make sure that checks for consensus sequence do not occur outside the array, therefore causing segfaults.
                        {
                            // Consensus sequence Trp-X-X-Trp || Trp-Ser/Thr-X-Cys according to https://www.uniprot.org/help/carbohyd
                            bool firstConsensus = false;
                            bool secondConsensus = false;
                            if ( linked[j].second.monomer()-3 > 0 )
                                if ( mmol[linked[j].second.polymer()][linked[j].second.monomer()+3].type().trim() == "TRP" ||
                                    mmol[linked[j].second.polymer()][linked[j].second.monomer()-3].type().trim() == "TRP" )
                                        firstConsensus = true;
                            if ( mmol[linked[j].second.polymer()][linked[j].second.monomer()+1].type().trim() == "SER" ||
                                mmol[linked[j].second.polymer()][linked[j].second.monomer()+1].type().trim() == "THR" &&
                                mmol[linked[j].second.polymer()][linked[j].second.monomer()+3].type().trim() == "CYS" )
                                    secondConsensus = true;

                            if (!firstConsensus && !secondConsensus) mg.add_root_annotation ( " Warning: this glycosylation point does not follow the Trp-X-X-Trp or Trp-Ser/Thr-X-Cys consensus sequence. ");
                        }

                        list_of_glycans.push_back ( mg );
                        break;
                    }
                }
            }
        }
    }

    for ( int i = 0 ; i < potential_rootless_polysaccharides.size() ; i++ )  // create ligand glycan roots with first sugar
    {
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_rootless_polysaccharides[i].first ) ;

        if(linked.empty())
        {
            if ( clipper::MSugar::search_database(potential_rootless_polysaccharides[i].first.type().c_str()) )
            {
                clipper::MSugar rootSugar (mmol, potential_rootless_polysaccharides[i].second, potential_rootless_polysaccharides[i].first, manb, debug_output);
                list_of_sugars.push_back ( rootSugar );
                clipper::String root_sugar_chain_id = potential_rootless_polysaccharides[i].second.substr(0,1);
                
                if(debug_output)
                {
                    DBG << "Created the rootSugar object in linked.empty(): " << root_sugar_chain_id << "/" << rootSugar.type().trim() << "-" << rootSugar.id().trim() << std::endl;
                    DBG << "potential rootles roots is " << potential_rootless_polysaccharides[i].first.type() << std::endl;
                    DBG << "rootSugar is " << rootSugar.type() << std::endl;
                    DBG << "id is " << root_sugar_chain_id << std::endl;
                }
                
                clipper::MGlycan mg (   potential_rootless_polysaccharides[i].second,
                                        rootSugar,
                                        root_sugar_chain_id,
                                        debug_output,
                                        this->expression_system );
                mg.set_kind_of_glycan ( "ligand" );
                list_of_glycans.push_back ( mg );
            }
        }
        else
        {
            for ( int j = 0 ; j < linked.size() ; j++ )
            {
                const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()];


                // Might need to revert this change. To test with running Privateer again on a local pdb_mirror.
                if ( clipper::MSugar::search_database(potential_rootless_polysaccharides[i].first.type().c_str()) )
                {
                    clipper::MSugar rootSugar (mmol, potential_rootless_polysaccharides[i].second, potential_rootless_polysaccharides[i].first, manb, debug_output);
                    list_of_sugars.push_back ( rootSugar );
                    clipper::String root_sugar_chain_id = potential_rootless_polysaccharides[i].second.substr(0,1);
                    
                    if(debug_output)
                    {
                        DBG << "Created the rootSugar object where !linked.empty(): " << root_sugar_chain_id << "/" << rootSugar.type().trim() << "-" << rootSugar.id().trim() << std::endl;
                        DBG << "potential rootles roots is " << potential_rootless_polysaccharides[i].first.type() << std::endl;
                        DBG << "rootSugar is " << rootSugar.type() << std::endl;
                        DBG << "id is " << mmol[linked[j].second.polymer()].id().trim() << std::endl;
                    }
                    clipper::MGlycan mg (   potential_rootless_polysaccharides[i].second,
                                            rootSugar,
                                            root_sugar_chain_id,
                                            debug_output,
                                            this->expression_system );

                    mg.set_kind_of_glycan ( "ligand" );
                    list_of_glycans.push_back ( mg );
                    break;
                }
            }
        }
    }



    for ( int i = 0 ; i < list_of_glycans.size() ; i++ )
    {
        std::vector < clipper::MSugar >& sugar_list = list_of_glycans[i].get_sugars();
        clipper::MSugar first_sugar = clipper::MSugar ( sugar_list.front() );

        if(debug_output)
        {
            DBG << "Extending tree of glycan from root of " << list_of_glycans[i].get_root_for_filename() << "\twith detected sugar = " << first_sugar.type() << "-" << first_sugar.id() << "/" <<first_sugar.get_seqnum() << std::endl;
        }
        
        if(first_sugar.ring_members().size() == 5 || first_sugar.ring_members().size() == 6)
            extend_tree ( list_of_glycans[i] , first_sugar );

        list_of_glycans[i].set_annotations( this->expression_system );
    }
}

void MGlycology::pyinit ( const clipper::MiniMol& mmol, const clipper::MAtomNonBond& manb, bool debug_output, std::string expression_system )
{
    this->debug_output = debug_output;
    this->manb = &manb;
    this->mmol = &mmol;

    this->expression_system = expression_system;

    std::vector<std::pair<clipper::MMonomer, clipper::String>> potential_n_roots;
    std::vector<std::pair<clipper::MMonomer, clipper::String>> potential_o_roots;
    std::vector<std::pair<clipper::MMonomer, clipper::String>> potential_s_roots;
    std::vector<std::pair<clipper::MMonomer, clipper::String>> potential_c_roots;
    std::vector<std::pair<clipper::MMonomer, clipper::String>> potential_rootless_polysaccharides;
    std::vector<std::pair<clipper::MMonomer, clipper::String>> special_rootless_polysaccharides;
    

    for ( int pol = 0; pol < mmol.size() ; pol++ )
        for ( int mon = 0 ; mon < mmol[pol].size() ; mon++ )
        {
            std::vector<clipper::MAtom> KDOlikerecipe;
            for (int atom = 0; atom < mmol[pol][mon].size(); atom++)
            {
                clipper::MAtom tmpAtom = mmol[pol][mon][atom];
                if(get_altconf(tmpAtom) != ' ')
                {
                    const char altconf = get_altconf(tmpAtom);
                    
                    std::string altConfSymbol(1, altconf);
                    std::string O1_with_altconf = "O1 :" + altConfSymbol;
                    std::string C1_with_altconf = "C1 :" + altConfSymbol;
                    std::string OA_with_altconf = "O1A :" + altConfSymbol;
                    std::string OB_with_altconf = "O1B :" + altConfSymbol;

                    if (tmpAtom.id().trim() == O1_with_altconf)
                        potential_rootless_polysaccharides.push_back( std::make_pair(mmol[pol][mon], mmol[pol].id()) );
                    
                    if (tmpAtom.id().trim() == C1_with_altconf || tmpAtom.id().trim() == OA_with_altconf || tmpAtom.id().trim() == OB_with_altconf)
                        KDOlikerecipe.push_back(tmpAtom);
                }
                else
                    if (tmpAtom.id().trim() == "O1") 
                        potential_rootless_polysaccharides.push_back( std::make_pair(mmol[pol][mon], mmol[pol].id()) );
                    
                    if (tmpAtom.id().trim() == "C1" || tmpAtom.id().trim() == "O1A" || tmpAtom.id().trim() == "O1B")
                        KDOlikerecipe.push_back(tmpAtom);
            }
            if(KDOlikerecipe.size() == 3)
                special_rootless_polysaccharides.push_back( std::make_pair(mmol[pol][mon], mmol[pol].id()) );

            // Will need to keep this list up to date with the latest discoveries
            // To do: include check on other ligands, such as lipids (e.g. ceramide O-glycosylation)
            if ( mmol[pol][mon].type() == "ASN" ) potential_n_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // n-linked GlcNAc ?
            else if ( mmol[pol][mon].type() == "ARG" ) potential_n_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // Arginine rhamnosylation?
            else if ( mmol[pol][mon].type() == "THR" ) potential_o_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // o-linked stuff ?
            else if ( mmol[pol][mon].type() == "SER" ) potential_o_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) );
            else if ( mmol[pol][mon].type() == "LYS" ) potential_o_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) );
            else if ( mmol[pol][mon].type() == "TYR" ) potential_o_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) );
            else if ( mmol[pol][mon].type() == "CYS" ) potential_s_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // s-linked stuff ?
            else if ( mmol[pol][mon].type() == "TRP" ) potential_c_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // C-linked stuff for C/TRP-mannosylation
        }
    for(int i = 0; i < special_rootless_polysaccharides.size(); i++)
        potential_rootless_polysaccharides.push_back(special_rootless_polysaccharides[i]);
    for ( int i = 0 ; i < potential_n_roots.size() ; i++ )  // create n-glycan roots with first sugar
    {
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_n_roots[i].first ) ;

        for ( int j = 0 ; j < linked.size() ; j++ )
        {
            const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()];

            if ( tmpmon.type().trim() == "NAG" )
            {
                if ( clipper::MSugar::search_database( tmpmon.type().c_str() ) )
                {
                    clipper::MSugar sugar( mmol, mmol[linked[j].second.polymer()].id().trim(), tmpmon, manb, debug_output );
                    if(sugar.ring_members().size() == 5 || sugar.ring_members().size() == 6)
                    {
                        list_of_sugars.push_back ( sugar );

                        if(debug_output)
                        {
                            DBG << "Created the MSugar object" << std::endl;
                        }

                        if(debug_output)
                        {
                            DBG << "potential n roots is " << potential_n_roots[i].first.type() << std::endl;
                            DBG << "sugar is " << sugar.type() << std::endl;
                            DBG << "id is " << mmol[linked[j].second.polymer()].id().trim() << std::endl;
                        }

                        clipper::String root_sugar_chain_id = mmol[linked[j].second.polymer()].id().trim().substr(0,1);
                        clipper::MGlycan mg (   potential_n_roots[i].second,
                                                potential_n_roots[i].first,
                                                list_of_sugars.back(),
                                                root_sugar_chain_id,
                                                debug_output,
                                                this->expression_system );

                        if(debug_output)
                        {
                            DBG << "Exited the glycan constructor!" << std::endl;
                        }

                        mg.set_kind_of_glycan ( "n-glycan" );

                        if ( linked[j].second.monomer()+2 < mmol[linked[j].second.polymer()].size() )
                        {
                            if ( mmol[linked[j].second.polymer()][linked[j].second.monomer()+2].type().trim() != "THR" &&
                                mmol[linked[j].second.polymer()][linked[j].second.monomer()+2].type().trim() != "SER" )
                                // this is not a consensus ASN-GlcNAc glycosylation point
                                mg.add_root_annotation ( " Warning: this glycosylation point does not follow the Asn-X-Thr/Ser consensus sequence. ");
                        }

                        clipper::MAtom o5 = sugar.ring_members()[0];              // O5
                        clipper::MAtom c1 = sugar.ring_members()[1];              // C1
                        clipper::MAtom nd2= sugar.anomeric_substituent();         // ND2 usually
                        clipper::MAtom cg = potential_n_roots[i].first.find("CG", clipper::MM::ANY);      // CG
                        clipper::MAtom cb = potential_n_roots[i].first.find("CB", clipper::MM::ANY);      // CB
                        clipper::ftype phi, psi;

                        phi   = clipper::Coord_orth::torsion (  o5.coord_orth(),
                                                                c1.coord_orth(),
                                                                nd2.coord_orth(),
                                                                cg.coord_orth() );

                        psi   = clipper::Coord_orth::torsion (  c1.coord_orth(),
                                                                nd2.coord_orth(),
                                                                cg.coord_orth(),
                                                                cb.coord_orth() );

                        if ( psi < 0 )
                            psi = clipper::Util::twopi() + psi;

                        mg.set_glycosylation_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi) );
                        mg.add_torsions_for_plots(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), potential_n_roots[i].first.type().trim(), nd2, sugar.type().trim(), c1);

                        list_of_glycans.push_back ( mg );
                        break;
                    }
                }
            }
            else if ( tmpmon.type().trim() == "NDG" )
            {
                if ( clipper::MSugar::search_database( tmpmon.type().c_str() ) )
                {
                    clipper::MSugar sugar( mmol, mmol[linked[j].second.polymer()].id().trim(), tmpmon, manb, debug_output );
                    if(sugar.ring_members().size() == 5 || sugar.ring_members().size() == 6)
                    {
                        list_of_sugars.push_back ( sugar );
                        clipper::String root_sugar_chain_id = mmol[linked[j].second.polymer()].id().trim().substr(0,1);
                        clipper::MGlycan mg (   potential_n_roots[i].second,
                                                potential_n_roots[i].first,
                                                list_of_sugars.back(),
                                                root_sugar_chain_id,
                                                debug_output,
                                                this->expression_system );

                        mg.set_kind_of_glycan ( "n-glycan" );

                        clipper::MAtom o5 = sugar.ring_members()[0];              // O5
                        clipper::MAtom c1 = sugar.ring_members()[1];              // C1
                        clipper::MAtom nd2= sugar.anomeric_substituent();         // ND2 usually
                        clipper::MAtom cg = potential_n_roots[i].first.find("CG", clipper::MM::ANY);      // CG
                        clipper::MAtom cb = potential_n_roots[i].first.find("CB", clipper::MM::ANY);      // CB
                        clipper::ftype phi, psi;

                        phi   = clipper::Coord_orth::torsion (  o5.coord_orth(),
                                                                c1.coord_orth(),
                                                                nd2.coord_orth(),
                                                                cg.coord_orth() );

                        psi   = clipper::Coord_orth::torsion (  c1.coord_orth(),
                                                                nd2.coord_orth(),
                                                                cg.coord_orth(),
                                                                cb.coord_orth() );

                        mg.set_glycosylation_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi) );
                        mg.add_torsions_for_plots(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), potential_n_roots[i].first.type().trim(), nd2, sugar.type().trim(), c1);


                        if ( linked[j].second.monomer()+2 < mmol[linked[j].second.polymer()].size() )
                        {
                            if ( mmol[linked[j].second.polymer()][linked[j].second.monomer()+2].type().trim() != "THR" &&
                                mmol[linked[j].second.polymer()][linked[j].second.monomer()+2].type().trim() != "SER" )
                                // this is not a consensus ASN-GlcNAc glycosylation point
                                mg.add_root_annotation ( " Warning: this glycosylation point does not follow the Asn-X-Thr/Ser consensus sequence. ");
                        }

                        list_of_glycans.push_back ( mg );

                        break;
                    }
                }
            }
        }
    }


    for ( int i = 0 ; i < potential_o_roots.size() ; i++ )  // create o-glycan roots with first sugar
    {
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_o_roots[i].first ) ;

        for ( int j = 0 ; j < linked.size() ; j++ )
        {
            const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()];

            if (( tmpmon.type().trim() == "NGA" ) || ( tmpmon.type().trim() == "A2G" ) || (tmpmon.type().trim() == "FUC" ) || (tmpmon.type().trim() == "RAM" ) || (tmpmon.type().trim() == "BGC" ) || (tmpmon.type().trim() == "BMA" ) || (tmpmon.type().trim() == "NAG" )
                || (tmpmon.type().trim() == "NDG" ))
            {
                if ( clipper::MSugar::search_database( tmpmon.type().c_str() ) )
                {
                    clipper::MSugar sugar( mmol, mmol[linked[j].second.polymer()].id().trim(), tmpmon, manb, debug_output );
                    if(sugar.ring_members().size() == 5 || sugar.ring_members().size() == 6)
                    {
                        list_of_sugars.push_back ( sugar );
                        clipper::String root_sugar_chain_id = mmol[linked[j].second.polymer()].id().trim().substr(0,1);
                        clipper::MGlycan mg (   potential_o_roots[i].second,
                                                potential_o_roots[i].first,
                                                list_of_sugars.back(),
                                                root_sugar_chain_id,
                                                debug_output,
                                                this->expression_system );

                        mg.set_kind_of_glycan ( "o-glycan" );

                        clipper::MAtom o5 = sugar.ring_members()[0];              // O5
                        clipper::MAtom c1 = sugar.ring_members()[1];              // C1
                        clipper::MAtom og1= sugar.anomeric_substituent();         // OG/OG1 SER/THR
                        clipper::MAtom cg = potential_o_roots[i].first.find("CB", clipper::MM::ANY);      // CB
                        clipper::MAtom cb = potential_o_roots[i].first.find("CA", clipper::MM::ANY);      // CA
                        clipper::ftype phi, psi;

                        phi   = clipper::Coord_orth::torsion (  o5.coord_orth(),
                                                                c1.coord_orth(),
                                                                og1.coord_orth(),
                                                                cg.coord_orth() );

                        psi   = clipper::Coord_orth::torsion (  c1.coord_orth(),
                                                                og1.coord_orth(),
                                                                cg.coord_orth(),
                                                                cb.coord_orth() );

                        mg.set_glycosylation_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi) );
                        mg.add_torsions_for_plots(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), potential_o_roots[i].first.type().trim(), og1, sugar.type().trim(), c1);


                        list_of_glycans.push_back ( mg );
                        break;
                    }
                }
            }
        }
    }


    for ( int i = 0 ; i < potential_s_roots.size() ; i++ )  // create o-glycan roots with first sugar
    {
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_s_roots[i].first ) ;

        for ( int j = 0 ; j < linked.size() ; j++ )
        {
            const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()];

            if (( tmpmon.type().trim() == "NGA" ) || ( tmpmon.type().trim() == "A2G" ) || (tmpmon.type().trim() == "FUC" ) || (tmpmon.type().trim() == "RAM" ) || (tmpmon.type().trim() == "BGC" ) || (tmpmon.type().trim() == "BMA" ) || (tmpmon.type().trim() == "NAG" )
                || (tmpmon.type().trim() == "NDG" ))
            {
                if ( clipper::MSugar::search_database( tmpmon.type().c_str() ) )
                {
                    clipper::MSugar sugar( mmol, mmol[linked[j].second.polymer()].id().trim(), tmpmon, manb, debug_output );
                    if(sugar.ring_members().size() == 5 || sugar.ring_members().size() == 6)
                    {
                        list_of_sugars.push_back ( sugar );
                        clipper::String root_sugar_chain_id = mmol[linked[j].second.polymer()].id().trim().substr(0,1);
                        clipper::MGlycan mg (   potential_s_roots[i].second,
                                                potential_s_roots[i].first,
                                                list_of_sugars.back(),
                                                root_sugar_chain_id,
                                                debug_output,
                                                this->expression_system );
                        mg.set_kind_of_glycan ( "s-glycan" );
                        list_of_glycans.push_back ( mg );
                        break;
                    }
                }
            }
        }
    }

    for ( int i = 0 ; i < potential_c_roots.size() ; i++ )  // create c-glycan roots with first sugar
    {
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_c_roots[i].first ) ;

        for ( int j = 0 ; j < linked.size() ; j++ )
        {
            const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()];

            if (( tmpmon.type().trim() == "MAN") || (tmpmon.type().trim() == "BMA") || (tmpmon.type().trim() == "7D1") || (tmpmon.type().trim() == "M1P") || (tmpmon.type().trim() == "M6P") || (tmpmon.type().trim() == "MBF") || (tmpmon.type().trim() == "MMA") || (tmpmon.type().trim() == "OPM") || (tmpmon.type().trim() == "M6D"))
            {
                if ( clipper::MSugar::search_database( tmpmon.type().c_str() ) )
                {
                    clipper::MSugar sugar( mmol, mmol[linked[j].second.polymer()].id().trim(), tmpmon, manb, debug_output );
                    if(sugar.ring_members().size() == 5 || sugar.ring_members().size() == 6)
                    {
                        list_of_sugars.push_back ( sugar );
                        clipper::String root_sugar_chain_id = mmol[linked[j].second.polymer()].id().trim().substr(0,1);
                        clipper::MGlycan mg (   potential_c_roots[i].second,
                                                potential_c_roots[i].first,
                                                list_of_sugars.back(),
                                                root_sugar_chain_id,
                                                debug_output,
                                                this->expression_system );
                        mg.set_kind_of_glycan ( "c-glycan" );

                        clipper::MAtom o5 = sugar.ring_members()[0];              // O5
                        clipper::MAtom c1 = sugar.ring_members()[1];              // C1
                        clipper::MAtom cd1= sugar.anomeric_substituent();         // CD1 TRP
                        clipper::MAtom cg = potential_c_roots[i].first.find("CG", clipper::MM::ANY);      // CG
                        clipper::MAtom cb = potential_c_roots[i].first.find("CB", clipper::MM::ANY);      // CB
                        clipper::ftype phi, psi;

                        phi   = clipper::Coord_orth::torsion (  o5.coord_orth(),
                                                                c1.coord_orth(),
                                                                cd1.coord_orth(),
                                                                cg.coord_orth() );

                        psi   = clipper::Coord_orth::torsion (  c1.coord_orth(),
                                                                cd1.coord_orth(),
                                                                cg.coord_orth(),
                                                                cb.coord_orth() );

                        mg.set_glycosylation_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi) );
                        mg.add_torsions_for_plots(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), potential_c_roots[i].first.type().trim(), cd1, sugar.type().trim(), c1);

                        if ( linked[j].second.monomer()+3 < mmol[linked[j].second.polymer()].size() )
                        // Make sure that checks for consensus sequence do not occur outside the array, therefore causing segfaults. 
                        {
                            // Consensus sequence Trp-X-X-Trp || Trp-Ser/Thr-X-Cys according to https://www.uniprot.org/help/carbohyd
                            bool firstConsensus = false;
                            bool secondConsensus = false;
                            if ( linked[j].second.monomer()-3 > 0 ) 
                                if ( mmol[linked[j].second.polymer()][linked[j].second.monomer()+3].type().trim() == "TRP" ||
                                    mmol[linked[j].second.polymer()][linked[j].second.monomer()-3].type().trim() == "TRP" )
                                        firstConsensus = true;
                            if ( mmol[linked[j].second.polymer()][linked[j].second.monomer()+1].type().trim() == "SER" || 
                                mmol[linked[j].second.polymer()][linked[j].second.monomer()+1].type().trim() == "THR" &&
                                mmol[linked[j].second.polymer()][linked[j].second.monomer()+3].type().trim() == "CYS" )
                                    secondConsensus = true;
                                
                            if (!firstConsensus && !secondConsensus) mg.add_root_annotation ( " Warning: this glycosylation point does not follow the Trp-X-X-Trp or Trp-Ser/Thr-X-Cys consensus sequence. ");
                        }

                        list_of_glycans.push_back ( mg );
                        break;
                    }
                }
            }
        }
    }

    for ( int i = 0 ; i < potential_rootless_polysaccharides.size() ; i++ )  // create ligand glycan roots with first sugar
    {
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_rootless_polysaccharides[i].first ) ;

        if(linked.empty())
        {
            if ( clipper::MSugar::search_database(potential_rootless_polysaccharides[i].first.type().c_str()) )
            {
                clipper::MSugar rootSugar (mmol, potential_rootless_polysaccharides[i].second, potential_rootless_polysaccharides[i].first, manb, debug_output);
                list_of_sugars.push_back ( rootSugar );
                clipper::String root_sugar_chain_id = potential_rootless_polysaccharides[i].second.substr(0,1);
                
                if(debug_output)
                {
                    DBG << "Created the rootSugar object in linked.empty(): " << root_sugar_chain_id << "/" << rootSugar.type().trim() << "-" << rootSugar.id().trim() << std::endl;
                    DBG << "potential rootles roots is " << potential_rootless_polysaccharides[i].first.type() << std::endl;
                    DBG << "rootSugar is " << rootSugar.type() << std::endl;
                    DBG << "id is " << root_sugar_chain_id << std::endl;
                }

                clipper::MGlycan mg (   potential_rootless_polysaccharides[i].second,
                                        rootSugar,
                                        root_sugar_chain_id,
                                        debug_output,
                                        this->expression_system );
                mg.set_kind_of_glycan ( "ligand" );
                list_of_glycans.push_back ( mg );
            }
        }
        else
        {
            for ( int j = 0 ; j < linked.size() ; j++ )
            {
                const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()];


                // Might need to revert this change. To test with running Privateer again on a local pdb_mirror.
                if ( clipper::MSugar::search_database(potential_rootless_polysaccharides[i].first.type().c_str()) )
                {
                    clipper::MSugar rootSugar (mmol, potential_rootless_polysaccharides[i].second, potential_rootless_polysaccharides[i].first, manb, debug_output);
                    list_of_sugars.push_back ( rootSugar );
                    clipper::String root_sugar_chain_id = potential_rootless_polysaccharides[i].second.substr(0,1);

                    if(debug_output)
                    {
                        DBG << "Created the rootSugar object where !linked.empty(): " << root_sugar_chain_id << "/" << rootSugar.type().trim() << "-" << rootSugar.id().trim() << std::endl;
                        DBG << "potential rootles roots is " << potential_rootless_polysaccharides[i].first.type() << std::endl;
                        DBG << "rootSugar is " << rootSugar.type() << std::endl;
                        DBG << "id is " << root_sugar_chain_id << std::endl;
                    }

                    clipper::MGlycan mg (   potential_rootless_polysaccharides[i].second,
                                            rootSugar,
                                            root_sugar_chain_id,
                                            debug_output,
                                            this->expression_system );
                    mg.set_kind_of_glycan ( "ligand" );
                    list_of_glycans.push_back ( mg );
                    break;
                }
            }
        }
    }


    for ( int i = 0 ; i < list_of_glycans.size() ; i++ )
    {
        std::vector < clipper::MSugar >& sugar_list = list_of_glycans[i].get_sugars();
        clipper::MSugar first_sugar = clipper::MSugar ( sugar_list.front() );
       
        if(debug_output)
        {
            DBG << "Extending tree of glycan from root of " << list_of_glycans[i].get_root_for_filename() << "\twith detected sugar = " << first_sugar.type() << "-" << first_sugar.id() << "/" <<first_sugar.get_seqnum() << std::endl;
        }
        
        if(first_sugar.ring_members().size() == 5 || first_sugar.ring_members().size() == 6)
            extend_tree ( list_of_glycans[i] , first_sugar );
        
        list_of_glycans[i].set_annotations( this->expression_system );
    }
}

std::string MGlycology::write_external_restraints ( bool restrain_rings,
                                                    bool restrain_links,
                                                    float weight ) {

  std::string restraints = "# External restraints for glycan refinement with Refmac5\n";
  restraints += "# Produced by Privateer MKIV, Glycojones team, University of York, UK.\n";
  std::vector<clipper::MGlycan> glycan_list = this->get_list_of_glycans();

  for ( int i = 0; i < glycan_list.size(); i++ ) {
    if ( restrain_rings ) {
      restraints += "\n# Ring conformation restraints for " + glycan_list[i].get_type()
                 + " at " + glycan_list[i].get_root_description() + "\n";
      restraints += glycan_list[i].write_ring_ext_restraints ( weight );
    }
    if ( restrain_links ) {
      restraints += "\n# Glycosidic bond conformation restraints\n" ;
      restraints += glycan_list[i].write_link_ext_restraints ( weight );
    }
  }
  restraints += "\n\n################ EOF ################\n" ;
  return restraints;
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
}                                                       // We will return a blank character if there is no code present or if it is, but is blank

void MGlycology::extend_tree ( clipper::MGlycan& mg, clipper::MSugar& msug )
{
    if(debug_output)
    {
        DBG << "clipper::MGlycan mg.number_of_nodes() = " << mg.number_of_nodes() << "\tmsug.id() = " << msug.id() << std::endl;
    }
    std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > contacts = get_contacts ( msug );

    const clipper::MiniMol& tmpmol = *this->mmol;

    if(debug_output)
    {
        DBG << "contacts.size() = " << contacts.size() << std::endl;
        for (int i = 0 ; i < contacts.size() ; i++ )
        {
            DBG << "Detected contacts[" << i << "].first.id().trim() from msug = " << contacts[i].first.id().trim() << "\t\t\tcontacts[" << i << "].second....id().trim() = " << tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()].type().trim() << "-" << tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()].id() << "-" << tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()][contacts[i].second.atom()].id() << std::endl;
            DBG << "Distance of contacts[" << i << "].first and  = " << "contacts[" << i << "].second = " << clipper::Coord_orth::length(contacts[i].first.coord_orth(), tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()][contacts[i].second.atom()].coord_orth()) << std::endl;
        }
    }
    for (int i = 0 ; i < contacts.size() ; i++ )
    {
        if (clipper::data::found_in_database ( tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()].type() ))
        {

            const std::vector<clipper::MSugar> sugar_list = mg.get_sugars();
            if (std::find_if(sugar_list.begin(), sugar_list.end(),
                [&](const clipper::MSugar& element) { return element.id() == tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()].id(); }) == sugar_list.end()) // prevent wrong circular linkages in glycosylation
            {
                clipper::MSugar tmpsug = clipper::MSugar ( *this->mmol, tmpmol[contacts[i].second.polymer()].id().trim(), tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()], *this->manb, debug_output );
                if(tmpsug.ring_members().size() == 5 || tmpsug.ring_members().size() == 6)
                {
                    clipper::MAtom acceptorAtom = tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()][contacts[i].second.atom()];
                    if(get_altconf(contacts[i].first) != ' ' && get_altconf(acceptorAtom) != ' ')
                    {
                        if(altconf_compatible(get_altconf(contacts[i].first), get_altconf(acceptorAtom)))
                        {
                            if(debug_output)
                            {
                                DBG << "parse_order ( contacts[" << i << "].first.id()) = " << parse_order ( contacts[i].first.id() ) << std::endl;
                            }
                            mg.link_sugars ( parse_order ( contacts[i].first.id() ), msug, tmpsug, contacts[i].first, acceptorAtom);
                            extend_tree ( mg, tmpsug );
                        }
                    }
                    else
                    {
                        if(debug_output)
                        {
                            DBG << "parse_order ( contacts[" << i << "].first.id()) = " << parse_order ( contacts[i].first.id() ) << std::endl;
                        }
                        mg.link_sugars ( parse_order ( contacts[i].first.id() ), msug, tmpsug, contacts[i].first, acceptorAtom);
                        extend_tree ( mg, tmpsug );
                    }
                }
            }
        }
    }
}


const std::vector < std::pair< clipper::String, clipper::MMonomer > > MGlycology::get_overlapping_residues ( const clipper::MMonomer& mm )
{
    // mm can be aminoacid (typically ASN) or sugar: ND2, O2, O3, O4, O6

    std::vector < clipper::MAtom > candidates; // look for the possibilities
    std::vector < std::pair< clipper::String, clipper::MMonomer > > tmpresults;
    std::vector < std::pair< clipper::String, clipper::MMonomer > > finalOutput;

    if ( mm.type().trim() == "ASN" )
    {
        int id1 = mm.lookup ( "ND2", clipper::MM::ANY );
        int id2 = mm.lookup ( "OD1", clipper::MM::ANY );

        if ( id1 != -1 )
        {
            candidates.push_back ( mm[id1] );

            if ( id2 != -1 )
                candidates.push_back ( mm[id2] );
        }
        else if ( id2 != -1 )
            candidates.push_back ( mm[id2] );
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
    else if ( mm.type().trim() == "TRP" )
    {
        int id = mm.lookup ( "CD1", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );
        else return tmpresults; // empty result

    }
    else
    {
        int id = mm.lookup ( "O1", clipper::MM::ANY );
        
        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O2", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O3", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O4", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );
        
        id = mm.lookup ( "O5", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O6", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );
        
        id = mm.lookup ( "O7", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O8", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S1", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S2", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S3", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S4", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S5", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S6", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S7", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S8", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N1", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N2", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N3", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N4", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N5", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N6", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N7", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N8", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F1", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F2", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F3", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F4", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F5", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F6", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F7", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F8", clipper::MM::ANY );

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
            std::vector < clipper::MAtomIndexSymmetry > contacts = this->manb->atoms_near ( candidates[i].coord_orth(), 2.5 );
            for (int j = 0 ; j < contacts.size() ; j++ )
            {
                if ((tmpmol[contacts[j].polymer()][contacts[j].monomer()].id().trim() != mm.id().trim())
                && ( clipper::Coord_orth::length ( tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()].coord_orth(), candidates[i].coord_orth() ) < 0.25 ))  // Beware: will report contacts that are not physically in contact, but needed for visualisation
                {                            //         of crappy structures in MG
                    tmpresults.push_back(std::make_pair(tmpmol[contacts[j].polymer()].id().trim(), tmpmol[contacts[j].polymer()][contacts[j].monomer()]));
                }
            }        
        }

        for(int i = 0; i < tmpresults.size(); i++)
        {
            std::pair< clipper::String, clipper::MMonomer > currentResult = tmpresults[i];
            
            auto searchresult_finalOutput = std::find_if(std::begin(finalOutput), std::end(finalOutput), [&currentResult](const std::pair<clipper::String, clipper::MMonomer> &vectorElement) {
                return currentResult.first == vectorElement.first && currentResult.second.id().trim() == vectorElement.second.id().trim() && currentResult.second.type().trim() == vectorElement.second.type().trim() && currentResult.second.seqnum() == vectorElement.second.seqnum();
            });
            auto searchresult_tmpresults = std::find_if(std::begin(tmpresults) + (i+1), std::end(tmpresults), [&currentResult](const std::pair<clipper::String, clipper::MMonomer> &vectorElement) {
                return currentResult.first == vectorElement.first && currentResult.second.id().trim() == vectorElement.second.id().trim() && currentResult.second.type().trim() == vectorElement.second.type().trim() && currentResult.second.seqnum() == vectorElement.second.seqnum();
            });
            
            if(searchresult_finalOutput == std::end(finalOutput) && searchresult_tmpresults != std::end(tmpresults))
                finalOutput.push_back(currentResult);

        }
        return finalOutput;
    }
}


const std::vector < std::pair< clipper::MAtom, clipper::MAtomIndexSymmetry > > MGlycology::get_contacts ( const clipper::MMonomer& mm )
{
    // mm can be aminoacid (typically ASN) or sugar: ND2, O2, O3, O4, O6

    std::vector < clipper::MAtom > candidates; // look for the possibilities
    std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > tmpresults;

    if ( mm.type().trim() == "ASN" )
    {
        int id1 = mm.lookup ( "ND2", clipper::MM::ANY );
        int id2 = mm.lookup ( "OD1", clipper::MM::ANY );

        if ( id1 != -1 )
        {
            candidates.push_back ( mm[id1] );

            if ( id2 != -1 )
                candidates.push_back ( mm[id2] );
        }
        else if ( id2 != -1 )
            candidates.push_back ( mm[id2] );
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
    else if ( mm.type().trim() == "TRP" )
    {
        int id = mm.lookup ( "CD1", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );
        else return tmpresults; // empty result

    }
    else
    {
        int id = mm.lookup ( "O1", clipper::MM::ANY );
        
        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O2", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O3", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O4", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );
        
        id = mm.lookup ( "O5", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O6", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );
        
        id = mm.lookup ( "O7", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O8", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S1", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S2", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S3", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S4", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S5", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S6", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S7", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S8", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N1", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N2", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N3", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N4", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N5", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N6", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N7", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N8", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F1", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F2", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F3", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F4", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F5", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F6", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F7", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F8", clipper::MM::ANY );

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
            std::vector < clipper::MAtomIndexSymmetry > contacts = this->manb->atoms_near ( candidates[i].coord_orth(), 2.5 );

            for (int j = 0 ; j < contacts.size() ; j++ )
            {
                if ((tmpmol[contacts[j].polymer()][contacts[j].monomer()].id().trim() != mm.id().trim())
                && ( clipper::Coord_orth::length ( tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()].coord_orth(), candidates[i].coord_orth() ) <= 2.5 ))  // Beware: will report contacts that are not physically in contact, but needed for visualisation
                {                            //         of crappy structures in MG
                    if ( altconf_compatible(get_altconf(tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()]),
                                             get_altconf(tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()])))
                    {
                        clipper::MAtom tmpAtom = tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()];
                        
                        if(debug_output)
                            DBG << "tmpAtom.id().trim() = " << tmpAtom.id().trim() << std::endl;

                        if(get_altconf(tmpAtom) != ' ')
                        {
                            const char altconf = get_altconf(tmpAtom);
                            

                            std::string altConfSymbol(1, altconf);
                            std::string C1_with_altconf = "C1 :" + altConfSymbol;
                            std::string C2_with_altconf = "C2 :" + altConfSymbol;

                            if(debug_output)
                                DBG << "tmpAtom.id().trim(): " << tmpAtom.id().trim() << "\tC1_with_altconf:" << C1_with_altconf << std::endl;

                            if(tmpAtom.id().trim() == C1_with_altconf || tmpAtom.id().trim() == C2_with_altconf) 
                            {
                                std::pair < clipper::MAtom , clipper::MAtomIndexSymmetry > link_tmp;
                                link_tmp.first = candidates[i];
                                link_tmp.second = contacts[j];
                                tmpresults.push_back ( link_tmp );
                                if(debug_output)
                                {
                                    DBG << "Contact detected between atoms of altconf = " << get_altconf(tmpAtom) << std::endl;
                                    DBG << "link_tmp.first(MAtom).id() = " << candidates[i].id() << "\tlink_tmp.second(MAtomNonBond) aka tmpAtom = " << tmpAtom.id() << std::endl;
                                }
                            }
                        }
                        else 
                        {
                            if(tmpAtom.id().trim() == "C1" || tmpAtom.id().trim() == "C2") 
                            {
                                std::pair < clipper::MAtom , clipper::MAtomIndexSymmetry > link_tmp;
                                link_tmp.first = candidates[i];
                                link_tmp.second = contacts[j];
                                tmpresults.push_back ( link_tmp );
                                if(debug_output)
                                {
                                    DBG << "link_tmp.first(MAtom).id() = " << candidates[i].id() << "\tlink_tmp.second(MAtomNonBond) aka tmpAtom = " << tmpAtom.id() << std::endl;
                                }
                            }
                        }
                    }
                }
            }        
        }

        std::sort(tmpresults.begin(), tmpresults.end(), [&tmpmol](const std::pair<clipper::MAtom,clipper::MAtomIndexSymmetry> &left, const std::pair<clipper::MAtom,clipper::MAtomIndexSymmetry> &right) {
            if(tmpmol[left.second.polymer()][left.second.monomer()].type().trim() == tmpmol[right.second.polymer()][right.second.monomer()].type().trim() && tmpmol[left.second.polymer()][left.second.monomer()].id() == tmpmol[right.second.polymer()][right.second.monomer()].id())
            {
                clipper::ftype distanceLeft = clipper::Coord_orth::length(left.first.coord_orth(), tmpmol[left.second.polymer()][left.second.monomer()][left.second.atom()].coord_orth());
                clipper::ftype distanceRight = clipper::Coord_orth::length(right.first.coord_orth(), tmpmol[right.second.polymer()][right.second.monomer()][right.second.atom()].coord_orth());
                return distanceLeft < distanceRight;
            }
            else return false;
        });
        return tmpresults;
    }
}