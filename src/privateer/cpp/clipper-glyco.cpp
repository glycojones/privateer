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

    this->sugar_diagnostics = clipper::MSugar::Diagnostics(false);
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

	if (this->sugar_ring_elements.size() == 5 && clipper::data::sugar_database[sugar_index].ref_conformation != "Pln")
	{
		this->cremerPople_furanose(*this->sugar_parent_molecule, mm);
		this->sugar_conformation = conformationFuranose(this->sugar_cremer_pople_params[1]);

		if(debug_output)
        {
			DBG << "After checking the conformation..." << std::endl;
		}
	}
	else if (this->sugar_ring_elements.size() == 6 && clipper::data::sugar_database[sugar_index].ref_conformation != "Pln")
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

    if (!sugar_ring_elements.empty() && sugar_ring_elements[1].name().trim().find("C1") != std::string::npos)
    {
        this->sugar_type = "aldose";
        this->sugar_denomination = "aldo";
    }
    else if(!sugar_ring_elements.empty() && sugar_ring_elements[1].name().trim().find("C2") != std::string::npos)
    {
        this->sugar_type = "ketose";
        this->sugar_denomination = "keto";
    }
    else
    {
        this->sugar_type = "linear";
        this->sugar_denomination = "linear";
    }

	if (sugar_ring_elements.size() == 5)
        this->sugar_denomination = this->sugar_denomination + "furanose";
	else if(sugar_ring_elements.size() == 6)
        this->sugar_denomination = this->sugar_denomination + "pyranose";
    else
        this->sugar_denomination = this->sugar_denomination + "_unsupported";

	this->sugar_denomination = clipper::String( this->sugar_anomer + "-" + this->sugar_handedness + "-" + this->sugar_denomination );

	// sanity check:

	this->sugar_sane = false;

	if(debug_output)
    {
		DBG << "Just before examining the ring..." << std::endl;
	}


	if (!sugar_ring_elements.empty() && examine_ring() ) sugar_diag_ring = true; else sugar_diag_ring = false;

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
        {
	        sugar_diag_chirality = true;
            this->sugar_diagnostics.sugar_diag_chirality_diagnostics.sugar_found_db = this->sugar_found_db;
            this->sugar_diagnostics.sugar_diag_chirality_diagnostics.sugar_handedness = sugar_handedness;
            this->sugar_diagnostics.sugar_diag_chirality_diagnostics.database_sugar_handedness = clipper::data::sugar_database[sugar_index].handedness;
            this->sugar_diagnostics.sugar_diag_chirality_diagnostics.final_result = sugar_diag_chirality;
            this->sugar_diagnostics.sugar_diag_chirality_diagnostics.initialized = true;
        }
	    else
        {
            sugar_diag_chirality = false;
            this->sugar_diagnostics.sugar_diag_chirality_diagnostics.sugar_found_db = this->sugar_found_db;
            this->sugar_diagnostics.sugar_diag_chirality_diagnostics.sugar_handedness = sugar_handedness;
            this->sugar_diagnostics.sugar_diag_chirality_diagnostics.database_sugar_handedness = clipper::data::sugar_database[sugar_index].handedness;
            this->sugar_diagnostics.sugar_diag_chirality_diagnostics.final_result = sugar_diag_chirality;
            this->sugar_diagnostics.sugar_diag_chirality_diagnostics.initialized = true;
        }


	    if ( ( ( sugar_anomer == "alpha") && ( clipper::data::sugar_database[sugar_index].anomer != "B" ) )
            || ( ( sugar_anomer == "beta") && ( clipper::data::sugar_database[sugar_index].anomer != "A" ) ) )
        {
            sugar_diag_anomer = true;
            this->sugar_diagnostics.sugar_diag_anomer_diagnostics.sugar_found_db = this->sugar_found_db;
            this->sugar_diagnostics.sugar_diag_anomer_diagnostics.sugar_anomer = sugar_anomer;
            this->sugar_diagnostics.sugar_diag_anomer_diagnostics.database_sugar_anomer = clipper::data::sugar_database[sugar_index].anomer;
            this->sugar_diagnostics.sugar_diag_anomer_diagnostics.final_result = sugar_diag_anomer;
            this->sugar_diagnostics.sugar_diag_anomer_diagnostics.initialized = true;
        }
	    else
        {
            sugar_diag_anomer = false;
            this->sugar_diagnostics.sugar_diag_anomer_diagnostics.sugar_found_db = this->sugar_found_db;
            this->sugar_diagnostics.sugar_diag_anomer_diagnostics.sugar_anomer = sugar_anomer;
            this->sugar_diagnostics.sugar_diag_anomer_diagnostics.database_sugar_anomer = clipper::data::sugar_database[sugar_index].anomer;
            this->sugar_diagnostics.sugar_diag_anomer_diagnostics.final_result = sugar_diag_anomer;
            this->sugar_diagnostics.sugar_diag_anomer_diagnostics.initialized = true;
        }

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
            {
                sugar_diag_puckering = true;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_found_db = this->sugar_found_db;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_min = ref_puckering - 0.18;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_max = ref_puckering + 0.15;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.comparison_operator = "and";
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_diag_conformation = sugar_diag_conformation;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.reference_puckering = ref_puckering;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.calculated_puckering_amplitude = puckering_amplitude();
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.final_result = sugar_diag_puckering;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.initialized = true;
            }
            else
            {
                sugar_diag_puckering = false;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_found_db = this->sugar_found_db;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_min = ref_puckering - 0.18;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_max = ref_puckering + 0.15;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.comparison_operator = "and";
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_diag_conformation = sugar_diag_conformation;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.reference_puckering = ref_puckering;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.calculated_puckering_amplitude = puckering_amplitude();
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.final_result = sugar_diag_puckering;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.initialized = true;
            }

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
            {
                sugar_diag_puckering = false;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_found_db = this->sugar_found_db;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_min = 0.42;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_max = 0.9;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.comparison_operator = "or";
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_diag_conformation = sugar_diag_conformation;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.reference_puckering = 0;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.calculated_puckering_amplitude = puckering_amplitude();
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.final_result = sugar_diag_puckering;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.initialized = true;
            }
            else
            {
                sugar_diag_puckering = true;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_found_db = this->sugar_found_db;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_min = 0.42;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_max = 0.9;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.comparison_operator = "or";
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_diag_conformation = sugar_diag_conformation;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.reference_puckering = 0;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.calculated_puckering_amplitude = puckering_amplitude();
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.final_result = sugar_diag_puckering;
                this->sugar_diagnostics.sugar_diag_puckering_diagnostics.initialized = true;
            }

            sugar_diag_bonds_rmsd = sugar_diag_angles_rmsd = true;
        }

	    if ( sugar_diag_puckering && sugar_diag_anomer && sugar_diag_chirality && sugar_diag_ring )
            sugar_sane = true;

        sugar_diagnostics.update_diagnostic_status(true);
        sugar_diagnostics.update_sugar_sane_status(sugar_sane);
	}
	else
	{
	    sugar_diag_anomer=true;
	    sugar_diag_chirality=true;  // perform a generic test based on rough ideal values

        this->sugar_diagnostics.sugar_diag_anomer_diagnostics.sugar_found_db = this->sugar_found_db;
        this->sugar_diagnostics.sugar_diag_anomer_diagnostics.sugar_anomer = sugar_anomer;
        this->sugar_diagnostics.sugar_diag_anomer_diagnostics.database_sugar_anomer = "NotFound";
        this->sugar_diagnostics.sugar_diag_anomer_diagnostics.final_result = sugar_diag_anomer;
        this->sugar_diagnostics.sugar_diag_anomer_diagnostics.initialized = true;

        this->sugar_diagnostics.sugar_diag_chirality_diagnostics.sugar_found_db = this->sugar_found_db;
        this->sugar_diagnostics.sugar_diag_chirality_diagnostics.sugar_handedness = sugar_handedness;
        this->sugar_diagnostics.sugar_diag_chirality_diagnostics.database_sugar_handedness = "NotFound";
        this->sugar_diagnostics.sugar_diag_chirality_diagnostics.final_result = sugar_diag_chirality;
        this->sugar_diagnostics.sugar_diag_chirality_diagnostics.initialized = true;

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
        {
            sugar_diag_puckering = false;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_found_db = this->sugar_found_db;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_min = 0.42;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_max = 0.9;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.comparison_operator = "or";
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_diag_conformation = sugar_diag_conformation;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.reference_puckering = 0;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.calculated_puckering_amplitude = puckering_amplitude();
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.final_result = sugar_diag_puckering;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.initialized = true;
        }
        else
        {
            sugar_diag_puckering = true;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_found_db = this->sugar_found_db;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_min = 0.42;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_max = 0.9;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.comparison_operator = "or";
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_diag_conformation = sugar_diag_conformation;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.reference_puckering = 0;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.calculated_puckering_amplitude = puckering_amplitude();
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.final_result = sugar_diag_puckering;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.initialized = true;
        }

	    if ( sugar_diag_puckering && sugar_diag_anomer && sugar_diag_chirality && sugar_diag_ring )
            sugar_sane = true;

        sugar_diagnostics.update_diagnostic_status(true);
        sugar_diagnostics.update_sugar_sane_status(sugar_sane);
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

    this->sugar_diagnostics = clipper::MSugar::Diagnostics(false);
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

    if (this->sugar_ring_elements.size() == 5 && clipper::data::sugar_database[sugar_index].ref_conformation != "Pln")
    {
        this->cremerPople_furanose(*this->sugar_parent_molecule, mm);
        this->sugar_conformation = conformationFuranose(this->sugar_cremer_pople_params[1]);

        if(debug_output)
        {
            DBG << "After checking the conformation..." << std::endl;
        }
    }
    else if (this->sugar_ring_elements.size() == 6 && clipper::data::sugar_database[sugar_index].ref_conformation != "Pln")
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


    if (!sugar_ring_elements.empty() && sugar_ring_elements[1].name().trim().find("C1") != std::string::npos)
    {
        this->sugar_type = "aldose";
        this->sugar_denomination = "aldo";
    }
    else if(!sugar_ring_elements.empty() && sugar_ring_elements[1].name().trim().find("C2") != std::string::npos)
    {
        this->sugar_type = "ketose";
        this->sugar_denomination = "keto";
    }
    else
    {
        this->sugar_type = "linear";
        this->sugar_denomination = "linear";
    }

	if (sugar_ring_elements.size() == 5)
        this->sugar_denomination = this->sugar_denomination + "furanose";
	else if(sugar_ring_elements.size() == 6)
        this->sugar_denomination = this->sugar_denomination + "pyranose";
    else
        this->sugar_denomination = this->sugar_denomination + "_unsupported";

	this->sugar_denomination = clipper::String( this->sugar_anomer + "-" + this->sugar_handedness + "-" + this->sugar_denomination );

	// sanity check:

	this->sugar_sane = false;

	if(debug_output)
    {
		DBG << "Just before examining the ring..." << std::endl;
	}

	if (!sugar_ring_elements.empty() && examine_ring() ) sugar_diag_ring = true; else sugar_diag_ring = false;

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
    {
        sugar_diag_chirality = true;
        this->sugar_diagnostics.sugar_diag_chirality_diagnostics.sugar_found_db = this->sugar_found_db;
        this->sugar_diagnostics.sugar_diag_chirality_diagnostics.sugar_handedness = sugar_handedness;
        this->sugar_diagnostics.sugar_diag_chirality_diagnostics.database_sugar_handedness = validation_data.handedness;
        this->sugar_diagnostics.sugar_diag_chirality_diagnostics.final_result = sugar_diag_chirality;
        this->sugar_diagnostics.sugar_diag_chirality_diagnostics.initialized = true;
    }
    else
    {
        sugar_diag_chirality = false;
        this->sugar_diagnostics.sugar_diag_chirality_diagnostics.sugar_found_db = this->sugar_found_db;
        this->sugar_diagnostics.sugar_diag_chirality_diagnostics.sugar_handedness = sugar_handedness;
        this->sugar_diagnostics.sugar_diag_chirality_diagnostics.database_sugar_handedness = validation_data.handedness;
        this->sugar_diagnostics.sugar_diag_chirality_diagnostics.final_result = sugar_diag_chirality;
        this->sugar_diagnostics.sugar_diag_chirality_diagnostics.initialized = true;
    }


    if ( ( ( sugar_anomer == "alpha") && ( validation_data.anomer != "B" ) )
	    || ( ( sugar_anomer == "beta") && ( validation_data.anomer != "A" ) ) )
    {
        sugar_diag_anomer = true;
        this->sugar_diagnostics.sugar_diag_anomer_diagnostics.sugar_found_db = true;
        this->sugar_diagnostics.sugar_diag_anomer_diagnostics.sugar_anomer = sugar_anomer;
        this->sugar_diagnostics.sugar_diag_anomer_diagnostics.database_sugar_anomer = clipper::data::sugar_database[sugar_index].anomer;
        this->sugar_diagnostics.sugar_diag_anomer_diagnostics.final_result = sugar_diag_anomer;
        this->sugar_diagnostics.sugar_diag_anomer_diagnostics.initialized = true;
    }
    else
    {
        sugar_diag_anomer = false;
        this->sugar_diagnostics.sugar_diag_anomer_diagnostics.sugar_found_db = true;
        this->sugar_diagnostics.sugar_diag_anomer_diagnostics.sugar_anomer = sugar_anomer;
        this->sugar_diagnostics.sugar_diag_anomer_diagnostics.database_sugar_anomer = validation_data.anomer;
        this->sugar_diagnostics.sugar_diag_anomer_diagnostics.final_result = sugar_diag_anomer;
        this->sugar_diagnostics.sugar_diag_anomer_diagnostics.initialized = true;
    }

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
        {
            sugar_diag_puckering = true;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_found_db = true;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_min = ref_puckering - 0.18;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_max = ref_puckering + 0.15;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.comparison_operator = "and";
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_diag_conformation = sugar_diag_conformation;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.reference_puckering = ref_puckering;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.calculated_puckering_amplitude = puckering_amplitude();
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.final_result = sugar_diag_puckering;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.initialized = true;
        }
        else
        {
            sugar_diag_puckering = false;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_found_db = true;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_min = ref_puckering - 0.18;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_max = ref_puckering + 0.15;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.comparison_operator = "and";
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_diag_conformation = sugar_diag_conformation;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.reference_puckering = ref_puckering;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.calculated_puckering_amplitude = puckering_amplitude();
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.final_result = sugar_diag_puckering;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.initialized = true;
        }

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
        {
            sugar_diag_puckering = false;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_found_db = true;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_min = 0.42;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_max = 0.9;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.comparison_operator = "or";
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_diag_conformation = sugar_diag_conformation;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.reference_puckering = 0;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.calculated_puckering_amplitude = puckering_amplitude();
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.final_result = sugar_diag_puckering;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.initialized = true;

        }
        else
        {
            sugar_diag_puckering = true;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_found_db = true;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_min = 0.42;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.range_max = 0.9;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.comparison_operator = "or";
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.sugar_diag_conformation = sugar_diag_conformation;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.reference_puckering = 0;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.calculated_puckering_amplitude = puckering_amplitude();
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.final_result = sugar_diag_puckering;
            this->sugar_diagnostics.sugar_diag_puckering_diagnostics.initialized = true;
        }

        sugar_diag_bonds_rmsd = sugar_diag_angles_rmsd = true;
    }

    if ( sugar_diag_puckering && sugar_diag_anomer && sugar_diag_chirality && sugar_diag_ring )
        sugar_sane = true;

    sugar_diagnostics.update_diagnostic_status(true);
    sugar_diagnostics.update_sugar_sane_status(sugar_sane);

    if(debug_output)
    {
	    DBG << "Just after examining the ring, exiting the constructor, good job!" << std::endl;
    }

}

MSugar::MSugar(const clipper::MiniMol& ml, const clipper::MMonomer& mm, const clipper::MAtomNonBond& nb, char alt_conf )
{

    copy(mm,clipper::MM::COPY_MPC);	// import_data from MMonomer


    this->sugar_supported = true;
    this->sugar_parent_molecule = &ml;
    this->sugar_parent_molecule_nonbond = &nb; // store pointers
    this->sugar_index = db_not_checked;
    this->sugar_index = 9999; // default value for "not found in database".
    this->sugar_alternate_confcode = " "; // initially, we would like to suppose this
    this->sugar_context = "";

    #if DUMP
        std::cout << std::endl ;
        DBG << "looking for " << this->id() << " " << this->type().trim() << " on the database..." << std::endl;
        // alt_conf != ' ' ? DBG << "Alternate locator supplied: " << alt_conf << std::endl : true;
    #endif

    #if DUMP
        std::cout << "Size of (*this).size() " << (*this).size() << std::endl;
        for(int i = 0; i < (*this).size(); i++)
        {
            std::cout << "Atom ID: (*this)[" << i << "].id " << (*this)[i].id() << std::endl;
        }
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
        #if DUMP
            DBG << "found it! " << std::endl;
        #endif

        std::vector<clipper::String> buffer = clipper::data::sugar_database[sugar_index].ring_atoms.trim().split(" ");


        for (int i=0 ; i < buffer.size() ; i++)
        {
            int index_atom = 0;

            if ( alt_conf != ' ' ) // alternate conformations are present
            {
                alt_conf == 'A' ? sugar_alternate_confcode = " :A" : sugar_alternate_confcode = " :B";

                index_atom = this->lookup(buffer[i].trim()+sugar_alternate_confcode,clipper::MM::UNIQUE);

                #if DUMP
                    DBG << "index_atom in line 146" << index_atom << std::endl;
                #endif

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
            else
            {

                index_atom = this->lookup(buffer[i],clipper::MM::ANY);


                #if DUMP
                    DBG << "index_atom in line 165 = " << index_atom << " value of buffer[" << i << "] =" << buffer[i] << "test" << std::endl;
                #endif

                if (index_atom == -1)
                {
                    this->sugar_supported = false;
                    this->sugar_sane = false; // we don't support cyclic sugars with less or more than 5-6 ring atoms
                    this->sugar_denomination = "    unsupported    ";
                    this->sugar_anomer = "X";
                    this->sugar_handedness = "X";


                    return;
                }
            }
            #if DUMP
                DBG << "trying to push (*this)[index_atom] in line 181" << (*this)[index_atom].id() << std::endl;
            #endif
            sugar_ring_elements.push_back((*this)[index_atom]);
        }
    }
    else
    {
        this->sugar_ring_elements = this->ringMembers();
    }


	if (this->sugar_ring_elements.size() == 5)
	{
		this->cremerPople_furanose(*this->sugar_parent_molecule, mm);
		this->sugar_conformation = conformationFuranose(this->sugar_cremer_pople_params[1]);

		#if DUMP
			DBG << "After checking the conformation..." << std::endl;
		#endif
	}
	else if (this->sugar_ring_elements.size() == 6)
	{
		this->cremerPople_pyranose(*this->sugar_parent_molecule, mm);
		this->sugar_conformation = conformationPyranose(this->sugar_cremer_pople_params[1], this->sugar_cremer_pople_params[2]);

        #if DUMP
			DBG << "After checking the conformation..." << std::endl;
		#endif
	}
	else
	{
		this->sugar_supported = false;
		this->sugar_sane = false; // we don't support cyclic sugars with less or more than 5-6 ring atoms
		this->sugar_denomination = "    unsupported    ";
		this->sugar_anomer = "X";
		this->sugar_handedness = "X";
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

	#if DUMP
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

	#if DUMP
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
        clipper::ftype distance = clipper::Coord_orth::length( ring_atoms[5].coord_orth(), mmol.atom(neighbourhood[i]).coord_orth() );
        if(neighbourhood[i].symmetry() == 0)
        {
            if  (( mmol.atom( neighbourhood[i] ).element().trim() != "H" ) && (mmol.atom(neighbourhood[i]).name().trim() != ring_atoms[5].name().trim()))
            // the target substituent could be anything apart from H, in-ring oxygen or in-ring carbon
            {
                std::vector<clipper::MSugar::Diagnostics::sugar_diag_ring_diagnostics_atom_pair> dummy_diagnostic;
                if ( bonded ( mmol.atom(neighbourhood[i]), ring_atoms[5], dummy_diagnostic) ) // check link
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
            clipper::ftype distance = clipper::Coord_orth::length( ring_atoms[4].coord_orth(), mmol.atom(neighbourhood[i]).coord_orth() );
            if((neighbourhood[i].symmetry() == 0))
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
                    std::vector<clipper::MSugar::Diagnostics::sugar_diag_ring_diagnostics_atom_pair> dummy_diagnostic;
                    if ( bonded(mmol.atom(neighbourhood[i]), ring_atoms[4], dummy_diagnostic) ) // check link
                        if ( !is_part_of_ring ( mmol.atom(neighbourhood[i]), ring_atoms ) ) // eliminate ring neighbours
                        {

                            if ( mmol.atom(neighbourhood[i]).element().trim() != "C" )
                                nz5_substituent = mmol.atom ( neighbourhood[i] ); // don't grab a carbon as substituent unless there's really no other option
                            else if ( nz5_substituent.name().trim() == "XXX" )
                                nz5_substituent = mmol.atom ( neighbourhood[i] );
                        }
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
                clipper::ftype distance = clipper::Coord_orth::length( ring_atoms[1].coord_orth(), mmol.atom(neighbourhood[i]).coord_orth() );
                if(neighbourhood[i].symmetry() == 0)
                {
                    if(debug_output)
                        DBG << "Anomer computation, current neighbour atom: " << mmol.atom(neighbourhood[i]).id().trim() << std::endl;
                    if (( mmol.atom(neighbourhood[i]).element().trim() != "H" )
                    && ( mmol.atom(neighbourhood[i]).name().trim() != anomeric_carbon.name().trim()) )
                    {
                        if ( bonded(neighbourhood[i], ring_atoms[1] ) )
                        {
                            if(debug_output)
                                DBG << "Anomer computation, bonded = true " << mmol.atom(neighbourhood[i]).id().trim() << std::endl;
                            if ( !is_part_of_ring(mmol.atom(neighbourhood[i]), ring_atoms))
                            {
                                if(debug_output)
                                    DBG << "Anomer computation, not part of the ring " << mmol.atom(neighbourhood[i]).id().trim() << std::endl;
                                if ( altconf_compatible(get_altconf ( mmol.atom(neighbourhood[i]) ), get_altconf ( ring_atoms[1] ) ))
                                {
                                    if(debug_output)
                                        DBG << "Anomer computation, compatible altconf " << mmol.atom(neighbourhood[i]).id().trim() << std::endl;
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
                    clipper::ftype distance = clipper::Coord_orth::length( configurational_carbon.coord_orth(), mmol.atom(neighbourhood[i2]).coord_orth() );
                    if(neighbourhood[i2].symmetry() == 0)
                    {
                        if ( !is_part_of_ring(mmol.atom(neighbourhood[i2]),ring_atoms)
                            && (mmol.atom(neighbourhood[i2]).element().trim() != "H" )
                            && altconf_compatible(get_altconf(mmol.atom(neighbourhood[i2])), get_altconf(anomeric_carbon)))
                //			&& (( get_altconf(mmol.atom(neighbourhood[i2])) == ' ' )
                //			|| ( get_altconf(mmol.atom(neighbourhood[i2])) == 'A' ) ) ) // the target substituent could be anything apart from H, in-ring oxygen or in-ring carbon
                        {
                            std::vector<clipper::MSugar::Diagnostics::sugar_diag_ring_diagnostics_atom_pair> dummy_diagnostic;
                            if ( bonded (mmol.atom(neighbourhood[i2]), configurational_carbon, dummy_diagnostic) ) // check link
                                if (( mmol.atom(neighbourhood[i2]).name().trim() != ring_atoms[i-1].name().trim())
                                    && ( mmol.atom(neighbourhood[i2]).name().trim() != ring_atoms[0].name().trim() )
                                    && ( mmol.atom(neighbourhood[i2]).name().trim() != configurational_carbon.name().trim() ))
                                        // eliminate ring neighbours
                                        configurational_substituent = mmol.atom(neighbourhood[i2]);
                        }
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
            clipper::ftype distance = clipper::Coord_orth::length( configurational_carbon.coord_orth(), mmol.atom(neighbourhood[i2]).coord_orth() );
            if(neighbourhood[i2].symmetry() == 0)
            {
                if ( !is_part_of_ring(mmol.atom(neighbourhood[i2]),ring_atoms)
                    && (mmol.atom(neighbourhood[i2]).element().trim() != "H" )
                    && altconf_compatible(get_altconf(mmol.atom(neighbourhood[i2])), get_altconf(configurational_carbon)))
    //				&& (( get_altconf(mmol.atom(neighbourhood[i2])) == ' ' )
    //				|| ( get_altconf(mmol.atom(neighbourhood[i2])) == 'A' ) ) ) // the target substituent could be anything apart from H, in-ring oxygen or in-ring carbon
                {   // set next carbon and configurational (subs, carbon)
                    std::vector<clipper::MSugar::Diagnostics::sugar_diag_ring_diagnostics_atom_pair> dummy_diagnostic;
                    if ( bonded (mmol.atom(neighbourhood[i2]), configurational_carbon, dummy_diagnostic) ) // check link
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

    if(debug_output)
        std::cout << ": called from get_stereochemistry, distance = " << distance << "\tma_one.id() " << sugar_parent_molecule->atom(ma_one).id() << "\tma_two.id() " << ma_two.id().trim() << std::endl;

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
            if ((distance >  1.16 ) && ( distance < 1.61 ))
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

    for (i = 1 ; i < sugar_ring_elements.size() - 1; i++ )  // calculate bond angles (all), lengths (all) and torsions (minus last one)
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

    for (i = 0 ; i < sugar_ring_elements.size() - 1; i++)
    {
        if (!bonded(sugar_ring_elements[i], sugar_ring_elements[i+1], this->sugar_diagnostics.sugar_diag_ring_diagnostics.ring_atom_diagnostic))
        {
            if(debug_output)
            {
                DBG << std::endl << "Returning false - sugar_ring_elements[" << i << "].id().trim() "  << sugar_ring_elements[i].id().trim() << "\t\t\tsugar_ring_elements[" << i+1 << "].id().trim() " << sugar_ring_elements[i+1].id().trim() << std::endl;
            }
            this->sugar_diagnostics.sugar_diag_ring_diagnostics.initialized = true;
            this->sugar_diagnostics.sugar_diag_ring_diagnostics.final_result = false;
            return false;
        }
    }

    // comment for future person editing this code: i++ gets incremented after the loop does a cycle, so loop finishes when i = 4, but after that cycle i gets iterated to 5, thus the bottom segment makes sense.
    if (!bonded(sugar_ring_elements[i], sugar_ring_elements[0], this->sugar_diagnostics.sugar_diag_ring_diagnostics.ring_atom_diagnostic))
    {
        if(debug_output)
        {
            DBG << std::endl << "Returning false - sugar_ring_elements[" << i << "].id().trim() "  << sugar_ring_elements[i].id().trim() << "\t\t\tsugar_ring_elements[0_" << 0 << "].id().trim() " << sugar_ring_elements[0].id().trim() << std::endl;
        }
        this->sugar_diagnostics.sugar_diag_ring_diagnostics.initialized = true;
        this->sugar_diagnostics.sugar_diag_ring_diagnostics.final_result = false;
        return false;
    }
    else
    {
        this->sugar_diagnostics.sugar_diag_ring_diagnostics.initialized = true;
        this->sugar_diagnostics.sugar_diag_ring_diagnostics.final_result = true;
        return true;
    }

}

/*! Internal function for checking whether two atoms are bonded or not
 * 	\param ma_one A clipper::MAtom object
 * 	\param ma_two A clipper::MAtom object
 * 	\return A boolean value with the obvious answer
 */

bool MSugar::bonded(const clipper::MAtom& ma_one, const clipper::MAtom& ma_two, std::vector<clipper::MSugar::Diagnostics::sugar_diag_ring_diagnostics_atom_pair>& ring_atom_diagnostic) const
{
    clipper::ftype distance = clipper::Coord_orth::length( ma_one.coord_orth(), ma_two.coord_orth() );
    if(debug_output)
    {
        DBG << std::endl << "Received - ma_one.id().trim() " << ma_one.id().trim() << "\t\t\tma_two.id().trim() " << ma_two.id().trim() << "\tdistance = " << distance << std::endl;
    }
    clipper::MSugar::Diagnostics::sugar_diag_ring_diagnostics_atom_pair current_pair;
    current_pair.measured_distance = distance;
    if ( ma_one.element().trim() == "C" )
    {
        if ( ma_two.element().trim() == "C" )
        {
            if ((distance >  1.18 ) && ( distance < 1.62 ))
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 1.18;
                current_pair.distance_max = 1.62;
                current_pair.atom_result = true;
                ring_atom_diagnostic.push_back(current_pair);
                return true; // C-C or C=C
            }
            else
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 1.18;
                current_pair.distance_max = 1.62;
                current_pair.atom_result = false;
                ring_atom_diagnostic.push_back(current_pair);
                return false;
            }
        }
        else if ( ma_two.element().trim() == "N" )
        {
            if ((distance >  1.24 ) && ( distance < 1.62 ))
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 1.24;
                current_pair.distance_max = 1.62;
                current_pair.atom_result = true;
                ring_atom_diagnostic.push_back(current_pair);
                return true; // C-N or C=N
            }
            else
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 1.24;
                current_pair.distance_max = 1.62;
                current_pair.atom_result = false;
                ring_atom_diagnostic.push_back(current_pair);
                return false;
            }
        }
        else if ( ma_two.element().trim() == "O" )
        {
            if ((distance > 1.16 ) && ( distance < 1.60 ))
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 1.16;
                current_pair.distance_max = 1.60;
                current_pair.atom_result = true;
                ring_atom_diagnostic.push_back(current_pair);
                return true; // C-O or C=O
            }
            else
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 1.16;
                current_pair.distance_max = 1.60;
                current_pair.atom_result = false;
                ring_atom_diagnostic.push_back(current_pair);
                return false;
            }
        }
        else if ( ma_two.element().trim() == "S" )
        {
            if (( distance > 1.2 ) && ( distance < 2.0 ))
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 1.2;
                current_pair.distance_max = 2.0;
                current_pair.atom_result = true;
                ring_atom_diagnostic.push_back(current_pair);
                return true; // C-S or C=S
            }
            else
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 1.2;
                current_pair.distance_max = 2.0;
                current_pair.atom_result = false;
                ring_atom_diagnostic.push_back(current_pair);
                return false;
            }
        }
        else if ( ma_two.element().trim() == "H" )
        {
            if ((distance >  0.96 ) && ( distance < 1.14 ))
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 0.96;
                current_pair.distance_max = 1.14;
                current_pair.atom_result = true;
                ring_atom_diagnostic.push_back(current_pair);
                return true; // C-H
            }
            else
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 0.96;
                current_pair.distance_max = 1.14;
                current_pair.atom_result = false;
                ring_atom_diagnostic.push_back(current_pair);
                return false;
            }
        }
    }
    else if ( ma_one.element().trim() == "N" )
    {
        if ( ma_two.element().trim() == "C" )
        {
            if ((distance >  1.24 ) && ( distance < 1.62 ))
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 1.24;
                current_pair.distance_max = 1.62;
                current_pair.atom_result = true;
                ring_atom_diagnostic.push_back(current_pair);
                return true; // N-C or N=C
            }
            else
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 1.24;
                current_pair.distance_max = 1.62;
                current_pair.atom_result = false;
                ring_atom_diagnostic.push_back(current_pair);
                return false;
            }
        }
        else if ( ma_two.element().trim() == "H" )
        {
            if ((distance >  0.90 ) && ( distance < 1.10 ))
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 0.90;
                current_pair.distance_max = 1.10;
                current_pair.atom_result = true;
                ring_atom_diagnostic.push_back(current_pair);
                return true; // N-H
            }
            else
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 0.90;
                current_pair.distance_max = 1.10;
                current_pair.atom_result = false;
                ring_atom_diagnostic.push_back(current_pair);
                return false;
            }
        }
    }
    else if ( ma_one.element().trim() == "O" )
    {
        if ( ma_two.element().trim() == "C" )
        {
            if ((distance >  1.16 ) && ( distance < 1.60 ))
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 1.16;
                current_pair.distance_max = 1.60;
                current_pair.atom_result = true;
                ring_atom_diagnostic.push_back(current_pair);
                return true; // O-C or O=C
            }
            else
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 1.16;
                current_pair.distance_max = 1.60;
                current_pair.atom_result = false;
                ring_atom_diagnostic.push_back(current_pair);
                return false;
            }
        }
        else if ( ma_two.element().trim() == "H" )
        {
            if ((distance >  0.88 ) && ( distance < 1.04 ))
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 0.88;
                current_pair.distance_max = 1.04;
                current_pair.atom_result = true;
                ring_atom_diagnostic.push_back(current_pair);
                return true; // O-H
            }
            else
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 0.88;
                current_pair.distance_max = 1.04;
                current_pair.atom_result = false;
                ring_atom_diagnostic.push_back(current_pair);
                return false;
            }
        }
    }
    else if ( ma_one.element().trim() == "S" )
    {
        if ( ma_two.element().trim() == "C" )
        {
            if ((distance >  1.26 ) && ( distance < 2.00 ))
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 1.26;
                current_pair.distance_max = 2.00;
                current_pair.atom_result = true;
                ring_atom_diagnostic.push_back(current_pair);
                return true; // O-C or O=C
            }
            else
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 1.26;
                current_pair.distance_max = 2.00;
                current_pair.atom_result = false;
                ring_atom_diagnostic.push_back(current_pair);
                return false;
            }
        }
        else if ( ma_two.element().trim() == "H" )
        {
            if ((distance >  0.78 ) && ( distance < 1.24 ))
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 0.78;
                current_pair.distance_max = 1.24;
                current_pair.atom_result = true;
                ring_atom_diagnostic.push_back(current_pair);
                return true; // O-H
            }
            else
            {
                current_pair.ma_one_atom = ma_one.id().trim();
                current_pair.ma_one_element = ma_one.element().trim();
                current_pair.ma_two_atom = ma_two.id().trim();
                current_pair.ma_two_element = ma_two.element().trim();
                current_pair.distance_min = 0.78;
                current_pair.distance_max = 1.24;
                current_pair.atom_result = false;
                ring_atom_diagnostic.push_back(current_pair);
                return false;
            }
        }
    }
    else if (( distance > 1.2) && (distance < 2.0))
    {
        current_pair.ma_one_atom = ma_one.id().trim();
        current_pair.ma_one_element = ma_one.element().trim();
        current_pair.ma_two_atom = ma_two.id().trim();
        current_pair.ma_two_element = ma_two.element().trim();
        current_pair.distance_min = 1.2;
        current_pair.distance_max = 2.0;
        current_pair.atom_result = true;
        ring_atom_diagnostic.push_back(current_pair);
        return true; // unknown bond
    }
    else
    {
        current_pair.ma_one_atom = ma_one.id().trim();
        current_pair.ma_one_element = ma_one.element().trim();
        current_pair.ma_two_atom = ma_two.id().trim();
        current_pair.ma_two_element = ma_two.element().trim();
        current_pair.distance_min = 1.2;
        current_pair.distance_max = 2.0;
        current_pair.atom_result = false;
        ring_atom_diagnostic.push_back(current_pair);
        return false;
    }

    current_pair.ma_one_atom = ma_one.id().trim();
    current_pair.ma_one_element = ma_one.element().trim();
    current_pair.ma_two_atom = ma_two.id().trim();
    current_pair.ma_two_element = ma_two.element().trim();
    current_pair.distance_min = 0;
    current_pair.distance_max = 0;
    current_pair.atom_result = false;
    ring_atom_diagnostic.push_back(current_pair);
    return false; // in case we haven't found any match
}


/*std::vector < std::pair< clipper::MAtomIndexSymmetry, float > > MSugar::get_stacked_residues ( std::string algorithm,
                                                                                                        float distance,
                                                                                                        float theta,
                                                                                                        float phi ) const
{
    std::vector<clipper::MAtom> ch_hydrogens;
    clipper::MAtom ma;
    clipper::Coord_orth centre_apolar;
    std::vector<std::pair<clipper::MAtomIndexSymmetry, float > > results;

    std::cout << "WE'RE ACTUALLY HERE" << std::endl;

    if ( this->type_of_sugar() == "beta-D-aldopyranose" ) {
      ch_hydrogens.push_back ( this->find("H1",clipper::MM::ANY) );
      ch_hydrogens.push_back ( this->find("H3",clipper::MM::ANY) );
      ch_hydrogens.push_back ( this->find("H5",clipper::MM::ANY) );
    }
    else if ( this->type_of_sugar() == "alpha-D-aldopyranose" ) {
      ch_hydrogens.push_back ( this->find("H3",clipper::MM::ANY) );
      ch_hydrogens.push_back ( this->find("H5",clipper::MM::ANY) );
    }
    else if ( this->type_of_sugar() == "beta-L-aldopyranose" ) {
      ch_hydrogens.push_back ( this->find("H1",clipper::MM::ANY) );
      ch_hydrogens.push_back ( this->find("H3",clipper::MM::ANY) );
      ch_hydrogens.push_back ( this->find("H5",clipper::MM::ANY) );
    }
    else if ( this->type_of_sugar() == "alpha-L-aldopyranose" ) {
      ch_hydrogens.push_back ( this->find("H3",clipper::MM::ANY) );
      ch_hydrogens.push_back ( this->find("H5",clipper::MM::ANY) );
    }
    else if ( this->type_of_sugar() == "beta-L-ketopyranose" ) {
      ch_hydrogens.push_back ( this->find("H4",clipper::MM::ANY) );
      ch_hydrogens.push_back ( this->find("H6",clipper::MM::ANY) );
    }
    else if ( this->type_of_sugar() == "alpha-L-ketopyranose" ) {
      ch_hydrogens.push_back ( this->find("H4",clipper::MM::ANY) );
      ch_hydrogens.push_back ( this->find("H6",clipper::MM::ANY) );
    }
    else if ( this->type() == "XYP" ) {
      ch_hydrogens.push_back ( this->find("H1B",clipper::MM::ANY) );
      ch_hydrogens.push_back ( this->find("H3B",clipper::MM::ANY) );
      ch_hydrogens.push_back ( this->find("H5B2",clipper::MM::ANY) );
    }
    else if ( this->type() == "XYS" ) {
      ch_hydrogens.push_back ( this->find("H3",clipper::MM::ANY) );
      ch_hydrogens.push_back ( this->find("H51",clipper::MM::ANY) );
    }
    else // monosaccharide is unsupported, return empty results
      return results;

    const std::vector<clipper::MAtomIndexSymmetry> neighbourhood = this->sugar_parent_molecule_nonbond->atoms_near(centre_apolar, 5.0);

    const clipper::MiniMol& mmol  = *sugar_parent_molecule;

	for ( int k = 0 ; k < neighbourhood.size() ; k++ )
	{
        const clipper::MMonomer& mmon = mmol[neighbourhood[k].polymer()][neighbourhood[k].monomer()];

        if (( mmon.type() != "TRP" ) && // might be worth extending to cover GLU, ASP, GLN, ASN
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

    if ( algorithm == "hudson" ) {
  		if ( distance < 4.5 ) // check distance, "compliant" with Hudson et al., JACS 2015
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
    else {
      if ( distance < 4.3 ) { // compliant with modified Brandl-Weiss algorithm
        // Implement rest of algo here
      }
    }
	}
    return results;
}
*/


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
bool MGlycan::link_sugars ( int link, clipper::MSugar& first_sugar, clipper::MSugar& next_sugar, clipper::MAtom& donorAtom, clipper::MAtom& acceptorAtom, bool noncircular, privateer::json::GlobalTorsionZScore& torsions_zscore_database)
{
    int index = 0;
    bool found = false;

    if(debug_output)
    {
        DBG << "Linking " << first_sugar.chain_id().trim() << "/" <<  first_sugar.id().trim() << first_sugar.type() << " with " << next_sugar.chain_id().trim() << "/" <<  next_sugar.id().trim() << next_sugar.type() << std::endl;
    }

    for ( int i = 0 ; i < node_list.size() ; i++ )
        if (    strcmp( node_list[i].get_sugar().id().c_str(), first_sugar.id().c_str()) == 0 &&
                strcmp( node_list[i].get_sugar().type().c_str(), first_sugar.type().c_str()) == 0 &&
                strcmp( node_list[i].get_sugar().chain_id().c_str(), first_sugar.chain_id().c_str()) == 0 &&
                node_list[i].get_sugar().seqnum() == first_sugar.seqnum())
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

    // Linkage new_connection ( link, next_sugar.anomer(), node_list.size()-1, noncircular );
    // if(noncircular)
    // {
    //     sugars.push_back ( next_sugar );
    //     node_list[index].add_connection ( new_connection ); // add the new connection to the previous node
    // }
    // else if(!noncircular)
    // {
    //     node_list.front().add_circular_connection(new_connection );
    // }

    // bool MGlycan::link_sugars ( int link, clipper::MSugar& first_sugar, clipper::MSugar& next_sugar, clipper::MAtom& donorAtom, clipper::MAtom& acceptorAtom, bool noncircular )
    // mg.link_sugars ( parse_order ( contacts[i].first, msug ), msug, tmpsug, contacts[i].first, acceptorAtom, true);
    // mg.link_sugars ( parse_order ( contacts[i].first, msug ), tmpsug, msug, acceptorAtom, contacts[i].first, false);

    // This is where I gave up when it came to cyclic glycans. Code behaviour is identical to the state before this whole quest, despite numerous line changes.
    // clipper::MSugar lastInsertedSugar = node_list.back().get_sugar();
    // if(debug_output)
    //     DBG << "lastInsertedSugar.id() = " << lastInsertedSugar.id().trim() << "\tnext_sugar.id() = " << next_sugar.id().trim() << "\tfirst_sugar.id() " << first_sugar.id().trim() << "\tnoncircual = " << std::boolalpha << noncircular << std::endl;
    // if( lastInsertedSugar.chain_id().trim() == next_sugar.chain_id().trim() &&
    //     lastInsertedSugar.id().trim() == next_sugar.id().trim() &&
    //     lastInsertedSugar.type().trim() == next_sugar.type().trim() &&
    //     lastInsertedSugar.seqnum() == next_sugar.seqnum() )
    // {
    //     if(debug_output)
    //     {
    //         DBG << "The linkage between the two sugars has already been added in previous iteration." << std::endl;
    //     }
    //     return true;
    // }

    if(noncircular)
    {
        Node new_node( next_sugar );      // create a new node with the next sugar
        node_list.push_back ( new_node ); // add the new sugar to the node list
    }

    if(debug_output)
    {
        DBG << "Creating new connection.\nlink: " << link << "\tnext_sugar.anomer(): " << next_sugar.anomer() << "\tnode_list.size()-1: " << node_list.size()-1 << std::endl;
        DBG << "Anomeric substituent: " << next_sugar.anomeric_substituent().id().trim() << " located at: " << next_sugar.anomeric_substituent().coord_orth().format() << std::endl;
    }

    Linkage new_connection ( link, next_sugar.anomer(), node_list.size()-1, noncircular );

    clipper::ftype  omega_nine  = 0.0,
                    omega_eight = 0.0,
                    omega_seven = 0.0,
                    omega_six   = 0.0,
                    omega       = 0.0,
                    psi         = 0.0,
                    phi         = 0.0,
                    phi_cone_ctwo_oeight_ceight = 0.0;

    clipper::MAtom actual_c1, first_sugar_c1, c1, o1, c2, o2, c3, o3, c4, o4, c5, o5_next_sugar, o5, next_sugar_ring_oxygen, c6, o6, first_sugar_ring_oxygen, c7, o7, c8, o8, c9, o9;
    if ( link == 9 )
    {
        // next_sugar 2 SIA
        // first_sugar 1 SLB
        if(next_sugar.ring_members().size() == 6)
        {
            o6 = next_sugar.ring_members()[0];              // O6
            c2 = next_sugar.ring_members()[1];              // C2
            c1 = next_sugar[next_sugar.lookup("C1",clipper::MM::ANY)];
            o9 = next_sugar.anomeric_substituent();         // O9 usually
        }
        else if(next_sugar.ring_members().size() == 5)
        {
            o6 = next_sugar.ring_members()[0];              // O5
            c1 = next_sugar[next_sugar.lookup("C1",clipper::MM::ANY)];
            c2 = next_sugar.ring_members()[1];              // C2
            o9 = next_sugar.anomeric_substituent();         // O8 usually
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
            o8 = first_sugar[first_sugar.lookup("C8",clipper::MM::ANY)];
            c9 = first_sugar[first_sugar.lookup("C9",clipper::MM::ANY)];

        }
        else if(first_sugar.ring_members().size() == 5)
        {
            c6 = first_sugar.ring_members()[4];
            first_sugar_ring_oxygen = first_sugar.ring_members()[0]; // O6 for omega_seven
            c7 = first_sugar[first_sugar.lookup("C7",clipper::MM::ANY)];
            o7 = first_sugar[first_sugar.lookup("O7",clipper::MM::ANY)];
            c8 = first_sugar[first_sugar.lookup("C8",clipper::MM::ANY)];
            o8 = first_sugar[first_sugar.lookup("C8",clipper::MM::ANY)];
            c9 = first_sugar[first_sugar.lookup("C9",clipper::MM::ANY)];
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
                                               o9.coord_orth(),
                                               c9.coord_orth() );

        phi_cone_ctwo_oeight_ceight   = clipper::Coord_orth::torsion (  c1.coord_orth(),
                                                                        c2.coord_orth(),
                                                                        o9.coord_orth(),
                                                                        c9.coord_orth() );

        psi   = clipper::Coord_orth::torsion ( c2.coord_orth(),
                                               o9.coord_orth(),
                                               c9.coord_orth(),
                                               c8.coord_orth() );

        omega_seven = clipper::Coord_orth::torsion (    o7.coord_orth(),
                                                        c7.coord_orth(),
                                                        c6.coord_orth(),
                                                        first_sugar_ring_oxygen.coord_orth() );

        omega_eight = clipper::Coord_orth::torsion (    o8.coord_orth(),
                                                        c8.coord_orth(),
                                                        c7.coord_orth(),
                                                        o7.coord_orth() );


        if(first_sugar.type().trim() == "SIA" || first_sugar.type().trim() == "SLB")
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
        add_torsions_for_detected_linkages(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, first_sugar.seqnum(), next_sugar.seqnum());
        if(!torsions_zscore_database.database_array.empty())
        {
            float Phi = clipper::Util::rad2d(phi);
            float Psi = clipper::Util::rad2d(psi);
            new_connection.calculate_and_set_zscore(Phi, Psi, first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, torsions_zscore_database);
        }
    }
    else if ( link == 8 )
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


        if(first_sugar.type().trim() == "SIA" || first_sugar.type().trim() == "SLB")
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
        add_torsions_for_detected_linkages(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, first_sugar.seqnum(), next_sugar.seqnum());
        if(!torsions_zscore_database.database_array.empty())
        {
            float Phi = clipper::Util::rad2d(phi);
            float Psi = clipper::Util::rad2d(psi);
            new_connection.calculate_and_set_zscore(Phi, Psi, first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, torsions_zscore_database);
        }
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
        add_torsions_for_detected_linkages(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, first_sugar.seqnum(), next_sugar.seqnum());
        if(!torsions_zscore_database.database_array.empty())
        {
            float Phi = clipper::Util::rad2d(phi);
            float Psi = clipper::Util::rad2d(psi);
            new_connection.calculate_and_set_zscore(Phi, Psi, first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, torsions_zscore_database);
        }
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
        add_torsions_for_detected_linkages(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, first_sugar.seqnum(), next_sugar.seqnum());
        if(!torsions_zscore_database.database_array.empty())
        {
            float Phi = clipper::Util::rad2d(phi);
            float Psi = clipper::Util::rad2d(psi);
            new_connection.calculate_and_set_zscore(Phi, Psi, first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, torsions_zscore_database);
        }
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
        add_torsions_for_detected_linkages(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, first_sugar.seqnum(), next_sugar.seqnum());
        if(!torsions_zscore_database.database_array.empty())
        {
            float Phi = clipper::Util::rad2d(phi);
            float Psi = clipper::Util::rad2d(psi);
            new_connection.calculate_and_set_zscore(Phi, Psi, first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, torsions_zscore_database);
        }
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
        add_torsions_for_detected_linkages(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, first_sugar.seqnum(), next_sugar.seqnum());
        if(!torsions_zscore_database.database_array.empty())
        {
            float Phi = clipper::Util::rad2d(phi);
            float Psi = clipper::Util::rad2d(psi);
            new_connection.calculate_and_set_zscore(Phi, Psi, first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, torsions_zscore_database);
        }
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
        add_torsions_for_detected_linkages(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, first_sugar.seqnum(), next_sugar.seqnum());
        if(!torsions_zscore_database.database_array.empty())
        {
            float Phi = clipper::Util::rad2d(phi);
            float Psi = clipper::Util::rad2d(psi);
            new_connection.calculate_and_set_zscore(Phi, Psi, first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, torsions_zscore_database);
        }
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
        add_torsions_for_detected_linkages(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, first_sugar.seqnum(), next_sugar.seqnum());
        if(!torsions_zscore_database.database_array.empty())
        {
            float Phi = clipper::Util::rad2d(phi);
            float Psi = clipper::Util::rad2d(psi);
            new_connection.calculate_and_set_zscore(Phi, Psi, first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, torsions_zscore_database);
        }
    }
    else if ( link == 1 )
    {
        // { "K5B" ,	 "N", 	 "D", 	 "4,7-anhydro-3-deoxy-D-gluco-oct-2-ulosonic acid" ,"O7 C4 C5 C6 C7", 0.621, "2ev", 0.003, 2.048 },
        // { "KD5" ,	 "N", 	 "D", 	 "4,7-anhydro-3-deoxy-D-manno-oct-2-ulosonic acid", "O7 C4 C5 C6 C7", 0.500, "Oh5", 0.095, 8.063 },
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
            if(next_sugar.type().trim() == "K5B" || next_sugar.type().trim() == "KD5")
            {
                c4 = next_sugar.ring_members()[1];
                c5 = next_sugar.ring_members()[2];
                o1 = first_sugar.anomeric_substituent();
            }
            else
            {
                o5 = first_sugar.ring_members()[0];              // O5
                o5_next_sugar = next_sugar.ring_members()[0];
                // c1 = next_sugar[next_sugar.lookup("C1",clipper::MM::ANY)]; // not going to be in the ring, but going to form the glycosidic bond regardless
                c1 = next_sugar.ring_members()[1];              // C2
                actual_c1 = next_sugar[next_sugar.lookup("C1", clipper::MM::ANY)];
                o1 = next_sugar.anomeric_substituent();         // O1 usually
                c3 = next_sugar.ring_members()[2];             // C3
            }
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
            add_torsions_for_detected_linkages(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, first_sugar.seqnum(), next_sugar.seqnum());
            if(!torsions_zscore_database.database_array.empty())
            {
                float Phi = clipper::Util::rad2d(phi);
                float Psi = clipper::Util::rad2d(psi);
                new_connection.calculate_and_set_zscore(Phi, Psi, first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, torsions_zscore_database);
            }
        }
        else if(next_sugar.ring_members().size() == 5)
        {
            if(next_sugar.type().trim() == "K5B" || next_sugar.type().trim() == "KD5")
            {
                phi   = clipper::Coord_orth::torsion (  c4.coord_orth(),
                                                        c5.coord_orth(),
                                                        o1.coord_orth(),
                                                        first_sugar_c1.coord_orth() );

                // Can't figure this shit out, to come back later.
                psi   = clipper::Coord_orth::torsion (  c5.coord_orth(),
                                                        o1.coord_orth(),
                                                        first_sugar_c1.coord_orth(),
                                                        c2.coord_orth() );
            }
            else
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
            }

            new_connection.set_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), clipper::Util::rad2d(omega), clipper::Util::rad2d(omega_six), clipper::Util::rad2d(omega_seven), clipper::Util::rad2d(omega_eight), clipper::Util::rad2d(omega_nine),  clipper::Util::rad2d(phi_cone_ctwo_oeight_ceight) );
            if(debug_output)
            {
                std::vector<clipper::ftype32> torsions = new_connection.get_torsions();
                DBG << "Torsions for when next_sugar has 5 members, link = " << link << ", phi = " << torsions[0] << "\t\tpsi = " << torsions[1] << "\t\tomega = " << torsions[2] << std::endl;
            }
            new_connection.set_linkage_atoms(donorAtom, acceptorAtom);
            add_torsions_for_detected_linkages(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, first_sugar.seqnum(), next_sugar.seqnum());
            if(!torsions_zscore_database.database_array.empty())
            {
                float Phi = clipper::Util::rad2d(phi);
                float Psi = clipper::Util::rad2d(psi);
                new_connection.calculate_and_set_zscore(Phi, Psi, first_sugar.type().trim(), donorAtom, next_sugar.type().trim(), acceptorAtom, torsions_zscore_database);
            }
        }
    }

    if(noncircular)
    {
        sugars.push_back ( next_sugar );
        node_list[index].add_connection ( new_connection ); // add the new connection to the previous node
    }
    else if(!noncircular)
    {
        node_list.front().add_circular_connection(new_connection );
    }

    return false;
}

void MGlycan::add_torsions_for_detected_linkages(float Phi, float Psi, clipper::String first_residue_name, clipper::MAtom first_atom, clipper::String second_residue_name, clipper::MAtom second_atom, int first_residue_seqnum, int second_residue_seqnum)
{
    if(!all_torsions_within_mglycan.empty())
    {
        auto search_result = std::find_if(all_torsions_within_mglycan.begin(), all_torsions_within_mglycan.end(), [first_residue_name, second_residue_name](MGlycanTorsionSummary& element)
        {
            return 
            first_residue_name == element.first_residue_name && 
            second_residue_name == element.second_residue_name 
            // NEED TO ADD CHECK FOR CORRECT LINKAGE DESIGNATION HERE
            ;
        });

        if(search_result != std::end(all_torsions_within_mglycan))
        {
            MGlycanTorsionSummary& found_torsion_description = *search_result;
            std::string donorPosition = std::regex_replace(first_atom.name().trim(), std::regex(R"([^\d])"), "");
            std::string acceptorPosition = std::regex_replace(second_atom.name().trim(), std::regex(R"([^\d])"), "");
            found_torsion_description.linkage_descriptors.push_back(std::make_pair(donorPosition, acceptorPosition));
            found_torsion_description.atoms.push_back(std::make_pair(first_atom, second_atom));
            found_torsion_description.torsions.push_back(std::make_pair(Phi, Psi));

            auto torsion_seach_result = std::find_if(found_torsion_description.combined_torsions.begin(), found_torsion_description.combined_torsions.end(), [donorPosition, acceptorPosition](std::pair<std::pair<std::string, std::string>, std::vector<std::pair<float,float>>>& element)
                {
                    return donorPosition == element.first.first && 
                    acceptorPosition == element.first.second;
                });

            if (torsion_seach_result != std::end(found_torsion_description.combined_torsions)) { 
                found_torsion_description.combined_torsions[torsion_seach_result-found_torsion_description.combined_torsions.begin()].second.push_back(std::make_pair(Phi, Psi));
            }
            else { 
                std::pair<float,float> tmp_torsions = std::make_pair(Phi, Psi);
                std::vector<std::pair<float,float>> tmp_torsion_vector; 
                tmp_torsion_vector.push_back(tmp_torsions);
                found_torsion_description.combined_torsions.push_back(std::make_pair(std::make_pair(donorPosition, acceptorPosition),tmp_torsion_vector));
            }

        }
        else
        {
            MGlycanTorsionSummary new_torsion;
            std::string type;
            if(clipper::data::is_amino_acid(first_residue_name) || clipper::data::is_amino_acid(second_residue_name))
                type = "protein-sugar";
            else
                type = "sugar-sugar";

            std::string donorPosition = std::regex_replace(first_atom.name().trim(), std::regex(R"([^\d])"), "");
            std::string acceptorPosition = std::regex_replace(second_atom.name().trim(), std::regex(R"([^\d])"), "");
            new_torsion.type = type;
            new_torsion.first_residue_name = first_residue_name;
            new_torsion.second_residue_name = second_residue_name;
            new_torsion.atoms.push_back(std::make_pair(first_atom, second_atom));
            // If donorPosition adn acceptorPositon are in linkkage_ddesc
            // Go to vector vector and add to that one not the top level vector 


            std::vector<std::pair<float,float>> tmp_vector = {std::make_pair(Phi, Psi)};
            new_torsion.combined_torsions.push_back(std::make_pair(std::make_pair(donorPosition, acceptorPosition), tmp_vector));
            
            new_torsion.linkage_descriptors.push_back(std::make_pair(donorPosition, acceptorPosition));
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

        std::string donorPosition = std::regex_replace(first_atom.name().trim(), std::regex(R"([^\d])"), "");
        std::string acceptorPosition = std::regex_replace(second_atom.name().trim(), std::regex(R"([^\d])"), "");
        // std::cout << first_atom.name().trim() << " " << second_atom.name().trim() << std::endl;
        first_torsion.type = type;
        first_torsion.first_residue_name = first_residue_name;
        first_torsion.second_residue_name = second_residue_name;
        first_torsion.atoms.push_back(std::make_pair(first_atom, second_atom));

        std::vector<std::pair<float,float>> tmp_vector = {std::make_pair(Phi, Psi)};
        first_torsion.combined_torsions.push_back(std::make_pair(std::make_pair(donorPosition, acceptorPosition), tmp_vector));

        first_torsion.linkage_descriptors.push_back(std::make_pair(donorPosition, acceptorPosition));
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

    if(node_list.size() > 0)
    {
        bool glycanHasCircularLinkage = false;
        int rootCircularLinkageAttachedTo;
        for(int i = 0; i < node_list[0].number_of_connections(); i++)
        {
            if(!node_list[0].get_connection(i).connection_is_non_circular())
                glycanHasCircularLinkage = true;
                rootCircularLinkageAttachedTo = node_list[0].get_connection(i).get_linked_node_id();
                break;
        }

        if(glycanHasCircularLinkage)
        {

            clipper::MSugar lastResidue = node_list[rootCircularLinkageAttachedTo].get_sugar();
            std::string lastResidueUniqueRES = clipper::data::convert_to_wurcs_residue_code ( lastResidue.type().trim() );
            if (std::find(uniqueResidues.begin(), uniqueResidues.end(), lastResidueUniqueRES) == uniqueResidues.end()) {
                    uniqueResidues.push_back(lastResidueUniqueRES);
            }

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
        }
        else
        {
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

        bool glycanHasCircularLinkage = false;
        int rootCircularLinkageAttachedTo;
        for(int i = 0; i < node_list[0].number_of_connections(); i++)
        {
            if(!node_list[0].get_connection(i).connection_is_non_circular())
                glycanHasCircularLinkage = true;
                rootCircularLinkageAttachedTo = node_list[0].get_connection(i).get_linked_node_id();
                break;
        }

        if(glycanHasCircularLinkage)
        {
            std::string msug_wurcs_string;
            msug = node_list[0].get_sugar();
            msug_wurcs_string = clipper::data::convert_to_wurcs_residue_code(msug.type().trim());

            std::vector<std::string>::iterator residueAssigner = std::find(uniqueResidueList.begin(), uniqueResidueList.end(), msug_wurcs_string);

            int residueID = std::distance(uniqueResidueList.begin(), residueAssigner);

            wurcs_string += std::to_string(residueID + 1);
            wurcs_string += "-";

            for (int i = 1; i < node_list.size(); i++)
            {
                if(i != rootCircularLinkageAttachedTo)
                {
                    msug = node_list[i].get_sugar();
                    msug_wurcs_string = clipper::data::convert_to_wurcs_residue_code(msug.type().trim());

                    std::vector<std::string>::iterator residueAssigner = std::find(uniqueResidueList.begin(), uniqueResidueList.end(), msug_wurcs_string);

                    int residueID = std::distance(uniqueResidueList.begin(), residueAssigner);

                    wurcs_string += std::to_string(residueID + 1);
                    if ( i < (node_list.size() - 1) && ( ((i+1) != rootCircularLinkageAttachedTo || i < (node_list.size() - 2))) )
                        wurcs_string += "-";
                }
            }
        }
        else
        {
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
                if(node_list[0].get_connection(j).connection_is_non_circular())
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

                    clipper::String acceptorAtomID = linkage_atoms.second.id().trim();
                    const char *acceptorLinkagePositionChar = acceptorAtomID.c_str();

                    int acceptorLinkagePosition = std::atoi(&acceptorLinkagePositionChar[1]);

                    wurcs_string += std::to_string(acceptorLinkagePosition);
                    if(linkage_atoms.first.element().trim() == "S" || linkage_atoms.first.element().trim() == "F" || linkage_atoms.first.element().trim() == "N")
                        wurcs_string += "*" + linkage_atoms.first.element().trim() + "*";


                    wurcs_string += "_";
                }
                else
                {
                        std::ostringstream linkagePosition;
                        connectedToNodeID = node_list[0].get_connection(j).get_linked_node_id();
                        msug = node_list[connectedToNodeID].get_sugar();
                        linkagePosition << node_list[0].get_connection(j).get_order();
                        std::pair<clipper::MAtom, clipper::MAtom> linkage_atoms = node_list[0].get_connection(j).get_linkage_atoms();

                        clipper::String acceptorAtomID = linkage_atoms.first.id().trim();
                        const char *acceptorLinkagePositionChar = acceptorAtomID.c_str();

                        int acceptorLinkagePosition = std::atoi(&acceptorLinkagePositionChar[1]);

                        wurcs_string += convertNumberToLetter(0);
                        wurcs_string += std::to_string(acceptorLinkagePosition);

                        wurcs_string += "-";

                        wurcs_string += convertNumberToLetter(connectedToNodeID);

                        wurcs_string += linkagePosition.str();
                        if(linkage_atoms.first.element().trim() == "S" || linkage_atoms.first.element().trim() == "F" || linkage_atoms.first.element().trim() == "N")
                            wurcs_string += "*" + linkage_atoms.first.element().trim() + "*";


                        wurcs_string += "_";
                }
            }
        }

        // Describe the rest of the linkages.
        for (int i = 1; i < node_list.size(); i++)
        {
            if(debug_output)
                DBG << "Getting sugar at index " << i << "/" << node_list.size() << std::endl;

            msug = node_list[i].get_sugar();

            if(debug_output)
            {
                DBG << std::endl << "Type of sugar via ::MSugar.full_type() = " << msug.full_type() << std::endl;
                DBG << "Number of connections for msug/node_list[" << i << "]: " << node_list[i].number_of_connections() << std::endl;
                DBG << "Residue code via ::MSugar.type().trim() = " << msug.type().trim() << std::endl;
                DBG << "number of connections " << node_list[i].number_of_connections() << std::endl;
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

                    clipper::String acceptorAtomID = linkage_atoms.second.id().trim();
                    const char *acceptorLinkagePositionChar = acceptorAtomID.c_str();

                    int acceptorLinkagePosition = std::atoi(&acceptorLinkagePositionChar[1]);

                    wurcs_string += std::to_string(acceptorLinkagePosition);
                    if(linkage_atoms.first.element().trim() == "S" || linkage_atoms.first.element().trim() == "F" || linkage_atoms.first.element().trim() == "N")
                        wurcs_string += "*" + linkage_atoms.first.element().trim() + "*";

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

float MGlycan::calculate_zscore(float phi, float psi, privateer::json::TorsionsZScoreDatabase& matched_linkage)
{
    float count_mean = matched_linkage.summary.first;
    float count_stddev = matched_linkage.summary.second;
    std::vector<std::unordered_map<std::string, int>> bin_data = matched_linkage.bin_data;

    int count = 0;
    for(const auto& bin: bin_data) {
            if ((bin.at("lower_phi") <= phi)) {
                if (phi < bin.at("higher_phi")) {
                    if (bin.at("lower_psi") <= psi){
                        if  (psi < bin.at("higher_psi")) {
                            count = bin.at("count");
                        }
                    }
                }
            }
        }

    if (count < 0) {
        fail("Something has gone wrong with the bin count. If the result is -1, this is likely due to presumed inclusion of the current model in torsion linkage dataset, check whether the name of the file is already a PDB code.");
    }

    float z_score = (count - count_mean) / count_stddev;
    return z_score;
}


void MGlycan::Linkage::calculate_and_set_zscore(float Phi, float Psi, clipper::String first_residue_name, clipper::MAtom first_atom, clipper::String second_residue_name, clipper::MAtom second_atom, privateer::json::GlobalTorsionZScore& torsions_zscore_database)
{
    std::string donor_position = std::regex_replace(first_atom.name().trim(), std::regex(R"([^\d])"), "");
    std::string acceptor_position = std::regex_replace(second_atom.name().trim(), std::regex(R"([^\d])"), "");
    std::string donor_sugar = first_residue_name;
    std::string acceptor_sugar = second_residue_name;

    std::vector<std::string> list_of_linkages_with_enough_datapoints{
        "ASN-1,2-NAG",
        "NAG-1,4-NAG",
        "NAG-1,4-BMA",
        "BMA-1,3-MAN",
        "BMA-1,6-MAN",
        "MAN-1,2-MAN",
        "MAN-1,3-MAN",
        "MAN-1,6-MAN", 
        "NAG-1,6-FUC", 
        "NAG-1,3-FUC", 
        "MAN-1,2-NAG", 
        "NAG-1,4-GAL",
    };
    std::string linkage_name = donor_sugar + "-" + acceptor_position + "," + donor_position + "-" + acceptor_sugar;

    auto search_result_in_torsions_zscore_db = std::find_if(torsions_zscore_database.database_array.begin(), torsions_zscore_database.database_array.end(), [donor_sugar, donor_position, acceptor_position, acceptor_sugar](privateer::json::TorsionsZScoreDatabase& element)
    {
        return donor_sugar == element.donor_sugar && donor_position == element.donor_end && acceptor_position == element.acceptor_end && acceptor_sugar == element.acceptor_sugar;
    });

    auto search_result_linkage_enough_datapoints = std::find(list_of_linkages_with_enough_datapoints.begin(), list_of_linkages_with_enough_datapoints.end(), linkage_name);
    bool linkage_enough_datapoints_result = false;
    if (search_result_linkage_enough_datapoints != list_of_linkages_with_enough_datapoints.end())
        linkage_enough_datapoints_result = true;


    if(search_result_in_torsions_zscore_db != std::end(torsions_zscore_database.database_array))
    {
        privateer::json::TorsionsZScoreDatabase& found_torsion_description = *search_result_in_torsions_zscore_db;
        float linkage_score = MGlycan::Linkage::calculate_zscore(Phi, Psi, found_torsion_description);
        if(isfinite(linkage_score) && ( !isnan(linkage_score) | !isinf(linkage_score)) )
        {
            this->set_linkage_zscore(linkage_score);
            this->set_linkage_enough_datapoints(linkage_enough_datapoints_result);
        }
    }
}

float MGlycan::Linkage::calculate_zscore(float phi, float psi, privateer::json::TorsionsZScoreDatabase& matched_linkage)
{
    float count_mean = matched_linkage.summary.first;
    float count_stddev = matched_linkage.summary.second;
    std::vector<std::unordered_map<std::string, int>> bin_data = matched_linkage.bin_data;

    int count = 0;
    for(const auto& bin: bin_data) {
            if ((bin.at("lower_phi") <= phi)) {
                if (phi < bin.at("higher_phi")) {
                    if (bin.at("lower_psi") <= psi){
                        if  (psi < bin.at("higher_psi")) {
                            count = bin.at("count");
                        }
                    }
                }
            }
        }

    if (count < 0) {
        fail("Something has gone wrong with the bin count. If the result is -1, this is likely due to presumed inclusion of the current model in torsion linkage dataset, check whether the name of the file is already a PDB code.");
    }

    float z_score = (count - count_mean) / count_stddev;
    return z_score;
}


///////////////////////// MGlycology ///////////////////////////////



MGlycology::MGlycology ( const clipper::MiniMol& mmol, bool debug_output, std::string expression_system )
{
    const clipper::MAtomNonBond nb = MAtomNonBond ( mmol, 1.0 );
    privateer::json::GlobalTorsionZScore torsions_zscore_database;

    this->debug_output = debug_output;
    this->init( mmol, nb, torsions_zscore_database, debug_output, expression_system );
}


MGlycology::MGlycology ( const clipper::MiniMol& mmol, const clipper::MAtomNonBond& manb, privateer::json::GlobalTorsionZScore& torsions_zscore_database, bool debug_output, std::string expression_system )
{
    this->init(mmol, manb, torsions_zscore_database, debug_output, expression_system);
}

void MGlycology::init ( const clipper::MiniMol& mmol, const clipper::MAtomNonBond& manb, privateer::json::GlobalTorsionZScore& torsions_zscore_database, bool debug_output, std::string expression_system )
{
    this->debug_output = debug_output;
    this->manb = &manb;
    this->mmol = &mmol;

    this->expression_system = expression_system;

    std::vector<std::pair<clipper::MMonomer, clipper::String>> potential_n_roots;
    std::vector<std::pair<clipper::MMonomer, clipper::String>> potential_o_roots;
    std::vector<std::pair<clipper::MMonomer, clipper::String>> potential_s_roots;
    std::vector<std::pair<clipper::MMonomer, clipper::String>> potential_c_roots;
    std::vector<std::pair<clipper::MMonomer, clipper::String>> potential_p_roots;
    std::vector<std::pair<clipper::MMonomer, clipper::String>> potential_rootless_polysaccharides;

    std::vector < clipper::MGlycan > list_of_glycans_modelled_as_glycosylation;
    std::vector < clipper::MGlycan > list_of_glycans_modelled_as_ligands;
    std::vector < clipper::MGlycan > list_of_glycans_modelled_as_single_residues;
    std::vector < clipper::MSugar > accounted_for_sugars;

    for ( int pol = 0; pol < mmol.size() ; pol++ )
        for ( int mon = 0 ; mon < mmol[pol].size() ; mon++ )
        {
            if ( mmol[pol][mon].type().trim() == "ASN" )        potential_n_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // n-linked GlcNAc ?
            else if ( mmol[pol][mon].type().trim() == "ARG" )   potential_n_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // Arginine rhamnosylation?
            else if ( mmol[pol][mon].type().trim() == "LYS" )   potential_n_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) );
            else if ( mmol[pol][mon].type().trim() == "THR" )   potential_o_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // o-linked stuff ?
            else if ( mmol[pol][mon].type().trim() == "SER" )   potential_o_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) );
            else if ( mmol[pol][mon].type().trim() == "TYR" )   potential_o_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) );
            else if ( mmol[pol][mon].type().trim() == "ASP" )   potential_o_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) );
            else if ( mmol[pol][mon].type().trim() == "GLU" )   potential_o_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) );
            else if ( mmol[pol][mon].type().trim() == "HYP" )   potential_o_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // hydroxyproline
            else if ( mmol[pol][mon].type().trim() == "LYZ" )   potential_o_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // hydroxylysine
            else if ( mmol[pol][mon].type().trim() == "CYS" )   potential_s_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // s-linked stuff ?
            else if ( mmol[pol][mon].type().trim() == "TRP" )   potential_c_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // C-linked stuff for C/TRP-mannosylation
            else if ( mmol[pol][mon].type().trim() == "SEP" )   potential_p_roots.push_back ( std::make_pair(mmol[pol][mon], mmol[pol].id()) ); // Phosphoglycosylation - no recorded cases on PDB yet.
            else if ( clipper::data::found_in_database(mmol[pol][mon].type().trim()) && !clipper::data::is_nucleic_acid(mmol[pol][mon].type().trim()) )      potential_rootless_polysaccharides.push_back( std::make_pair(mmol[pol][mon], mmol[pol].id()) );
        }
    for ( int i = 0 ; i < potential_n_roots.size() ; i++ )  // create n-glycan roots with first sugar
    {
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_n_roots[i].first, potential_n_roots[i].second ) ;

        for ( int j = 0 ; j < linked.size() ; j++ )
        {
            const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()];

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

                    clipper::MAtom aa_atom_alpha;
                    clipper::MAtom aa_atom_bravo;

                    if(potential_n_roots[i].first.type().trim() == "ASN")
                    {
                        aa_atom_alpha = potential_n_roots[i].first.find("CG", clipper::MM::ANY);      // CG
                        aa_atom_bravo = potential_n_roots[i].first.find("CB", clipper::MM::ANY);      // CB
                        if(nd2.name().trim() == "XXX") nd2 = potential_n_roots[i].first.find("ND2", clipper::MM::ANY);
                    }
                    else if(potential_n_roots[i].first.type().trim() == "ARG")
                    {
                        aa_atom_alpha = potential_n_roots[i].first.find("CZ", clipper::MM::ANY);      // CZ
                        aa_atom_bravo = potential_n_roots[i].first.find("NE", clipper::MM::ANY);      // NE
                        if(nd2.name().trim() == "XXX") nd2 = potential_n_roots[i].first.find("NH1", clipper::MM::ANY);
                    }
                    else if(potential_n_roots[i].first.type().trim() == "LYS")
                    {
                        aa_atom_alpha = potential_n_roots[i].first.find("CE", clipper::MM::ANY);      // CE
                        aa_atom_bravo = potential_n_roots[i].first.find("CD", clipper::MM::ANY);      // CD
                        if(nd2.name().trim() == "XXX") nd2 = potential_n_roots[i].first.find("NZ", clipper::MM::ANY);
                    }
                    else
                    {
                        continue;
                    }
                    clipper::ftype phi, psi;

                    phi   = clipper::Coord_orth::torsion (  o5.coord_orth(),
                                                            c1.coord_orth(),
                                                            nd2.coord_orth(),
                                                            aa_atom_alpha.coord_orth() );

                    psi   = clipper::Coord_orth::torsion (  c1.coord_orth(),
                                                            nd2.coord_orth(),
                                                            aa_atom_alpha.coord_orth(),
                                                            aa_atom_bravo.coord_orth() );

                    if ( psi < 0 )
                        psi = clipper::Util::twopi() + psi;

                    mg.set_glycosylation_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi) );
                    mg.add_torsions_for_detected_linkages(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), potential_n_roots[i].first.type().trim(), nd2, sugar.type().trim(), c1, potential_n_roots[i].first.seqnum(), sugar.seqnum()); 

                    // This is hella cursed. A really cursed hacky implementation just to support linkage highlights in SNFG diagrams for ASN-NAG linkage.
                    // Ideally clipper::MGlycan::Linkage should have been reimplemented, but that would have taken too much time.
                    if(!torsions_zscore_database.database_array.empty())
                    {
                        mg.set_protein_sugar_linkage_zscore_attempt_to_calculate(true);
                        std::string amino_acid = potential_n_roots[i].first.type().trim();
                        std::string donor_position = std::regex_replace(nd2.name().trim(), std::regex(R"([^\d])"), "");
                        std::string first_sugar = sugar.type().trim();
                        std::string acceptor_position = std::regex_replace(c1.name().trim(), std::regex(R"([^\d])"), "");
                        auto search_result_in_torsions_zscore_db = std::find_if(torsions_zscore_database.database_array.begin(), torsions_zscore_database.database_array.end(), [amino_acid, donor_position, acceptor_position, first_sugar](privateer::json::TorsionsZScoreDatabase& element)
                        {
                            return amino_acid == element.donor_sugar && donor_position == element.donor_end && acceptor_position == element.acceptor_end && first_sugar == element.acceptor_sugar;
                        });

                        if(search_result_in_torsions_zscore_db != std::end(torsions_zscore_database.database_array))
                        {

                            privateer::json::TorsionsZScoreDatabase& found_torsion_description = *search_result_in_torsions_zscore_db;
                            float linkage_score = mg.calculate_zscore(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), found_torsion_description);
                            mg.set_protein_sugar_linkage_zscore(linkage_score);
                        }
                    }

                    if ( linked[j].second.monomer()+2 < mmol[linked[j].second.polymer()].size() )
                    {
                        if ( mmol[linked[j].second.polymer()][linked[j].second.monomer()+2].type().trim() != "THR" &&
                            mmol[linked[j].second.polymer()][linked[j].second.monomer()+2].type().trim() != "SER" )
                            // this is not a consensus ASN-GlcNAc glycosylation point
                            mg.add_root_annotation ( " Warning: this glycosylation point does not follow the Asn-X-Thr/Ser consensus sequence. ");
                    }

                    list_of_glycans_modelled_as_glycosylation.push_back ( mg );

                    break;
                }
            }
        }
    }


    for ( int i = 0 ; i < potential_o_roots.size() ; i++ )  // create o-glycan roots with first sugar
    {
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_o_roots[i].first, potential_o_roots[i].second ) ;

        for ( int j = 0 ; j < linked.size() ; j++ )
        {
            const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()];

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

                    clipper::MAtom aa_atom_alpha;
                    clipper::MAtom aa_atom_bravo;

                    if(potential_o_roots[i].first.type().trim() == "THR")
                    {
                        aa_atom_alpha = potential_o_roots[i].first.find("CB", clipper::MM::ANY);      // CB
                        aa_atom_bravo = potential_o_roots[i].first.find("CA", clipper::MM::ANY);      // CA
                        if(og1.name().trim() == "XXX") og1 = potential_o_roots[i].first.find("OG1", clipper::MM::ANY);
                    }
                    else if(potential_o_roots[i].first.type().trim() == "SER")
                    {
                        aa_atom_alpha = potential_o_roots[i].first.find("CB", clipper::MM::ANY);      // CB
                        aa_atom_bravo = potential_o_roots[i].first.find("CA", clipper::MM::ANY);      // CA
                        if(og1.name().trim() == "XXX") og1 = potential_o_roots[i].first.find("OG", clipper::MM::ANY);
                    }
                    else if(potential_o_roots[i].first.type().trim() == "TYR")
                    {
                        aa_atom_alpha = potential_o_roots[i].first.find("CZ", clipper::MM::ANY);      // CZ
                        aa_atom_bravo = potential_o_roots[i].first.find("CE1", clipper::MM::ANY);     // CE1 - come back to this after you figure out what to do about CE2
                        if(og1.name().trim() == "XXX") og1 = potential_o_roots[i].first.find("OH", clipper::MM::ANY);
                    }
                    else if(potential_o_roots[i].first.type().trim() == "ASP")
                    {
                        aa_atom_alpha = potential_o_roots[i].first.find("CG", clipper::MM::ANY);      // CG
                        aa_atom_bravo = potential_o_roots[i].first.find("CB", clipper::MM::ANY);     // CB
                        if(og1.name().trim() == "XXX") og1 = potential_o_roots[i].first.find("OD2", clipper::MM::ANY);
                    }
                    else if(potential_o_roots[i].first.type().trim() == "GLU")
                    {
                        aa_atom_alpha = potential_o_roots[i].first.find("CD", clipper::MM::ANY);      // CD
                        aa_atom_bravo = potential_o_roots[i].first.find("CG", clipper::MM::ANY);     // CG
                        if(og1.name().trim() == "XXX") og1 = potential_o_roots[i].first.find("OE2", clipper::MM::ANY);
                    }
                    else if(potential_o_roots[i].first.type().trim() == "HYP")
                    {
                        aa_atom_alpha = potential_o_roots[i].first.find("CG", clipper::MM::ANY);      // CG
                        aa_atom_bravo = potential_o_roots[i].first.find("CB", clipper::MM::ANY);     // CB - come back to this after you figure out what to do about CD
                        if(og1.name().trim() == "XXX") og1 = potential_o_roots[i].first.find("OD1", clipper::MM::ANY);
                    }
                    else if(potential_o_roots[i].first.type().trim() == "LYZ")
                    {
                        aa_atom_alpha = potential_o_roots[i].first.find("CD", clipper::MM::ANY);      // CD
                        aa_atom_bravo = potential_o_roots[i].first.find("CG", clipper::MM::ANY);     // CB - come back to this after you figure out what to do about CE
                        if(og1.name().trim() == "XXX") og1 = potential_o_roots[i].first.find("OH", clipper::MM::ANY);
                    }
                    else
                    {
                        continue;
                    }


                    clipper::ftype phi, psi;
                    // potential_o_roots[i].first.find("CB", clipper::MM::ANY);      // CB


                    phi   = clipper::Coord_orth::torsion (  o5.coord_orth(),
                                                            c1.coord_orth(),
                                                            og1.coord_orth(),
                                                            aa_atom_alpha.coord_orth() );

                    psi   = clipper::Coord_orth::torsion (  c1.coord_orth(),
                                                            og1.coord_orth(),
                                                            aa_atom_alpha.coord_orth(),
                                                            aa_atom_bravo.coord_orth() );

                    if ( psi < 0 )
                        psi = clipper::Util::twopi() + psi;

                    mg.set_glycosylation_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi) );
                    mg.add_torsions_for_detected_linkages(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), potential_o_roots[i].first.type().trim(), og1, sugar.type().trim(), c1, potential_o_roots[i].first.seqnum(), sugar.seqnum());


                    list_of_glycans_modelled_as_glycosylation.push_back ( mg );
                    break;
                }
            }
        }
    }


    for ( int i = 0 ; i < potential_s_roots.size() ; i++ )  // create s-glycan roots with first sugar
    {
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_s_roots[i].first, potential_s_roots[i].second ) ;

        for ( int j = 0 ; j < linked.size() ; j++ )
        {
            const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()];

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

                    clipper::MAtom o5 = sugar.ring_members()[0];              // O5
                    clipper::MAtom c1 = sugar.ring_members()[1];              // C1
                    clipper::MAtom sg = sugar.anomeric_substituent();         // SG CYS

                    clipper::MAtom aa_atom_alpha;
                    clipper::MAtom aa_atom_bravo;

                    if(potential_s_roots[i].first.type().trim() == "CYS")
                    {
                        aa_atom_alpha = potential_s_roots[i].first.find("CB", clipper::MM::ANY);      // CB
                        aa_atom_bravo = potential_s_roots[i].first.find("CA", clipper::MM::ANY);      // CA
                        if(sg.name().trim() == "XXX") sg = potential_s_roots[i].first.find("SG", clipper::MM::ANY);
                    }
                    else
                    {
                        continue;
                    }


                    clipper::ftype phi, psi;


                    phi   = clipper::Coord_orth::torsion (  o5.coord_orth(),
                                                            c1.coord_orth(),
                                                            sg.coord_orth(),
                                                            aa_atom_alpha.coord_orth() );

                    psi   = clipper::Coord_orth::torsion (  c1.coord_orth(),
                                                            sg.coord_orth(),
                                                            aa_atom_alpha.coord_orth(),
                                                            aa_atom_bravo.coord_orth() );

                    if ( psi < 0 )
                        psi = clipper::Util::twopi() + psi;

                    mg.set_glycosylation_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi) );
                    mg.add_torsions_for_detected_linkages(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), potential_s_roots[i].first.type().trim(), sg, sugar.type().trim(), c1, potential_s_roots[i].first.seqnum(), sugar.seqnum());


                    list_of_glycans_modelled_as_glycosylation.push_back ( mg );
                    break;
                }
            }
        }
    }

    for ( int i = 0 ; i < potential_c_roots.size() ; i++ )  // create c-glycan roots with first sugar
    {
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_c_roots[i].first, potential_c_roots[i].second ) ;

        for ( int j = 0 ; j < linked.size() ; j++ )
        {
            const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()];

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
                        DBG << "potential c roots is " << potential_c_roots[i].first.type() << std::endl;
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
                    if ( sugar.conformation_name() == "1c4" )
                    {
                        sugar.override_conformation_diag ( true );
                    }

                    clipper::MAtom o5 = sugar.ring_members()[0];              // O5
                    clipper::MAtom c1 = sugar.ring_members()[1];              // C1
                    clipper::MAtom cd1 = sugar.anomeric_substituent();         // CD1 TRP

                    clipper::MAtom aa_atom_alpha;
                    clipper::MAtom aa_atom_bravo;

                    if(potential_c_roots[i].first.type().trim() == "TRP")
                    {
                        aa_atom_alpha = potential_c_roots[i].first.find("CG", clipper::MM::ANY);      // CB
                        aa_atom_bravo = potential_c_roots[i].first.find("CB", clipper::MM::ANY);      // CA
                        if(cd1.name().trim() == "XXX") cd1 = potential_c_roots[i].first.find("CD1", clipper::MM::ANY);
                    }
                    else
                    {
                        continue;
                    }

                    clipper::ftype phi, psi;

                    phi   = clipper::Coord_orth::torsion (  o5.coord_orth(),
                                                            c1.coord_orth(),
                                                            cd1.coord_orth(),
                                                            aa_atom_alpha.coord_orth() );

                    psi   = clipper::Coord_orth::torsion (  c1.coord_orth(),
                                                            cd1.coord_orth(),
                                                            aa_atom_alpha.coord_orth(),
                                                            aa_atom_bravo.coord_orth() );


                    if ( psi < 0 )
                        psi = clipper::Util::twopi() + psi;


                    mg.set_glycosylation_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi) );
                    mg.add_torsions_for_detected_linkages(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), potential_c_roots[i].first.type().trim(), cd1, sugar.type().trim(), c1, potential_c_roots[i].first.seqnum(), sugar.seqnum());

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

                    list_of_glycans_modelled_as_glycosylation.push_back ( mg );
                    break;
                }
            }
        }
    }

    for ( int i = 0 ; i < potential_p_roots.size() ; i++ )  // create c-glycan roots with first sugar
    {
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_p_roots[i].first, potential_p_roots[i].second ) ;

        for ( int j = 0 ; j < linked.size() ; j++ )
        {
            const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()];

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
                        DBG << "potential p roots is " << potential_p_roots[i].first.type() << std::endl;
                        DBG << "sugar is " << sugar.type() << std::endl;
                        DBG << "id is " << mmol[linked[j].second.polymer()].id().trim() << std::endl;
                    }

                    clipper::String root_sugar_chain_id = mmol[linked[j].second.polymer()].id().trim().substr(0,1);
                    clipper::MGlycan mg (   potential_p_roots[i].second,
                                            potential_p_roots[i].first,
                                            list_of_sugars.back(),
                                            root_sugar_chain_id,
                                            debug_output,
                                            this->expression_system );
                    mg.set_kind_of_glycan ( "p-glycan" );

                    clipper::MAtom o5 = sugar.ring_members()[0];              // O5
                    clipper::MAtom c1 = sugar.ring_members()[1];              // C1
                    clipper::MAtom op = sugar.anomeric_substituent();         // O3P/O2P

                    clipper::MAtom aa_atom_alpha;
                    clipper::MAtom aa_atom_bravo;

                    if(potential_p_roots[i].first.type().trim() == "SEP")
                    {
                        aa_atom_alpha = potential_p_roots[i].first.find("P", clipper::MM::ANY);      // P
                        aa_atom_bravo = potential_p_roots[i].first.find("OG", clipper::MM::ANY);      // OG
                    }
                    else
                    {
                        continue;
                    }

                    clipper::ftype phi, psi;

                    phi   = clipper::Coord_orth::torsion (  o5.coord_orth(),
                                                            c1.coord_orth(),
                                                            op.coord_orth(),
                                                            aa_atom_alpha.coord_orth() );

                    psi   = clipper::Coord_orth::torsion (  c1.coord_orth(),
                                                            op.coord_orth(),
                                                            aa_atom_alpha.coord_orth(),
                                                            aa_atom_bravo.coord_orth() );

                    if ( psi < 0 )
                        psi = clipper::Util::twopi() + psi;

                    mg.set_glycosylation_torsions ( clipper::Util::rad2d(phi), clipper::Util::rad2d(psi) );
                    mg.add_torsions_for_detected_linkages(clipper::Util::rad2d(phi), clipper::Util::rad2d(psi), potential_p_roots[i].first.type().trim(), op, sugar.type().trim(), c1, potential_p_roots[i].first.seqnum(), sugar.seqnum());

                    list_of_glycans_modelled_as_glycosylation.push_back ( mg );
                    break;
                }
            }
        }
    }

    // Might have to split extension process into two, first traditional, then ligands only. Delete roots from reserved_sugars that have been assigned to other glycans.
    for ( int i = 0 ; i < list_of_glycans_modelled_as_glycosylation.size() ; i++ )
    {
        std::vector < clipper::MSugar >& sugar_list = list_of_glycans_modelled_as_glycosylation[i].get_sugars();
        clipper::MSugar first_sugar = clipper::MSugar ( sugar_list.front() );

        if(debug_output)
        {
            DBG << "Extending tree of glycan from root of " << list_of_glycans_modelled_as_glycosylation[i].get_root_for_filename() << "\twith detected sugar = " << first_sugar.type() << "-" << first_sugar.id() << "/" <<first_sugar.get_seqnum() << std::endl;
        }

        if(first_sugar.ring_members().size() == 5 || first_sugar.ring_members().size() == 6)
            extend_tree ( list_of_glycans_modelled_as_glycosylation[i] , first_sugar, accounted_for_sugars, torsions_zscore_database );

        list_of_glycans_modelled_as_glycosylation[i].set_annotations( this->expression_system );

        sugar_list = list_of_glycans_modelled_as_glycosylation[i].get_sugars(); // update the sugars
        accounted_for_sugars.insert(accounted_for_sugars.end(), sugar_list.begin(), sugar_list.end());
    }

    // pretty cursed STL function to remove sugars from vectors that have already been included in other clipper::MGlycan objects contained within list_of_glycans_modelled_as_glycosylation
    potential_rootless_polysaccharides.erase(std::remove_if(potential_rootless_polysaccharides.begin(), potential_rootless_polysaccharides.end(),
                                                            [&accounted_for_sugars](std::pair<clipper::MMonomer, clipper::String>& rootless_glycan) {
                                                                if (std::find_if(accounted_for_sugars.begin(), accounted_for_sugars.end(),
                                                                        [&rootless_glycan](clipper::MSugar& accounted_for_sugar){
                                                                        return  accounted_for_sugar.chain_id().trim() == rootless_glycan.second.trim() &&
                                                                                accounted_for_sugar.id().trim() == rootless_glycan.first.id().trim() &&
                                                                                accounted_for_sugar.type().trim() == rootless_glycan.first.type().trim() &&
                                                                                accounted_for_sugar.seqnum() == rootless_glycan.first.seqnum(); }) != accounted_for_sugars.end())
                                                                {
                                                                    return true;
                                                                }
                                                                else
                                                                {
                                                                    return false;
                                                                }
                                                            }), potential_rootless_polysaccharides.end());


    for ( int i = 0 ; i < potential_rootless_polysaccharides.size(); i++ )  // create ligand glycan roots with first sugar
    {
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_rootless_polysaccharides[i].first, potential_rootless_polysaccharides[i].second ) ;

        if(!linked.empty())
        {
            for ( int j = 0 ; j < linked.size() ; j++ )
            {
                const clipper::MMonomer& tmpmon = mmol[linked[j].second.polymer()][linked[j].second.monomer()];


                // Might need to revert this change. To test with running Privateer again on a local pdb_mirror.
                if ( clipper::MSugar::search_database(potential_rootless_polysaccharides[i].first.type().c_str()) && ( !clipper::data::is_nucleic_acid(potential_rootless_polysaccharides[i].first.type().trim()) && !clipper::data::is_nucleic_acid(tmpmon.type().trim()) ) )
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
                    list_of_glycans_modelled_as_ligands.push_back ( mg );
                    break;
                }
            }
        }
    }

    std::vector<int> vector_indices_to_remove;
    for ( int i = 0 ; i < list_of_glycans_modelled_as_ligands.size() ; i++ )
    {
        std::vector < clipper::MSugar >& sugar_list = list_of_glycans_modelled_as_ligands[i].get_sugars();
        clipper::MSugar first_sugar = clipper::MSugar ( sugar_list.front() );

        if(debug_output)
        {
            DBG << "Extending tree of glycan from root of " << list_of_glycans_modelled_as_ligands[i].get_root_for_filename() << "\twith detected sugar = " << first_sugar.type() << "-" << first_sugar.id() << "/" <<first_sugar.get_seqnum() << std::endl;
        }

        int initialLength = list_of_glycans_modelled_as_ligands[i].get_sugars().size();
        if(first_sugar.ring_members().size() == 5 || first_sugar.ring_members().size() == 6)
            extend_tree ( list_of_glycans_modelled_as_ligands[i] , first_sugar, accounted_for_sugars, torsions_zscore_database );
        int lengthAfterTreeExtension = list_of_glycans_modelled_as_ligands[i].get_sugars().size();

        if(lengthAfterTreeExtension == initialLength)
            vector_indices_to_remove.push_back(i);
        else
        {
            list_of_glycans_modelled_as_ligands[i].set_annotations( this->expression_system );

            sugar_list = list_of_glycans_modelled_as_ligands[i].get_sugars(); // update the sugars
            accounted_for_sugars.insert(accounted_for_sugars.end(), sugar_list.begin(), sugar_list.end());
        }
    }

    for (int i = vector_indices_to_remove.size() - 1; i >= 0; i--)
    {
        list_of_glycans_modelled_as_ligands.erase(list_of_glycans_modelled_as_ligands.begin()+vector_indices_to_remove[i]);
    }

    // pretty cursed STL function to remove sugars from vectors that have already been included in other clipper::MGlycan objects contained within list_of_glycans_modelled_as_glycosylation
    potential_rootless_polysaccharides.erase(std::remove_if(potential_rootless_polysaccharides.begin(), potential_rootless_polysaccharides.end(),
                                                            [&accounted_for_sugars](std::pair<clipper::MMonomer, clipper::String>& rootless_glycan) {
                                                                if (std::find_if(accounted_for_sugars.begin(), accounted_for_sugars.end(),
                                                                        [&rootless_glycan](clipper::MSugar& accounted_for_sugar){
                                                                        return  accounted_for_sugar.chain_id().trim() == rootless_glycan.second.trim() &&
                                                                                accounted_for_sugar.id().trim() == rootless_glycan.first.id().trim() &&
                                                                                accounted_for_sugar.type().trim() == rootless_glycan.first.type().trim() &&
                                                                                accounted_for_sugar.seqnum() == rootless_glycan.first.seqnum(); }) != accounted_for_sugars.end())
                                                                {
                                                                    return true;
                                                                }
                                                                else
                                                                {
                                                                    return false;
                                                                }
                                                            }), potential_rootless_polysaccharides.end());


    for ( int i = 0 ; i < potential_rootless_polysaccharides.size(); i++ )  // create ligand glycan roots with first sugar
    {
        std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > linked = get_contacts ( potential_rootless_polysaccharides[i].first, potential_rootless_polysaccharides[i].second ) ;

        if(linked.empty() && !clipper::data::is_nucleic_acid(potential_rootless_polysaccharides[i].first.type().trim()))
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
                list_of_glycans_modelled_as_single_residues.push_back ( mg );
            }
        }
    }

    // accounted_for_sugars.insert(accounted_for_sugars.end(), sugar_list.begin(), sugar_list.end());
    this->list_of_glycans.insert(this->list_of_glycans.end(), list_of_glycans_modelled_as_glycosylation.begin(), list_of_glycans_modelled_as_glycosylation.end());
    this->list_of_glycans.insert(this->list_of_glycans.end(), list_of_glycans_modelled_as_ligands.begin(), list_of_glycans_modelled_as_ligands.end());
    this->list_of_glycans.insert(this->list_of_glycans.end(), list_of_glycans_modelled_as_single_residues.begin(), list_of_glycans_modelled_as_single_residues.end());
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

void MGlycology::extend_tree ( clipper::MGlycan& mg, clipper::MSugar& msug, std::vector<clipper::MSugar>& accounted_for_sugars, privateer::json::GlobalTorsionZScore& torsions_zscore_database )
{
    if(debug_output)
    {
        DBG << "clipper::MGlycan mg.number_of_nodes() = " << mg.number_of_nodes() << "\tmsug.id() = " << msug.id() << std::endl;
    }
    std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > contacts = get_contacts ( msug, msug.chain_id() );

    const clipper::MiniMol& tmpmol = *this->mmol;

    if(debug_output)
    {
        DBG << "contacts.size() = " << contacts.size() << std::endl;
        for (int i = 0 ; i < contacts.size() ; i++ )
        {
            DBG << "Detected contacts[" << i << "].first.id().trim() from msug = " << msug.id().trim() << msug.type().trim() << "- " << contacts[i].first.id().trim() << "\t\t\tcontacts[" << i << "].second....id().trim() = " << tmpmol[contacts[i].second.polymer()].id().trim() << "/" << tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()].type().trim() << "-" << tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()].id() << "-" << tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()][contacts[i].second.atom()].id() << std::endl;
            DBG << "Distance of contacts[" << i << "].first and  = " << "contacts[" << i << "].second = " << clipper::Coord_orth::length(contacts[i].first.coord_orth(), tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()][contacts[i].second.atom()].coord_orth()) << std::endl;
        }
    }

    bool currentTerminalSugarAlreadyConnectedToRoot = false;

    for (int i = 0 ; i < contacts.size() ; i++ )
    {
        if (clipper::data::found_in_database ( tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()].type() ) && !clipper::data::is_nucleic_acid(tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()].type().trim()))
        {
            const std::vector<clipper::MSugar> sugar_list = mg.get_sugars();
            const clipper::MSugar root_sugar = sugar_list.front();
            auto search_result_already_part_of_same_glycan = std::find_if(sugar_list.begin(), sugar_list.end(), [&](const clipper::MSugar& element)
            {
                return      element.chain_id().trim() == tmpmol[contacts[i].second.polymer()].id().trim() &&
                            element.id().trim() == tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()].id().trim() &&
                            element.type().trim() == tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()].type().trim() &&
                            element.seqnum() == tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()].seqnum();
            });

            auto search_result_msug_part_of_same_glycan = std::find_if(sugar_list.begin(), sugar_list.end(), [&](const clipper::MSugar& element)
            {
                return      element.chain_id().trim() == msug.chain_id().trim() &&
                            element.id().trim() == msug.id().trim() &&
                            element.type().trim() == msug.type().trim() &&
                            element.seqnum() == msug.seqnum();
            });

            bool terminal_sugar_connected_to_root =     ( root_sugar.chain_id().trim() == tmpmol[contacts[i].second.polymer()].id().trim() &&
                                                        root_sugar.id() == tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()].id() &&
                                                        root_sugar.type().trim() == tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()].type().trim() &&
                                                        root_sugar.seqnum() == tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()].seqnum() &&
                                                        sugar_list.size() > 3 && std::distance(sugar_list.begin(), search_result_msug_part_of_same_glycan) > 2 ) &&
                                                        ( msug.type().trim() != "FUC" && msug.type().trim() != "FUL" );

            if(debug_output)
            {
                DBG << std::boolalpha << "terminal_sugar_connected_to_root = " << terminal_sugar_connected_to_root << "\tstd::distance = " << std::distance(sugar_list.begin(), search_result_msug_part_of_same_glycan) << std::endl;
            }

            if( ( search_result_already_part_of_same_glycan == std::end(sugar_list) ) && !terminal_sugar_connected_to_root)
            {
                clipper::MSugar tmpsug = clipper::MSugar ( *this->mmol, tmpmol[contacts[i].second.polymer()].id().trim(), tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()], *this->manb, debug_output );
                if(tmpsug.ring_members().size() == 5 || tmpsug.ring_members().size() == 6)
                {
                    // check if sugar is not a part of another glycan
                    auto search_result_part_of_another_glycan = std::find_if(accounted_for_sugars.begin(), accounted_for_sugars.end(), [&tmpsug](clipper::MSugar& reserved_sugar)
                    {
                        return  tmpsug.chain_id() == reserved_sugar.chain_id() &&
                                tmpsug.id().trim() == reserved_sugar.id().trim() &&
                                tmpsug.type().trim() == reserved_sugar.type().trim() &&
                                tmpsug.seqnum() == reserved_sugar.seqnum();

                    });
                    if(search_result_part_of_another_glycan == std::end(accounted_for_sugars))
                    {
                        clipper::MAtom acceptorAtom = tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()][contacts[i].second.atom()];
                        if(get_altconf(contacts[i].first) != ' ' && get_altconf(acceptorAtom) != ' ')
                        {
                            if(altconf_compatible(get_altconf(contacts[i].first), get_altconf(acceptorAtom)))
                            {
                                if(debug_output)
                                {
                                    DBG << "parse_order ( contacts[" << i << "].first.id()) = " << parse_order ( contacts[i].first, msug ) << std::endl;
                                }
                                mg.link_sugars ( parse_order ( contacts[i].first, msug ), msug, tmpsug, contacts[i].first, acceptorAtom, true, torsions_zscore_database);
                                extend_tree ( mg, tmpsug, accounted_for_sugars, torsions_zscore_database );
                            }
                        }
                        else
                        {
                            if(debug_output)
                            {
                                DBG << "parse_order ( contacts[" << i << "].first.id()) = " << parse_order ( contacts[i].first, msug ) << std::endl;
                            }
                            mg.link_sugars ( parse_order ( contacts[i].first, msug ), msug, tmpsug, contacts[i].first, acceptorAtom, true, torsions_zscore_database);
                            extend_tree ( mg, tmpsug, accounted_for_sugars, torsions_zscore_database );
                        }
                    }
                }
            }
            else if(terminal_sugar_connected_to_root && !currentTerminalSugarAlreadyConnectedToRoot)
            {
                clipper::MGlycan::Node rootNode = mg.get_node(0);
                bool currentPairConnected = false;
                for(int j = 0; j < rootNode.number_of_connections(); j++)
                {
                    int connectedToNodeID = rootNode.get_connection(j).get_linked_node_id();
                    clipper::MSugar rootConnectedToOtherMSugar = mg.get_node(connectedToNodeID).get_sugar();

                    if( msug.chain_id().trim() == rootConnectedToOtherMSugar.chain_id().trim() &&
                        msug.id().trim() == rootConnectedToOtherMSugar.id().trim() &&
                        msug.type().trim() == rootConnectedToOtherMSugar.type().trim() &&
                        msug.seqnum() == rootConnectedToOtherMSugar.seqnum() )
                    {
                        currentPairConnected = true;
                        break;
                    }
                }

                if(!currentPairConnected)
                {
                    clipper::MSugar tmpsug = sugar_list.front();
                    clipper::MAtom acceptorAtom = tmpmol[contacts[i].second.polymer()][contacts[i].second.monomer()][contacts[i].second.atom()];
                    if(get_altconf(contacts[i].first) != ' ' && get_altconf(acceptorAtom) != ' ')
                    {
                        if(altconf_compatible(get_altconf(contacts[i].first), get_altconf(acceptorAtom)))
                        {
                            if(debug_output)
                            {
                                DBG << "parse_order - terminal/root connection ( contacts[" << i << "].first.id()) = " << parse_order ( contacts[i].first, msug ) << std::endl;
                            }
                            mg.link_sugars ( parse_order ( contacts[i].first, msug ), tmpsug, msug, acceptorAtom, contacts[i].first, false, torsions_zscore_database);
                        }
                    }
                    else
                    {
                        if(debug_output)
                        {
                            DBG << "parse_order - terminal/root connection ( contacts[" << i << "].first.id()) = " << parse_order ( contacts[i].first, msug ) << std::endl;
                        }
                        mg.link_sugars ( parse_order ( contacts[i].first, msug ), tmpsug, msug, acceptorAtom, contacts[i].first, false, torsions_zscore_database);
                    }
                    currentTerminalSugarAlreadyConnectedToRoot = true;
                }
            }
        }
    }
}

int MGlycology::parse_order(clipper::MAtom& atom_in_sugar, clipper::MSugar& sugar)
{
    const clipper::MiniMol& tmpmol = *mmol;

    std::vector < clipper::MAtomIndexSymmetry > contacts = this->manb->atoms_near ( atom_in_sugar.coord_orth(), 2.5 );
    std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > tmpresults;
    for (int j = 0 ; j < contacts.size() ; j++ )
    {
        if ( (  (tmpmol[contacts[j].polymer()].id().trim() == sugar.chain_id().trim())
            &&  (tmpmol[contacts[j].polymer()][contacts[j].monomer()].id().trim() == sugar.id().trim())
            &&  (tmpmol[contacts[j].polymer()][contacts[j].monomer()].type().trim() == sugar.type().trim())
            &&  (tmpmol[contacts[j].polymer()][contacts[j].monomer()].seqnum() == sugar.seqnum())
            &&  (tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()].element().trim() == "C")
            )   &&  (clipper::Coord_orth::length ( tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()].coord_orth(), atom_in_sugar.coord_orth() ) <= 2.5 )
                &&  (contacts[j].symmetry() == 0))  // Beware: will report contacts that are not physically in contact, but needed for visualisation
        {                            //         of crappy structures in MG
            if ( altconf_compatible(get_altconf(tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()]), get_altconf(atom_in_sugar) ) )
            {
                std::pair < clipper::MAtom , clipper::MAtomIndexSymmetry > link_tmp;
                link_tmp.first = atom_in_sugar;
                link_tmp.second = contacts[j];
                tmpresults.push_back ( link_tmp );
                if(debug_output)
                {
                    DBG << "atom_in_sugar = " << atom_in_sugar.id().trim() << "\tlink_tmp.second(MAtomNonBond) aka tmpAtom = " << tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()].id().trim() << std::endl;
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

    if(!tmpresults.empty())
    {
        std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > closest_pair = tmpresults.front();

        clipper::MAtom atom_to_parse = tmpmol[closest_pair.second.polymer()][closest_pair.second.monomer()][closest_pair.second.atom()];
        clipper::String atom_to_parse_id = atom_to_parse.id().trim();
        const char *s = atom_to_parse_id.c_str();

        int result = std::atoi(&s[1]);

        if(debug_output)
        {
            DBG << "Returning parsed linkage identifier: " << result << "\t from atom_to_parse_id " << atom_to_parse_id << std::endl;
        }

        return result;
    }
    else
    {
        if(debug_output)
        {
            DBG << "Uh oh, unable to find closest atom within same sugar" << std::endl;
        }

        clipper::MAtom atom_to_parse = atom_in_sugar;
        clipper::String atom_to_parse_id = atom_to_parse.id().trim();
        const char *s = atom_to_parse_id.c_str();

        int result = std::atoi(&s[1]);

        if(debug_output)
        {
            DBG << "Returning parsed linkage identifier: " << result << "\t from atom_to_parse_id " << atom_to_parse_id << std::endl;
        }

        return result;
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
        int id = mm.lookup ( "O", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O1", clipper::MM::ANY );

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

        id = mm.lookup ( "O9", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S", clipper::MM::ANY );

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

        id = mm.lookup ( "S9", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N", clipper::MM::ANY );

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

        id = mm.lookup ( "N9", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F", clipper::MM::ANY );

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

        id = mm.lookup ( "F9", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C1", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C2", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C3", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C4", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C5", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C6", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C7", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C8", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C9", clipper::MM::ANY );

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
                if(contacts[j].symmetry() == 0)
                {
                    if ((tmpmol[contacts[j].polymer()][contacts[j].monomer()].id().trim() != mm.id().trim())
                    && ( clipper::Coord_orth::length ( tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()].coord_orth(), candidates[i].coord_orth() ) < 0.25 ))  // Beware: will report contacts that are not physically in contact, but needed for visualisation
                    {                            //         of crappy structures in MG
                        tmpresults.push_back(std::make_pair(tmpmol[contacts[j].polymer()].id().trim(), tmpmol[contacts[j].polymer()][contacts[j].monomer()]));
                    }
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


const std::vector < std::pair< clipper::MAtom, clipper::MAtomIndexSymmetry > > MGlycology::get_contacts ( const clipper::MMonomer& mm, const clipper::String monomer_chain_id )
{
    // mm can be aminoacid (typically ASN) or sugar: ND2, O2, O3, O4, O6

    std::vector < clipper::MAtom > candidates; // look for the possibilities
    std::vector < std::pair < clipper::MAtom, clipper::MAtomIndexSymmetry > > tmpresults;

    if(debug_output)
        DBG << "Input: " << monomer_chain_id << "/" << mm.id().trim() << "-" << mm.type().trim() << std::endl;

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
    else if ( mm.type().trim() == "ARG" )
    {
        int id1 = mm.lookup ( "NH2", clipper::MM::ANY );
        int id2 = mm.lookup ( "NH1", clipper::MM::ANY );

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
    else if ( mm.type().trim() == "LYS" )
    {
        int id = mm.lookup ( "NZ", clipper::MM::ANY );

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
    else if ( mm.type().trim() == "SER" )
    {
        int id = mm.lookup ( "OG", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );
        else return tmpresults; // empty result

    }
    else if ( mm.type().trim() == "TYR" )
    {
        int id = mm.lookup ( "OH", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );
        else return tmpresults; // empty result

    }
    else if ( mm.type().trim() == "ASP" )
    {
        int id1 = mm.lookup ( "OD2", clipper::MM::ANY );
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
    else if ( mm.type().trim() == "GLU" )
    {
        int id1 = mm.lookup ( "OE2", clipper::MM::ANY );
        int id2 = mm.lookup ( "OE1", clipper::MM::ANY );

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
    else if ( mm.type().trim() == "HYP" )
    {
        int id = mm.lookup ( "OD1", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );
        else return tmpresults; // empty result

    }
    else if ( mm.type().trim() == "LYZ" )
    {
        int id = mm.lookup ( "OH", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );
        else return tmpresults; // empty result

    }
    else if ( mm.type().trim() == "CYS" )
    {
        int id = mm.lookup ( "SG", clipper::MM::ANY );

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
    else if ( mm.type().trim() == "SEP" )
    {
        int id1 = mm.lookup ( "O2P", clipper::MM::ANY );
        int id2 = mm.lookup ( "O3P", clipper::MM::ANY );

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
    else
    {
        int id = mm.lookup ( "O", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O1", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O1A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O1B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O2", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O2A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O2B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O3", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O3A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O3B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O4", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O4A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O4B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O5", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O5A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O5B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O6", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O6A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O6B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O7", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O7A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O7B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O8", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O8A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O8B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O9", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O9A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "O9B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S1", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S1A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S1B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S2", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S2A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S2B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S3", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S3A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S3B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S4", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S4A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S4B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S5", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S5A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S5B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S6", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S6A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S6B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S7", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S7A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S7B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S8", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S8A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S8B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S9", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S9A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "S9B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N1", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N1A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N1B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N2", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N2A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N2B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N3", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N3A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N3B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N4", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N4A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N4B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N5", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N5A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N5B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N6", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N6A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N6B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N7", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N7A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N7B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N8", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N8A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N8B", clipper::MM::ANY );

        id = mm.lookup ( "N9", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N9A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "N9B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F1", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F1A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F1B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F2", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F2A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F2B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F3", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F3A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F3B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F4", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F4A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F4B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F5", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F5A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F5B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F6", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F6A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F6B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F7", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F7A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F7B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F8", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F8A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F8B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F9", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F9A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "F9B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C1", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C1A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C1B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C2", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C2A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C2B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C3", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C3A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C3B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C4", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C4A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C4B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C5", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C5A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C5B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C6", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C6A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C6B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C7", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C7A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C7B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C8", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C8A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C8B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C9", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C9A", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

        id = mm.lookup ( "C9B", clipper::MM::ANY );

        if ( id != -1 )
            candidates.push_back ( mm[id] );

    }

    if ( candidates.size() == 0 )
        return tmpresults;  // empty result
    else
    {
        const clipper::MiniMol& tmpmol = *mmol;

        char prefered_altconf = ' ';
        for ( int i = 0 ; i < candidates.size() ; i++ )
        {
            if(candidates[i].element().trim() != "C")
            {
                std::vector < clipper::MAtomIndexSymmetry > contacts = this->manb->atoms_near ( candidates[i].coord_orth(), 2.0 );

                bool altConfDetected = false;
                bool managedToFindContactInSameAltConf = false;
                for (int j = 0 ; j < contacts.size() ; j++ )
                {
                    if ( (  (tmpmol[contacts[j].polymer()].id().trim() != monomer_chain_id.trim())
                        ||  (tmpmol[contacts[j].polymer()][contacts[j].monomer()].id().trim() != mm.id().trim())
                        ||  (tmpmol[contacts[j].polymer()][contacts[j].monomer()].type().trim() != mm.type().trim())
                        ||  (tmpmol[contacts[j].polymer()][contacts[j].monomer()].seqnum() != mm.seqnum())
                        )   &&  (clipper::Coord_orth::length ( tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()].coord_orth(), candidates[i].coord_orth() ) <= 2.0 )
                            &&  (contacts[j].symmetry() == 0) )  // Beware: will report contacts that are not physically in contact, but needed for visualisation
                    //  if (   (tmpmol[contacts[j].polymer()][contacts[j].monomer()].id().trim() != mm.id().trim())
                    //     &&  (clipper::Coord_orth::length ( tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()].coord_orth(), candidates[i].coord_orth() ) <= 2.5 ))  // Beware: will report contacts that are not physically in contact, but needed for visualisation
                    {                            //         of crappy structures in MG
                        if ( altconf_compatible( get_altconf(tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()]), get_altconf(candidates[i]) ) )
                        {
                            clipper::MAtom tmpAtom = tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()];

                            if(debug_output)
                                DBG << "tmpAtom.id().trim() = " << tmpAtom.id().trim() << "\tfrom: " << tmpmol[contacts[j].polymer()].id().trim() << "/" << tmpmol[contacts[j].polymer()][contacts[j].monomer()].id().trim() << "-" << tmpmol[contacts[j].polymer()][contacts[j].monomer()].type().trim() << std::endl;

                            if(get_altconf(tmpAtom) != ' ')
                            {

                                if(tmpAtom.element().trim() == "C")
                                {
                                    std::pair < clipper::MAtom , clipper::MAtomIndexSymmetry > link_tmp;
                                    link_tmp.first = candidates[i];
                                    link_tmp.second = contacts[j];
                                    tmpresults.push_back ( link_tmp );
                                    managedToFindContactInSameAltConf = true;
                                    if(debug_output)
                                    {
                                        DBG << "Contact detected between atoms of altconf = " << get_altconf(tmpAtom) << std::endl;
                                        DBG << "link_tmp.first(MAtom).id() = " << candidates[i].id() << "\tlink_tmp.second(MAtomNonBond) aka tmpAtom = " << tmpAtom.id() << std::endl;
                                    }
                                }
                            }
                            else
                            {
                                if(tmpAtom.element().trim() == "C")
                                {
                                    std::pair < clipper::MAtom , clipper::MAtomIndexSymmetry > link_tmp;
                                    link_tmp.first = candidates[i];
                                    link_tmp.second = contacts[j];
                                    tmpresults.push_back ( link_tmp );
                                    managedToFindContactInSameAltConf = true;
                                    if(debug_output)
                                    {
                                        DBG << "link_tmp.first(MAtom).id() = " << candidates[i].id() << "\tlink_tmp.second(MAtomNonBond) aka tmpAtom = " << tmpAtom.id() << std::endl;
                                    }
                                }
                            }
                        }
                        else if( !altconf_compatible( get_altconf(tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()]), get_altconf(candidates[i]) ) )
                        {
                            altConfDetected = true;
                        }
                    }
                }

                if (altConfDetected && !managedToFindContactInSameAltConf)
                {
                    for (int j = 0 ; j < contacts.size() ; j++ )
                    {
                        if ( (  (tmpmol[contacts[j].polymer()].id().trim() != monomer_chain_id.trim())
                            ||  (tmpmol[contacts[j].polymer()][contacts[j].monomer()].id().trim() != mm.id().trim())
                            ||  (tmpmol[contacts[j].polymer()][contacts[j].monomer()].type().trim() != mm.type().trim())
                            ||  (tmpmol[contacts[j].polymer()][contacts[j].monomer()].seqnum() != mm.seqnum())
                            )  &&  (clipper::Coord_orth::length ( tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()].coord_orth(), candidates[i].coord_orth() ) <= 2.0 ))  // Beware: will report contacts that are not physically in contact, but needed for visualisation
                        //  if (   (tmpmol[contacts[j].polymer()][contacts[j].monomer()].id().trim() != mm.id().trim())
                        //     &&  (clipper::Coord_orth::length ( tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()].coord_orth(), candidates[i].coord_orth() ) <= 2.5 ))  // Beware: will report contacts that are not physically in contact, but needed for visualisation
                        {                            //         of crappy structures in MG

                            clipper::MAtom tmpAtom = tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()];
                            if(prefered_altconf == ' ' && get_altconf(tmpAtom) != ' ')
                                prefered_altconf = get_altconf(tmpAtom);

                            if(debug_output)
                                DBG << "tmpAtom.id().trim(), but this time with altconfs ayyy lmao = " << tmpAtom.id().trim() << "\tfrom: " << tmpmol[contacts[j].polymer()].id().trim() << "/" << tmpmol[contacts[j].polymer()][contacts[j].monomer()].id().trim() << "-" << tmpmol[contacts[j].polymer()][contacts[j].monomer()].type().trim() << "\tcurrent prefered_altConf = " << prefered_altconf << std::endl;

                            if(get_altconf(tmpAtom) == prefered_altconf && get_altconf(tmpAtom) != ' ')
                            {
                                const char altconf = get_altconf(tmpAtom);

                                if(tmpAtom.element().trim() == "C")
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
                            else if(get_altconf(tmpAtom) == ' ')
                            {
                                if(tmpAtom.element().trim() == "C")
                                {
                                    std::pair < clipper::MAtom , clipper::MAtomIndexSymmetry > link_tmp;
                                    link_tmp.first = candidates[i];
                                    link_tmp.second = contacts[j];
                                    tmpresults.push_back ( link_tmp );
                                    managedToFindContactInSameAltConf = true;
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
            else if(candidates[i].element().trim() == "C")
            {
                std::vector < clipper::MAtomIndexSymmetry > contacts = this->manb->atoms_near ( candidates[i].coord_orth(), 2.5 );

                bool altConfDetected = false;
                bool managedToFindContactInSameAltConf = false;
                for (int j = 0 ; j < contacts.size() ; j++ )
                {
                    if ( (  (tmpmol[contacts[j].polymer()].id().trim() != monomer_chain_id.trim())
                        ||  (tmpmol[contacts[j].polymer()][contacts[j].monomer()].id().trim() != mm.id().trim())
                        ||  (tmpmol[contacts[j].polymer()][contacts[j].monomer()].type().trim() != mm.type().trim())
                        ||  (tmpmol[contacts[j].polymer()][contacts[j].monomer()].seqnum() != mm.seqnum())
                        )   &&  (clipper::Coord_orth::length ( tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()].coord_orth(), candidates[i].coord_orth() ) <= 2.0 )
                            &&  (contacts[j].symmetry() == 0) )  // Beware: will report contacts that are not physically in contact, but needed for visualisation
                    //  if (   (tmpmol[contacts[j].polymer()][contacts[j].monomer()].id().trim() != mm.id().trim())
                    //     &&  (clipper::Coord_orth::length ( tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()].coord_orth(), candidates[i].coord_orth() ) <= 2.5 ))  // Beware: will report contacts that are not physically in contact, but needed for visualisation
                    {                            //         of crappy structures in MG
                        if ( altconf_compatible( get_altconf(tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()]), get_altconf(candidates[i]) ) )
                        {
                            clipper::MAtom tmpAtom = tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()];

                            if(debug_output)
                                DBG << "tmpAtom.id().trim() = " << tmpAtom.id().trim() << "\tfrom: " << tmpmol[contacts[j].polymer()].id().trim() << "/" << tmpmol[contacts[j].polymer()][contacts[j].monomer()].id().trim() << "-" << tmpmol[contacts[j].polymer()][contacts[j].monomer()].type().trim() << std::endl;

                            if(get_altconf(tmpAtom) != ' ')
                            {

                                if(tmpAtom.element().trim() == "O" || tmpAtom.element().trim() == "S" || tmpAtom.element().trim() == "N" || tmpAtom.element().trim() == "F" || mm.type().trim() == "TRP")
                                {
                                    std::pair < clipper::MAtom , clipper::MAtomIndexSymmetry > link_tmp;
                                    link_tmp.first = candidates[i];
                                    link_tmp.second = contacts[j];
                                    tmpresults.push_back ( link_tmp );
                                    managedToFindContactInSameAltConf = true;
                                    if(debug_output)
                                    {
                                        DBG << "Contact detected between atoms of altconf = " << get_altconf(tmpAtom) << std::endl;
                                        DBG << "link_tmp.first(MAtom).id() = " << candidates[i].id() << "\tlink_tmp.second(MAtomNonBond) aka tmpAtom = " << tmpAtom.id() << std::endl;
                                    }
                                }
                            }
                            else
                            {
                                if(tmpAtom.element().trim() == "O" || tmpAtom.element().trim() == "S" || tmpAtom.element().trim() == "N" || tmpAtom.element().trim() == "F" || mm.type().trim() == "TRP")
                                {
                                    std::pair < clipper::MAtom , clipper::MAtomIndexSymmetry > link_tmp;
                                    link_tmp.first = candidates[i];
                                    link_tmp.second = contacts[j];
                                    tmpresults.push_back ( link_tmp );
                                    managedToFindContactInSameAltConf = true;
                                    if(debug_output)
                                    {
                                        DBG << "link_tmp.first(MAtom).id() = " << candidates[i].id() << "\tlink_tmp.second(MAtomNonBond) aka tmpAtom = " << tmpAtom.id() << std::endl;
                                    }
                                }
                            }
                        }
                        else if( !altconf_compatible( get_altconf(tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()]), get_altconf(candidates[i]) ) )
                        {
                            altConfDetected = true;
                        }
                    }
                }

                if (altConfDetected && !managedToFindContactInSameAltConf)
                {
                    for (int j = 0 ; j < contacts.size() ; j++ )
                    {
                        if ( (  (tmpmol[contacts[j].polymer()].id().trim() != monomer_chain_id.trim())
                            ||  (tmpmol[contacts[j].polymer()][contacts[j].monomer()].id().trim() != mm.id().trim())
                            ||  (tmpmol[contacts[j].polymer()][contacts[j].monomer()].type().trim() != mm.type().trim())
                            ||  (tmpmol[contacts[j].polymer()][contacts[j].monomer()].seqnum() != mm.seqnum())
                            )  &&  (clipper::Coord_orth::length ( tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()].coord_orth(), candidates[i].coord_orth() ) <= 2.0 ))  // Beware: will report contacts that are not physically in contact, but needed for visualisation
                        //  if (   (tmpmol[contacts[j].polymer()][contacts[j].monomer()].id().trim() != mm.id().trim())
                        //     &&  (clipper::Coord_orth::length ( tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()].coord_orth(), candidates[i].coord_orth() ) <= 2.5 ))  // Beware: will report contacts that are not physically in contact, but needed for visualisation
                        {                            //         of crappy structures in MG

                            clipper::MAtom tmpAtom = tmpmol[contacts[j].polymer()][contacts[j].monomer()][contacts[j].atom()];
                            if(prefered_altconf == ' ' && get_altconf(tmpAtom) != ' ')
                                prefered_altconf = get_altconf(tmpAtom);

                            if(debug_output)
                                DBG << "tmpAtom.id().trim(), but this time with altconfs ayyy lmao = " << tmpAtom.id().trim() << "\tfrom: " << tmpmol[contacts[j].polymer()].id().trim() << "/" << tmpmol[contacts[j].polymer()][contacts[j].monomer()].id().trim() << "-" << tmpmol[contacts[j].polymer()][contacts[j].monomer()].type().trim() << "\tcurrent prefered_altConf = " << prefered_altconf << std::endl;

                            if(get_altconf(tmpAtom) == prefered_altconf && get_altconf(tmpAtom) != ' ')
                            {
                                const char altconf = get_altconf(tmpAtom);

                                if(tmpAtom.element().trim() == "O" || tmpAtom.element().trim() == "S" || tmpAtom.element().trim() == "N" || tmpAtom.element().trim() == "F")
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
                            else if(get_altconf(tmpAtom) == ' ')
                            {
                                if(tmpAtom.element().trim() == "O" || tmpAtom.element().trim() == "S" || tmpAtom.element().trim() == "N" || tmpAtom.element().trim() == "F")
                                {
                                    std::pair < clipper::MAtom , clipper::MAtomIndexSymmetry > link_tmp;
                                    link_tmp.first = candidates[i];
                                    link_tmp.second = contacts[j];
                                    tmpresults.push_back ( link_tmp );
                                    managedToFindContactInSameAltConf = true;
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