/*! \file clipper-glyco.h
  Header file for handling sugar data */

// clipper-glyco.h: a set of tools for handling sugars
// version  0.6.0
// 2013 Jon Agirre & Kevin Cowtan, The University of York
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk
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

#ifndef CLIPPER_GLYCO_H_INCLUDED
#define CLIPPER_GLYCO_H_INCLUDED

#include <stdlib.h>
#include <clipper/clipper.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-minimol.h>
#include "clipper-glyco_data.h"


namespace clipper
{

    //! MiniMol Sugar object for handling cyclic carbohydrate data
    /*! The MiniMol Sugar object is a derivation of clipper::MMonomer,
        and holds information specific to 5- and 6-membered cyclic sugars

        Upon creation, the input MMonomer is subjected to a sanity check
        against a database built into this library, which contains information
        like anomer, handedness, complete name and a list of the atoms integrating
        the sugar ring. Cremer-Pople parameters and conformation codes are provided
        for analysis and validation purposes. Information about connectivity is
        readily available right after creation too.

        It is strongly suggested that an MAtomNonBond object is supplied to the
        constructor, in order to avoid the creation (time expensive) of an object
        for each MSugar object.
    */
 
    class MSugar : public clipper::MMonomer
    {

        public:

            MSugar ();
            MSugar ( const clipper::MiniMol& mmol, const clipper::MMonomer& mmon ); //!< default constructor
            MSugar ( const clipper::MiniMol& mmol, const clipper::MMonomer& mmon, const clipper::MAtomNonBond& manb ); 
            //!< provide pre-calculated (time expensive) MAtomNonBond object. This object will tipically be re-used for many MSugar objects

	    clipper::String conformation_name() { return clipper::data::conformational_landscape[sugar_conformation]; } 
            //!< get a fixed three-character code for the conformation
			
            clipper::String conformation_name_iupac() { return clipper::data::iupac_conformational_landscape[sugar_conformation]; } 
            //!< get HTML-formatted, iupac-compliant codes describing the conformation of the sugar ring
			
            clipper::ftype puckering_amplitude() { return sugar_cremer_pople_params[0]; } 
            //!< convenience function for getting the puckering amplitude (in Angstroems) of the sugar ring
	    
            int conformation_code() { return sugar_conformation; } 
            //!< function for working with a numerical denomination of sugar conformation
	    
            const clipper::MAtom& anomeric_carbon() { return this->sugar_anomeric_carbon; } 
            //!< returns the anomeric carbon
		    
            const clipper::MAtom& anomeric_substituent() { return this->sugar_anomeric_substituent; } 
            //!< returns the anomeric substituent
			
            const clipper::MAtom& configurational_carbon() { return this->sugar_configurational_carbon; } 
            //!< returns the configurational carbon
	    
            const clipper::MAtom& configurational_substituent() { return this->sugar_configurational_substituent; } 
            //!< returns the configurational substituent
	    
            clipper::String anomer() { return sugar_anomer; } 
            //!< convenience function to get the anomer type, e.g. "alpha"
	    
            clipper::String handedness() { return sugar_handedness; } 
            //!< get the sugar's handedness: "D" for dextro, "L" for laevo, "N" for undetermined and "X" for unsupported (missing from the clipper::data sugar database)
	    
            clipper::String type_of_sugar() { return sugar_denomination; } 
            //!< returns a string containing an anomer-handedness-type denomination, e.g. "beta-D-aldopyranose"
	    
            int ring_cardinality() { return sugar_ring_elements.size(); } 
            //!< get the number of atoms forming the main sugar ring
			
            std::vector<clipper::MAtom> ring_members() { return sugar_ring_elements; } 
            //!< returns an standard vector containing those MAtom's forming the ring
			
            const clipper::Coord_orth& ring_centre() { return this->sugar_centre; } 
            //!< get the ring's centre (original coordinates)
			
            const clipper::Vec3<ftype>& ring_mean_plane() { return sugar_mean_plane; } 
            //!< get the vector normal to the sugar ring's mean plane, with origin in ring_centre()
	    
            std::vector<ftype> cremer_pople_params() { return sugar_cremer_pople_params; } 
            //!< returns Cremer-Pople parameters (Cremer and Pople, 1975) in the form { Q , Angle1, ... }
	    
            int potential_linkages() { return sugar_linked_to.size(); } 
            //!<  returns the number of potential linkages to other sugars
	    
            bool is_sane() { return sugar_sane; }
            //!< checks it against the internal clipper::data sugar database for correct anomer, handedness and ring members

            clipper::String full_name() { return sugar_name; }
			
            clipper::String full_type() { return sugar_type; }
	    
            std::vector<clipper::ftype> ring_angles() { return sugar_ring_angles; }
			
            std::vector<clipper::ftype> ring_bonds() { return sugar_ring_bonds; }
			
            std::vector<clipper::ftype> ring_torsions() { return sugar_ring_torsion; }
			
            clipper::ftype ring_bond_rmsd() { return sugar_ring_bond_rmsd; }
			
            clipper::ftype ring_angle_rmsd() { return sugar_ring_angle_rmsd; }
                        
            clipper::ftype get_bfactor() { return sugar_bfactor; }
                        
            bool in_database(clipper::String name) { return sugar_found_db; }
			
            bool is_supported() { return sugar_supported; }

	    static bool search_database(clipper::String name) 
            {
	        for (int i = 0; i < clipper::data::sugar_database_size ; i++) 
		    if (name.trim() == clipper::data::sugar_database[i].name_short.trim()) 
		        return true; 
					
		return false;
	    } //!< returns true if found
			
            bool ok_with_ring() { return sugar_diag_ring; }
			
            bool ok_with_bonds_rmsd() { return sugar_diag_bonds_rmsd; }  // compare with database, if not found report dev from ideal  
			
            bool ok_with_angles_rmsd() { return sugar_diag_angles_rmsd; } // same
			
            bool ok_with_anomer() { return sugar_diag_anomer; }
			
            bool ok_with_chirality() { return sugar_diag_chirality; }
            
            bool ok_with_conformation() { return sugar_diag_conformation; }

            bool ok_with_puckering() { return sugar_diag_puckering; }

            clipper::ftype get_rscc() { return sugar_rscc; }
                        
            void set_rscc ( clipper::ftype rscc_in ) { sugar_rscc = rscc_in; }
                        
            clipper::String get_diagnostic() { return sugar_diagnostic; }
                        
            void set_diagnostic ( clipper::String message ) { sugar_diagnostic = message; }                       


	private:

	    typedef std::vector< std::pair< clipper::MAtom, clipper::MAtom > > visited_arcs;
	    typedef std::pair< std::pair<clipper::MAtom, clipper::MAtom>, std::pair<clipper::MAtom, clipper::MAtom > > stereochemistry_pairs;

	    // int codes for conformations, types and anomers

	    static const int conf_pyranose_4C1 = 1;  static const int conf_pyranose_1C4 = 2;  static const int conf_pyranose_3OB = 3;  static const int conf_pyranose_B25 = 4;  
	    static const int conf_pyranose_14B = 5;  static const int conf_pyranose_B3O = 6;  static const int conf_pyranose_25B = 7;  static const int conf_pyranose_B14 = 8;  
	    static const int conf_pyranose_OE  = 9;  static const int conf_pyranose_E5 = 10;  static const int conf_pyranose_4E = 11;  static const int conf_pyranose_E3 = 12;
	    static const int conf_pyranose_2E = 13;  static const int conf_pyranose_E1 = 14;  static const int conf_pyranose_3E = 15;  static const int conf_pyranose_E2 = 16;  
	    static const int conf_pyranose_1E = 17;  static const int conf_pyranose_EO = 18;  static const int conf_pyranose_5E = 19;  static const int conf_pyranose_E4 = 20;  
	    static const int conf_pyranose_OH5 =21;  static const int conf_pyranose_4H5 = 22; static const int conf_pyranose_4H3 = 23; static const int conf_pyranose_2H3 = 24;
	    static const int conf_pyranose_2H1 =25;  static const int conf_pyranose_OH1 = 26; static const int conf_pyranose_3H2 = 27; static const int conf_pyranose_1H2 = 28; 
	    static const int conf_pyranose_1HO =29;  static const int conf_pyranose_5HO = 30; static const int conf_pyranose_5H4 = 31; static const int conf_pyranose_3H4 = 32; 
	    static const int conf_pyranose_OS2 = 33; static const int conf_pyranose_1S5 = 34; static const int conf_pyranose_1S3 = 35; static const int conf_pyranose_2SO = 36;
	    static const int conf_pyranose_5S1 = 37; static const int conf_pyranose_3S1 = 38; 
			
	    static const int conf_furanose_3T2 = 39; static const int conf_furanose_3EV = 40; static const int conf_furanose_3T4 = 41; static const int conf_furanose_EV4 = 42;  
	    static const int conf_furanose_OT4 = 43; static const int conf_furanose_OEV = 44; static const int conf_furanose_OT1 = 45; static const int conf_furanose_EV1 = 46;  
	    static const int conf_furanose_2T1 = 47; static const int conf_furanose_2EV = 48; static const int conf_furanose_2T3 = 49; static const int conf_furanose_EV3 = 50;  
	    static const int conf_furanose_4T3 = 51; static const int conf_furanose_4EV = 52; static const int conf_furanose_4TO = 53; static const int conf_furanose_EVO = 54;  
            static const int conf_furanose_1TO = 55; static const int conf_furanose_1EV = 56; static const int conf_furanose_1T2 = 57; static const int conf_furanose_EV2 = 58;  
			
	    static const int db_not_checked = 9999;  static const int db_not_found = 101010;

	    static const int anomer_alpha = 100;     static const int anomer_beta = 101; 	  static const int handedness_d = 200;     static const int handedness_l = 201;     
	    static const int type_aldose = 1000;     static const int type_ketose = 2000;


	    // We'll use the sugar_ prefix throughout the class for private members

	    const MiniMol*                  sugar_parent_molecule;
	    const MAtomNonBond*             sugar_parent_molecule_nonbond;
	    Coord_orth                      sugar_centre;
	    clipper::Vec3<clipper::ftype>   sugar_mean_plane;
	    std::vector<MAtom>              sugar_ring_elements;
	    std::vector<ftype>              sugar_cremer_pople_params;          // (1) puckering amplitude, (2..) angles
	    clipper::MAtom                  sugar_anomeric_carbon;
	    clipper::MAtom                  sugar_configurational_carbon;
	    clipper::MAtom                  sugar_anomeric_substituent;
	    clipper::MAtom                  sugar_configurational_substituent;
	    bool                            sugar_found_db;                     // true if the sugar's code is present in the reference data structure
	    bool                            sugar_sane;                         // true if passed sanity checks
	    bool                            sugar_supported;                    // false if the sugar has no connectivity, missing atoms, etc.
	    std::vector<clipper::ftype>     sugar_ring_bonds;                   // bond lengths among ring. [0] is O to anomeric carbon
	    std::vector<clipper::ftype>     sugar_ring_angles;                  // ring angles, starting with last_carbon-O-anomeric carbon
	    std::vector<clipper::ftype>     sugar_ring_torsion;                 // torsion angles, starting with C5-O5-C1-C2
	    int                             sugar_index;                        // 9999 if unchecked, 101010 if absent from the database, 0-400 if found
	    String                          sugar_denomination;                 // e.g. alpha-aldopyranose
	    String                          sugar_anomer;                       // "alpha", "beta" or "undetermined"
	    String                          sugar_handedness;                   // "D", "L" or "undetermined"
	    String                          sugar_type;                         // i.e. "aldose" or "ketose"
	    String                          sugar_name;
	    bool                            sugar_diag_ring;
	    bool                            sugar_diag_bonds_rmsd;
	    bool                            sugar_diag_angles_rmsd;
	    bool                            sugar_diag_anomer;
	    bool                            sugar_diag_chirality;
            bool                            sugar_diag_conformation;
            bool                            sugar_diag_puckering;
            clipper::String                 sugar_diagnostic;                   // full diagnostic to be used in Coot and ccp4i2
            clipper::ftype                  sugar_rscc;                         // RSCC to be used in Coot and ccp4i2
	    clipper::ftype                  sugar_ring_bond_rmsd;
	    clipper::ftype                  sugar_ring_angle_rmsd;
            clipper::ftype                  sugar_bfactor;
            int                             sugar_conformation;
	    clipper::String                 sugar_alternate_confcode;
	    std::vector < MSugar >          sugar_linked_to;		        // size: number of carbon atoms - 1 ;

	    // private methods

	    std::vector<clipper::ftype>                     cremerPople_pyranose ( const clipper::MiniMol& , clipper::MMonomer); // modifies object
	    std::vector<clipper::ftype>                     cremerPople_furanose ( const clipper::MiniMol& , clipper::MMonomer); // modifies object
	    std::vector<clipper::MAtom&>                    find_bonded ( const clipper::MAtom& ) const; // accesses object
	    int                                             conformationPyranose ( const clipper::ftype&, const clipper::ftype& ) const;
	    int                                             conformationFuranose ( const clipper::ftype& ) const;

	    std::vector<clipper::MAtom> 		    ringMembers ( ) const; // internal function, accesses object
	    const std::vector<clipper::MAtom> 		    findPath ( const clipper::MMonomer&, int, MSugar::visited_arcs& ) const;
	    bool 					    closes_ring ( const clipper::MAtom&, MSugar::visited_arcs& ) const;
	    const std::vector<clipper::MAtom>               findBonded ( const clipper::MAtom&, MSugar::visited_arcs& ) const; // accesses object
	    bool                                            lookup_visited ( const MSugar::visited_arcs&, const std::pair<clipper::MAtom, clipper::MAtom> ) const;
	    stereochemistry_pairs                           get_stereochemistry ( const clipper::MiniMol& ); // accesses object
	    bool                                            is_stereocentre ( const clipper::MAtom&, const clipper::MiniMol& ); // accesses object
	    bool                                            is_part_of_ring ( const clipper::MAtom&, const std::vector<clipper::MAtom> ) const;
	    bool                                            bonded ( const clipper::MAtomIndexSymmetry&, const clipper::MAtom& ) const;
	    bool                                            bonded ( const clipper::MAtom&, const clipper::MAtom& ) const; 
	    bool                                            lookup_database ( clipper::String );
	    bool                                            examine_ring();
    
            const char get_altconf ( const clipper::MAtom& ) const;
	    	
    }; // class MSugar


    class MGlycan
    {
        public:

            MGlycan () { } //!< null constructor
            MGlycan ( clipper::String chain, const clipper::MMonomer& root_aa, clipper::MSugar& root_sugar );
            
            bool link_sugars ( int link, clipper::MSugar& first_sugar, clipper::MSugar& next_sugar ); // true if there's been any problem

            clipper::MMonomer get_root () { return root.first; }
            clipper::String get_type () { return kind_of_glycan; } // n-glycan or o-glycan
            clipper::String print_linear ( const bool print_info, const bool html_format, const bool translate );
            clipper::String print_SVG ( bool vertical, bool print_info, bool colour_gradient );
            std::vector < clipper::MSugar > get_sugars () { return sugars; }
            clipper::String get_chain () { return chain; }
            void set_kind_of_glycan ( clipper::String input ) { kind_of_glycan = input; }
            clipper::String get_kind_of_glycan ( ) { return kind_of_glycan; }

        private:
            
            struct t_connection;

            typedef struct t_node
            {
                std::vector < t_connection > connections;
                clipper::MSugar sugar;
            } t_node;

            typedef struct t_connection
            {
                int index;            // carbon to which this is connected 
                clipper::String type; // anomer
                t_node* node;         // sugar connected to the present
            } t_connection;

            // private stuff
            std::vector < t_node > node_list; 
            clipper::String kind_of_glycan;                 // can be n-glycan or o-glycan
            std::pair < clipper::MMonomer, t_node > root;   // this is the root, should be null if this is a ligand saccharide
            std::vector < clipper::MSugar > sugars;         // vector of sugars    
            clipper::String chain;     
            clipper::String carbname_of(clipper::String name) const;
    }; // class MGlycan

    class MGlycology
    {
        public: 

            MGlycology () { } //!< null constructor
            MGlycology ( const clipper::MiniMol& );
            MGlycology ( const clipper::MiniMol&, const clipper::MAtomNonBond& );
            
            clipper::MGlycan get_glycan_by_id ( int id );
            clipper::MGlycan get_glycan_by_root ( clipper::MMonomer& root ) 
                { for (int i=0;i<list_of_glycans.size();i++) if (list_of_glycans[i].get_root().id()==root.id()) return list_of_glycans[i]; }
            std::vector < clipper::MGlycan > get_list_of_glycans () { return list_of_glycans; }
            std::vector < clipper::MSugar > get_sugar_list() { return list_of_sugars; } 

        private:
            
            // private properties
            std::vector < clipper::MSugar >  list_of_sugars; 
            std::vector < clipper::MGlycan > list_of_glycans;
            const clipper::MAtomNonBond* manb;
            const clipper::MiniMol * mmol;
            const std::vector < std::pair< clipper::MAtom, clipper::MAtomIndexSymmetry > > get_contacts ( const clipper::MMonomer& mm_one );
            int parse_order ( clipper::String str ) { const char *s = str.c_str(); return atoi(&s[2]); }
            void extend_tree ( clipper::MGlycan& mg, clipper::MSugar& msug );   
            const char get_altconf ( const clipper::MAtom& ) const;


            // private methods
    
    }; // class MGlycology

} // namespace clipper

#endif
