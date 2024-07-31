/*! \file clipper-glyco.h
  Header file for handling sugar data */

// clipper-glyco.h: a set of tools for handling sugars
// version  0.9.1
// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York


#ifndef CLIPPER_GLYCO_H_INCLUDED
#define CLIPPER_GLYCO_H_INCLUDED

#include <fstream>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <regex>
#include <clipper/clipper.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-minimol.h>
#include "clipper-glyco_data.h"
#include <clipper/minimol/minimol_utils.h>
#include "privateer-json.h"

inline bool altconf_compatible ( char m1, char m2 )
{
    // OLD VERSION:
        // if (( m1 == 'A' && m2 == 'B') || ( m1 == 'B' && m2 == 'A'))
        //     return false;
        // else
        //     return true;
    // would cause an endless while loop in is_stereocentre() when there would be altconfs like 'C'
    if(( m1 == m2))
        return true;
    else if((m1 == ' ' && m2 != ' ') || (m1 != ' ' && m2 == ' '))
        return true;
    else
        return false;

}
//!< false if one is from conf A and the other is from conf B


inline clipper::Vec3<clipper::ftype> find_aromatic_plane ( clipper::MMonomer mmon )
{
    clipper::Vec3<clipper::ftype> result (0.0, 0.0, 0.0);

    if ( mmon.type() == "TRP" )
    {
        bool foundAtomAlpha = false, foundAtomBravo = false, foundAtomCharlie = false;
        for(int atom = 0; atom < mmon.size(); atom++ )
        {
            clipper::MAtom currentAtom = mmon[atom];
            if(currentAtom.id().trim() == "CD1")
                foundAtomAlpha = true;
            if(currentAtom.id().trim() == "CD2")
                foundAtomBravo = true;
            if(currentAtom.id().trim() == "CE2")
                foundAtomCharlie = true;
        }

        if(foundAtomAlpha == false || foundAtomBravo == false || foundAtomCharlie == false)
            return result;
        clipper::Vec3<clipper::ftype> vec1 ( mmon.find ( "CD1", clipper::MM::ANY ).coord_orth().x() - mmon.find ( "CD2", clipper::MM::ANY ).coord_orth().x(),
                                            mmon.find ( "CD1", clipper::MM::ANY ).coord_orth().y() - mmon.find ( "CD2", clipper::MM::ANY ).coord_orth().y(),
                                            mmon.find ( "CD1", clipper::MM::ANY ).coord_orth().z() - mmon.find ( "CD2", clipper::MM::ANY ).coord_orth().z());

        clipper::Vec3<clipper::ftype> vec2 ( mmon.find ( "CE2", clipper::MM::ANY ).coord_orth().x() - mmon.find ( "CD2", clipper::MM::ANY ).coord_orth().x(),
                                            mmon.find ( "CE2", clipper::MM::ANY ).coord_orth().y() - mmon.find ( "CD2", clipper::MM::ANY ).coord_orth().y(),
                                            mmon.find ( "CE2", clipper::MM::ANY ).coord_orth().z() - mmon.find ( "CD2", clipper::MM::ANY ).coord_orth().z());
        result = clipper::Vec3<clipper::ftype>::cross ( vec1, vec2 );
        return result.unit();

    }
    else if ( mmon.type() == "TYR" )
    {
        bool foundAtomAlpha = false, foundAtomBravo = false, foundAtomCharlie = false;
        for(int atom = 0; atom < mmon.size(); atom++ )
        {
            clipper::MAtom currentAtom = mmon[atom];
            if(currentAtom.id().trim() == "CE1")
                foundAtomAlpha = true;
            if(currentAtom.id().trim() == "CG")
                foundAtomBravo = true;
            if(currentAtom.id().trim() == "CE2")
                foundAtomCharlie = true;
        }

        if(foundAtomAlpha == false || foundAtomBravo == false || foundAtomCharlie == false)
            return result;
        clipper::Vec3<clipper::ftype> vec2 ( mmon.find ( "CE1", clipper::MM::ANY ).coord_orth().x() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().x(),
                                            mmon.find ( "CE1", clipper::MM::ANY ).coord_orth().y() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().y(),
                                            mmon.find ( "CE1", clipper::MM::ANY ).coord_orth().z() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().z());

        clipper::Vec3<clipper::ftype> vec1 ( mmon.find ( "CE2", clipper::MM::ANY ).coord_orth().x() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().x(),
                                            mmon.find ( "CE2", clipper::MM::ANY ).coord_orth().y() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().y(),
                                            mmon.find ( "CE2", clipper::MM::ANY ).coord_orth().z() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().z());
        result = clipper::Vec3<clipper::ftype>::cross ( vec1, vec2 );
        return result.unit();
    }
    else if ( mmon.type() == "PHE" )
    {
        bool foundAtomAlpha = false, foundAtomBravo = false, foundAtomCharlie = false;
        for(int atom = 0; atom < mmon.size(); atom++ )
        {
            clipper::MAtom currentAtom = mmon[atom];
            if(currentAtom.id().trim() == "CE1")
                foundAtomAlpha = true;
            if(currentAtom.id().trim() == "CG")
                foundAtomBravo = true;
            if(currentAtom.id().trim() == "CE2")
                foundAtomCharlie = true;
        }

        if(foundAtomAlpha == false || foundAtomBravo == false || foundAtomCharlie == false)
            return result;
        clipper::Vec3<clipper::ftype> vec2 ( mmon.find ( "CE1", clipper::MM::ANY ).coord_orth().x() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().x(),
                                            mmon.find ( "CE1", clipper::MM::ANY ).coord_orth().y() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().y(),
                                            mmon.find ( "CE1", clipper::MM::ANY ).coord_orth().z() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().z());

        clipper::Vec3<clipper::ftype> vec1 ( mmon.find ( "CE2", clipper::MM::ANY ).coord_orth().x() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().x(),
                                            mmon.find ( "CE2", clipper::MM::ANY ).coord_orth().y() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().y(),
                                            mmon.find ( "CE2", clipper::MM::ANY ).coord_orth().z() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().z());
        result = clipper::Vec3<clipper::ftype>::cross ( vec1, vec2 );
        return result.unit();
    }
    else if ( mmon.type() == "HIS" )
    {
        bool foundAtomAlpha = false, foundAtomBravo = false, foundAtomCharlie = false;
        for(int atom = 0; atom < mmon.size(); atom++ )
        {
            clipper::MAtom currentAtom = mmon[atom];
            if(currentAtom.id().trim() == "CE1")
                foundAtomAlpha = true;
            if(currentAtom.id().trim() == "CG")
                foundAtomBravo = true;
            if(currentAtom.id().trim() == "NE2")
                foundAtomCharlie = true;
        }

        if(foundAtomAlpha == false || foundAtomBravo == false || foundAtomCharlie == false)
            return result;
        clipper::Vec3<clipper::ftype> vec2 ( mmon.find ( "CE1", clipper::MM::ANY ).coord_orth().x() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().x(),
                                            mmon.find ( "CE1", clipper::MM::ANY ).coord_orth().y() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().y(),
                                            mmon.find ( "CE1", clipper::MM::ANY ).coord_orth().z() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().z());

        clipper::Vec3<clipper::ftype> vec1 ( mmon.find ( "NE2", clipper::MM::ANY ).coord_orth().x() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().x(),
                                            mmon.find ( "NE2", clipper::MM::ANY ).coord_orth().y() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().y(),
                                            mmon.find ( "NE2", clipper::MM::ANY ).coord_orth().z() - mmon.find ( "CG", clipper::MM::ANY ).coord_orth().z());
        result = clipper::Vec3<clipper::ftype>::cross ( vec1, vec2 );
        return result.unit();
    }
    return result;
}


inline clipper::ftype get_angle ( clipper::Vec3<clipper::ftype> vec1, clipper::Vec3<clipper::ftype> vec2 )
{
    clipper::ftype angle = acos ( clipper::Vec3<clipper::ftype>::dot ( vec1, vec2 ) /
                                  (sqrt ( pow(vec1[0],2) + pow(vec1[1],2) + pow(vec1[2],2)) *
                                   sqrt ( pow(vec2[0],2) + pow(vec2[1],2) + pow(vec2[2],2))) ) ;

    return angle;
}

//! Clipper extensions for handling sugars. All code within this namespace shall respect Clipper conventions
namespace clipper
{

    //! A class for handling monosaccharides in the pyranose or furanose forms
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
            MSugar ( const clipper::MiniMol& mmol, const clipper::String chainID, const clipper::MMonomer& mmon, bool& debug_output, char alt_conf = ' ' ); //!< default constructor
            MSugar ( const clipper::MiniMol& mmol, const clipper::String chainID, const clipper::MMonomer& mmon, const clipper::MAtomNonBond& manb, bool& debug_output, char alt_conf = ' ' );
            MSugar ( const clipper::MiniMol& mmol, const clipper::String chainID, const clipper::MMonomer& mmon, const clipper::MAtomNonBond& manb, clipper::data::sugar_database_entry& validation_data, bool& debug_output, char alt_conf = ' '  );
            MSugar ( const clipper::MiniMol& mmol, const clipper::MMonomer& mmon, const clipper::MAtomNonBond& manb, char alt_conf = ' ' );
            //!< provide pre-calculated (time expensive) MAtomNonBond object. This object will tipically be re-used for many MSugar objects

            class Diagnostics
            {
                public:

                    struct sugar_diag_puckering_diagnostics_template
                    {
                        bool sugar_found_db;
                        float range_min;
                        float range_max;
                        std::string comparison_operator;
                        bool sugar_diag_conformation;
                        float reference_puckering;
                        float calculated_puckering_amplitude;
                        bool initialized = false;
                        bool final_result;
                    };

                    struct sugar_diag_anomer_diagnostics_template
                    {
                        bool sugar_found_db;
                        std::string sugar_anomer;
                        std::string database_sugar_anomer;
                        bool initialized = false;
                        bool final_result;
                    };

                    struct sugar_diag_chirality_diagnostics_template
                    {
                        bool sugar_found_db;
                        std::string sugar_handedness;
                        std::string database_sugar_handedness;
                        bool initialized = false;
                        bool final_result;
                    };

                    struct sugar_diag_ring_diagnostics_atom_pair
                    {
                        std::string ma_one_atom;
                        std::string ma_two_atom;
                        std::string ma_one_element;
                        std::string ma_two_element;
                        float distance_min;
                        float distance_max;
                        float measured_distance;
                        bool atom_result;
                    };

                    struct sugar_diag_ring_diagnostics_template
                    {
                        bool initialized = true;
                        bool final_result;
                        std::vector<sugar_diag_ring_diagnostics_atom_pair> ring_atom_diagnostic;
                    };

                    Diagnostics () { this->diagnostic_complete = false; };
                    Diagnostics ( bool sugar_diagnostic_complete ) { this->diagnostic_complete = sugar_diagnostic_complete; };
                    void update_diagnostic_status ( bool sugar_diagnostic_complete ) { this->diagnostic_complete = sugar_diagnostic_complete; }
                    void update_sugar_sane_status ( bool sugar_sane_status ) { this->sugar_sane = sugar_sane_status; }


                    sugar_diag_puckering_diagnostics_template get_puckering_diagnostics() { return sugar_diag_puckering_diagnostics; };
                    sugar_diag_anomer_diagnostics_template get_anomer_diagnostics() { return sugar_diag_anomer_diagnostics; };
                    sugar_diag_chirality_diagnostics_template get_chirality_diagnostics() { return sugar_diag_chirality_diagnostics; };
                    sugar_diag_ring_diagnostics_template get_ring_diagnostics() { return sugar_diag_ring_diagnostics; };

                    sugar_diag_puckering_diagnostics_template sugar_diag_puckering_diagnostics;
                    sugar_diag_anomer_diagnostics_template sugar_diag_anomer_diagnostics;
                    sugar_diag_chirality_diagnostics_template sugar_diag_chirality_diagnostics;
                    sugar_diag_ring_diagnostics_template sugar_diag_ring_diagnostics;
                private:
                    bool diagnostic_complete;
                    bool sugar_sane;

            };

            const bool operator== ( const clipper::MSugar& m2 ) const { return ( this->id() == m2.id() ); }

            const clipper::String conformation_name() const { return clipper::data::conformational_landscape[sugar_conformation]; }
            //!< get a fixed three-character code for the conformation

            const clipper::String conformation_name_iupac() const { return clipper::data::iupac_conformational_landscape[sugar_conformation]; }
            //!< get HTML-formatted, iupac-compliant codes describing the conformation of the sugar ring

            const clipper::ftype puckering_amplitude() const
            {
                if(sugar_cremer_pople_params.empty())
                    return -1;
                else
                    return sugar_cremer_pople_params[0];
            }
            //!< convenience function for getting the puckering amplitude (in Angstroems) of the sugar ring

            int conformation_code() const { return sugar_conformation; }
            //!< function for working with a numerical denomination of sugar conformation

            const clipper::MAtom& anomeric_carbon() const { return this->sugar_anomeric_carbon; }


            const clipper::MAtom& anomeric_substituent() const { return this->sugar_anomeric_substituent; }
            //!< returns the anomeric substituent

            const clipper::MAtom& configurational_carbon() const { return this->sugar_configurational_carbon; }
            //!< returns the configurational carbon

            const clipper::MAtom& configurational_substituent() const { return this->sugar_configurational_substituent; }
            //!< returns the configurational substituent

            clipper::String anomer() const { return sugar_anomer; }
            //!< convenience function to get the anomer type, e.g. "alpha"

            clipper::String handedness() const { return sugar_handedness; }
            //!< get the sugar's handedness: "D" for dextro, "L" for laevo, "N" for undetermined and "X" for unsupported (missing from the clipper::data sugar database)

            clipper::String chain_id() const { return sugar_chain_id; };

            clipper::String type_of_sugar() const { return sugar_denomination; }
            //!< returns a string containing an anomer-handedness-type denomination, e.g. "beta-D-aldopyranose"

            int ring_cardinality() const { return sugar_ring_elements.size(); }
            //!< get the number of atoms forming the main sugar ring

            std::vector<clipper::MAtom> ring_members() const { return sugar_ring_elements; }
            //!< returns an standard vector containing those MAtom's forming the ring

            const clipper::Coord_orth& ring_centre() const { return this->sugar_centre; }
            //!< get the ring's centre (original coordinates)

            const clipper::Vec3<ftype>& ring_mean_plane() const { return sugar_mean_plane; }
            //!< get the vector normal to the sugar ring's mean plane, with origin in ring_centre()

            std::vector<ftype> cremer_pople_params() const { return sugar_cremer_pople_params; }
            //!< returns Cremer-Pople parameters (Cremer and Pople, 1975) in the form { Q , Angle1, ... }

            int potential_linkages() const { return sugar_linked_to.size(); }
            //!<  returns the number of potential linkages to other sugars

            // std::vector < std::pair< clipper::MAtomIndexSymmetry, float > > get_stacked_residues ( std::string = "hudson",
            //                                                                                              float = 4.5,
            //                                                                                              float = 40.0,
            //                                                                                              float = 0.0 ) const ;
            //!< returns chain and monomer for stacked residues (restricted to TRP, TYR, PHE, HIS)
            /*!
              \sa get_stacked_residues()
              \param algorithm Specify "hudson" (default) or "plevin".
              \param distance Distance between the carbon and the centre of the ring.
              \param theta An angle, different depending on the algorithm.
              \param phi An angle used in the Plevin calculation, has to be greater than 120 deg.
            */

            bool is_sane() const { return sugar_sane; }
            //!< checks it against the internal clipper::data sugar database for correct anomer, handedness and ring members

            clipper::String full_name() const { return sugar_name_full; }

            clipper::String short_name() const { return sugar_name_short; }

            clipper::String pdb_id() const { return sugar_pdb_id; }

            int get_seqnum() const { return sugar_seqnum; }

            clipper::String full_type() const { return sugar_type; }

            std::vector<clipper::ftype> ring_angles() const { return sugar_ring_angles; }

            std::vector<clipper::ftype> ring_bonds() const { return sugar_ring_bonds; }

            std::vector<clipper::ftype> ring_torsions() const { return sugar_ring_torsion; }

            clipper::ftype ring_bond_rmsd() const { return sugar_ring_bond_rmsd; }

            clipper::ftype ring_angle_rmsd() const { return sugar_ring_angle_rmsd; }

            clipper::ftype get_bfactor() const { return sugar_bfactor; }

            bool in_database(clipper::String name) const { return sugar_found_db; }

            bool is_supported() const { return sugar_supported; }

            static bool search_database(clipper::String name)
            {
                for (int i = 0; i < clipper::data::sugar_database_size ; i++)
                    if (name.trim() == clipper::data::sugar_database[i].name_short.trim())
                        return true;

                return false;
            } //!< returns true if found

            bool ok_with_ring() const { return sugar_diag_ring; }

            bool ok_with_bonds_rmsd() const { return sugar_diag_bonds_rmsd; }  // compare with database, if not found report dev from ideal

            bool ok_with_angles_rmsd() const { return sugar_diag_angles_rmsd; } // same

            bool ok_with_anomer() const { return sugar_diag_anomer; }

            bool ok_with_chirality() const { return sugar_diag_chirality; }

            bool ok_with_conformation() const { return sugar_diag_conformation; }

            bool ok_with_puckering() const { return sugar_diag_puckering; }

            clipper::ftype get_rscc() const { return sugar_rscc; }

            clipper::ftype get_accum() const { return sugar_accum; }

            bool get_occupancy_check() const { return sugar_occupancy_check; }

            void set_rscc ( clipper::ftype rscc_in ) { sugar_rscc = rscc_in; }

            void set_accum_score ( clipper::ftype accum_in ) { sugar_accum = accum_in; }

            void override_conformation_diag ( bool is_it_ok ) { sugar_diag_conformation = is_it_ok; }

            void override_anomer_diag ( bool is_it_ok ) { sugar_diag_anomer = is_it_ok; }

            void set_occupancy_check ( bool occupancy_check_in ) { sugar_occupancy_check = occupancy_check_in; }

            clipper::String get_diagnostic() const { return sugar_diagnostic; }

            void set_diagnostic ( clipper::String message ) { sugar_diagnostic = message; }

            clipper::String get_context ( ) const { return sugar_context; }

            void set_context ( clipper::String context ) { this->sugar_context = context; }

            Diagnostics get_detailed_diagnostics () const { return sugar_diagnostics; };

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

            static const int anomer_alpha = 100;     static const int anomer_beta = 101;          static const int handedness_d = 200;     static const int handedness_l = 201;
            static const int type_aldose = 1000;     static const int type_ketose = 2000;


            // We'll use the sugar_ prefix throughout the class for private members
            Diagnostics                     sugar_diagnostics;
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
            String                          sugar_chain_id;
            String                          sugar_denomination;                 // e.g. alpha-aldopyranose
            String                          sugar_anomer;                       // "alpha", "beta" or "undetermined"
            String                          sugar_handedness;                   // "D", "L" or "undetermined"
            String                          sugar_type;                         // i.e. "aldose" or "ketose"
            String                          sugar_name_full;
            String                          sugar_name_short;
            String                          sugar_pdb_id;
            int                             sugar_seqnum;
            bool                            sugar_diag_ring;
            bool                            sugar_diag_bonds_rmsd;
            bool                            sugar_diag_angles_rmsd;
            bool                            sugar_diag_anomer;
            bool                            sugar_diag_chirality;
            bool                            sugar_diag_conformation;
            bool                            sugar_diag_puckering;
            bool                            sugar_occupancy_check;
            clipper::String                 sugar_diagnostic;               // full diagnostic to be used in Coot and ccp4i2
            clipper::ftype                  sugar_rscc;                     // RSCC to be used in Coot and ccp4i2
            clipper::ftype                  sugar_accum;                     // accum score, <Fo> in Cryo-EM, <mFo> in Xray.
            clipper::ftype                  sugar_ring_bond_rmsd;
            clipper::ftype                  sugar_ring_angle_rmsd;
            clipper::ftype                  sugar_bfactor;
            int                             sugar_conformation;
            clipper::String                 sugar_alternate_confcode;
            std::vector < MSugar >          sugar_linked_to;                    // size: number of carbon atoms - 1
            clipper::String                 sugar_context;                  // n-glycan, o-glycan, s-glycan, c-glycan, p-glycan or ligand
            int                             glycan_index;
            bool                            debug_output;

            // private methods

            std::vector<clipper::ftype>                     cremerPople_pyranose ( const clipper::MiniMol& , clipper::MMonomer); // modifies object
            std::vector<clipper::ftype>                     cremerPople_furanose ( const clipper::MiniMol& , clipper::MMonomer); // modifies object
            std::vector<clipper::MAtom&>                    find_bonded ( const clipper::MAtom& ) const; // accesses object
            int                                             conformationPyranose ( const clipper::ftype&, const clipper::ftype& ) const;
            int                                             conformationFuranose ( const clipper::ftype& ) const;

            std::vector<clipper::MAtom>                     ringMembers ( ) const; // internal function, accesses object
            const std::vector<clipper::MAtom>               findPath ( const clipper::MMonomer&, int, MSugar::visited_arcs& ) const;
            bool                                            closes_ring ( const clipper::MAtom&, MSugar::visited_arcs& ) const;
            const std::vector<clipper::MAtom>               findBonded ( const clipper::MAtom&, MSugar::visited_arcs& ) const; // accesses object
            bool                                            lookup_visited ( const MSugar::visited_arcs&, const std::pair<clipper::MAtom, clipper::MAtom> ) const;
            stereochemistry_pairs                           get_stereochemistry ( const clipper::MiniMol& ); // accesses object
            bool                                            is_stereocentre ( const clipper::MAtom&, const clipper::MiniMol& ); // accesses object
            bool                                            is_part_of_ring ( const clipper::MAtom&, const std::vector<clipper::MAtom> ) const;
            bool                                            bonded ( const clipper::MAtomIndexSymmetry&, const clipper::MAtom& ) const;
            bool                                            bonded ( const clipper::MAtom&, const clipper::MAtom&, std::vector<clipper::MSugar::Diagnostics::sugar_diag_ring_diagnostics_atom_pair>& ring_atom_diagnostic ) const;
            bool                                            lookup_database ( clipper::String );
            bool                                            examine_ring();

            const char get_altconf ( const clipper::MAtom& ) const;

    }; // class MSugar

    //! A class that holds two MSugar objects
    class MDisaccharide
    {
        public:

            MDisaccharide () {} //!< null constructor
            MDisaccharide ( clipper::MSugar& sugar_one, clipper::MSugar& sugar_two, bool& debug_output) { this->sugar_one = sugar_one; this->sugar_two = sugar_two; this->debug_output = debug_output; }
            MDisaccharide ( clipper::MiniMol& mmol, const clipper::MAtomNonBond& manb, const clipper::String chainID, clipper::MMonomer& mm, bool& debug_output );
            clipper::MSugar& get_first_sugar () { return sugar_one; }
            clipper::MSugar& get_second_sugar () { return sugar_two; }

            static int search_disaccharides(clipper::String name)
            {
                for (int i = 0; i < clipper::data::disaccharide_database_size ; i++)
                    if (name.trim() == clipper::data::disaccharide_database[i].name_short.trim())
                        return i;

                return -1;
            } //!< returns -1 if not found

        private:
            bool debug_output;
            clipper::MSugar sugar_one;
            clipper::MSugar sugar_two;

    }; // class MDisaccharide

    //! Stores a tree built with Node(s) and Linkage(s), which contain references to MSugar(s)
    // TO DO: Make this class a prototype that would be inherited both by MGlycan and MLigandGlycan. Currently, I'm just gonna add hack'ish way to work with
    //        ligand glycans via the introduction of an alternative constructor.
    class MGlycan
    {
        public:

            MGlycan () { } //!< null constructor
            MGlycan ( clipper::String chain, clipper::MMonomer& root_aa, clipper::MSugar& root_sugar, clipper::String& root_sugar_chain_id, bool& debug_output, std::string expression_system = "undefined" );
            MGlycan ( clipper::String chain, clipper::MSugar& root_sugar, clipper::String& root_sugar_chain_id, bool& debug_output, std::string expression_system = "undefined" );

            struct MGlycanTorsionSummary
            {
                std::string type;
                clipper::String first_residue_name; // donorResidue
                clipper::String second_residue_name; // acceptorResidue
                std::vector<std::pair<clipper::MAtom, clipper::MAtom>> atoms; // .first = donorAtom, .second = acceptorAtom
                std::vector<std::pair<std::string, std::string>> linkage_descriptors; // .first = donorPosition, .second = acceptorPosition
                std::vector<std::pair<std::pair<std::string, std::string>, std::vector<std::pair<float,float>>>> combined_torsions; // This needs to be changed 
                std::vector<std::pair<float, float>> torsions; // .first = Phi, .second = Psi
            };

            class Node;

            class Linkage
            {
                public:

                    /* we need a copy constructor, as we are holding a pointer
                    Linkage ( const Linkage& linkage ) : node  ( linkage.node ),
                                                         index ( linkage.index ),
                                                         type  ( linkage.type ),
                                                         torsion_phi ( linkage.torsion_phi ),
                                                         torsion_psi ( linkage.torsion_psi ) { }
                    */
                    // Linkage () {};
                    Linkage ( int index, std::string anomericity, int connect_to_id, bool noncircular )
                    {
                        node_id = connect_to_id;
                        this->index = index;
                        this->noncircular = noncircular;
                        this->type = anomericity;
                        torsion_phi = torsion_psi = torsion_omega = torsion_omega_six = torsion_omega_seven = torsion_omega_eight = torsion_omega_nine = 0.0;
                        donorAtom = acceptorAtom = clipper::MAtom();
                        this->linkage_zscore_calculated = false;
                        this->linkage_zscore = -42069;
                    }

                    // Linkage& operator= ( const Linkage& link ) { if ( this != &link ) node = link.node; return *this;  }

                    int get_order () const { return index; }
                    void set_order ( int order ) { index = order; }

                    int get_linked_node_id ( ) const { return node_id; }
                    bool connection_is_non_circular ( ) const { return noncircular; };
                    void modify_linked_node_id ( int modified_connect_to_id ) { node_id = modified_connect_to_id; }

                    std::string get_anomericity ( ) const { return type; }
                    //!< alpha or beta
                    void set_anomericity ( std::string anomericity ) { type = anomericity; }
                    //!< alpha or beta

                    std::string get_description ( bool is_ketose = false )
                    {
                        std::ostringstream message;
                        std::vector< float > torsions = get_torsions();

                        if ( is_ketose ) // ketoses have C1 outside of the ring. C2 is the anomeric centre
                            message << "[ " << get_anomericity() << " 2-" << get_order()
                                    << " ] with phi: " << torsions[0] << "; psi: " << torsions[1];
                        else
                            message << "[ " << get_anomericity() << " 1-" << get_order()
                                    << " ] with phi: " << torsions[0] << "; psi: " << torsions[1];

                        if ( torsions.size() == 3 )
                            message << "; omega: " << torsions[2] << ";";
                        else
                            message << ";";

                        message << annotation;

                        return message.str();

                    }

                    void add_annotation ( std::string message )
                    {
                        annotation = annotation + " " + message;
                    }

                    std::string get_annotation ( ) { return annotation; }

                    std::vector<float> get_torsions () const
                    {
                        std::vector<float> result;

                        result.push_back ( torsion_phi );
                        result.push_back ( torsion_psi );

                        if (index == 6 || (index == 1 && torsion_omega != 0.0))
                        {
                            result.push_back ( torsion_omega );
                        }

                        if (index == 7)
                        {
                            result.push_back ( torsion_omega_six );
                            result.push_back ( torsion_omega_seven );
                        }

                        if (index == 8)
                        {
                            result.push_back ( torsion_omega_seven );
                            result.push_back ( torsion_omega_eight );
                            result.push_back ( torsion_omega_nine );
                            result.push_back ( torsion_phi_cone_ctwo_oeight_ceight );
                        }

                        return result;

                    } //!< [0]=phi, [1]=psi, { [2]=omega if 1-6 linkage }

                    std::pair<clipper::MAtom, clipper::MAtom> get_linkage_atoms () const
                    {
                        return std::make_pair(donorAtom, acceptorAtom);
                    } //!< .first=donorAtom, .second=acceptorAtom

                    bool check_if_linkage_zscore_calculated () const
                    {
                        return linkage_zscore_calculated;
                    }

                    float get_linkage_zscore () const
                    {
                        return linkage_zscore;
                    }

                    bool get_linkage_enough_datapoints () const
                    {
                        return linkage_enough_datapoints;
                    }

                    void set_torsions ( float phi, float psi, float omega, float omega_six, float omega_seven, float omega_eight, float omega_nine, float phi_cone_ctwo_oeight_ceight )
                    {
                        torsion_phi         =   phi;
                        torsion_psi         =   psi;
                        torsion_omega       =   omega;
                        torsion_omega_six   =   omega_six;
                        torsion_omega_seven =   omega_seven;
                        torsion_omega_eight =   omega_eight;
                        torsion_omega_nine  =   omega_nine;
                        torsion_phi_cone_ctwo_oeight_ceight = phi_cone_ctwo_oeight_ceight;
                    }

                    void set_linkage_atoms( clipper::MAtom& inputDonorAtom, clipper::MAtom& inputAcceptorAtom) { donorAtom = inputDonorAtom; acceptorAtom = inputAcceptorAtom; };

                    void set_linkage_zscore( float input_zscore ) { this->linkage_zscore = input_zscore; this->linkage_zscore_calculated = true; };

                    void set_linkage_enough_datapoints( bool input_enough_datapoints ) { this->linkage_enough_datapoints = input_enough_datapoints; };

                    std::string format() const
                    {
                        std::stringstream s;
                        s << type << " 1-" << index << " to "
                          << node_id;

                        return s.str();
                    }

                    void calculate_and_set_zscore(float Phi, float Psi, clipper::String first_residue_name, clipper::MAtom first_atom, clipper::String second_residue_name, clipper::MAtom second_atom, privateer::json::GlobalTorsionZScore& torsions_zscore_database);

                    float calculate_zscore(float phi, float psi, privateer::json::TorsionsZScoreDatabase& matched_linkage);
                private:
                    float torsion_phi;
                    float torsion_psi;
                    float torsion_omega;        // for 1-6 linkages
                    float torsion_omega_six;    // for 1-7 linkages
                    float torsion_omega_seven;  // for 1-7 and 2-8 linkages
                    float torsion_omega_eight;  // for 2-8 linkages
                    float torsion_omega_nine;   // for 2-8 linkages
                    float torsion_phi_cone_ctwo_oeight_ceight; // for 2-8 linkages
                    float linkage_zscore;
                    bool linkage_enough_datapoints;
                    int index;                  // carbon to which this is connected
                    int node_id;                // sugar connected to by this linkage
                    bool linkage_zscore_calculated;
                    bool noncircular;
                    std::string type;           // anomer
                    std::string annotation;     // include validation information
                    clipper::MAtom donorAtom;
                    clipper::MAtom acceptorAtom;

            }; // class Linkage


            class Node
            {
                public:

                    Node () { initialised = false; }
                    Node ( clipper::MSugar& new_sugar ) { sugar = clipper::MSugar(new_sugar); initialised = true; }
                    Node ( const Node& node )  { this->sugar = clipper::MSugar(node.sugar); this->connections = node.connections; initialised = true; }

                    Node& operator= ( const Node& node ) // this operator seems unnecessary - will keep it in case implementation changes
                    {
                        if ( this != &node )
                        {
                            connections = node.connections;
                            sugar = node.sugar;
                        }
                        initialised = true;
                        return *this;
                    }

                    void add_annotation ( std::string message )
                    {
                        annotation = annotation + " " + message;
                    }

                    std::string get_annotation ( ) { return annotation; }

                    const clipper::MSugar& get_sugar () const { return sugar; }
                    void set_sugar ( clipper::MSugar& new_sugar ) { sugar = clipper::MSugar(new_sugar); }

                    //const std::vector< Linkage >& get_connections () const { return connections; }
                    int add_connection ( const Linkage& connection ) { connections.push_back ( connection ); return connections.size()-1; }
                    int add_circular_connection ( const Linkage& connection ) { connections.insert(connections.begin(), connection); return 0; }

                    const int number_of_connections ( ) const { return connections.size(); }

                    void remove_connection ( const int vectorIndex )
                    {
                        if ( vectorIndex > connections.size() - 1 )
                            connections.pop_back();
                        else
                            connections.erase(connections.begin() + vectorIndex);
                    }


                    Linkage& get_connection ( const int index )
                    {
                        if ( index > connections.size() -1 )
                            return connections.back();
                        else
                            return connections[index];
                    }

                    const std::string format() const
                    {
                        std::stringstream s;

                        s << "Address: " << this
                          << " Sugar: " << sugar.type()
                          << "; Connections: \n";

                        for ( int i = 0; i < number_of_connections() ; i++ )
                            s << connections[i].format() << "\n";

                        return s.str();
                    }

                    const bool is_initialised () { return this->initialised ; }

                private:
                    bool initialised;
                    std::vector < Linkage > connections;
                    clipper::MSugar sugar;
                    std::string annotation;

            }; // class Node

            bool link_sugars  ( int link, clipper::MSugar& first_sugar, clipper::MSugar& next_sugar, clipper::MAtom& donorAtom, clipper::MAtom& acceptorAtom, bool noncircular, privateer::json::GlobalTorsionZScore& torsions_zscore_database ); // true if there's been any problem
            void add_torsions_for_detected_linkages(float Phi, float Psi, clipper::String first_residue_name, clipper::MAtom first_atom, clipper::String second_residue_name, clipper::MAtom second_atom, int first_residue_seqnum, int second_residue_seqnum);
            std::vector<MGlycanTorsionSummary> return_torsion_summary_within_glycan() { return all_torsions_within_mglycan; };
            const std::pair < clipper::MMonomer, clipper::MSugar >& get_root () const { return this->root; }
            const clipper::String& get_type () const { return kind_of_glycan; } // n-glycan, o-glycan, s-glycan, c-glycan, p-glycan or ligand
            // std::string get_root_by_name () const { return get_root().first.type().trim() + "-" + get_root().first.id().trim() + "/" + get_chain().substr(0,1); }
            std::string get_root_by_name () const { return get_root().second.type().trim() + "-" + get_root().second.id().trim() + "/" + get_root_sugar_chainID().trim().substr(0,1) + "_" + get_root().first.type().trim() + "-" + get_root().first.id().trim() + "/" + get_chain().substr(0,1); }
            std::string get_root_for_filename () { return  get_root().second.type().trim() + get_root().second.id().trim() + "-[" + get_root_sugar_chainID().trim().substr(0,1) + "]_[" + get_chain().trim().substr(0,1) + "]-" + get_root().first.type().trim() + get_root().first.id().trim(); }

            std::string write_ring_ext_restraints ( float sigma );
            std::string write_link_ext_restraints ( float sigma );

            clipper::String print_linear ( const bool print_info, const bool html_format, const bool translate );
            clipper::String print_SVG ( bool vertical, bool print_info, bool colour_gradient );

            // NEW FUNCTIONS INTRODUCED DUE TO WURCS IMPLEMENTATION BEGIN //
            char convertNumberToLetter(int number); // need to be relocated, doesn't really belong under ::MGlycan.
            std::vector < std::string > obtain_unique_WURCS_residues();
            std::vector < std::string > obtain_unique_residue_codes();
            const int obtain_total_number_of_glycosidic_bonds();
            clipper::String generate_wurcs ();
            // NEW FUNCTIONS INTRODUCED DUE TO WURCS IMPLEMENTATION END //


            std::vector < clipper::MSugar >& get_sugars () { return sugars; }

            const Node& get_node ( int index ) const { if (index>node_list.size()-1) return node_list.back(); else return node_list[index]; }
            std::vector<Node> get_node_list_vector() { return node_list; };

            const clipper::String& get_chain () const { return chain; }
            const clipper::String& get_root_sugar_chainID () const { return chain_root_sugar; }
            // const std::string get_wurcs() const { return wurcs; }
            int number_of_nodes ( ) const { return node_list.size(); }

            void remove_node_at_index ( int index );
            void replace_sugar_at_index ( int index, clipper::MSugar& donor );
            void update_msugar_in_root ( clipper::MSugar& newmsug );

            int get_number_of_connections_at_index ( int index ) const { if (index>node_list.size()-1) return node_list.back().number_of_connections(); else return node_list[index].number_of_connections(); }

            void set_kind_of_glycan ( clipper::String input ) { kind_of_glycan = input; }

            std::vector<float> get_glycosylation_torsions () const
            {
                std::vector<float> result;
                result.push_back ( torsion_phi );
                result.push_back ( torsion_psi );
                return result;
            } //!< [0]=phi, [1]=psi

            std::string print_torsions () const
                {
                    std::string output;

                    output = "Phi: " + std::to_string(torsion_phi) + " Psi: " + std::to_string(torsion_psi);

                    return output;
                }

            void set_glycosylation_torsions ( float phi, float psi ) { torsion_phi = phi; torsion_psi = psi; }

            std::string get_link_description ( )
            {
                std::ostringstream message;

                message << root.second.anomer() << " link. "
                        << "Phi: " << torsion_phi << "; Psi: " << torsion_psi << ";"
                        << link_annotation ;

                return message.str();

            }

            std::string get_root_description ( )
            {
                std::ostringstream message;

                message << get_root_by_name()
                        << root_annotation;

                return message.str();
            }

            void add_root_annotation ( std::string annotation )
            {
                root_annotation = root_annotation + annotation;
            }

            void add_link_annotation ( std::string annotation )
            {
                link_annotation = link_annotation + annotation;
            }

            void set_annotations ( std::string expression_system );  // function that annotates glycobiologic validation

            void set_protein_sugar_linkage_zscore_attempt_to_calculate (bool input_value ) {this->protein_sugar_linkage_zscore_attempt_to_calculate = input_value; };
            bool get_protein_sugar_linkage_zscore_attempt_to_calculate ( ) { return protein_sugar_linkage_zscore; };
            void set_protein_sugar_linkage_zscore (float input_value ) {this->protein_sugar_linkage_zscore = input_value; };
            float get_protein_sugar_linkage_zscore ( ) { return protein_sugar_linkage_zscore; };
            float calculate_zscore(float phi, float psi, privateer::json::TorsionsZScoreDatabase& matched_linkage);

        private:
            bool debug_output;
            clipper::String kind_of_glycan;                 // can be n-glycan, o-glycan, s-glycan, c-glycan, p-glycan or ligand
            std::pair < clipper::MMonomer, clipper::MSugar > root;     // this is the root, should be null if this is a ligand saccharide
            std::vector < Node > node_list;                 // interlinked nodes
            clipper::String chain;                          // Chain ID for this glycan
            clipper::String chain_root_sugar;              // Chain ID for root sugar
            std::vector < clipper::MSugar > sugars;         // vector of sugars
            clipper::ftype torsion_psi, torsion_phi;        // Torsions of the protein-glycan link
            std::string root_annotation, link_annotation, expression_system;
            std::vector<MGlycanTorsionSummary> all_torsions_within_mglycan;
            float protein_sugar_linkage_zscore = 42069;    // crappy hack implemented only for ASN-NAG linkage. For anyone who is going to be working on this
                                                            // in near future, consider modifying clipper::MGlycan::Linkage class to also be associated
                                                            // with amino acid - sugar linkages, rather than sugar - sugar linkages.
            bool protein_sugar_linkage_zscore_attempt_to_calculate = false;

    }; // class MGlycan

    class MGlycology
    {
        public:

            MGlycology () { } //!< null constructor
            MGlycology ( const clipper::MiniMol&, bool debug_output, std::string expression_system = "undefined" );
            MGlycology ( const clipper::MiniMol&, const clipper::MAtomNonBond&, privateer::json::GlobalTorsionZScore&, bool debug_output, std::string expression_system = "undefined" );

            void init ( const clipper::MiniMol&, const clipper::MAtomNonBond&, privateer::json::GlobalTorsionZScore&,  bool debug_output, std::string expression_system );
            clipper::MGlycan get_glycan_by_id ( int id ) { return list_of_glycans[id]; };
            clipper::MGlycan get_glycan_by_root ( clipper::MMonomer& root )
            {
                for (int i=0;i<list_of_glycans.size();i++)
                    if (list_of_glycans[i].get_root().first.id()==root.id())
                        return list_of_glycans[i];
                return clipper::MGlycan();
            }
            std::vector < clipper::MGlycan > get_list_of_glycans () const { return list_of_glycans; }
            std::vector < clipper::MSugar > get_sugar_list() { return list_of_sugars; }
            std::string write_external_restraints ( bool restrain_rings,
                                                    bool restrain_links,
                                                    float sigma = 0.1 );

        private:

            // private properties
            bool debug_output;
            std::vector < clipper::MSugar >  list_of_sugars;
            std::vector < clipper::MGlycan > list_of_glycans;
            const clipper::MAtomNonBond* manb;
            const clipper::MiniMol * mmol;

            // private methods
            const std::vector < std::pair< clipper::String, clipper::MMonomer > > get_overlapping_residues ( const clipper::MMonomer& mm );
            const std::vector < std::pair< clipper::MAtom, clipper::MAtomIndexSymmetry > > get_contacts ( const clipper::MMonomer& mm, const clipper::String monomer_chain_id );
            int parse_order ( clipper::MAtom& atom_in_sugar, clipper::MSugar& sugar );
            void extend_tree ( clipper::MGlycan& mg, clipper::MSugar& msug, std::vector<clipper::MSugar>& accounted_for_sugars, privateer::json::GlobalTorsionZScore& torsions_zscore_database );
            const char get_altconf ( const clipper::MAtom& ) const;
            std::string expression_system;

    }; // class MGlycology

} // namespace clipper


#endif
