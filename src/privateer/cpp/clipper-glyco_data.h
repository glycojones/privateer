/*! \file clipper-glyco_data.h
  Header file for sugar data */

// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York



#ifndef CLIPPER_GLYCO_DATA
#define CLIPPER_GLYCO_DATA

#include <clipper/clipper.h>
#include <unordered_set>
#include <set>
#include <algorithm>


#ifdef _MSC_VER
#ifdef PRIVATEER_BUILD_DLL
#define PRIVATEER_DLL_IMPORT(type) __declspec(dllexport) type
#else
#define PRIVATEER_DLL_IMPORT(type) __declspec(dllimport) type
#endif
#else
#define PRIVATEER_DLL_IMPORT(type) type
#endif

namespace clipper
{
    namespace data
    {

        struct sugar_database_entry
        {
            String name_short;       // three-character name
            String anomer;           // A=alpha, B=beta, N=unknown
            String handedness;       // D, L or N for unknown
            String name_long;        // longer name
            String ring_atoms;       // e.g. ring_atoms = { "C1","C2","C3","C4","C5","O5" };
            ftype  ref_puckering;    // puckering amplitude calculated from CCP4 monomer library-idealised coordinates
            String ref_conformation; // conformation detected from the idealised coordinates
            ftype  ref_bonds_rmsd;   // rms bonds calculated from CCP4 monomer library-idealised coordinates
            ftype  ref_angles_rmsd;  // rms angles calculated from CCP4 monomer library-idealised coordinates
        };

        struct disaccharide                   // for those particular cases
        {
            String name_short;
            String name_long;
            sugar_database_entry sugar_one;
            sugar_database_entry sugar_two;
        };

        #ifdef _MSC_VER
            extern PRIVATEER_DLL_IMPORT(const String) iupac_conformational_landscape[];
            extern PRIVATEER_DLL_IMPORT(const String) conformational_landscape[];
            extern PRIVATEER_DLL_IMPORT(const sugar_database_entry) sugar_database[];
            extern PRIVATEER_DLL_IMPORT(const disaccharide) disaccharide_database[];
            extern PRIVATEER_DLL_IMPORT(const int) disaccharide_database_size;
            extern PRIVATEER_DLL_IMPORT(const int) sugar_database_size;
        #else
            extern const String iupac_conformational_landscape[];
            extern const String conformational_landscape[];
            extern const sugar_database_entry sugar_database[];
            extern const disaccharide disaccharide_database[];
            extern const int disaccharide_database_size;
            extern const int sugar_database_size;
        #endif  


        bool found_in_database ( clipper::String name );
        bool found_in_database ( std::string name );
        bool is_amino_acid ( std::string name );
        std::string carbname_of ( std::string name );
        std::string alternative_anomer ( std::string name );
        std::vector<std::string> alternative_monomer( std::string name );
        bool residue_has_alternate_anomer( std::string name );
        bool residue_has_alternate_monomer( std::string name, bool glucose_only );
        std::string convert_to_wurcs_residue_code( std::string name );
        std::string get_anomer( std::string name );

    } // namespace data

} // namespace clipper

#endif
