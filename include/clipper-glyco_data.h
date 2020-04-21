/*! \file clipper-glyco_data.h
  Header file for sugar data */

//C (C) 2014 Jon Agirre & Kevin Cowtan, The University of York
//L
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


#ifndef CLIPPER_GLYCO_DATA
#define CLIPPER_GLYCO_DATA

#include <clipper/clipper.h>

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

        extern const String iupac_conformational_landscape[];
        extern const String conformational_landscape[];
        extern const sugar_database_entry sugar_database[];
        extern const disaccharide disaccharide_database[];
        extern const int disaccharide_database_size;
        extern const int sugar_database_size;

        bool found_in_database ( clipper::String name );
        bool found_in_database ( std::string name );
        std::string carbname_of ( std::string name );
        std::string convert_to_wurcs_residue_code( std::string name );

    } // namespace data

} // namespace clipper

#endif
