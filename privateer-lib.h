/*! \file privateer-lib.h
  A collection of tools that help the Privateer do his job */
// version  MKII
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

#ifndef PRIVATEER_LIB_H_INCLUDED
#define PRIVATEER_LIB_H_INCLUDED

#include <fstream>
#include <algorithm>
#include <clipper/clipper.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-minimol.h>
#include "privateer-ss-find.h"


namespace privateer
{

    // Data for Privateer
    namespace data
    {
        struct dummy_atom
        {
            clipper::String atom_name;
            clipper::ftype x;
            clipper::ftype y;
            clipper::ftype z;
        };

        struct fingerprint
        {
            clipper::String name_short;
            clipper::String context;
            clipper::String anomer;
            clipper::String handedness;
            clipper::String name_long;
            int num_control_points;
            dummy_atom peaks[20];
            dummy_atom voids[20];
            dummy_atom atoms[20];
        };

        struct build_options
        {
            bool nglycans;
            bool oglycans;
            bool ligands;
        };

        extern const fingerprint fingerprint_list[];
        extern const int fingerprint_list_size;
    } // namespace data

    namespace coot
    {
        // Coot support, Scheme
        void insert_coot_prologue_scheme ( std::fstream& );
        void insert_coot_epilogue_scheme ( std::fstream& );
        void insert_coot_files_loadup_scheme ( std::fstream&, const clipper::String&, const clipper::String&, const clipper::String&, const clipper::String&, bool mode );
        void insert_coot_go_to_sugar_scheme ( std::fstream&, const clipper::Coord_orth& sugar_centre, const clipper::String& diagnostic );
        void insert_coot_statusbar_text_scheme ( std::fstream&, clipper::String& );

        // Coot support, Python
        void insert_coot_files_loadup_python ( std::fstream&, const clipper::String&, const clipper::String&, const clipper::String&, const clipper::String&, bool mode );
        void insert_coot_prologue_python ( std::fstream& );
        void insert_coot_epilogue_python ( std::fstream& );
        void insert_coot_go_to_sugar_python ( std::fstream&, const clipper::Coord_orth& sugar_centre, const clipper::String& diagnostic );
        void insert_coot_statusbar_text_python ( std::fstream&, clipper::String& );
    }
        
    // General Privateer functions
    clipper::ftype real_space_correlation ( const clipper::Xmap<float>&, const clipper::Xmap<float>& );

    clipper::MMonomer get_ideal_monomer ( const privateer::data::fingerprint& fp );
    clipper::MMonomer get_peak_monomer  ( const privateer::data::fingerprint& fp );
    clipper::MMonomer get_void_monomer  ( const privateer::data::fingerprint& fp );

    void process_building_options ( clipper::String building_options, privateer::data::build_options &flags );
    clipper::MiniMol build_sugars ( clipper::Xmap<float>& xwrk, privateer::data::build_options& options, double step, int nhit );

} // namespace privateer

#endif
