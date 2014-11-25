/* \file privateer-lib.cpp
	A set of utilities that helps the Privateer do his job */
// version  0.1.0
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

#include "privateer-lib.h"

void privateer::insert_coot_prologue_scheme ( std::fstream& output )
{
    output  << "; This script has been created by Privateer-validate (Agirre and Cowtan, 2013-15)\n"	
            << "(set-graphics-window-size 1873 968)\n"
	    << "(set-graphics-window-position 0 0)\n"
            << "(set-go-to-atom-window-position 0 19)\n"
            << "(set-display-control-dialog-position 366 20)\n"
            << "(vt-surface 2)\n"
            << "(set-clipping-front  0.00)\n"
            << "(set-clipping-back  0.00)\n"
            << "(set-map-radius 10.00)\n"
            << "(set-iso-level-increment  0.0500)\n"
            << "(set-diff-map-iso-level-increment  0.0050)\n"
            << "(set-colour-map-rotation-on-read-pdb 21.00)\n"
            << "(set-colour-map-rotation-on-read-pdb-flag 1)\n"
            << "(set-colour-map-rotation-on-read-pdb-c-only-flag 1)\n"
            << "(set-swap-difference-map-colours 0)\n"
            << "(set-background-colour  0.00  0.00  0.00)\n"
            << "(set-symmetry-size 13.00)\n"
            << "(set-symmetry-colour-merge  0.50)\n"
            << "(set-symmetry-colour  0.10  0.20  0.80)\n"
            << "(set-symmetry-atom-labels-expanded 0)\n"
            << "(set-active-map-drag-flag 1)\n"
            << "(set-show-aniso 0)\n"
            << "(set-aniso-probability 50.00)\n"
            << "(set-smooth-scroll-steps 40)\n"
            << "(set-smooth-scroll-limit 10.00)\n"
            << "(set-font-size 2)\n"
            << "(set-rotation-centre-size  0.10)\n"
            << "(set-default-bond-thickness 5)\n"
            << "(scale-zoom  0.20)\n"
            << "(set-nomenclature-errors-on-read \"auto-correct\")"
            << "(set-run-state-file-status 0)\n";
}

void privateer::insert_coot_files_loadup_scheme ( std::fstream& output, const clipper::String& pdb, const clipper::String& mapbest, const clipper::String& mapdiff, const clipper::String& mapomit, bool mode )
{
    if (!mode) output << "(handle-read-draw-molecule \"" << pdb << "\")\n";
    
    if ( mapbest == "" ) // no map output
    {			
        output << "(set-last-map-colour 1.00  0.13  0.89)\n"
	       << "(interesting-things-gui \"Validation report from Privateer\"\n\t(list\n\t\t";
    }
    else
    {		
        if (!mode) 
            output << "(handle-read-ccp4-map \"" << mapbest << "\" 0)\n" << "(handle-read-ccp4-map \"" << mapomit << "\" 1)\n";
	
        output << "(set-last-map-colour 1.00  0.13  0.89)\n"
	       << "(interesting-things-gui \"Validation report from Privateer\"\n\t(list\n\t\t";
    }
}
 
void privateer::insert_coot_files_loadup_python ( std::fstream& output, const clipper::String& pdb, const clipper::String& mapbest, const clipper::String& mapdiff, const clipper::String& mapomit, bool mode )
{
	if (!mode) output  << "handle_read_draw_molecule (\"" << pdb << "\")\n";
	
        if ( mapbest == "" ) // no map output
	{
		output << "set_last_map_colour  (1.00,  0.13,  0.89)\n"
		       << "interesting_things_gui (\"Validation report from Privateer\",[\n";
	}
	else
	{
	    if (!mode)
                output << "handle_read_ccp4_map (\"" << mapbest << "\", 0)\n" << "handle_read_ccp4_map (\"" << mapomit << "\", 1)\n";
            
            output << "set_last_map_colour  (1.00,  0.13,  0.89)\n"
		   << "interesting_things_gui (\"Validation report from Privateer\",[\n";
	}
} 

void privateer::insert_coot_epilogue_scheme ( std::fstream& output )
{	
	output  << "\n\n))\n(set-scroll-wheel-map 3)\n"
			<< "(set-matrix 60.00)\n"
                        << "(set-refine-with-torsion-restraints 1)\n"
			<< "(set-show-symmetry-master 0)\n";
}

void privateer::insert_coot_prologue_python ( std::fstream& output )
{

    output  << "# This script has been created by Privateer-validate (Agirre and Cowtan, 2013-15)\n"		
            << "set_graphics_window_size (1873, 968)\n"
	    << "set_graphics_window_position (0, 0)\n"
	    << "set_go_to_atom_window_position (0, 19)\n"
	    << "vt_surface (2)\n"
	    << "set_clipping_front  (0.00)\n"
	    << "set_clipping_back  (0.00)\n"
	    << "set_map_radius (10.00)\n"
	    << "set_iso_level_increment  (0.0500)\n"
	    << "set_diff_map_iso_level_increment  (0.0050)\n"
	    << "set_colour_map_rotation_on_read_pdb (21.00)\n"
	    << "set_colour_map_rotation_on_read_pdb_flag (1)\n"
	    << "set_colour_map_rotation_on_read_pdb_c_only_flag (1)\n"
	    << "set_swap_difference_map_colours (0)\n"
	    << "set_background_colour  (0.00,  0.00,  0.00)\n"
	    << "set_symmetry_size (13.00)\n"
	    << "set_symmetry_colour_merge  (0.50)\n"
	    << "set_symmetry_colour  (0.10,  0.20,  0.80)\n"
	    << "set_symmetry_atom_labels_expanded (0)\n"
	    << "set_active_map_drag_flag (1)\n"
	    << "set_show_aniso (0)\n"
	    << "set_aniso_probability (50.00)\n"
	    << "set_smooth_scroll_steps (40)\n"
	    << "set_smooth_scroll_limit (10.00)\n"
	    << "set_font_size (2)\n"
	    << "set_rotation_centre_size (0.10)\n"
	    << "set_default_bond_thickness (4)\n"
	    << "scale_zoom (0.20)\n"
	    << "set_nomenclature_errors_on_read (\"auto-correct\")\n"
	    << "set_run_state_file_status (0)\n"
            << "toggle_idle_spin_function\n";	
}

void privateer::insert_coot_epilogue_python ( std::fstream& output )
{
	output  << "\n\n])\nset_scroll_wheel_map (3)\n"
			<< "set_matrix (60.00)\n"
                        << "set_refine_with_torsion_restraints (1)\n"
			<< "set_show_symmetry_master (0)\n";
}

void privateer::insert_coot_go_to_sugar_scheme ( std::fstream& output, const clipper::Coord_orth& sugar_centre, const clipper::String& diagnostic )
{
	output  << "\t(list\t\"" << diagnostic << "\"\t" << sugar_centre.x() << "\t" << sugar_centre.y() << "\t" << sugar_centre.z() << ")\n";
}

void privateer::insert_coot_go_to_sugar_python ( std::fstream& output, const clipper::Coord_orth& sugar_centre, const clipper::String& diagnostic )
{	
	output  << "\t[\"" << diagnostic << "\",\t" << sugar_centre.x() << ",\t" << sugar_centre.y() << ",\t" << sugar_centre.z() << "],\n";
}

void privateer::insert_coot_statusbar_text_scheme ( std::fstream& output, clipper::String& text)
{
	output  << "(add-status-bar-text \"" << text << "\")" ; 
}

void privateer::insert_coot_statusbar_text_python ( std::fstream& output, clipper::String& text)
{
	output  << "add_status_bar_text (\"" << text << "\")" ; 
}

clipper::ftype privateer::real_space_correlation ( const clipper::Xmap<float>& map1, const clipper::Xmap<float>& map2 ) 
{
	return 0.0;
}
