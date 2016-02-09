
// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2016 Jon Agirre & Kevin Cowtan
// York Structural Biology Laboratory
// The University of York
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk
//

#include "privateer-lib.h"

void privateer::coot::insert_coot_prologue_scheme ( std::fstream& output )
{
    output  << "; This script has been created by Privateer (Agirre, Iglesias, Rovira, Davies, Wilson and Cowtan, 2013-16)\n"
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

void privateer::coot::insert_coot_files_loadup_scheme ( std::fstream& output, const clipper::String& pdb, const clipper::String& mapbest, const clipper::String& mapdiff, const clipper::String& mapomit, bool mode )
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

void privateer::coot::insert_coot_files_loadup_python ( std::fstream& output, const clipper::String& pdb, const clipper::String& mapbest, const clipper::String& mapdiff, const clipper::String& mapomit, bool mode )
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

void privateer::coot::insert_coot_epilogue_scheme ( std::fstream& output )
{
    output  << "\n\n))\n(set-scroll-wheel-map 3)\n"
            << "(set-matrix 60.00)\n"
                        << "(set-refine-with-torsion-restraints 1)\n"
            << "(set-show-symmetry-master 0)\n";
}

void privateer::coot::insert_coot_prologue_python ( std::fstream& output )
{

    output  << "# This script has been created by Privateer (Agirre, Iglesias, Rovira, Davies, Wilson and Cowtan, 2013-16)\n"
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

void privateer::coot::insert_coot_epilogue_python ( std::fstream& output )
{
    output  << "\n\n])\nset_scroll_wheel_map (3)\n"
            << "set_matrix (60.00)\n"
                        << "set_refine_with_torsion_restraints (1)\n"
            << "set_show_symmetry_master (0)\n";
}

void privateer::coot::insert_coot_command ( std::fstream& output, std::string command )
{
    output << command << "\n" ;
}

void privateer::coot::insert_coot_go_to_sugar_scheme ( std::fstream& output, const clipper::Coord_orth& sugar_centre, const clipper::String& diagnostic )
{
    output  << "\t(list\t\"" << diagnostic << "\"\t" << sugar_centre.x() << "\t" << sugar_centre.y() << "\t" << sugar_centre.z() << ")\n";
}

void privateer::coot::insert_coot_go_to_sugar_python ( std::fstream& output, const clipper::Coord_orth& sugar_centre, const clipper::String& diagnostic )
{
    output  << "\t[\"" << diagnostic << "\",\t" << sugar_centre.x() << ",\t" << sugar_centre.y() << ",\t" << sugar_centre.z() << "],\n";
}

void privateer::coot::insert_coot_statusbar_text_scheme ( std::fstream& output, clipper::String& text)
{
    output  << "(add-status-bar-text \"" << text << "\")" ;
}

void privateer::coot::insert_coot_statusbar_text_python ( std::fstream& output, clipper::String& text )
{
    output  << "add_status_bar_text (\"" << text << "\")" ;
}

float privateer::real_space_correlation ( const clipper::Xmap<float>& map1, const clipper::Xmap<float>& map2 )
{
    return 0.0; // to be implemented
}




///////// Privateer's glycoplot library /////////

std::string privateer::glycoplot::get_colour ( Colour colour, bool original_style )
{
    if ( original_style )    // Essentials of glycobiology, 3rd edition //
    {
        switch (colour)
        {
            case blue:
                return "#0000fa;";
            case green:
                return "#00c832;";
            case black:
                return "#000000;";
            case orange:
                return "#fa6400;";
            case yellow:
                return "#ffff00;";
            case tan:
                return "#966432;";
            case purple:
                return "#7d007d;";
            case red:
                return "#fa0000;";
            case cyan:
                return "#c8fafa;";
        }
        return "#ffffff;";
    }
    else                      // Use washed-out colours, less disruptive in a publication //
    {
        switch (colour)
        {
            case blue:
                return "#0000fa;";
            case green:
                return "#00c832;";
            case black:
                return "#000000;";
            case orange:
                return "#fa6400;";
            case yellow:
                return "#ffff00;";
            case tan:
                return "#966432;";
            case purple:
                return "#7d007d;";
            case red:
                return "#fa0000;";
            case cyan:
                return "#c8fafa;";
        }
        return "#ffffff;";
    }
}

void privateer::glycoplot::Plot::write_svg_header   ( std::fstream& of )
{
    
    of << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n\n"
       << "<!-- Generator: Privateer (YSBL, University of York, distributed by CCP4) -->\n"
       << "<!-- Please reference: Agirre, Iglesias, Rovira, Davies, Wilson & Cowtan (2015) Nat Struct & Mol Biol 22(11), 833-834 -->\n\n"
       << "<svg xmlns:dc=\"http://purl.org/dc/elements/1.1/\"\n"
       << "  xmlns:cc=\"http://creativecommons.org/ns#\"\n"
       << "  xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"\n"
       << "  xmlns:svg=\"http://www.w3.org/2000/svg\"\n"
       << "  xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n"
       << "  xmlns=\"http://www.w3.org/2000/svg\"\n"
       << "  version=\"1.1\"\n"
       << "  width=\"" << get_width() << "\" height=\"" << get_height() << "\">\n\n" ;
    
}

void privateer::glycoplot::Plot::write_svg_definitions( std::fstream& of )
{
    of << "  <defs>\n"
       << "    <!-- GlcNAc --> "
       <<  "<rect width =\"50\" height=\"50\" id=\"glcnac\" style=\" stroke:"
       << get_colour ( black, original_colour_scheme ) << " fill:" << get_colour ( blue, original_colour_scheme )
       << "stroke-width:0.6;\" />\n"
       << "  </defs>\n\n" ;
}

void privateer::glycoplot::Plot::write_svg_contents ( std::fstream& of )
{
    for (int i = 0; i < list_of_shapes.size() ; i ++)
    {
        of << list_of_shapes[i]->get_XML();
    }
}


void privateer::glycoplot::Plot::write_svg_footer ( std::fstream& of )
{
    of << "</svg>" ;
}


std::string privateer::glycoplot::Plot::get_svg_string_header   ( )
{
    std::ostringstream of;
    
    of << "<svg version=\"1.1\" width=\"" << get_width() << "\" height=\"" << get_height() << "\">\n" ;
    
    return of.str();
}


std::string privateer::glycoplot::Plot::get_svg_string_contents ( )
{
    std::ostringstream of;
    
    for (int i = 0; i < list_of_shapes.size() ; i ++)
    {
        of << list_of_shapes[i]->get_XML();
    }
    
    return of.str();
}


std::string privateer::glycoplot::Plot::get_svg_string_footer ( )
{
    std::ostringstream of;
    
    of << "</svg>" ;
    
    return of.str();
}


bool privateer::glycoplot::Plot::write_to_file  ( std::string file_path )
{
    std::fstream out;
    
    out.open( file_path, std::fstream::out);
    
    write_svg_header      ( out );
    write_svg_definitions ( out );
    write_svg_contents    ( out );
    write_svg_footer      ( out );
    
    out.close();
    
    return false;
}


std::string privateer::glycoplot::Plot::get_XML  ( )
{
    return get_svg_string_header() + get_svg_string_contents() + get_svg_string_footer();
}


bool privateer::glycoplot::Plot::plot_glycan ( clipper::MGlycan glycan )
{
    
    privateer::glycoplot::GlcNAc *glcnac = new privateer::glycoplot::GlcNAc(2, 3, "This GlcNAc is quite all right" ); // ignored for now
    //privateer::glycoplot::Man *man = new privateer::glycoplot::Man(2,3);
    //privateer::glycoplot::Fuc *fuc = new privateer::glycoplot::Fuc(2,3);
    add_shape(glcnac);
    //add_shape(man);
    //add_shape(fuc);
    std::cout << "adding one GlcNAc" << std::endl;
    
    return false;
}


std::string privateer::glycoplot::GlcNAc::get_XML ()
{
    std::ostringstream tmp;
    
    tmp   <<  "  <use xlink:href=\"#glcnac\" x=\"" << get_x() << "\""
          <<  " y=\"" << get_y() << "\" id=\"" << get_id() << "\" " << "title=\"" << get_tooltip() << "\">\n"
          <<  "    <title>" << get_tooltip() << "</title>\n"
          <<  "  </use>";
    
    return tmp.str();
}


std::string privateer::glycoplot::Man::get_XML ()
{
    std::ostringstream tmp;
    
    tmp   <<  "<circle \n" <<
              "      r =\""             << get_radius()      << "\"\n" <<
              "      cx=\""             << get_x()           << "\"\n" <<
              "      cy=\""             << get_y()           << "\"\n" <<
              "      id=\""             << get_id()          << "\"\n" <<
              "      title=\""          << get_tooltip()     << "\"\n" <<
              "/>\n";
    
    return tmp.str();
}


std::string privateer::glycoplot::Fuc::get_XML ()
{
    std::ostringstream tmp;
    
    tmp   <<  "<polygon \n" <<
              "       points='210 100, 210 200, 270 150'"  <<
              "      rx=\"0\"\n"        <<
              "      ry=\"0\"\n"        <<
              "      id=\""             << get_id()          << "\"\n"  <<
              "      title=\""          << get_tooltip()     << "\"\n"  <<
              "/>\n";
    
    return tmp.str();
}