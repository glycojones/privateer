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

void privateer::insert_coot_statusbar_text_python ( std::fstream& output, clipper::String& text )
{
    output  << "add_status_bar_text (\"" << text << "\")" ;
}

clipper::ftype privateer::real_space_correlation ( const clipper::Xmap<float>& map1, const clipper::Xmap<float>& map2 )
{
    return 0.0;
}

clipper::MMonomer privateer::get_ideal_monomer ( const privateer::data::fingerprint& fp )
{
    clipper::MMonomer tmp_mon;
    tmp_mon.set_type ( fp.name_short );
    tmp_mon.set_id ( 0 );

    for ( int index = 0; index < fp.num_control_points ; index++ )
    {
        clipper::MAtom tmp_atm;
        clipper::Coord_orth coords(fp.atoms[index].x, fp.atoms[index].y, fp.atoms[index].z );
        tmp_atm.set_coord_orth ( coords );
        tmp_atm.set_name ( fp.atoms[index].atom_name );
        tmp_atm.set_element ( fp.atoms[index].atom_name[0] );
        tmp_atm.set_id ( 1 );
        tmp_atm.set_occupancy ( 1.0 );
        tmp_atm.set_u_iso ( 0.25 ); // aiming for an iso-B of ~20.0
        tmp_mon.insert ( tmp_atm );
    }
    return tmp_mon;
}

clipper::MMonomer privateer::get_peak_monomer ( const privateer::data::fingerprint& fp )
{
    clipper::MMonomer tmp_mon;
    tmp_mon.set_type ( fp.name_short );
    tmp_mon.set_id ( 0 );

    for ( int index = 0; index < fp.num_control_points ; index++ )
    {
        clipper::MAtom tmp_atm;
        clipper::Coord_orth coords(fp.peaks[index].x, fp.peaks[index].y, fp.peaks[index].z );
        tmp_atm.set_coord_orth ( coords );
        tmp_atm.set_name ( fp.peaks[index].atom_name );
        tmp_atm.set_element ( fp.peaks[index].atom_name[0] );
        tmp_atm.set_id ( 1 );
        tmp_atm.set_occupancy ( 1.0 );
        tmp_atm.set_u_iso ( 0.25 ); // aiming for an iso-B of ~20.0
        tmp_mon.insert ( tmp_atm );
    }
    return tmp_mon;
}

clipper::MMonomer privateer::get_void_monomer ( const privateer::data::fingerprint& fp )
{
    clipper::MMonomer tmp_mon;
    tmp_mon.set_type ( fp.name_short );
    tmp_mon.set_id ( 0 );

    for ( int index = 0; index < fp.num_control_points ; index++ )
    {
        clipper::MAtom tmp_atm;
        clipper::Coord_orth coords(fp.voids[index].x, fp.voids[index].y, fp.voids[index].z );
        tmp_atm.set_coord_orth ( coords );
        tmp_atm.set_name ( fp.voids[index].atom_name );
        tmp_atm.set_element ( fp.voids[index].atom_name[0] );
        tmp_atm.set_id ( 1 );
        tmp_atm.set_occupancy ( 1.0 );
        tmp_atm.set_u_iso ( 0.25 ); // aiming for an iso-B of ~20.0
        tmp_mon.insert ( tmp_atm );
    }
    return tmp_mon;
}

void privateer::process_building_options ( clipper::String building_options, privateer::data::build_options &flags )
{
    std::vector < clipper::String > buffer = building_options.split( "," );

    for ( int i = 0; i < buffer.size(); i++ )
    {
        if ( buffer[i].trim() == "all" )
        {
            flags.nglycans = true;
            flags.oglycans = true;
            flags.ligands = true;
            return;
        }
        else if ( buffer[i].trim() == "nglycans" )
        {
            flags.nglycans = true;
        }
        else if ( buffer[i].trim() == "oglycans" )
        {
            flags.oglycans = true;
        }
        else if ( buffer[i].trim() == "ligands" )
        {
            flags.ligands = true;
        }
    }
    return;
}

clipper::MiniMol privateer::build_sugars ( clipper::Xmap<float>& xwrk, privateer::data::build_options& options, double step, int nhit )
{
    double rad, sigcut;
    clipper::MPolymer mprep;
    typedef std::pair<clipper::Coord_orth,clipper::Coord_orth> Pair_coord;
    std::vector<Pair_coord> all_co;


    // get cutoff (for optimisation)
    clipper::Map_stats stats( xwrk );
    sigcut = stats.mean() + 0.5*stats.std_dev();

    clipper::MAtom peak_atom, void_atom;
    clipper::MMonomer peaks, voids, ideal;
    
    ideal = privateer::get_ideal_monomer ( privateer::data::fingerprint_list[0] );
    peaks = privateer::get_peak_monomer  ( privateer::data::fingerprint_list[0] );
    voids = privateer::get_void_monomer  ( privateer::data::fingerprint_list[0] );

    mprep.insert ( ideal, -1 );

    // set up targets
    for ( int r = 0; r < peaks.size(); r++ )
        all_co.push_back( Pair_coord( peaks[r].coord_orth(), voids[r].coord_orth() ) );

    // get map radius
    clipper::Atom_list atoms = voids.atom_list();
    double r2 = 0.0;

    for ( int a = 0; a < atoms.size(); a++ )
    {
        double d2 = atoms[a].coord_orth().lengthsq();
        if ( d2 > r2 )
            r2 = d2;
    }

    rad = sqrt( r2 ) + 1.0;

    // make a list of rotations
    std::vector<clipper::RTop_orth> rots;

    // make a list of rotation ops to try
    float glim = 360.0;  // gamma
    float blim = 180.0;  // beta
    float alim = 360.0;  // alpha

    // do a uniformly sampled search of orientation space
    float anglim = clipper::Util::min( alim, glim );

    for ( float bdeg=step/2; bdeg < 180.0; bdeg += step )
    {
        float beta = clipper::Util::d2rad(bdeg);
        float spl = anglim/clipper::Util::intf(cos(0.5*beta)*anglim/step+1);
        float smi = anglim/clipper::Util::intf(sin(0.5*beta)*anglim/step+1);
        for ( float thpl=spl/2; thpl < 720.0; thpl += spl )
            for ( float thmi=smi/2; thmi < 360.0; thmi += smi )
            {
                float adeg = clipper::Util::mod(0.5*(thpl+thmi),360.0);
                float gdeg = clipper::Util::mod(0.5*(thpl-thmi),360.0);

                if ( adeg <= alim && bdeg <= blim && gdeg <= glim )
                {
                    float alpha = clipper::Util::d2rad(adeg);
                    float gamma = clipper::Util::d2rad(gdeg);
                    clipper::Euler_ccp4 euler( alpha, beta, gamma );
                    rots.push_back(clipper::RTop_orth(clipper::Rotation(euler).matrix()));
                }
            }
    }

    // feature search
    SSfind ssfind;
    ssfind.prep_xmap( xwrk, rad );
    ssfind.prep_search( xwrk );
    std::vector<SearchResult> results =

    ssfind.search( all_co, rots, sigcut, 0.0 );

    std::sort( results.begin(), results.end() );
    std::reverse( results.begin(), results.end() );
    std::cout << results.size() << std::endl;

    for ( int i = 0; i < clipper::Util::min( int(results.size()), nhit ); i++ )
        std::cout << i << " " << results[i].score << " " << results[i].rot << " " << results[i].trn << std::endl;

    // output
    clipper::MiniMol mol_new( xwrk.spacegroup(), xwrk.cell() );
    const clipper::String chainid1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const clipper::String chainid2 = "abcdefghijklmnopqrstuvwxyz";

    clipper::Grid_sampling grid = xwrk.grid_sampling();

    for ( int i = 0; i < clipper::Util::min( int(results.size()), nhit ); i++ )
    {
        clipper::String id; int roff;

        if ( i < 26 )
        {
            id = chainid1.substr( i, 1 );
            roff = 0;
        }
        else
        {
            id = chainid2.substr( (i-26)/100, 1 );
            roff = 10*(i%100);
        }

        clipper::MPolymer mprot = mprep;
        int ir = results[i].rot;
        int it = results[i].trn;
        clipper::RTop_orth rtop( rots[ir].rot(),
                                 xwrk.coord_orth( grid.deindex(it).coord_map() ) );
        mprot.transform( rtop );
        mprot.set_id( id );

        for ( int r = 0; r < mprot.size(); r++ )
            mprot[r].set_seqnum( roff+mprot[r].seqnum() );

        mol_new.insert( mprot );

    }

    return mol_new;
}

const privateer::data::fingerprint privateer::data::fingerprint_list[] =
{
    {
        "ARA", "ligand", "A", "L", "alpha-L-Arabinose", 10,
        { { "AO5", 2.826, 0.443, 0.013 }, { "AC5", 2.182, 0.000,-1.193 }, { "AC4", 0.772, 0.418,-1.248 },
          { "AO4", 0.631, 1.605,-1.628 }, { "AC1", 2.206, 0.000, 1.183 }, { "AO1", 2.875, 0.588, 2.298 },
          { "AC2", 0.731, 0.512, 1.249 }, { "AO2", 0.131, 0.300, 2.188 }, { "AC3", 0.000, 0.000, 0.000 }, { "AO3",-1.125, 0.610, 0.003 } } ,
        { { "V01", 2.758, 2.629, 3.794 }, { "V02", 2.490, 5.743,-5.634 }, { "V03", 0.870, 1.227, 5.015 },
          { "V04",-5.376, 2.851, 3.476 }, { "V05",-3.243, 3.161, 4.321 }, { "V06",-1.485,-2.860, 2.153 },
          { "V07",-3.001,-0.086,-6.081 }, { "V08", 2.412,-2.704,-0.215 }, { "V09", 3.461, 0.901,-3.155 }, { "V10",-0.204, 4.505,-2.072 } } ,
        { { "C1",  2.178, 0.000, 1.206 }, { "C2",  0.733, 0.508, 1.236 }, { "C3",  0.000, 0.000, 0.000 },
          { "C4",  0.728, 0.462,-1.271 }, { "C5",  2.172, 0.000,-1.203 }, { "O1",  2.874, 0.596, 2.298 },
          { "O2",  0.088, 0.115, 2.439 }, { "O3", -1.357, 0.457,-0.031 }, { "O4",  0.635, 1.868,-1.478 }, { "O5", 2.817, 0.460, 0.003 } }
    },
    {
        "NAG", "nglycan", "B", "D", "N-Acetyl-D-Glucosamine", 3,
        { { "P", 0.0, 0.0, 0.0 }, { "P", 0.0, 0.0, 0.0 }, { "P", 0.0, 0.0, 0.0 } } ,
        { { "V", 0.0, 0.0, 0.0 }, { "V", 0.0, 0.0, 0.0 }, { "V", 0.0, 0.0, 0.0 } } ,
        { { "C", 0.0, 0.0, 0.0 }, { "C", 0.0, 0.0, 0.0 }, { "C", 0.0, 0.0, 0.0 } }
    },

    {
        "BMA", "nglycan", "B", "D", "beta-D-mannopyranose", 3,
        { { "P", 0.0, 0.0, 0.0 }, { "P", 0.0, 0.0, 0.0 }, { "P", 0.0, 0.0, 0.0 } } ,
        { { "V", 0.0, 0.0, 0.0 }, { "V", 0.0, 0.0, 0.0 }, { "V", 0.0, 0.0, 0.0 } } ,
        { { "C", 0.0, 0.0, 0.0 }, { "C", 0.0, 0.0, 0.0 }, { "C", 0.0, 0.0, 0.0 } }
    }
};



