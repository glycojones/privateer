// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

#ifndef PRIVATEER_LIB_H_INCLUDED
#define PRIVATEER_LIB_H_INCLUDED

#include <fstream>
#include <iomanip>
#include <algorithm>
#include <math.h>
#include "clipper-glyco.h"
#include <clipper/clipper.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include <clipper/clipper-minimol.h>
#include <ccp4srs/ccp4srs_manager.h>
#include <ccp4srs/ccp4srs_defs.h>
#include "privateer-restraints.h"

typedef clipper::HKL_data_base::HKL_reference_index HRI;

inline const std::string b2s ( bool b )
{
  return b ? "yes" : "no";
}

namespace privateer
{
    namespace coot
    {
        // Coot support, Scheme
        void insert_coot_prologue_scheme ( std::fstream& );
        void insert_coot_epilogue_scheme ( std::fstream& );
        void insert_coot_files_loadup_scheme ( std::fstream&, const clipper::String&, const clipper::String&, const clipper::String&, const clipper::String&, bool mode, const clipper::String& pdbblobs, bool blobsoutput);
        void insert_coot_go_to_blob_scheme ( std::fstream& output, const clipper::Coord_orth& blob_centre, const clipper::String& diagnostic );
        void insert_coot_go_to_sugar_scheme ( std::fstream&, const clipper::Coord_orth& sugar_centre, const clipper::String& diagnostic );
        void insert_coot_statusbar_text_scheme ( std::fstream&, clipper::String& );
        void insert_coot_command ( std::fstream& output, std::string command );

        // Coot support, Python
        void insert_coot_files_loadup_python ( std::fstream&, const clipper::String&, const clipper::String&, const clipper::String&, const clipper::String&, bool mode, const clipper::String& pdbblobs, bool blobsoutput);
        void insert_coot_prologue_python ( std::fstream& );
        void insert_coot_epilogue_python ( std::fstream& );
        void insert_coot_go_to_blob_python ( std::fstream& output, const clipper::Coord_orth& blob_centre, const clipper::String& diagnostic );
        void insert_coot_go_to_sugar_python ( std::fstream&, const clipper::Coord_orth& sugar_centre, const clipper::String& diagnostic );
        void insert_coot_statusbar_text_python ( std::fstream&, clipper::String& );
    }

    namespace util
    {
        bool calculate_sigmaa_maps (const clipper::Atom_list& list_of_atoms,
                                    const clipper::HKL_data<clipper::data32::F_sigF>& reflection_data,
                                    const clipper::HKL_data<clipper::data32::F_phi>& simulated_cryoem_reflection_data,
                                    clipper::Xmap<float>& best_map,
                                    clipper::Xmap<float>& difference_map,
                                    bool ignore_set_null,
                                    bool useMTZ, 
                                    int n_refln = 1000,
                                    int n_param = 20);
        bool check_if_matom_initialized(clipper::MAtom& input);
        void print_usage();
        void print_supported_code_list ();
        char get_altconformation(clipper::MAtom ma);
        void write_refmac_keywords ( std::vector < std::string > code_list );
        bool write_libraries ( std::vector < std::string > code_list, float esd = 5.0 );
        bool compute_and_print_external_validation ( const std::vector<clipper::String> validation_options,
                                                     clipper::data::sugar_database_entry& external_validation );
        void print_XML ( std::vector < std::pair < clipper::String, clipper::MSugar > > sugarList,
                         std::vector < clipper::MGlycan > list_of_glycans,
                         std::vector<std::vector<std::pair<std::pair<clipper::MGlycan, std::vector<int>>,float>>>& list_of_glycans_associated_to_permutations,
                         clipper::String pdbname,
                         std::vector<privateer::json::GlycomicsDatabase>& glycomics_database );
        bool read_coordinate_file_mtz (clipper::MMDBfile& mfile, clipper::MiniMol& mmol, clipper::String& ippdb, bool batch);
        bool read_coordinate_file_mrc (clipper::MMDBfile& mfile, clipper::MiniMol& mmol, clipper::String& ippdb, clipper::Xmap<double>& input_map, bool batch);
        clipper::Xmap<float> read_map_file ( std::string path );
        int find_index_of_value_from_wurcs ( std::vector<privateer::json::GlycomicsDatabase>& glycomics_database, std::string inputwurcs );
        std::vector<char> number_of_conformers ( clipper::MMonomer& mmon );
        void print_monosaccharide_summary (bool batch, bool showGeom, int pos_slash, bool useMRC, std::vector<std::pair<clipper::String, clipper::MSugar>>& ligandList, FILE *output, clipper::HKL_info& hklinfo, clipper::String input_model);
        void print_monosaccharide_summary_python (bool batch, bool showGeom, int pos_slash, bool useMRC, std::vector<std::pair<clipper::String, clipper::MSugar>>& ligandList, clipper::HKL_info& hklinfo, clipper::String input_model);
        float calculate_linkage_zscore(float phi, float psi, privateer::json::TorsionsZScoreDatabase& matched_linkage);
        std::string retrieve_input_PDB_code(clipper::String input_model_path);
    }

    namespace glycanbuilderplot
    {
        enum Colour { blue, rootblue, red, rootred, yellow, rootyellow, orange, green, purple, pink, cyan, tan, black, brown, white };
        enum Link_type { up, up_side, branch_side, side, down_side, down };

        std::string get_colour ( Colour colour, bool original_style, bool inverted = false  );

        inline const std::string get_svg_tooltip ( clipper::MSugar sugar, bool validation )
        {
            std::ostringstream str;
            str << std::setprecision(2) << std::fixed
                << "<tspan>" << sugar.type() << " "+sugar.id().trim() << " in "
                << sugar.conformation_name() << " conformation. </tspan>"
                << "<tspan>Mean B-factor: " << sugar.get_bfactor() << ". </tspan>"
                << "<tspan>Detected type: " << sugar.type_of_sugar() << ". </tspan>";
            if ( validation )
            {
                if ( sugar.ok_with_conformation() && sugar.ok_with_anomer() &&
                     sugar.ok_with_chirality() && sugar.ok_with_puckering() )
                    str << "<tspan>No issues have been detected.</tspan>";
                else
                {
                    str << "Detected issues: ";

                    if ( ! sugar.ok_with_anomer() )
                        str << "<tspan>Anomer mismatch; </tspan>";
                    if ( ! sugar.ok_with_chirality() )
                        str << "<tspan>Absolute configuration mismatch; </tspan>";
                    if ( ! sugar.ok_with_puckering() )
                        str << "<tspan>Unlikely ring puckering amplitude; </tspan>";
                    if ( ! sugar.ok_with_conformation() )
                        str << "<tspan>The sugar is in a high-energy conformation; </tspan>";
                }
             }

            return str.str();
        }

        class Shape
        {
            public:
                Shape() { }
                Shape( int x, int y ) { pos_x = x; pos_y = y; }
                virtual ~Shape() { };
                void set_pos( int x, int y ) { pos_x = x; pos_y =y ; }
                int  get_y  ( ) { return pos_y; }
                int  get_x  ( ) { return pos_x; }
                std::string get_id() { return svg_id; }
                virtual std::string get_XML ( ) = 0;
                void set_tooltip ( std::string tooltip ) { this->tooltip = tooltip;  }
                std::string get_tooltip ( ) { return this->tooltip; }
                void set_mmdbsel ( std::string mmdbsel ) { this->mmdbsel = mmdbsel;  }
                std::string get_mmdbsel ( ) { return this->mmdbsel; }

            protected:
                int pos_x;
                int pos_y;
                std::string svg_id;
                std::string tooltip;
                std::string mmdbsel;
        };


        class Plot
        {
            public:
                Plot()
                {
                    this->width  = 800;
                    this->height = 300;

                    viewbox.push_back(0);
                    viewbox.push_back(0);
                    viewbox.push_back(width);
                    viewbox.push_back(height);

                    this->vertical=false;
                    this->original_colour_scheme = true;
                    this->title = "Figure generated with Privateer";
                    this->inverted_background = false;
                    this->validation = true;
                    this->add_links = true;
                } //!< default constructor
                

                Plot( bool vertical, bool original_style, std::string title, bool inverted_background = false, bool smaller=false, bool validation = true, bool add_links = true )
                {
                    vertical ? width  = 300 : width  = 800;
                    vertical ? height = 800 : height = 300;

                    viewbox.push_back(0);
                    viewbox.push_back(0);
                    viewbox.push_back(width);
                    viewbox.push_back(height);

                    this->smaller = smaller;
                    this->vertical=vertical;
                    this->original_colour_scheme = original_style;
                    this->title = title;
                    this->inverted_background = inverted_background;
                    this->validation = validation;
                    this->add_links = add_links;
                }

                /*! Set the SVG viewport
                 * 	\param viewbox A std::vector<int> containing, in this order: ( x0, y0, x1, y1 )
                 * 	\return False if there haven't been any problems
                 */

                bool set_viewbox ( std::vector<int>& viewbox )
                {
                    if ( viewbox.size() != 4 )
                        return true;

                    this->viewbox = viewbox;

                    return false;
                }


                /*! Swap blacks and whites - only useful if you're Jon Agirre?
                 * 	\return true if blacks and whites have been swapped
                 */

                bool invert_background ( )
                {
                    inverted_background ? inverted_background = false : inverted_background = true;
                    return inverted_background;
                } //!< swap blacks and whites - only useful if you're Jon Agirre?

                /*! Get the SVG viewbox
                 * 	\return A std::vector<int> containing, in this order: ( x0, y0, x1, y1 )
                 */

                const std::string get_viewbox ( )
                {
                    std::stringstream str;
                    str << viewbox[0] << " "
                        << viewbox[1] << " "
                        << viewbox[2] << " "
                        << viewbox[3];
                    return str.str();
                }

                void tighten_viewbox ( );

                int get_width  ( ) { return width;  }
                int get_height ( ) { return height; }


                /*! Set the SVG's width and height.
                 * 	\param width Image width
                 * 	\param height Image height
                 */

                void set_size ( int width, int height )
                {
                    this->width = width; this->height = height;
                }

                /*! Set the SVG's width, height, and adjust viewbox accordingly.
                 * 	\param width Image width
                 * 	\param height Image height
                 */

                void set_size_and_viewbox ( int width, int height )
                {
                    this->width = width; this->height = height;
                    viewbox.clear();
                    viewbox.push_back(0);
                    viewbox.push_back(0);
                    viewbox.push_back(width);
                    viewbox.push_back(height);
                }

                void set_title ( std::string title ) { this->title = title; }
                std::string get_title () { return this->title; }
                bool plot_glycan ( clipper::MGlycan glycan );
                bool plot_demo ( ); //!< creates a demo plot with all the blocks, links and roots Privateer can generate
                bool write_to_file  ( std::string file_path ); //!< returns true if there have been any problems
                std::string write_to_string ( );
                std::string get_XML ();
                void delete_shapes ( ) { for (int i=0; i<list_of_shapes.size();i++) delete list_of_shapes[i]; }

                std::string get_svg_string_header      ( );
                std::string get_svg_string_contents    ( );
                std::string get_svg_string_footer      ( );

            private:
                int width;
                int height;
                std::string title;
                bool vertical;
                bool smaller;
                bool original_colour_scheme;
                bool inverted_background;
                bool validation;
                bool add_links;

                void write_svg_header        ( std::fstream& of );
                void write_svg_definitions   ( std::fstream& of );
                void write_svg_contents      ( std::fstream& of );
                void write_svg_footer        ( std::fstream& of );
                // Need to rework this bit so that ostringstream would brand into fstream etc.
                void write_svg_header_ostringstream        ( std::ostringstream& of );
                void write_svg_definitions_ostringstream   ( std::ostringstream& of );
                void write_svg_contents_ostringstream      ( std::ostringstream& of );
                void write_svg_footer_ostringstream        ( std::ostringstream& of );

                void recursive_paint ( clipper::MGlycan mg, clipper::MGlycan::Node node, int x, int y, bool oxford_angles = false );

                // we want to draw linkages first, so that the blocks are then drawn on top of them

                void add_block ( Shape * block ) { list_of_shapes.push_back ( block );}
                void add_link  ( Shape * link  ) { list_of_shapes.insert ( list_of_shapes.begin(), link ); }
                std::vector<int> viewbox;
                std::vector < Shape * > list_of_shapes;

        };


        class Square : public Shape
        {
            public:
                Square() { }
                Square( int x, int y, int width, int height ) { set_pos(x, y); this->width=width; this->height=height; }
                virtual ~Square() {};
                int  get_width  ( ) { return width;  }
                int  get_height ( ) { return height; }
                void set_size ( int w, int h ) { width=w; height=h; }
                virtual std::string get_XML ( ) = 0;

            protected:
                int width;
                int height;
                int width_border;
        };


        class Diamond : public Shape
        {
            public:
                Diamond() { } //!< null constructor
                Diamond( int x, int y, int width, int height ) { set_pos(x, y); this->width=width; this->height=height; } //!< constructor
                virtual ~Diamond() {};
                int  get_width  ( ) { return width;  }
                int  get_height ( ) { return height; }
                void set_size ( int w, int h ) { width=w; height=h; }
                virtual std::string get_XML ( ) = 0;

            protected:
                int width;
                int height;
                int width_border;
        };


        class Star : public Shape
        {
            public:
                Star() { } //!< null constructor
                Star( int x, int y, int width, int height ) { set_pos(x, y); this->width=width; this->height=height; } //!< constructor
                virtual ~Star() {};
                int  get_width  ( ) { return width;  }
                int  get_height ( ) { return height; }
                void set_size ( int w, int h ) { width=w; height=h; }
                virtual std::string get_XML ( ) = 0;

            protected:
                int width;
                int height;
                int width_border;
        };


        class Hexagon : public Shape
        {
            public:
                Hexagon() { } //!< null constructor
                Hexagon( int x, int y, int width, int height ) { set_pos(x, y); this->width=width; this->height=height; } //!< constructor
                virtual ~Hexagon() {};
                int  get_width  ( ) { return width;  }
                int  get_height ( ) { return height; }
                void set_size ( int w, int h ) { width=w; height=h; }
                virtual std::string get_XML ( ) = 0;

            protected:
                int width;
                int height;
                int width_border;
        };


        class Triangle : public Shape
        {
            public:
                Triangle() { } //!< null constructor
                Triangle( int x, int y, int side ) { set_pos(x, y); this->side=side; } //!< constructor
                virtual ~Triangle() {};
                void set_side ( int s ) { side=s; }
                int get_side  ( ) { return side;  }
                virtual std::string get_XML ( )=0;

            protected:
                int side;
                int width_border;
        };



        class Circle : public Shape
        {
            public:
                Circle() { } //!< null constructor
                Circle( int x, int y, int radius ) { set_pos(x, y); this->radius=radius; } //!< constructor
                virtual ~Circle() {};
                void set_radius ( int r ) { radius=r; }
                int  get_radius  ( ) { return radius;  }
                virtual std::string get_XML ( ) = 0;

            protected:
                int radius;
                int width_border;
        };


        // Standard conventions for Glycobiology

        // hexoses

        class Glc : public virtual Circle
        {
            public:
                Glc() { } //!< null constructor
                Glc( int x, int y, std::string message, std::string mmdbsel = "" ) { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel);}
                std::string get_XML ( );

        };

        class Gal : public virtual Circle
        {
            public:
                Gal() { } //!< null constructor
                Gal( int x, int y, std::string message, std::string mmdbsel = "" ) { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel); }
                std::string get_XML ( );

        };

        class Man : public virtual Circle
        {
            public:
                Man() { } //!< null constructor
                Man( int x, int y, std::string message, std::string mmdbsel = "" ) { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel);}
                std::string get_XML ( );

        };

        class Fuc : public virtual Triangle
        {
            public:
                Fuc() { } //!< null constructor
                Fuc( int x, int y, std::string message, std::string mmdbsel = "" ) { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel);}
                std::string get_XML ( );

        };

        class Xyl : public virtual Star
        {
            public:
                Xyl() { } //!< null constructor
                Xyl( int x, int y, std::string message,  std::string mmdbsel = "") { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel); }
                std::string get_XML ( );

        };


        // hexosamines

        class GlcN : public virtual Square
        {
            public:
                GlcN() { } //!< null constructor
                GlcN( int x, int y, std::string message, std::string mmdbsel = "") { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel); }
                std::string get_XML ( );

        };

        class GalN : public virtual Square
        {
            public:
                GalN() { } //!< null constructor
                GalN( int x, int y, std::string message, std::string mmdbsel = "") { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel); }
                std::string get_XML ( );

        };

        class ManN : public virtual Square
        {
            public:
                ManN() { } //!< null constructor
                ManN( int x, int y, std::string message, std::string mmdbsel = "") { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel); }
                std::string get_XML ( );

        };


        // N-acetyl hexosamines

        class GlcNAc : public virtual Square
        {
            public:
                GlcNAc() { } //!< null constructor
                GlcNAc( int x, int y, std::string message, std::string mmdbsel = "") { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel); }
                std::string get_XML ( );

        };

        class GalNAc : public virtual Square
        {
            public:
                GalNAc() { } //!< null constructor
                GalNAc( int x, int y, std::string message, std::string mmdbsel = "") { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel); }
                std::string get_XML ( );

        };

        class ManNAc : public virtual Square
        {
            public:
                ManNAc() { } //!< null constructor
                ManNAc( int x, int y, std::string message, std::string mmdbsel = "") { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel); }
                std::string get_XML ( );

        };


        // acidic sugars

        class Neu5Ac : public virtual Diamond
        {
            public:
                Neu5Ac() { } //!< null constructor
                Neu5Ac( int x, int y, std::string message, std::string mmdbsel = "") { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel); }
                std::string get_XML ( );

        };

        class Neu5Gc : public virtual Diamond
        {
            public:
                Neu5Gc() { } //!< null constructor
                Neu5Gc( int x, int y, std::string message, std::string mmdbsel = "") { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel); }
                std::string get_XML ( );

        };

        class KDN : public virtual Diamond
        {
            public:
                KDN() { } //!< null constructor
                KDN( int x, int y, std::string message, std::string mmdbsel = "") { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel); }
                std::string get_XML ( );

        };

        class GlcA : public virtual Diamond
        {
            public:
                GlcA() { } //!< null constructor
                GlcA( int x, int y, std::string message, std::string mmdbsel = "") { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel); }
                std::string get_XML ( );

        };

        class IdoA : public virtual Diamond
        {
            public:
                IdoA() { } //!< null constructor
                IdoA( int x, int y, std::string message, std::string mmdbsel = "") { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel); }
                std::string get_XML ( );

        };

        class GalA : public virtual Diamond
        {
            public:
                GalA() { } //!< null constructor
                GalA( int x, int y, std::string message, std::string mmdbsel = "") { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel); }
                std::string get_XML ( );

        };

        class ManA : public virtual Diamond
        {
            public:
                ManA() { } //!< null constructor
                ManA( int x, int y, std::string message, std::string mmdbsel = "") { set_pos(x, y); set_tooltip ( message ); set_mmdbsel(mmdbsel); }
                std::string get_XML ( );

        };

        class Unk : public virtual Hexagon
        {
            public:
                Unk() { } //!< null constructor
                Unk( int x, int y, const char letter, std::string message, std::string mmdbsel = "") { set_pos(x, y); set_tooltip ( message ); code += letter;  set_mmdbsel(mmdbsel);}
                std::string get_XML ( );
            private:
                std::string code;

        };

        // bond types

        // class AlphaBond : public virtual Shape
        // {
        //     public:
        //         AlphaBond() { } //!< null constructor
        //         AlphaBond( int x, int y, Link_type bond, std::string message, std::string mmdbsel = "" ) { set_pos(x, y); set_tooltip ( message ); set_bond_type (bond); this->set_mmdbsel ( mmdbsel );}
        //         std::string get_XML ( );

        //     private:
        //         Link_type bond_type;
        //         void set_bond_type ( Link_type bond ) { this->bond_type = bond; }
        //         Link_type get_bond_type () { return bond_type; }
        // };

        class Bond : public virtual Shape
        {
            public:
                Bond() { } //!< null constructor
                Bond( int x, int y, Link_type bond, std::string message, std::string mmdbsel = "" ) { set_pos(x, y); set_tooltip ( message ); set_bond_type (bond); this->set_mmdbsel ( mmdbsel ); }
                Bond( int x, int y, std::string anomerSymbol, Link_type bond, std::string message, std::string mmdbsel = "" ) { set_pos(x, y); set_tooltip ( message ); set_bond_type (bond); this->anomerSymbol = anomerSymbol; this->set_mmdbsel ( mmdbsel ); }
                Bond( int x, int y, Link_type bond, std::string anomerSymbol, std::string linkagePosition, std::string message, std::string mmdbsel = "" ) { set_pos(x, y); set_tooltip ( message ); set_bond_type (bond); this->anomerSymbol = anomerSymbol; this->linkagePosition = linkagePosition; this->set_mmdbsel ( mmdbsel ); }
                std::string get_XML ( );

            private:
                std::string anomerSymbol;
                std::string linkagePosition;
                Link_type bond_type;
                void set_bond_type ( Link_type bond ) { this->bond_type = bond; }
                Link_type get_bond_type () { return bond_type; }
        };


        // glycan root

        class GlycanRoot : public virtual Shape
        {
            public:
                GlycanRoot() {} //!< null constructor

                GlycanRoot( int x, int y, std::string link_atom, std::string root_name, std::string root_id, std::string message, std::string mmdbsel="" )
                {
                    set_pos(x, y);
                    set_tooltip ( message );
                    set_link_atom ( link_atom );
                    set_root_name ( root_name );
                    set_root_id   ( root_id   );
                    set_mmdbsel ( mmdbsel );
                }

                std::string get_XML ( );

                void set_link_atom ( std::string name ) { this->link_atom = name; }
                std::string get_link_atom ( ) { return this->link_atom; }

                void set_root_name ( std::string name ) { this->root_name = name; }
                std::string get_root_name ( ) { return this->root_name; }

                void set_root_id ( std::string id ) { this->root_id = id; }
                std::string get_root_id ( ) { return this->root_id; }

            private:
                std::string link_atom;
                std::string root_name;
                std::string root_id;
        };
    }

    namespace scripting
    {
        inline std::string carbname_of ( std::string name ) { return clipper::data::carbname_of ( name ); }
        inline bool found_in_database ( std::string name ) { return clipper::data::found_in_database ( name ); }
        std::string get_annotated_glycans ( std::string pdb_filename, bool original_colour_scheme = false, std::string expression_system = "undefined" );
        std::string get_annotated_glycans_hierarchical ( std::string pdb_filename, bool original_colour_scheme = false, std::string expression_system = "undefined"  );
        std::string print_wurcs ( std::string pdb_filename, std::string expression_system = "undefined");
        std::string print_node ( const clipper::MiniMol& mmol, const clipper::MGlycan& mg, const clipper::MGlycan::Node& node, const std::string chain, const clipper::MGlycan::Linkage& connection );
        void svg_graphics_demo ( bool original_colour_scheme, bool inverted_background = false );
        void compute_linkage_torsion_zscores_for_glycan(std::vector<privateer::json::TorsionsZScoreDatabase>& torsions_zscore_database, std::vector<clipper::MGlycan::MGlycanTorsionSummary>& glycan_torsions); 
        inline void write_refmac_keywords ( std::vector < std::string > code_list ) { return privateer::util::write_refmac_keywords(code_list); }
        inline bool write_libraries ( std::vector < std::string > code_list, float esd ) { return privateer::util::write_libraries(code_list, esd); }
    }

} // namespace privateer

#endif
