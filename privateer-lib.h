
// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2016 Jon Agirre & Kevin Cowtan
// York Structural Biology Laboratory
// The University of York
// mailto: jon.agirre@york.ac.uk
// mailto: kevin.cowtan@york.ac.uk
//

#ifndef PRIVATEER_LIB_H_INCLUDED
#define PRIVATEER_LIB_H_INCLUDED

#include <fstream>
#include <algorithm>
#include "clipper-glyco.h"
#include <clipper/clipper.h>
#include <clipper/clipper-mmdb.h>
#include <clipper/clipper-minimol.h>


namespace privateer
{
    namespace coot
    {
        // Coot support, Scheme
        void insert_coot_prologue_scheme ( std::fstream& );
        void insert_coot_epilogue_scheme ( std::fstream& );
        void insert_coot_files_loadup_scheme ( std::fstream&, const clipper::String&, const clipper::String&, const clipper::String&, const clipper::String&, bool mode );
        void insert_coot_go_to_sugar_scheme ( std::fstream&, const clipper::Coord_orth& sugar_centre, const clipper::String& diagnostic );
        void insert_coot_statusbar_text_scheme ( std::fstream&, clipper::String& );
        void insert_coot_command ( std::fstream& output, std::string command );

        // Coot support, Python
        void insert_coot_files_loadup_python ( std::fstream&, const clipper::String&, const clipper::String&, const clipper::String&, const clipper::String&, bool mode );
        void insert_coot_prologue_python ( std::fstream& );
        void insert_coot_epilogue_python ( std::fstream& );
        void insert_coot_go_to_sugar_python ( std::fstream&, const clipper::Coord_orth& sugar_centre, const clipper::String& diagnostic );
        void insert_coot_statusbar_text_python ( std::fstream&, clipper::String& );
    }
    
    float real_space_correlation ( const clipper::Xmap<float>&, const clipper::Xmap<float>& );

    
    namespace glycoplot
    {
        enum Colour { blue, red, yellow, orange, green, purple, cyan, tan, black };
        enum Link_type { up, up_side, side, down_side, down };
        
        static std::string get_colour ( Colour colour, bool original_style );
        
        class Shape
        {
            public:
                Shape() { } //!< null constructor
                Shape( int x, int y ) { pos_x = x; pos_y = y; } //!< constructor
                void set_pos( int x, int y ) { pos_x = x; pos_y =y ; }
                int  get_y  ( ) { return pos_y; }
                int  get_x  ( ) { return pos_x; }
                std::string get_id() { return svg_id; }
                virtual std::string get_XML ( ) = 0;
                void set_tooltip ( std::string tooltip ) { this->tooltip = tooltip;  }
                std::string get_tooltip ( ) { return this->tooltip; }
            
            protected:
                int pos_x;
                int pos_y;
                std::string svg_id;
                std::string tooltip;
        };
        
        
        class Plot
        {
            public:
                Plot()
                {
                    this->width  = 800;
                    this->height = 300;
                    this->vertical=false;
                    this->original_colour_scheme = true;
                    this->title = "Figure generated with Privateer";
                } //!< default constructor
            
                Plot( bool vertical, bool bigger, bool original_style, std::string title )
                {
                    bigger ? width  = 800 : width  = 266;
                    bigger ? height = 300 : height = 100;
                    this->vertical=vertical;
                    this->original_colour_scheme = original_style;
                    this->title = title;
                }
            
                int get_width  ( ) { return width;  }
                int get_height ( ) { return height; }
                void set_title ( std::string title ) { this->title = title; }
                std::string get_title () { return this->title; }
                bool plot_glycan ( clipper::MGlycan glycan );
                bool write_to_file  ( std::string file_path ); //!< returns true if there have been any problems
                std::string get_XML ();
            
            private:
                int width;
                int height;
                std::string title;
                bool vertical;
                bool original_colour_scheme;
                void write_svg_header        ( std::fstream& of );
                void write_svg_definitions   ( std::fstream& of );
                void write_svg_contents      ( std::fstream& of );
                void write_svg_footer        ( std::fstream& of );
            
                std::string get_svg_string_header      ( );
                std::string get_svg_string_definitions ( );
                std::string get_svg_string_contents    ( );
                std::string get_svg_string_footer      ( );
            
                void add_shape ( Shape * shape ) { list_of_shapes.push_back ( shape );}
                std::vector < Shape * > list_of_shapes;
                
        };
        
        
        class Square : public Shape
        {
            public:
                Square() { } //!< null constructor
                Square( int x, int y, int width, int height ) { set_pos(x, y); this->width=width; this->height=height; } //!< constructor
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
                void set_radius ( int r ) { radius=r; }
                int  get_radius  ( ) { return radius;  }
                virtual std::string get_XML ( ) = 0;
            
            protected:
                int radius;
                int width_border;
        };
        
        
        // Standard conventions for Glycobiology
        
        class GlcNAc : public virtual Square
        {
            public:
                GlcNAc() { } //!< null constructor
                GlcNAc( int x, int y, std::string message) { set_pos(x, y); set_tooltip ( message ); }
                std::string get_XML ( );
            
        };
        
        class Man : public virtual Circle
        {
            public:
                Man() { } //!< null constructor
                Man( int x, int y, std::string message ) { set_pos(x, y); set_tooltip ( message ); }
                std::string get_XML ( );
                
        };
        
        class Fuc : public virtual Triangle
        {
            public:
                Fuc() { } //!< null constructor
                Fuc( int x, int y, std::string message ) { set_pos(x, y); set_tooltip ( message ); }
                std::string get_XML ( );
                
        };
        
        class AlphaBond : public virtual Shape
        {
            public:
                AlphaBond() { } //!< null constructor
                AlphaBond( int x, int y, Link_type bond, std::string message ) { set_pos(x, y); set_tooltip ( message ); set_bond_type (bond);}
                std::string get_XML ( );
            
            private:
                Link_type bond_type;
                void set_bond_type ( Link_type bond ) { this->bond_type = bond; }
                Link_type get_bond_type () { return bond_type; }
        };
        
        class BetaBond : public virtual Shape
        {
            public:
                BetaBond() { } //!< null constructor
                BetaBond( int x, int y, Link_type bond, std::string message ) { set_pos(x, y); set_tooltip ( message ); set_bond_type (bond); }
                std::string get_XML ( );
            
            private:
                Link_type bond_type;
                void set_bond_type ( Link_type bond ) { this->bond_type = bond; }
                Link_type get_bond_type () { return bond_type; }
        };
    }
    
} // namespace privateer

#endif
