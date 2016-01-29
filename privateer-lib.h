
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
    
    clipper::ftype real_space_correlation ( const clipper::Xmap<float>&, const clipper::Xmap<float>& );

    
    namespace glycoplot
    {
        
        class Shape
        {
            public:
                Shape() { } //!< null constructor
                Shape( int x, int y ) { pos_x = x; pos_y = y; } //!< constructor
                void set_pos( int x, int y ) { pos_x = x; pos_y =y ; }
                int  get_y  ( ) { return pos_y; }
                int  get_x  ( ) { return pos_x; }
                void set_tooltip ( std::string tooltip ) { this->tooltip = tooltip; }
            
            private:
                int pos_x;
                int pos_y;
                std::string tooltip;
        };
        
        class Plot
        {
            public:
                Plot() {} //!< null constructor
                Plot( int width, int height, bool vertical ) { this->width=width; this->height=height; this->vertical=vertical; }
                bool plot_to_file ( std::fstream& ); //!< returns true if there have been any problems
            private:
                int width;
                int height;
                bool vertical;
                bool write_svg_header   ( std::fstream& );
                bool write_svg_contents ( std::fstream& );
                bool write_svg_footer   ( std::fstream& );
                
                std::vector < Shape > list_of_shapes;
                
        };
        
        class Square : public Shape
        {
            public:
                Square() { } //!< null constructor
                Square( int x, int y, int width, int height ) { set_pos(x, y); this->width=width; this->height=height; } //!< constructor
                int  get_width  ( ) { return width;  }
                int  get_height ( ) { return height; }
                std::string paint ( );
            
            private:
                int width;
                int height;
                std::string colour_fill;
                std::string colour_border;
                int width_border;
        };
        
        class Circle : public Shape
        {
            public:
                Circle() { } //!< null constructor
                Circle( int x, int y, int radius ) { set_pos(x, y); this->radius=radius; } //!< constructor
                int  get_radius  ( ) { return radius;  }
                std::string paint ( );
                
            private:
                int radius;
                std::string colour_fill;
                std::string colour_border;
                int width_border;
        };
        
        class Triangle : public Shape
        {
            public:
                Triangle() { } //!< null constructor
                Triangle( int x, int y, int side ) { set_pos(x, y); this->side=side; } //!< constructor
                int get_side  ( ) { return side;  }
                std::string paint ( );
                
            private:
                int side;
                std::string colour_fill;
                std::string colour_border;
                int width_border;
        };
    }
    
} // namespace privateer

#endif
