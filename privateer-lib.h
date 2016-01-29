
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

} // namespace privateer

#endif
