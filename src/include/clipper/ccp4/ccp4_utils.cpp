/* mtz_utils.cpp: ccp4 utils for the clipper libraries */
//C Copyright (C) 2000-2006 Kevin Cowtan and University of York
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


#include "ccp4_utils.h"

#include <ccp4/ccp4_general.h>
#include <ccp4/ccp4_program.h>

#include <iostream>
#include <cstdlib>


CCP4CommandInput::CCP4CommandInput( int argc, char** argv, bool echo )
{
  for ( int arg = 0; arg < argc; arg++ ) {
    std::string thisarg( argv[arg] );
    if ( thisarg.length() > 2 )
      if ( thisarg[0] == '-' && thisarg[1] == '-' )
	thisarg = thisarg.substr(1);
    if ( thisarg == "-stdin" ) {
      std::string line;
      while ( !std::getline(std::cin,line).eof() )
	if ( line.length() > 0 )
	  if ( line[0] != '#' ) {
	    int i = line.find_first_not_of( " \t" );
	    if ( i != std::string::npos ) {
	      line = line.substr( i );
	      i = line.find_first_of( " \t" );
	      std::string word = line.substr( 0, i );
	      if ( word[0] != '-' ) word = "-" + word;
	      push_back( word );
	      if ( i != std::string::npos ) {
		line = line.substr( i );
		i = line.find_first_not_of( " \t" );
		if ( i != std::string::npos ) {
		  line = line.substr( i );
		  i = line.find_last_not_of( " \t" );
		  push_back( line.substr( 0, i+1 ) );
		}
	      }
	    }
	  }
    } else {
      push_back( thisarg );
    }
  }
  // output version number
  if ( size() == 2 ) {
    if ( (*this)[1] == "-i" ) {
      CCP4::ccp4_prog_info();
      std::exit(0);
    }
  }
  // echo output
  if ( echo ) {
    for ( int i = 1; i < size(); i++ ) {
      char c1 = ' ', c2 = ' ';  // spot keywords vs negative numbers
      if ( (*this)[i].length() > 0 ) c1 = (*this)[i][0];
      if ( (*this)[i].length() > 1 ) c2 = (*this)[i][1];
      if ( c1 == '-' && ( c2 < '0' || c2 > '9' ) )
	std::cout << std::endl << (*this)[i].substr(1);  // keywords on newline
      else
	std::cout << " \t"     << (*this)[i];            // args on sameline
    }
    std::cout << std::endl;
  }
}


CCP4Program::CCP4Program( const char* name, const char* vers, const char* rcsdate )
{
  name_ = name;
  html = ( getenv( "CCP_SUPPRESS_HTML" ) == NULL ); 
  summ = ( getenv( "CCP_SUPPRESS_SUMMARY" ) == NULL ); 
  CCP4::ccp4ProgramName( (char*)name );
  CCP4::ccp4_prog_vers( (char*)vers );
  CCP4::ccp4RCSDate( (char*)rcsdate );
  summary_beg();
  if ( html ) std::cout << "<html> <!-- CCP4 HTML LOGFILE -->" << std::endl
			<< "<hr>" << std::endl << "<pre>" << std::endl;
  CCP4::ccp4_banner();
  summary_end();
  CCP4::ccp4ProgramTime(1);
}


CCP4Program::~CCP4Program()
{
  std::cout << std::endl;
  summary_beg();
  std::cout << name_ << ": " << msg_ << std::endl;
  CCP4::ccp4ProgramTime(0);
  if ( html ) std::cout << "</pre>" << std::endl << "</html>" << std::endl;
  summary_end();
}


void CCP4Program::summary_beg() const
{
  if ( summ ) {
    if ( html )
      std::cout << "<B><FONT COLOR='#FF0000'><!--SUMMARY_BEGIN-->" << std::endl;
    else
      std::cout << "<!--SUMMARY_BEGIN-->" << std::endl;
  }
}


void CCP4Program::summary_end() const
{
  if ( summ ) {
    if ( html )
      std::cout << "<!--SUMMARY_END--></FONT></B>" << std::endl;
    else
      std::cout << "<!--SUMMARY_END-->" << std::endl;
  }
}
