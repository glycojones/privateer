/*
     overview.h: overview of CINCH - Crystallography IN C Headers
     Copyright (C) 2003  CCLRC, Martyn Winn

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
*/
/** @mainpage CINCH - Crystallography IN C Headers
 *
 * @note I'm fed up of the uninspiring "new library" and here's
 *       my first attempt at a funky name. Alternatives welcome!
 *
 * @verbatim

<!-- ::INDEX_INFO::C libraries::Library::::CCP4 C/C++ Software Library:::::::: -->

   @endverbatim
 *
 * @section aims Aims

The CCP4 software suite is based around a library of routines which
cover common tasks, such as file opening, parsing keyworded input,
reading and writing of standard data formats, applying symmetry
operations, etc.  Programs in the suite call these routines which, as
well as saving the programmer some effort, ensure that the varied
programs in the suite have a similar look-and-feel. 
<p>
Since 2002, there has been a major effort to re-write
much of the CCP4 library into C/C++. The aims are:

<ul>
<li>To implement a better representation of the underlying data model.
For example, Eugene Krissinel's MMDB library acts on a data structure
which represents the various levels of structure of a protein model.
The new MTZ library encapsulates the crystal/dataset hierarchy that
is increasingly being used by programs.
<li>To maintain support for existing programs. In particular, the 
existing Fortran APIs will be maintained, although they will now often
be only wrappers to functions in the new library. It is hoped that many
existing programs will be migrated to using the new library directly.
<li>To provide support for scripting. It is possible to generate APIs 
for Python, Tcl and Perl automatically from the core C code. Thus, much
of the standard CCP4 functionality wil be available to scripts used
e.g. in ccp4i or the molecular graphics project.
</ul>

This incremental approach, maintaining the existing suite while
improving the underlying code, puts constraints on what is possible, but
is considered more appropriate for a collaborative project like CCP4.

 * @section start This documentation

<p>
This documentation is generated automatically by 
<a href="http://www.doxygen.org/">Doxygen</a> from
comment sections in the code. It is therefore detailed and extensive. 
The library divides roughly into the following sections:
<dl>
<dt>CMTZ library
<dd>See the @ref cmtz_page page for C/C++ programmers, and the
@ref cmtz_f_page page for Fortran programmers.
<dt>CMAP library
<dd>See the @ref cmap_page page for C/C++ programmers, and the
@ref cmap_f_page page for Fortran programmers.
<dt>MMDB library
<dd>See Eugene's <a href="http://www.ebi.ac.uk/~keb/cldoc">documentation</a>.
<dt>CSYM library
<dd>See the @ref csym_page page for C/C++ programmers, and the
@ref csym_f_page page for Fortran programmers.
<dt>CCP4 utility library
<dd>See the @ref utilities_page page for C/C++ programmers.
<dt>Low level disk i/o
<dd>See the @ref diskio_f_page page for Fortran programmers.
</dl>

 */

