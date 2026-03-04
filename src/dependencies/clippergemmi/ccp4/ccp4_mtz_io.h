/*! \file lib/mtz_io.h
    Header file for reflection data mtz importer
*/
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


#ifndef CLIPPER_CCP4_MTZ_IO
#define CLIPPER_CCP4_MTZ_IO


#include "ccp4_mtz_types.h"


namespace clipper
{

  //! CCP4MTZ import/export type registry
  /*! This class acts as a registry of import/export data types, which
    provide a translation from a data element name to an MTZ column
    type.

    The registry is instantiated statically as \c mtz_type_registry
    and initialised with a list of built-in datatypes.
  */

  class CCP4MTZ_type_registry
  {
  public:
    //! constructor: initialise the registry with some built-in types
    CCP4MTZ_type_registry();
    //! add a new type to the registry
    static void add_type( const String& name, const String& type, const ftype32& scale );
    //! return MTZ column type
    static String type( const String& name );
    //! return MTZ column scale
    static ftype32 scale( const String& name );
    //! add a new group to the registry
    static void add_group( const String& name, const String& type );
    //! return MTZ group type
    static String group_type( const String& name );
  private:
    static char names[200][12];
    static char types[200][4];
    static ftype32 scales[200];
    static char groupnames[200][20];
    static char grouptypes[200][8];
  };




  //! MTZ import/export parent class for clipper objects
  /*! This is the import/export class which can be linked to an mtz
    file and be used to transfer data into or out of a Clipper
    data structure.

    Note that to access the MTZ file efficiently, data reads and
    writes are deferred until the file is closed.

    \anchor MTZpaths \par MTZpaths: MTZ column specification

    Note that the specification of the MTZ column names is quite
    versatile. The MTZ crystal and dataset must be specified, although
    the wildcard '*' may replace a complete name. Several MTZ columns
    will correspond to a single datalist. This may be handled in two
    ways:
    -# A simple name. The corresponding MTZ columns will be named
    after the datalist name, a dot, the datalist type, a dot, and a
    type name for the indivudal column,
    i.e. /crystal/dataset/datalist.listtype.coltype This is the
    default Clipper naming convention for MTZ data.
    -# A comma separated list of MTZ column names enclosed in square
    brackets.  This allows MTZ data from legacy applications to be
    accessed. \n
    Thus, for example, an MTZPATH of
    \code
    native/CuKa/fsigfdata
    \endcode
    expands to MTZ column names of
    \code
    fsigfdata.F_sigF.F
    fsigfdata.F_sigF.sigF
    \endcode
    with a crystal called \c native and a dataset called \c CuKa. An MTZPATH of
    \code
    native/CuKa/[FP,SIGFP]
    \endcode
    expands to MTZ column names of
    \code
    FP
    SIGFP
    \endcode
    with a crystal called \c native and a dataset called \c CuKa.

    \archor MTZ_iotypes \par MTZ_iotypes: Import/export types

    For an HKL_data object to be imported or exported, an MTZ_iotype
    for that datatype must exist in the
    MTZ_iotypes_registry. MTZ_iotypes are defined for all the built-in
    datatypes. If you need to store a user defined type in an MTZ
    file, then register that type with the MTZ_iotypes_registry. */
  class CCP4MTZfile
  {
   public:
    enum CCP4MTZcolumn_label_mode { Default, Normal, Legacy };

    //! Constructor: does nothing
    CCP4MTZfile();
    //! Destructor: close any file that was left open
    ~CCP4MTZfile();

    //! Open a file for read access
    void open_read( const String filename_in );
    //! Close a file after reading
    void close_read();
    //! Open a file for read access
    void open_append( const String filename_in, const String filename_out );
    //! Close a file after reading
    void close_append();
    //! Open a file for read access
    void open_write( const String filename_out );
    //! Close a file after reading
    void close_write();

    //! get file spacegroup
    const Spacegroup& spacegroup() const;
    //! get file cell
    const Cell& cell() const;
    //! get file resolution
    const Resolution& resolution() const;
    //! get file HKL sampling
    const HKL_sampling& hkl_sampling() const;

    //! read the reflection list from the MTZ
    void import_hkl_list( HKL_info& target );
    //! import parameters of HKL_info object from the MTZ
    void import_hkl_info( HKL_info& target, const bool generate = true );
    //! import crystal info from the MTZ
    void import_crystal( MTZcrystal& cxtl, const String mtzpath );
    //! import dataset info from the MTZ
    void import_dataset( MTZdataset& cset, const String mtzpath );
    //! mark a hkl_data for import from MTZ
    void import_hkl_data( HKL_data_base& cdata, const String mtzpath );

    //! write the reflection list to the MTZ (along with cell, spacegroup)
    void export_hkl_info( const HKL_info& target );
    //! export crystal info from the MTZ
    void export_crystal( const MTZcrystal& cxtl, const String mtzpath );
    //! export dataset info from the MTZ
    void export_dataset( const MTZdataset& cset, const String mtzpath );
    //! mark a hkl_data for export to MTZ
    void export_hkl_data( const HKL_data_base& cdata, const String mtzpath );

    //! mark a chkl_data container for import from MTZ
    void import_chkl_data( Container& target, const String mtzpath, const String path = "" );
    //! mark a chkl_data container for export to MTZ
    void export_chkl_data( Container& target, const String mtzpath );

    //! return a vector of strings of the file column names/path
    std::vector<String> column_paths() const;
    //! return a vector of strings of column names/paths just allocated
    const std::vector<String>& assigned_paths() const;

    //! get title for MTZ file
    String title() const;
    //! set title for MTZ file
    void set_title( const String& title );
    //! get history for MTZ file
    std::vector<String> history() const;
    //! add to history for MTZ file
    void set_history( const std::vector<String>& history );

    //! return number of reflections in file
    int num_reflections() const;
    //! return sort order
    std::vector<int> sort_order() const;
    //! return low resolution limits from file (A)
    ftype32 low_res_limit() const;
    //! return high resolution limits from file (A)
    ftype32 high_res_limit() const;
    //! CCP4 space group number
    int ccp4_spacegroup_number() const;
    //! space group confidence
    char spacegroup_confidence() const;
    //! set space group confidence
    void set_spacegroup_confidence(const char& spg_confidence);

    //! set default column label interpretation mode
    void set_column_label_mode( CCP4MTZcolumn_label_mode mode );
    //! set header verbosity level
    void set_verbose( int verbose );

    //! \deprecated
    std::vector<String> column_labels() const;
    //! \deprecated
    void import_hkl_data( HKL_data_base& cdata, MTZdataset& cset, MTZcrystal& cxtl, const String mtzpath );
    //! \deprecated
    void export_hkl_data( const HKL_data_base& cdata, const MTZdataset& cset, const MTZcrystal& cxtl, const String mtzpath );


    // index from clipper data lists to mtz columns (by name)
    struct datacolinf{String label,type,source,grpname,grptype; int grpposn; };
    struct datasetinf{MTZdataset dataset;std::vector<datacolinf> columns;};
    struct crystalinf{MTZcrystal crystal;std::vector<datasetinf> datasets;};
    struct hkldatacol{String path; ftype scale;};
  private:
    enum MTZmode { NONE, READ, WRITE, APPEND };

    //! mtz object
    String filename_in_, filename_out_;
    std::vector<crystalinf> crystals;
    HKL_info*       hkl_info_i;
    const HKL_info* hkl_info_o;
    std::vector<HKL_data_base*      > hkl_data_i;
    std::vector<const HKL_data_base*> hkl_data_o;
    std::vector<std::vector<hkldatacol> > hkl_data_cols;
    std::vector<String> assigned_paths_;
    //! mode
    MTZmode mode;
    CCP4MTZcolumn_label_mode colmode;
    int verbose_;

    //! File spacegroup, cell, resolution
    Spacegroup spacegroup_;
    Cell cell_;
    Resolution resolution_;
    HKL_sampling hkl_sampling_;
    //! title, history
    String title_;
    std::vector<String> history_;
    //! Number of reflections
    int num_reflections_;
    //! Sort order
    std::vector<int> sort_order_;
    //! resolution limits from file header
    ftype32 reslim_lo, reslim_hi;
    //! CCP4 space group number
    int ccp4_spacegroup_number_;
    //! space group confidence
    char spacegroup_confidence_;

    // generic methods
    bool match_path( const String& path, int& x, int& s, int& c );
  };


} // namespace clipper

#endif
