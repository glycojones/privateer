/* mtz_io.cpp: class file for reflection data  mtz importer */
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

#include "ccp4_mtz_io.h"

#include "../core/hkl_datatypes.h"

#include <ccp4/cmtzlib.h>
#include <ccp4/csymlib.h>
#include <ccp4/ccp4_parser.h>


namespace clipper {


char CCP4MTZ_type_registry::names[200][12];
char CCP4MTZ_type_registry::types[200][4];
ftype32 CCP4MTZ_type_registry::scales[200];
char CCP4MTZ_type_registry::groupnames[200][20];
char CCP4MTZ_type_registry::grouptypes[200][8];
CCP4MTZ_type_registry mtz_type_registry;


CCP4MTZ_type_registry::CCP4MTZ_type_registry()
{
  for ( int j = 0; j < 200; j++ ) names[j][0] = '\0';
  add_type( "I", "J", 1.0 );
  add_type( "sigI", "Q", 1.0 );
  add_type( "I+", "K", 1.0 );
  add_type( "sigI+", "M", 1.0 );
  add_type( "I-", "K", 1.0 );
  add_type( "sigI-", "M", 1.0 );
  add_type( "covI+-", "C", 1.0 );
  add_type( "F", "F", 1.0 );
  add_type( "sigF", "Q", 1.0 );
  add_type( "F+", "G", 1.0 );
  add_type( "sigF+", "L", 1.0 );
  add_type( "F-", "G", 1.0 );
  add_type( "sigF-", "L", 1.0 );
  add_type( "covF+-", "C", 1.0 );
  add_type( "E", "F", 1.0 );
  add_type( "sigE", "Q", 1.0 );
  add_type( "E+", "G", 1.0 );
  add_type( "sigE+", "L", 1.0 );
  add_type( "E-", "G", 1.0 );
  add_type( "sigE-", "L", 1.0 );
  add_type( "covE+-", "C", 1.0 );
  add_type( "A", "A", 1.0 );
  add_type( "B", "A", 1.0 );
  add_type( "C", "A", 1.0 );
  add_type( "D", "A", 1.0 );
  add_type( "phi", "P", Util::rad2d(1.0) );
  add_type( "fom", "W", 1.0 );
  add_type( "flag", "I", 1.0 );
  add_type( "Dano", "D", 1.0 );
  add_type( "sigDano", "Q", 1.0 );
  for ( int j = 0; j < 200; j++ ) groupnames[j][0] = '\0';
  add_group( "I_sigI", "ISIG" );
  add_group( "I_sigI_ano", "IANO" );
  add_group( "F_sigF", "FSIG" );
  add_group( "F_sigF_ano", "FANO" );
  add_group( "E_sigE", "ESIG" );
  add_group( "E_sigE_ano", "EANO" );
  add_group( "F_phi", "FPHI" );
  add_group( "Phi_fom", "PHIW" );
  add_group( "ABCD", "ABCD" );
  add_group( "D_sigD", "DSIG" );
  add_group( "Flag", "FREE" );
}

void CCP4MTZ_type_registry::add_type( const String& name, const String& type, const ftype32& scale )
{
  int i, j;
  for ( j = 0; j < 200; j++ )
    if ( names[j][0] == '\0' ) break;
  if ( j == 200 )
    Message::message( Message_fatal( "CCP4MTZ_type_registry: registry full." ) );
  for ( i = 0; i < Util::min( int(name.length()), 11 ); i++ )
    names[j][i] = name[i];
  names[j][i] = '\0';
  for ( i = 0; i < Util::min( int(type.length()), 3 ); i++ )
    types[j][i] = type[i];
  types[j][i] = '\0';
  scales[j] = scale;
}

String CCP4MTZ_type_registry::type( const String& name )
{
  int j;
  for ( j = 0; j < 200; j++ )
    if ( String( names[j] ) == name ) break;
  if ( j == 200 )
    Message::message( Message_fatal( "CCP4MTZ_type_registry: name not found." ) );
  return String( types[j] );
}

ftype32 CCP4MTZ_type_registry::scale( const String& name )
{
  int j;
  for ( j = 0; j < 200; j++ )
    if ( String( names[j] ) == name ) break;
  if ( j == 200 )
    Message::message( Message_fatal( "CCP4MTZ_type_registry: name not found." ) );
  return scales[j];
}

void CCP4MTZ_type_registry::add_group( const String& name, const String& type )
{
  int i, j;
  for ( j = 0; j < 200; j++ )
    if ( groupnames[j][0] == '\0' ) break;
  if ( j == 200 )
    Message::message( Message_fatal( "CCP4MTZ_type_registry: registry full." ) );
  for ( i = 0; i < Util::min( int(name.length()), 19 ); i++ )
    groupnames[j][i] = name[i];
  groupnames[j][i] = '\0';
  for ( i = 0; i < Util::min( int(type.length()), 7 ); i++ )
    grouptypes[j][i] = type[i];
  grouptypes[j][i] = '\0';
}

String CCP4MTZ_type_registry::group_type( const String& name )
{
  int j;
  for ( j = 0; j < 200; j++ )
    if ( String( groupnames[j] ) == name ) break;
  if ( j == 200 )
    Message::message( Message_fatal( "CCP4MTZ_type_registry: name not found." ) );
  return String( grouptypes[j] );
}



/* return the mtz path of a column as a string */
String mtz_col_path( CMtz::MTZ* mtz, CMtz::MTZCOL* col )
{
  char* resultchr = CMtz::MtzColPath( mtz, col );
  String resultstr( resultchr );
  free(resultchr);
  return resultstr;
}

/* build a lookup of columns in a crystal */
std::vector<CMtz::MTZCOL*> build_lookup( CMtz::MTZ* mtzin )
{
  std::vector<CMtz::MTZCOL*> result;
  for (int x=0; x < CMtz::MtzNxtal(mtzin); x++) {
    CMtz::MTZXTAL* xtl = CMtz::MtzIxtal(mtzin,x);
    for (int s=0; s < CMtz::MtzNsetsInXtal(xtl); s++) {
      CMtz::MTZSET* set = CMtz::MtzIsetInXtal(xtl,s);
      for (int c=0; c < CMtz::MtzNcolsInSet(set); c++)
        result.push_back( CMtz::MtzIcolInSet(set,c) );
    }
  }
  return result;
}

/* Turn a clipper path name into a list of MTZ column paths. */
const std::vector<String> mtz_assign( const String& assign, const String& type, const String& ftype, const int& f_size )
{
  // Translate aassignment to mtz column names
  // ... loop through each list in turn
  std::vector<String> mtz_names( f_size, "MNF" );
  // interpret list name in terms of columns
  if ( assign.find( "[" ) == String::npos ) {
    // name is a single path: add type extenders
    std::vector<String> list = ftype.split(" ");
    for (int i=0; i < list.size(); i++) mtz_names[i] = assign + "." + type + "." + list[i];
  } else {
    // name is a list of mtz columns: extract column names from list
    String pref = assign.substr( 0, assign.find("[") );
    String post = assign.substr( assign.find("[") + 1,
				 assign.find("]") - assign.find("[") - 1 );
    std::vector<String> list = post.split(", ");
    for (int i=0; i < list.size(); i++) mtz_names[i] = pref + list[i];
  }
  return mtz_names;
}

/* Check whether a dummy input column name has been given to suppress
  input or output of some of the data. */
bool is_virtual_col( const String& path )
{
  // test for virtual column name
  String name = path.tail();
  return ( name=="MNF" || name=="NAN" || name=="mnf" || name=="nan" );
}

/* Do glob matching on two strings. The second string may contain '*'. */
bool match_glob( const String& s1, const String& s2 )
{
  if ( s2.find('*') == String::npos ) return s1 == s2;
  String s2a = s2.substr( 0, s2.find('*') );
  String s2b = s2.substr( s2.find('*') + 1 );
  return ( s2a == s1.substr(                          0, s2a.length() ) &&
	   s2b == s1.substr( s1.length() - s2b.length(), s2b.length() ) );
}

/* Read spacegroup info from mtzin */
void read_spacegroup( CMtz::MTZ* mtzin, Spacegroup& sg )
{
  ftype rsym[4][4];
  String symops;
  for ( int i = 0; i < mtzin->mtzsymm.nsym; i++ ) {
    for ( int j = 0; j < 4; j++ )
      for ( int k = 0; k < 4; k++ )
        rsym[j][k] = mtzin->mtzsymm.sym[i][j][k];
    symops += Symop(rsym).format() + ";";
  }
  sg.init( Spgr_descr( symops, Spgr_descr::Symops ) );
}

/* Write spacegroup info to mtzout */
void write_spacegroup( CMtz::MTZ* mtzout, const Spacegroup& sg,
		       const char& spacegroup_confidence )
{
  // tables of MTZ symbols
  char mtzlauetab[231][8]={"?","1","-1","2","2","2","m","m","m","m","2/m","2/m","2/m","2/m","2/m","2/m","222","222","222","222","222","222","222","222","222","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mm2","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","mmm","4","4","4","4","4","4","-4","-4","4/m","4/m","4/m","4/m","4/m","4/m","422","422","422","422","422","422","422","422","422","422","4mm","4mm","4mm","4mm","4mm","4mm","4mm","4mm","4mm","4mm","4mm","4mm","-4m2","-4m2","-4m2","-4m2","-42m","-42m","-42m","-42m","-42m","-42m","-4m2","-4m2","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","4/mmm","3","3","3","3","-3","-3","312","321","312","321","312","321","32","3m1","31m","3m1","31m","3m","3m","-31m","-31m","-3m1","-3m1","-3m","-3m","6","6","6","6","6","6","-6","6/m","6/m","622","622","622","622","622","622","6mm","6mm","6mm","6mm","-62m","-62m","-6m2","-6m2","6/mmm","6/mmm","6/mmm","6/mmm","23","23","23","23","23","m-3","m-3","m-3","m-3","m-3","m-3","m-3","432","432","432","432","432","432","432","432","-43m","-43m","-43m","-43m","-43m","-43m","m-3m","m-3m","m-3m","m-3m","m-3m","m-3m","m-3m","m-3m","m-3m","m-3m"};
  unsigned int ccp4_code[][2] = { {0x90a34743,1003},{0x7851c0d1,1004},{0x458c0b8c,1005},{0x939daf66,1006},{0xd9a7d85a,1007},{0xec8b9a56,1008},{0x5b2a5e1b,1009},{0x9b92779f,1010},{0x31194672,1011},{0x97e84b5c,1012},{0xc9e50c74,1013},{0xdc48ad2f,1014},{0xb0d2ac85,1015},{0xd756eb1b,1017},{0x852190f0,1018},{0x3879e795,1020},{0xcf830dc8,1021},{0x2d514a9f,1022},{0x56ae76fe,1023},{0xadb55f32,1059},{0xb5c8e336,1094},{0xd5a0aa2d,1146},{0xd9a29bac,1148},{0xa20b8591,1155},{0xb951b4f7,1160},{0x219be015,1161},{0x1c80e47a,1166},{0xbb691c91,1167},{0xf4722972,1197},{0x1e4dbca6,2005},{0x15463be2,2014},{0x010e6701,2017},{0x82570fbb,2018},{0x8aa408da,3005},{0xee8326bb,3018},{0x6042abb6,4005},{0xd7e36ffb,5005} };
  String mtzlaue = String("PG")+mtzlauetab[sg.spacegroup_number()];
  String mtzsymb = sg.symbol_hm();
  if ( sg.symbol_hm_ext() == "H" || sg.symbol_hm_ext() == "R" )
    mtzsymb = sg.symbol_hm_ext() + mtzsymb.substr(1);
  int mtzspgn = sg.spacegroup_number();
  for ( int s = 0; s < sizeof(ccp4_code)/sizeof(ccp4_code[0]); s++ )
    if ( ccp4_code[s][0] == sg.hash() ) mtzspgn = ccp4_code[s][1];
  // now write the records
  mtzout->mtzsymm.spcgrp = mtzspgn;
  mtzout->mtzsymm.nsym   = sg.num_symops();
  mtzout->mtzsymm.nsymp  = sg.num_primops();
  mtzout->mtzsymm.symtyp = mtzsymb[0];
  mtzout->mtzsymm.spg_confidence = spacegroup_confidence;
  strncpy( mtzout->mtzsymm.spcgrpname, mtzsymb.c_str(), 11 );
  strncpy( mtzout->mtzsymm.pgname, mtzlaue.c_str(), 11 );
  for ( int i = 0; i < sg.num_symops(); i++ ) {
    for ( int j = 0; j < 3; j++ )
      for ( int k = 0; k < 3; k++ )
	mtzout->mtzsymm.sym[i][j][k] = sg.symop(i).rot()(j,k);
    for ( int j = 0; j < 3; j++ )
      mtzout->mtzsymm.sym[i][j][3] = sg.symop(i).trn()[j];
  }

  /*
  CSym::ccp4_symop symops[sg.num_symops()];
  for ( int i = 0; i < sg.num_symops(); i++ ) {
    for ( int j = 0; j < 3; j++ )
      for ( int k = 0; k < 3; k++ )
	symops[i].rot[j][k] = sg.symop(i).rot()(j,k);
    for ( int j = 0; j < 3; j++ )
      symops[i].trn[j] = sg.symop(i).trn()[j];
  }
  CSym::CCP4SPG* cspg =
    CSym::ccp4_spgrp_reverse_lookup( sg.num_symops(), symops );
  if ( cspg != NULL ) {
    float msymops[192][4][4] = {{{}}};
    for ( int i = 0; i < cspg->nsymop; i++ )
      CCP4::rotandtrn_to_mat4( msymops[i], cspg->symop[i] );
    char symbol[20];
    strncpy( symbol, cspg->symbol_xHM, 20 );
    CSym::ccp4spg_name_de_colon( symbol );
    char lattice[2] = { symbol[0], '\0' };
    int err = ccp4_lwsymm( mtzout,
			   cspg->nsymop, cspg->nsymop_prim, msymops,
			   lattice,
			   cspg->spg_ccp4_num,
			   symbol, cspg->point_group );
    if ( err == 0 ) Message::message( Message_fatal( "CCP4MTZfile::close_write: Unable to write spacegroup info. lwsymm failed." ) );
    CSym::ccp4spg_free( &cspg );
  } else {
    Message::message( Message_fatal( "CCP4MTZfile::close_write: Unable to write spacegroup info.\nCCP4 not configured or unsupported spacegroup." ) );
  }
  */
}

/* Read hierarchy from mtzin */
void read_hierarchy( CMtz::MTZ* mtzin, std::vector<CCP4MTZfile::crystalinf>& crystals )
{
  crystals.clear();
  CCP4MTZfile::crystalinf newxtl;
  CCP4MTZfile::datasetinf newset;
  CCP4MTZfile::datacolinf newcol;
  for (int x=0; x < CMtz::MtzNxtal(mtzin); x++) {
    CMtz::MTZXTAL* xtl = CMtz::MtzIxtal(mtzin,x);
    newxtl.crystal = MTZcrystal( xtl->xname, xtl->pname, Cell( Cell_descr( xtl->cell[0], xtl->cell[1], xtl->cell[2], xtl->cell[3], xtl->cell[4], xtl->cell[5] ) ) );
    crystals.push_back( newxtl );
    for (int s=0; s < CMtz::MtzNsetsInXtal(xtl); s++) {
      CMtz::MTZSET* set = CMtz::MtzIsetInXtal(xtl,s);
      newset.dataset = MTZdataset( set->dname, set->wavelength );
      crystals.back().datasets.push_back( newset );
      for (int c=0; c < CMtz::MtzNcolsInSet(set); c++) {
        const CMtz::MTZCOL* col = CMtz::MtzIcolInSet(set,c);
	newcol.label = col->label;
	newcol.type  = col->type;
	newcol.source  = col->colsource;
	newcol.grpname = col->grpname;
	newcol.grptype = col->grptype;
	newcol.grpposn = col->grpposn;
	crystals.back().datasets.back().columns.push_back( newcol );
      }
    }
  }
}

/* Write columns to mtzout */
void write_hierarchy( CMtz::MTZ* mtzout, std::vector<CCP4MTZfile::crystalinf>& crystals )
{
  for (int x = 0; x < crystals.size(); x ++ ) {
    const MTZcrystal& cxtl = crystals[x].crystal;
    float cp[6];
    cp[0] = cxtl.a(); cp[3] = cxtl.alpha_deg();
    cp[1] = cxtl.b(); cp[4] = cxtl.beta_deg();
    cp[2] = cxtl.c(); cp[5] = cxtl.gamma_deg();
    CMtz::MTZXTAL* mxtl = CMtz::MtzAddXtal( mtzout, cxtl.crystal_name().c_str(),
					    cxtl.project_name().c_str(), cp );
    for (int s = 0; s < crystals[x].datasets.size(); s ++ ) {
      const MTZdataset& cset = crystals[x].datasets[s].dataset;
      CMtz::MTZSET* mset = CMtz::MtzAddDataset( mtzout, mxtl,
						cset.dataset_name().c_str(),
						cset.wavelength() );
      for (int c = 0; c < crystals[x].datasets[s].columns.size(); c ++ ) {
	CCP4MTZfile::datacolinf& ccol = crystals[x].datasets[s].columns[c];
	CMtz::MTZCOL* mcol = CMtz::MtzAddColumn( mtzout, mset,
						 ccol.label.c_str(),
						 ccol.type.c_str() );
	mcol->source = (ccol.source==" !NEW! ") ? 0 : 1;
	strcpy( mcol->colsource, ccol.source.substr(0,36).c_str() );
	strcpy( mcol->grpname, ccol.grpname.substr(0,30).c_str() );
	strcpy( mcol->grptype, ccol.grptype.substr(0,4).c_str() );
	mcol->grpposn = ccol.grpposn;
      }
    }
  }
}

/* Build a column cross-reference */
void reference_cols( CMtz::MTZ* mtz, std::vector<CMtz::MTZCOL*> lookup, std::vector<std::vector<CCP4MTZfile::hkldatacol> > hkl_data_cols, std::vector<std::vector<int> >& cols, std::vector<std::vector<ftype> >& scls )
{
  // Lookup the input columns
  int nlst = hkl_data_cols.size();
  cols.clear(); scls.clear();
  cols.resize( nlst ); scls.resize( nlst );
  // handle each data list in turn. Each may have several mtz cols
  for ( int lst = 0; lst < nlst; lst++ ) {
    int ncols = hkl_data_cols[lst].size();
    cols[lst].resize( ncols, -1 );
    scls[lst].resize( ncols, 1.0 );
    // assign the columns to mtz indexes
    for ( int col = 0; col < ncols; col++ )  // loop over columns in list
      if ( !is_virtual_col( hkl_data_cols[lst][col].path ) ) {
	int mcol;  // loop over cols in MTZ
	for ( mcol = 0; mcol < lookup.size(); mcol++ )
	  if ( hkl_data_cols[lst][col].path ==
	       mtz_col_path( mtz, lookup[mcol] ) ) break;
	if ( mcol == lookup.size() )
	  Message::message( Message_fatal( "CCP4MTZfile - internal error" ) );
	cols[lst][col] = mcol;  // store column index
	scls[lst][col] = hkl_data_cols[lst][col].scale;
      }
  }
}

/*! Look up the crystal, dataset and column corresponding to the given
  mtz path. True if exact column match found. */
bool CCP4MTZfile::match_path( const String& path, int& x, int& s, int& c )
{
  std::vector<String> names = path.split("/");
  for ( x = 0; x < crystals.size(); x ++ )
   if ( match_glob(crystals[x].crystal.crystal_name(),names[0]) )
    for ( s = 0; s < crystals[x].datasets.size(); s ++ )
     if ( match_glob(crystals[x].datasets[s].dataset.dataset_name(),names[1] ) )
      for ( c = 0; c < crystals[x].datasets[s].columns.size(); c ++ )
       if ( match_glob(crystals[x].datasets[s].columns[c].label,names[2]) )
	return true;
  c = -1;
  for ( x = 0; x < crystals.size(); x ++ )
   if ( match_glob(crystals[x].crystal.crystal_name(),names[0]) )
    for ( s = 0; s < crystals[x].datasets.size(); s ++ )
     if ( match_glob(crystals[x].datasets[s].dataset.dataset_name(),names[1] ) )
      return false;
  s = -1;
  for ( x = 0; x < crystals.size(); x ++ )
   if ( match_glob(crystals[x].crystal.crystal_name(),names[0]) )
    return false;
  x = -1;
  return false;
}

/*! Constructing an CCP4MTZfile does nothing except flag the object as not
  attached to any file for either input or output */
CCP4MTZfile::CCP4MTZfile()
{
  mode = NONE;
  colmode = Default;
  verbose_ = 0;
}

/*! Close any files which were left open. This is particularly
  important since to access the MTZ file efficiently, data reads and
  writes are deferred until the file is closed. */
CCP4MTZfile::~CCP4MTZfile()
{
  switch ( mode ) {
  case READ:
    close_read(); break;
  case WRITE:
    close_write(); break;
  case APPEND:
    close_append(); break;
  case NONE:
    break;
  }
}

/*! The file is opened for reading. This CCP4MTZfile object will remain
  attached to this file until it is closed. Until that occurs, no
  other file may be opened with this object, however another CCP4MTZfile
  object could be used to access another file.
  \param filename_in The input filename or pathname. */
void CCP4MTZfile::open_read( const String filename_in )
{
  if ( mode != NONE )
    Message::message( Message_fatal( "CCP4MTZfile: open_read - File already open" ) );

  // read the mtz
  filename_in_ = filename_in;

  // open file
  CMtz::MTZ* mtzin = CMtz::MtzGet( filename_in_.c_str(), 0 );
  if ( mtzin == NULL ) Message::message( Message_fatal( "CCP4MTZfile: open_read - File missing or corrupted: "+filename_in_ ) );
  CMtz::MtzAssignHKLtoBase( mtzin );
  // get the list of column names/paths
  read_hierarchy( mtzin, crystals );
  int x, s, c;
  match_path( "/HKL_base/HKL_base/H", x, s, c );
  if ( x < 0 ) Message::message( Message_fatal( "CCP4MTZfile: No cell! " ) );
  cell_ = crystals[x].crystal;
  // get spacegroup by decoding symops
  read_spacegroup( mtzin, spacegroup_ );
  // get resolution:
  float minres, maxres;
  CMtz::MtzResLimits( mtzin, &minres, &maxres );
  reslim_hi = 0.9999 / sqrt( Util::max(minres,maxres) );
  if ( Util::min(minres,maxres) > 0.0 )
    reslim_lo =  ( 1.0001 / sqrt( Util::min(minres,maxres) ) );
  else
    reslim_lo = 1000000.0;
  resolution_.init( reslim_hi );
  hkl_sampling_ = HKL_sampling( cell_, resolution_ );
  // get title , history
  title_ = String( mtzin->title, 70 );
  history_ = std::vector<String>( mtzin->histlines );
  for ( int h = 0; h < history_.size(); h++ )
    history_[h] = String( mtzin->hist + MTZRECORDLENGTH*h, MTZRECORDLENGTH );
  num_reflections_ = CMtz::MtzNref( mtzin );
  // Sort order
  int isort[5];
  CMtz::ccp4_lrsort( mtzin, isort );
  sort_order_.assign(isort, isort+5);
  ccp4_spacegroup_number_ = mtzin->mtzsymm.spcgrp;
  spacegroup_confidence_  = mtzin->mtzsymm.spg_confidence;

  // output
  if ( verbose_ >= 1 ) {
    printf( "MTZ file: %s\n", filename_in_.c_str() );
  }
  if ( verbose_ >= 2 ) {
    printf( "  Number of reflections:%10i        Number of columns:%4i\n",
	    CMtz::MtzNref( mtzin ), CMtz::MtzNcol( mtzin ) );
  }
  if ( verbose_ >= 3 ) {
    printf( "  Resolution range: %10.3f - %10.3f\n",
	    reslim_lo, reslim_hi );
    printf( "  Spacegroup: %-24s  (number %i) %c\n", mtzin->mtzsymm.spcgrpname,
	    mtzin->mtzsymm.spcgrp, mtzin->mtzsymm.spg_confidence );
  }

  // close file
  CMtz::MtzFree( mtzin );

  mode = READ;
}

/*! Close the file after reading. This command also actually fills in
  the data in any HKL_data structures which have been marked for
  import. */
void CCP4MTZfile::close_read()
{
  if ( mode != READ )
    Message::message( Message_fatal( "CCP4MTZfile: no file open for read" ) );

  // fast path for reading header only
  if ( hkl_data_i.size() == 0 ) {
    // these should be redundent
    hkl_data_i.clear(); hkl_data_cols.clear(); crystals.clear(); assigned_paths_.clear();
    // tidy up
    mode = NONE;
    filename_in_ = filename_out_ = "";
    return;
  }

  // close the mtz:
  // this actually imports the data from the mtz file
  // - we save up all the work to do it on a single pass

  HKL hkl;
  float fdata_i[1000]; // file data: allow up to 1000 columns on file
  int   idata_i[1000]; // flag data: allow up to 1000 columns on file
  xtype pdata[100];  // program data: allow up to 100 floats per group type
  float res;
  int ref, lst, col;

  // open file
  CMtz::MTZ* mtzin = CMtz::MtzGet( filename_in_.c_str(), 0 );
  if ( mtzin == NULL ) Message::message( Message_fatal( "CCP4MTZfile: close_read - File missing or corrupted: "+filename_in_ ) );
  CMtz::MtzAssignHKLtoBase( mtzin );
  std::vector<CMtz::MTZCOL*> lookup_i = build_lookup( mtzin );

  // lookup hkl cols
  int hklcol_i[3] = {0,1,2};
  for ( int i = 0; i < lookup_i.size(); i++ ) {
    if ( mtz_col_path( mtzin, lookup_i[i] ).tail() == "H" ) hklcol_i[0] = i;
    if ( mtz_col_path( mtzin, lookup_i[i] ).tail() == "K" ) hklcol_i[1] = i;
    if ( mtz_col_path( mtzin, lookup_i[i] ).tail() == "L" ) hklcol_i[2] = i;
  }

  // Lookup the input columns
  std::vector<std::vector<int> > cols;
  std::vector<std::vector<ftype> > scls;
  reference_cols( mtzin, lookup_i, hkl_data_cols, cols, scls );

  // update the data lists to ensure size consistency
  for ( lst = 0; lst < hkl_data_i.size(); lst++ ) hkl_data_i[lst]->update();

  // Import the data
  for ( ref = 0; ref < CMtz::MtzNref( mtzin ); ref++ ) {
    // read reflection
    CMtz::ccp4_lrreff( mtzin, &res, fdata_i, idata_i, (const CMtz::MTZCOL**)
		       &lookup_i[0], lookup_i.size(), ref+1 );
    hkl = HKL( Util::intr(fdata_i[hklcol_i[0]]),
	       Util::intr(fdata_i[hklcol_i[1]]),
	       Util::intr(fdata_i[hklcol_i[2]]) );
    // do each data list
    for ( lst = 0; lst < cols.size(); lst++ ) {
      // for each list, grab the relevent columns
      for ( col = 0; col < cols[lst].size(); col++ ) {
	// set to NaN unless readable and present
	Util::set_null( pdata[col] );
	//pdata[col] = Util::nanf();
	// try and read the data
	if ( cols[lst][col] >= 0 )
	  if ( !idata_i[cols[lst][col]] )
	    pdata[col] = xtype( fdata_i[cols[lst][col]] / scls[lst][col] );
      }
      // now set the data
      hkl_data_i[lst]->data_import( hkl, pdata );
    }
  }

  // tidy up
  CMtz::MtzFree( mtzin );

  hkl_data_i.clear();
  hkl_data_cols.clear();
  crystals.clear();
  assigned_paths_.clear();
  filename_in_ = filename_out_ = "";

  mode = NONE;
}

/*! The file is opened for writing. This will be a new file, created
  entirely from data from within the program, rather than by extending
  an existing file. Similar restrictions apply as for open_read().

  In practice the open_append() method is usually preferred.
  \param filename_out The output filename or pathname. */
void CCP4MTZfile::open_write( const String filename_out )
{
  if ( mode != NONE )
    Message::message( Message_fatal( "CCP4MTZfile: open_write - File already open" ) );

  // make the h,k,l columns to the local list
  crystalinf newxtl; datasetinf newset; datacolinf newcol;
  newxtl.crystal = MTZcrystal( "HKL_base", "HKL_base", Cell() );
  newset.dataset = MTZdataset( "HKL_base", 9.999 );
  newcol.label = newcol.type = "H";
  // add columns
  crystals.push_back(newxtl);
  crystals[0].datasets.push_back(newset);
  crystals[0].datasets[0].columns.push_back(newcol);
  newcol.label = "K";
  crystals[0].datasets[0].columns.push_back(newcol);
  newcol.label = "L";
  crystals[0].datasets[0].columns.push_back(newcol);
  // title, history
  title_ = "From Clipper CCP4MTZfile";
  history_.clear();
  // Default spacegroup_confidence = unknown
  set_spacegroup_confidence('X');

  filename_out_ = filename_out;
  mode = WRITE;
}

/*! Close the file after writing. This command also actually writes
  the data reflection list from the HKL_info object and the data from
  any HKL_data objects which have been marked for import. */
void CCP4MTZfile::close_write()
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "CCP4MTZfile: no file open for write" ) );

  // close the mtz:
  // this actually exports the data from the mtz file
  // - we save up all the work to do it on a single pass

  HKL hkl;
  float fdata[1000]; // file data: allow up to 1000 columns on file
  int   idata[1000]; // flag data: allow up to 1000 columns on file
  xtype pdata[100];  // program data: allow up to 100 floats per group type
  int ref, lst, col, c;

  // set the cell for the HKL crystal
  crystals[0].crystal = MTZcrystal( "HKL_base", "HKL_base", hkl_info_o->cell());

  // NOW THE ACTUAL MTZ STUFF
  // create the mtz object
  CMtz::MTZ* mtzout = CMtz::MtzMalloc( 0, 0 );

  // write general info
  mtzout->refs_in_memory = 0;
  write_spacegroup( mtzout, spacegroup(), spacegroup_confidence() );
  CMtz::ccp4_lwtitl( mtzout, (char*)title_.c_str(), 0 );
  std::vector<char> newhist( history_.size()*MTZRECORDLENGTH, ' ' );
  if ( history_.size() > 0 ) {
    for ( int h = 0; h < history_.size(); h++ )
      for ( int c = 0; c < history_[h].size(); c++ )
	newhist[h*MTZRECORDLENGTH+c] = history_[h][c];
    CMtz::MtzAddHistory( mtzout, (const char (*)[MTZRECORDLENGTH])&(newhist[0]),
			 history_.size() );
  }

  // add all the stored columns to the mtz list
  write_hierarchy( mtzout, crystals );
  std::vector<CMtz::MTZCOL*> lookup_o = build_lookup( mtzout );

  // Lookup the output columns
  std::vector<std::vector<int> > cols;
  std::vector<std::vector<ftype> > scls;
  reference_cols( mtzout, lookup_o, hkl_data_cols, cols, scls );

  // write the data
  int nref = hkl_info_o->num_reflections();
  mtzout->fileout = CMtz::MtzOpenForWrite( filename_out_.c_str() );
  for ( ref = 0; ref < nref; ref++ ) {
    // read/append reflection
    hkl = hkl_info_o->hkl_of(ref);
    fdata[0] = float( hkl.h() );
    fdata[1] = float( hkl.k() );
    fdata[2] = float( hkl.l() );
    idata[0] = idata[1] = idata[2] = 0;
    // output reflection
    for ( lst = 0; lst < cols.size(); lst++ ) {
      // for each list, fetch the data
      hkl_data_o[lst]->data_export( hkl, pdata );
      // and copy to the relevent columns
      for ( col = 0; col < cols[lst].size(); col++ ) {
	// set to mtz->mnf, unless we have a value for it
	if ( !Util::is_nan( pdata[col] ) && cols[lst][col] >= 0 ) {
	  idata[ cols[lst][col] ] = 0;
	  fdata[ cols[lst][col] ] = float( pdata[col] * scls[lst][col] );
	} else {
	  idata[ cols[lst][col] ] = 1;
	  fdata[ cols[lst][col] ] = 0.0;
	}
      }
    }
    // write appended record
    // help: no equivalent to lrrefl! Need idata and lookup
    for ( c = 0; c < lookup_o.size(); c++ )
      if (idata[c]) fdata[c] = CCP4::ccp4_nan().f;
    CMtz::ccp4_lwrefl( mtzout, fdata, &lookup_o[0], lookup_o.size(), ref+1 );
  }

  // tidy up
  CMtz::MtzPut( mtzout, filename_out_.c_str() );
  CMtz::MtzFree( mtzout );

  hkl_data_o.clear();
  hkl_data_cols.clear();
  crystals.clear();
  assigned_paths_.clear();
  filename_in_ = filename_out_ = "";

  mode = NONE;
}

/*! A file is opened for appending. One file is opened for reading,
  and a second is opened for writing. The second file will contain all
  the information from the first, plus any additional columns exported
  from HKL_data objects.
  \param filename_in The input filename or pathname.
  \param filename_out The output filename or pathname. */
void CCP4MTZfile::open_append( const String filename_in, const String filename_out )
{
  if ( mode != NONE )
    Message::message( Message_fatal( "CCP4MTZfile: open_append - File already open" ) );

  // read the mtz
  filename_in_ = filename_in;
  filename_out_ = filename_out;

  // open file
  CMtz::MTZ* mtzin = CMtz::MtzGet( filename_in_.c_str(), 0 );
  if ( mtzin == NULL ) Message::message( Message_fatal( "CCP4MTZfile: open_append - File missing or corrupted: "+filename_in_ ) );
  CMtz::MtzAssignHKLtoBase( mtzin );
  // get the list of column names/paths
  read_hierarchy( mtzin, crystals );
  int x, s, c;
  match_path( "/HKL_base/HKL_base/H", x, s, c );
  if ( x < 0 ) Message::message( Message_fatal( "CCP4MTZfile: No cell! " ) );
  cell_ = crystals[x].crystal;
  // get spacegroup by decoding symops
  read_spacegroup( mtzin, spacegroup_ );
  // get resolution:
  float minres, maxres;
  CMtz::MtzResLimits( mtzin, &minres, &maxres );
  reslim_hi = 0.9999 / sqrt( Util::max(minres,maxres) );
  if ( Util::min(minres,maxres) > 0.0 )
    reslim_lo =  ( 1.0001 / sqrt( Util::min(minres,maxres) ) );
  else
    reslim_lo = 1000000.0;
  resolution_.init( reslim_hi );
  hkl_sampling_ = HKL_sampling( cell_, resolution_ );
  // get title , history
  title_ = String( mtzin->title, 70 );
  history_ = std::vector<String>( mtzin->histlines );
  for ( int h = 0; h < history_.size(); h++ )
    history_[h] = String( mtzin->hist + MTZRECORDLENGTH*h, MTZRECORDLENGTH );
  num_reflections_ = CMtz::MtzNref( mtzin );
  // Sort order
  int isort[5];
  CMtz::ccp4_lrsort( mtzin, isort );
  sort_order_.assign(isort, isort+5);
  ccp4_spacegroup_number_ = mtzin->mtzsymm.spcgrp;
  spacegroup_confidence_  = mtzin->mtzsymm.spg_confidence;

  // close file
  CMtz::MtzFree( mtzin );

  mode = APPEND;
}

/*! Close the files after appending. This command actually copies the
  input file to the output file, adding data from any HKL_data objects
  which have been marked for import. */
void CCP4MTZfile::close_append()
{
  if ( mode != APPEND )
    Message::message( Message_fatal( "CCP4MTZfile: no file open for append" ) );

  // close the mtz:
  // this actually exports the data from the mtz file
  // - we save up all the work to do it on a single pass

  HKL hkl;
  float fdata_i[1000]; // file data: allow up to 1000 columns on file
  int   idata_i[1000]; // flag data: allow up to 1000 columns on file
  float fdata[1000]; // file data: allow up to 1000 columns on file
  int   idata[1000]; // flag data: allow up to 1000 columns on file
  xtype pdata[100];  // program data: allow up to 100 floats per group type
  float res;
  int ref, lst, col, c;

  // NOW THE ACTUAL MTZ STUFF
  // open file
  CMtz::MTZ* mtzin = CMtz::MtzGet( filename_in_.c_str(), 0 );
  if ( mtzin == NULL ) Message::message( Message_fatal( "CCP4MTZfile: close_append - File missing or corrupted: "+filename_in_ ) );
  CMtz::MtzAssignHKLtoBase( mtzin );
  std::vector<CMtz::MTZCOL*> lookup_i = build_lookup( mtzin );

  // lookup hkl cols
  int hklcol_i[3] = {0,1,2};
  for ( int i = 0; i < lookup_i.size(); i++ ) {
    if ( mtz_col_path( mtzin, lookup_i[i] ).tail() == "H" ) hklcol_i[0] = i;
    if ( mtz_col_path( mtzin, lookup_i[i] ).tail() == "K" ) hklcol_i[1] = i;
    if ( mtz_col_path( mtzin, lookup_i[i] ).tail() == "L" ) hklcol_i[2] = i;
  }

  // create the mtz object
  CMtz::MTZ* mtzout = CMtz::MtzMalloc( 0, 0 );

  // write general info
  mtzout->refs_in_memory = 0;
  mtzout->mtzsymm.nsym   = mtzin->mtzsymm.nsym;
  mtzout->mtzsymm.nsymp  = mtzin->mtzsymm.nsymp;
  mtzout->mtzsymm.symtyp = mtzin->mtzsymm.symtyp;
  mtzout->mtzsymm.spcgrp = mtzin->mtzsymm.spcgrp;
  mtzout->mtzsymm.spg_confidence = mtzin->mtzsymm.spg_confidence;
  strncpy( mtzout->mtzsymm.spcgrpname,mtzin->mtzsymm.spcgrpname, 11 );
  strncpy( mtzout->mtzsymm.pgname,    mtzin->mtzsymm.pgname    , 11 );
  for ( int i = 0; i < mtzin->mtzsymm.nsym; i++ )
    for ( int j = 0; j < 4; j++ )
      for ( int k = 0; k < 4; k++ )
	mtzout->mtzsymm.sym[i][j][k] = mtzin->mtzsymm.sym[i][j][k];
  CMtz::ccp4_lwtitl( mtzout, (char*)title_.c_str(), 0 );
  std::vector<char> newhist( history_.size()*MTZRECORDLENGTH, ' ' );
  if ( history_.size() > 0 ) {
    for ( int h = 0; h < history_.size(); h++ )
      for ( int c = 0; c < history_[h].size(); c++ )
	newhist[h*MTZRECORDLENGTH+c] = history_[h][c];
    CMtz::MtzAddHistory( mtzout, (const char (*)[MTZRECORDLENGTH])&(newhist[0]),
			 history_.size() );
  }

  // add all the stored columns to the mtz list
  write_hierarchy( mtzout, crystals );
  std::vector<CMtz::MTZCOL*> lookup_o = build_lookup( mtzout );

  // Lookup the output columns
  std::vector<std::vector<int> > cols;
  std::vector<std::vector<ftype> > scls;
  reference_cols( mtzout, lookup_o, hkl_data_cols, cols, scls );

  // make list of data to append
  std::vector<int> cpcols( lookup_i.size() );
  for ( int i = 0; i < lookup_i.size(); i++ ) {
    String name_i = mtz_col_path( mtzin, lookup_i[i] );
    for ( c = 0; c < lookup_o.size(); c++ )
      if ( name_i == mtz_col_path( mtzout, lookup_o[c] ) ) break;
    if ( c == lookup_o.size() )
      Message::message( Message_fatal( "CCP4MTZfile - internal error" ) );
    cpcols[i] = c;
  }

  // write the data
  int nref = CMtz::MtzNref( mtzin );
  mtzout->fileout = CMtz::MtzOpenForWrite( filename_out_.c_str() );
  for ( ref = 0; ref < nref; ref++ ) {
    // read reflection
    CMtz::ccp4_lrreff( mtzin, &res, fdata_i, idata_i, (const CMtz::MTZCOL**)
		       &lookup_i[0], lookup_i.size(), ref+1 );
    hkl = HKL( Util::intr(fdata_i[hklcol_i[0]]),
	       Util::intr(fdata_i[hklcol_i[1]]),
	       Util::intr(fdata_i[hklcol_i[2]]) );
    // copy across all data from old file
    for ( int i = 0; i < cpcols.size(); i++ ) {
      fdata[cpcols[i]] = fdata_i[i];
      idata[cpcols[i]] = idata_i[i];
    }
    // output reflection
    for ( lst = 0; lst < cols.size(); lst++ ) {
      // for each list, fetch the data
      hkl_data_o[lst]->data_export( hkl, pdata );
      // and copy to the relevent columns
      for ( col = 0; col < cols[lst].size(); col++ ) {
	// set to mtz->mnf, unless we have a value for it
	if ( !Util::is_nan( pdata[col] ) && cols[lst][col] >= 0 ) {
	  idata[ cols[lst][col] ] = 0;
	  fdata[ cols[lst][col] ] = float( pdata[col] * scls[lst][col] );
	} else {
	  idata[ cols[lst][col] ] = 1;
	  fdata[ cols[lst][col] ] = 0.0;
	}
      }
    }
    // write appended record
    // help: no equivalent to lrrefl! Need idata and lookup
    for ( c = 0; c < lookup_o.size(); c++ )
      if (idata[c]) fdata[c] = CCP4::ccp4_nan().f;
    CMtz::ccp4_lwrefl( mtzout, fdata, &lookup_o[0], lookup_o.size(), ref+1 );
  }

  // pointers
  if ( mtzout->xml != NULL ) mtzout->xml = strdup( mtzin->xml );
  if ( mtzout->unknown_headers != NULL ) {
    mtzout->unknown_headers = (char *) malloc(mtzin->n_unknown_headers*MTZRECORDLENGTH*sizeof(char) );
    strncpy( mtzout->unknown_headers, mtzin->unknown_headers, mtzin->n_unknown_headers*MTZRECORDLENGTH );
  }
  mtzout->n_unknown_headers = mtzin->n_unknown_headers;

  // tidy up
  CMtz::MtzPut( mtzout, filename_out_.c_str() );
  CMtz::MtzFree( mtzin );
  CMtz::MtzFree( mtzout );

  hkl_data_i.clear();
  hkl_data_o.clear();
  hkl_data_cols.clear();
  crystals.clear();
  assigned_paths_.clear();
  filename_in_ = filename_out_ = "";

  mode = NONE;
}

/*! Get the spacegroup from the MTZ file. \return The spacegroup. */
const Spacegroup& CCP4MTZfile::spacegroup() const
{ return spacegroup_; }

/*! Get the base cell from the MTZ file. \return The cell. */
const Cell& CCP4MTZfile::cell() const
{ return cell_; }

/*! Get the resolution limit from the MTZ file. \return The resolution. */
const Resolution& CCP4MTZfile::resolution() const
{ return resolution_; }

/*! Get the HKL sampling from the MTZ file. \return The hkl_sampling. */
const HKL_sampling& CCP4MTZfile::hkl_sampling() const
{ return hkl_sampling_; }

/*! Import the list of reflection HKLs from an MTZ file into an
  HKL_info object. If the resolution limit of the HKL_info object is
  lower than the limit of the file, any excess reflections will be
  rejected, as will any systematic absences or duplicates.
  \param target The HKL_info object to be initialised. */
void CCP4MTZfile::import_hkl_list( HKL_info& target )
{
  if ( mode != READ )
    Message::message( Message_fatal( "CCP4MTZfile: no file open for read" ) );

  float fdata[1000], s; // file data: allow up to 1000 columns on file
  int   idata[1000];    // flag data: allow up to 1000 columns on file
  HKL hkl;
  std::vector<HKL> hkls;

  // open file
  CMtz::MTZ* mtzin = CMtz::MtzGet( filename_in_.c_str(), 0 );
  if ( mtzin == NULL ) Message::message( Message_fatal( "CCP4MTZfile: import_hkl_list - File missing or corrupted: "+filename_in_ ) );
  CMtz::MtzAssignHKLtoBase( mtzin );
  std::vector<CMtz::MTZCOL*> lookup_i = build_lookup( mtzin );

  // lookup hkl cols
  int hklcol_i[3] = {0,1,2};
  for ( int i = 0; i < lookup_i.size(); i++ ) {
    if ( mtz_col_path( mtzin, lookup_i[i] ).tail() == "H" ) hklcol_i[0] = i;
    if ( mtz_col_path( mtzin, lookup_i[i] ).tail() == "K" ) hklcol_i[1] = i;
    if ( mtz_col_path( mtzin, lookup_i[i] ).tail() == "L" ) hklcol_i[2] = i;
  }

  // read the reflections from the mtz
  ftype slim = target.resolution().invresolsq_limit();

  // read the reflections
  for ( int ref = 0; ref < CMtz::MtzNref( mtzin ); ref++ ) {
    // read reflection
    CMtz::ccp4_lrrefl( mtzin, &s, fdata, idata, ref );
    hkl = HKL( Util::intr(fdata[hklcol_i[0]]),
	       Util::intr(fdata[hklcol_i[1]]),
	       Util::intr(fdata[hklcol_i[2]]) );
    // check the resolution against the master cell
    bool in_res = hkl.invresolsq(target.cell()) < slim;
    // also check against any hkl_datas which have been assigned
    //for (int lst=0; lst < list_map_i.size(); lst++) in_res = in_res ||
    //        (hkl.invresolsq(list_map_i[lst].list->base_cell()) < slim);
    // store reflection
    if ( in_res ) hkls.push_back( hkl );
  }
  // close file
  CMtz::MtzFree( mtzin );

  target.add_hkl_list( hkls );
}

/*! Import a complete HKL_info object. The supplied HKL_info object is
  examined, and if any of the parameters (spacegroup, cell, or
  resolution) are unset, then they will be set using values from the
  file. The reflections list will then be generated (the default), or
  imported from the file.

  This method is a shortcut which can generally replace most common
  combinations of calls to import_spacegroup(), import_cell(),
  import_resolution() and import_hkl_list().

  \param target The HKL_info object to be initialised.
  \param generate Generate the list of HKLs rather than importing it from the file. */
void CCP4MTZfile::import_hkl_info( HKL_info& target, const bool generate )
{
  // first check if the HKL_info params are already set
  Spacegroup s = target.spacegroup();
  Cell       c = target.cell();
  Resolution r = target.resolution();
  // import any missing params
  if ( s.is_null() ) s = spacegroup_;
  if ( c.is_null() ) c = cell_;
  if ( r.is_null() ) r = resolution_;
  target.init( s, c, r );
  // now make the HKL list
  if ( generate )
    target.generate_hkl_list();
  else
    import_hkl_list( target );
}

/*! \param cxtl The crystal to import.
    \param mtzpath The mtz path of the crystal. */
void CCP4MTZfile::import_crystal( MTZcrystal& cxtl, const String mtzpath )
{
  int x, s, c;
  String colpath = mtzpath;
  if ( colmode == Legacy )
    if ( colpath.find("/")==String::npos && colpath.find("[")==String::npos )
      colpath = "/*/*/["+colpath+"]";
  if ( colpath.find( "[" ) != String::npos ) {
    std::vector<String> list = colpath.split( "[], " );
    colpath = list[0]+list[1];
  }
  match_path( colpath, x, s, c );
  if ( x >= 0 )
    cxtl = crystals[x].crystal;
  else
    Message::message( Message_fatal("CCP4MTZfile: No such crystal: "+colpath) );
  return;
}

/*! \param cset The dataset to import.
    \param mtzpath The mtz path of the dataset. */
void CCP4MTZfile::import_dataset( MTZdataset& cset, const String mtzpath )
{
  int x, s, c;
  String colpath = mtzpath;
  if ( colmode == Legacy )
    if ( colpath.find("/")==String::npos && colpath.find("[")==String::npos )
      colpath = "/*/*/["+colpath+"]";
  if ( colpath.find( "[" ) != String::npos ) {
    std::vector<String> list = colpath.split( "[], " );
    colpath = list[0]+list[1];
  }
  match_path( colpath, x, s, c );
  if ( x >= 0 && s >= 0 )
    cset = crystals[x].datasets[s].dataset;
  else
    Message::message( Message_fatal("CCP4MTZfile: No such dataset: "+colpath) );
  return;
}

void CCP4MTZfile::export_crystal( const MTZcrystal& cxtl, const String mtzpath )
{
  int x, s, c;
  match_path( mtzpath, x, s, c );
  if ( x < 0 ) {
    crystalinf newxtl;
    newxtl.crystal = MTZcrystal( mtzpath.split("/")[0], cxtl.project_name(), cxtl );
    crystals.push_back( newxtl );
  }
}

void CCP4MTZfile::export_dataset( const MTZdataset& cset, const String mtzpath ){
  int x, s, c;
  match_path( mtzpath, x, s, c );
  if ( x < 0 ) Message::message( Message_fatal( "CCP4MTZfile: export_dataset - Missing crystal: "+mtzpath ) );
  if ( s < 0 ) {
    datasetinf newset;
    newset.dataset = MTZdataset( mtzpath.split("/")[1], cset.wavelength() );
    crystals[x].datasets.push_back( newset );
  }
}

/*! Import data from an MTZ file into an HKL_data object. The dataset
  and crystal information from the first corresponding MTZ column are
  also returned.

  An MTZ column type must be present in the MTZ_type_registry for the
  HKL_data type element name concerned.

  This routine does not actually read any data, but rather marks the
  data to be read when the file is closed.

  For container objects import_chkl_data() is preferred.
  \param cdata The HKL_data object into which data is to be imported.
  \param mtzpath The MTZ column names, as a path. See \ref MTZpaths
  for details. */
void CCP4MTZfile::import_hkl_data( HKL_data_base& cdata, const String mtzpath )
{
  if ( mode != READ )
    Message::message( Message_fatal( "CCP4MTZfile: import_hkl_data - no file open for read" ) );

  // legacy column label handling
  String colpath = mtzpath;
  if ( colmode == Legacy )
    if ( colpath.find("/")==String::npos && colpath.find("[")==String::npos )
      colpath = "/*/*/["+colpath+"]";
  // add the exported data columns to the local list
  int ncols = cdata.data_size();
  std::vector<String> col_names = mtz_assign( colpath, cdata.type(),
					      cdata.data_names(), ncols );
  std::vector<String> dat_names = cdata.data_names().split(" ");
  std::vector<hkldatacol> newcols(ncols);
  // assign the columns to mtz indexes
  for ( int col=0; col < ncols; col++ )  // loop over columns in list
    if ( !is_virtual_col( col_names[col] ) ) {
      int x, s, c;
      match_path( col_names[col], x, s, c );
      if ( c < 0 ) Message::message( Message_fatal( "CCP4MTZfile: import_hkl_data - Missing column, crystal or dataset: "+col_names[col] ) );
      String ctype = CCP4MTZ_type_registry::type( dat_names[col] );
      String mtype = crystals[x].datasets[s].columns[c].type;
      if ( ctype != mtype ) Message::message( Message_warn( "CCP4MTZfile: Mtz column type mismatch: "+crystals[x].datasets[s].columns[c].label+" "+mtype+"-"+ctype ) );
      newcols[col].path = "/" + crystals[x].crystal.crystal_name()
	                + "/" + crystals[x].datasets[s].dataset.dataset_name()
                        + "/" + crystals[x].datasets[s].columns[c].label;
      newcols[col].scale = CCP4MTZ_type_registry::scale( dat_names[col] );
      assigned_paths_.push_back(  // store names for user query
        "/" + crystals[x].crystal.crystal_name() +
	"/" + crystals[x].datasets[s].dataset.dataset_name() +
	"/" + crystals[x].datasets[s].columns[c].label +
	" " + crystals[x].datasets[s].columns[c].type );

      // initialise any uninitialised datasets (CLIPPER2)
      if ( cdata.is_null() ) {
	HKL_data<datatypes::Flag> flag( spacegroup_, cell_, hkl_sampling_ );
	cdata.init( spacegroup_, crystals[x].crystal, hkl_sampling_ );
      }
    } else {
      newcols[col].path = "NAN";
      assigned_paths_.push_back( "/*/*/NAN -" );  // store names for user query
    }
  hkl_data_i.push_back( &cdata );
  hkl_data_cols.push_back( newcols );
}



/*! Export a complete HKL_info object, including spacegroup, cell, and
  list of reflection HKLs from an HKL_info object to an MTZ file. This
  is compulsory when writing an MTZ file, but forbidden when
  appending, since the HKLs will then come from the input MTZ.
  \param target The HKL_info object to supply the parameters. */
void CCP4MTZfile::export_hkl_info( const HKL_info& target )
{
  if ( mode != WRITE )
    Message::message( Message_fatal( "CCP4MTZfile: export_hkl_info - no file open for write" ) );

  spacegroup_ = target.spacegroup();
  cell_       = target.cell();
  resolution_ = target.resolution();

  hkl_info_o = &target;
}


/*! Export data from an HKL_data object into an MTZ file. MTZdataset and
  crystal information must be supplied, and will be applied to all
  columns in the output MTZ.

  An MTZ column type must be present in the MTZ_type_registry for the
  HKL_data type element name concerned.

  This routine does not actually write any data, but rather marks the
  data to be written when the file is closed.

  Normally export_chkl_data() is preferred.
  \param cdata The HKL_data object from which data is to be exported.
  \param mtzpath The MTZ column names, as a path. See \ref MTZpaths
  for details. */
void CCP4MTZfile::export_hkl_data( const HKL_data_base& cdata, const String mtzpath )
{
  if ( mode != WRITE && mode != APPEND )
    Message::message( Message_fatal( "CCP4MTZfile: export_hkl_data - no file open for write/append" ) );

  // add the exported data columns to the local list
  int ncols = cdata.data_size();
  std::vector<String> col_names = mtz_assign( mtzpath, cdata.type(),
					      cdata.data_names(), ncols );
  std::vector<String> dat_names = cdata.data_names().split(" ");
  std::vector<hkldatacol> newcols(ncols);
  clipper::String grpname = mtzpath.tail().split("[],")[0];
  clipper::String grptype = CCP4MTZ_type_registry::group_type( cdata.type() );
  // assign the columns to mtz indexes
  for ( int col=0; col < ncols; col++ ) {  // loop over columns in list
    int x, s, c;
    match_path( col_names[col], x, s, c );
    if ( x < 0 || s < 0 ) Message::message( Message_fatal( "CCP4MTZfile: export_hkl_data - Missing crystal or dataset: "+col_names[col] ) );
    datacolinf newcol;
    newcol.label = col_names[col].tail();
    newcol.type = CCP4MTZ_type_registry::type( dat_names[col] );
    newcol.source = " !NEW! ";
    newcol.grpname = grpname;
    newcol.grptype = grptype;
    newcol.grpposn = col;
    c = crystals[x].datasets[s].columns.size();
    crystals[x].datasets[s].columns.push_back(newcol);
    newcols[col].path = "/" + crystals[x].crystal.crystal_name()
                      + "/" + crystals[x].datasets[s].dataset.dataset_name()
                      + "/" + crystals[x].datasets[s].columns[c].label;
    newcols[col].scale = CCP4MTZ_type_registry::scale( dat_names[col] );
    assigned_paths_.push_back(  // store names for user query
      "/" + crystals[x].crystal.crystal_name() +
      "/" + crystals[x].datasets[s].dataset.dataset_name() +
      "/" + crystals[x].datasets[s].columns[c].label +
      " " + crystals[x].datasets[s].columns[c].type );
  }
  hkl_data_o.push_back( &cdata );
  hkl_data_cols.push_back( newcols );
}


/*! Import data from an MTZ into a CHKL_data object. If they don't
  already exist, then CMTZcrystal and CMTZdataset objects will be created to
  match the MTZ crystal and dataset information for the first MTZ
  column used. These new objects will be children of the parent
  CHKL_info object, and this CHKL_data will be moved to become a child
  of the CMTZdataset object.

  Thus, to import data into a CHKL_data, you must first create a
  CHKL_data anywhere below a parent HKL_info. Then call this method,
  and the object will be moved to a position below the HKL_info
  corresponding to its position in the data hierarchy in the MTZ
  file. The following code imports data, dataset, and crystal from an
  MTZ file:

  \code
  CHKL_info myhkl; // must be given cell, spacegroup, and HKL list.
  ...
  CHKL_data<F_sigF> mydata
  CCP4MTZfile file = open_read("in.mtz");
  file.import_chkl_data( mydata, "native/CuKa/[FP,SIGFP]" );
  file.close_read();
  \endcode

  An MTZ column type must be present in the MTZ_type_registry for the
  HKL_data type element name concerned.

  This routine does not actually read any data, but rather marks the
  data to be read when the file is closed.
  \param target The HKL_data object into which data is to be imported.
  \param mtzpath The MTZ column names, as a path. See \ref MTZpaths
  for details.
  \param path Where to put this in the data hierarchy, overriding the
  MTZ crystal and dataset. */
void CCP4MTZfile::import_chkl_data( Container& target, const String mtzpath, const String path )
{
  if ( mode != READ )
    Message::message( Message_fatal( "CCP4MTZfile: no file open for read" ) );

  // get this object
  HKL_data_base* hp = dynamic_cast<HKL_data_base*>(&target);
  if ( hp == NULL )
    Message::message( Message_fatal( "CCP4MTZfile: import object not HKL_data" ) );

  // and the parent hkl_info
  CHKL_info* chkl = target.parent_of_type_ptr<CHKL_info>();
  if ( chkl == NULL )
    Message::message( Message_fatal( "CCP4MTZfile: import HKL_data has no HKL_info" ) );

  // now import the data
  MTZcrystal xtl;
  MTZdataset set;
  import_crystal ( xtl, mtzpath );
  import_dataset ( set, mtzpath );
  import_hkl_data( *hp, mtzpath );

  Container *cxtl, *cset;
  // check for matching crystals
  String xtlname = path.notail().notail().tail();
  if ( xtlname == "" ) xtlname = xtl.crystal_name();
  cxtl = chkl->find_path_ptr( xtlname );
  if ( cxtl == NULL ) {
    cxtl = new CMTZcrystal( *chkl, xtlname, xtl );
    cxtl->set_destroyed_with_parent(); // (use garbage collection for this obj)
  }
  // check for matching datasets
  String setname = path.notail().tail();
  if ( setname == "" ) setname = set.dataset_name();
  cset = cxtl->find_path_ptr( setname );
  if ( cset == NULL ) {
    cset = new CMTZdataset( *cxtl, setname, set );
    cset->set_destroyed_with_parent(); // (use garbage collection for this obj)
  }
  // move the data to the new dataset
  String datname = path.tail();
  if ( datname == "" ) datname = mtzpath.tail();
  target.move( cset->path() + "/" + datname );
}


/*! Export data from a CHKL_data object to an MTZ. The object must
  have a parent Cdataset and CMTZcrystal to provide the MTZ crystal and
  dataset information. The MTZ file will be checked for names matching
  the names of these objects, and the new MTZ columns will be added to
  the corresponding dataset if it exists, otherwise it will be
  created.

  An MTZ column type must be present in the MTZ_type_registry for the
  HKL_data type element name concerned.

  This routine does not actually write any data, but rather marks the
  data to be written when the file is closed.
  \param target The HKL_data object from which data is to be exported.
  \param mtzpath The MTZ column names, as a path. See \ref MTZpaths
  for details. */
void CCP4MTZfile::export_chkl_data( Container& target, const String mtzpath )
{
  if ( mode != WRITE && mode != APPEND )
    Message::message( Message_fatal( "CCP4MTZfile: no file open for write/append" ) );

  // get this object
  HKL_data_base* hp = dynamic_cast<HKL_data_base*>( &target );
  if ( hp == NULL )
    Message::message( Message_fatal( "CCP4MTZfile: export object not HKL_data" ) );
  // get parent objects
  MTZdataset* dp = target.parent_of_type_ptr<MTZdataset>();
  if ( dp == NULL )
    Message::message( Message_fatal( "CCP4MTZfile: HKL_data has no parent MTZdataset" ) );
  MTZcrystal* xp = target.parent_of_type_ptr<MTZcrystal>();
  if ( xp == NULL )
    Message::message( Message_fatal( "CCP4MTZfile: HKL_data has no parent MTZcrystal" ) );

  export_crystal( *xp, mtzpath );
  export_dataset( *dp, mtzpath );
  export_hkl_data( *hp, mtzpath );
}


/*! Return a vector of all the columns in the MTZ file, including
  crystal, dataset and type information. The result is a vector of
  Strings. Each String contains text of the form
  '/CrystalName/DatasetName/ColumnName Type'. */
std::vector<String> CCP4MTZfile::column_paths() const
{
  std::vector<String> result;
  int x, s, c;
  for ( x = 0; x < crystals.size(); x ++ )
    for ( s = 0; s < crystals[x].datasets.size(); s ++ )
      for ( c = 0; c < crystals[x].datasets[s].columns.size(); c ++ )
	result.push_back( "/" + crystals[x].crystal.crystal_name() +
			  "/" + crystals[x].datasets[s].dataset.dataset_name() +
			  "/" + crystals[x].datasets[s].columns[c].label +
			  " " + crystals[x].datasets[s].columns[c].type );
  return result;
}

/*! Return a vector of all the columns in the MTZ file used for the
  previous call to import_hkl_data or export_hkl_data. The result is a
  vector of Strings. Each String contains text of the form
  '/CrystalName/DatasetName/ColumnName Type'.

  This function is commonly used to add an output column to the same
  crystal and dataset as a particular input column. e.g.
  \code
  String opcol = mtzfile.assigned_paths()[0].notail() + "/[FWT,PHIWT]";
  \endcode */
const std::vector<String>& CCP4MTZfile::assigned_paths() const
{ return assigned_paths_; }

String CCP4MTZfile::title() const
{ return title_; }

void CCP4MTZfile::set_title( const String& title )
{ title_ = title; }

std::vector<String> CCP4MTZfile::history() const
{ return history_; }

int CCP4MTZfile::num_reflections() const
{ return num_reflections_; }

std::vector<int> CCP4MTZfile::sort_order() const
{ return sort_order_; }

ftype32 CCP4MTZfile::low_res_limit() const
{ return reslim_lo; }

ftype32 CCP4MTZfile::high_res_limit() const
{ return reslim_hi; }

int CCP4MTZfile::ccp4_spacegroup_number() const
{ return ccp4_spacegroup_number_; }

char CCP4MTZfile::spacegroup_confidence() const
{ return spacegroup_confidence_; }

void CCP4MTZfile::set_spacegroup_confidence(const char& spg_confidence)
{ spacegroup_confidence_ = spg_confidence; }

void CCP4MTZfile::set_history( const std::vector<String>& history )
{
  std::vector<String> histnew;
  histnew.insert( histnew.end(), history.begin() , history.end()  );
  histnew.insert( histnew.end(), history_.begin(), history_.end() );
  history_ = histnew;
}

void CCP4MTZfile::set_column_label_mode( CCP4MTZcolumn_label_mode mode )
{ colmode = mode; }

void CCP4MTZfile::set_verbose( int verbose )
{ verbose_ = verbose; }

/* old form for compatibility */
void CCP4MTZfile::import_hkl_data( HKL_data_base& cdata, MTZdataset& cset, MTZcrystal& cxtl, const String mtzpath )
{
  import_crystal ( cxtl, mtzpath );
  import_dataset ( cset, mtzpath );
  import_hkl_data( cdata, mtzpath );
}

/* old form for compatibility */
void CCP4MTZfile::export_hkl_data( const HKL_data_base& cdata, const MTZdataset& cset, const MTZcrystal& cxtl, const String mtzpath )
{
  export_crystal( cxtl, mtzpath );
  export_dataset( cset, mtzpath );
  export_hkl_data( cdata, mtzpath );
}


/*! Get the column labels from the MTZ file. \return Array of labels. */
std::vector<String> CCP4MTZfile::column_labels() const
{
  std::vector<String> result;
  for ( int x = 0; x < crystals.size(); x++ )
    for ( int s = 0; s < crystals[x].datasets.size(); s++ )
      for ( int c = 0; c < crystals[x].datasets[s].columns.size(); c++ )
	result.push_back( "/"+crystals[x].crystal.crystal_name()+"/"+crystals[x].datasets[s].dataset.dataset_name()+"/" +crystals[x].datasets[s].columns[c].label+" "+crystals[x].datasets[s].columns[c].type);
  return result;
}


} // namespace clipper
