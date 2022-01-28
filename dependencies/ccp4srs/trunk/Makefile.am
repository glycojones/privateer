ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = ccp4srs/libccp4srs.la
bin_PROGRAMS = pdb2to3
noinst_PROGRAMS = srsgen srsview

ccp4srs_libccp4srs_la_SOURCES = \
	ccp4srs/ccp4srs_angle.cpp      ccp4srs/ccp4srs_manager.cpp \
	ccp4srs/ccp4srs_atom.cpp       ccp4srs/ccp4srs_monomer.cpp \
	ccp4srs/ccp4srs_base.cpp       ccp4srs/ccp4srs_plane.cpp   \
	ccp4srs/ccp4srs_bond.cpp       ccp4srs/ccp4srs_torsion.cpp \
	ccp4srs/ccp4srs_chem.cpp       ccp4srs/ccp4srs_tree.cpp    \
	ccp4srs/ccp4srs_chicenter.cpp  ccp4srs/ccp4srs_types.cpp   \
	ccp4srs/ccp4srs_container.cpp  ccp4srs/memio_.cpp          \
	ccp4srs/ccp4srs_entry.cpp      ccp4srs/mpfile_.cpp         \
	ccp4srs/ccp4srs_index.cpp      ccp4srs/ccp4srs_smiles.cpp

pkginclude_HEADERS = \
	ccp4srs/ccp4srs_angle.h        ccp4srs/ccp4srs_index.h    \
	ccp4srs/ccp4srs_atom.h         ccp4srs/ccp4srs_manager.h  \
	ccp4srs/ccp4srs_base.h         ccp4srs/ccp4srs_monomer.h  \
	ccp4srs/ccp4srs_bond.h         ccp4srs/ccp4srs_plane.h    \
	ccp4srs/ccp4srs_chem.h         ccp4srs/ccp4srs_torsion.h  \
	ccp4srs/ccp4srs_chicenter.h    ccp4srs/ccp4srs_tree.h     \
	ccp4srs/ccp4srs_container.h    ccp4srs/ccp4srs_types.h    \
	ccp4srs/ccp4srs_defs.h         ccp4srs/memio_.h           \
	ccp4srs/ccp4srs_entry.h        ccp4srs/mpfile_.h          \
        ccp4srs/ccp4srs_smiles.h

ccp4srs_libccp4srs_la_LIBADD = -lmmdb2 -lz
ccp4srs_libccp4srs_la_LDFLAGS = -no-undefined

srsgen_SOURCES = srsgen.cpp
srsgen_LDADD = ccp4srs/libccp4srs.la -lmmdb2
pdb2to3_SOURCES = srsrem.cpp
pdb2to3_LDADD = ccp4srs/libccp4srs.la -lmmdb2
srsview_SOURCES = srsview.cpp
srsview_LDADD = ccp4srs/libccp4srs.la -lmmdb2

