
/* privateer-frontend_scripting.i */

%module privateer
%include "std_string.i"
%include "std_vector.i"
%template() std::vector<std::string>;

%{

#include "clipper-glyco_data.h"
#include "privateer-lib.h"

%}

namespace privateer
{
    namespace scripting
    {
        std::string carbname_of ( std::string name );
        bool found_in_database ( std::string name );
        std::string get_annotated_glycans ( std::string pdb_filename, bool original_colour_scheme = false, std::string expression_system = "undefined" );
        std::string get_annotated_glycans_hierarchical ( std::string pdb_filename,  bool original_colour_scheme = false, std::string expression_system = "undefined" );
        void privateer::scripting::svg_graphics_demo ( bool original_colour_scheme, bool inverted_background );
    }

    namespace util
    {
        void write_refmac_keywords ( std::vector < std::string > code_list );
        bool write_libraries ( std::vector < std::string > code_list, float esd );
    }

    namespace glycoplot
    {
        std::string get_colour ( privateer::glycoplot::Colour colour, bool original_style );
        enum Colour { blue, red, yellow, orange, green, purple, cyan, tan, black, white };
    }
}

%pythoncode
%{

def print_glycosidic_torsions ( pdb_filename = "", first="ASN", second="NAG" ) :

    import privateer
    xml = privateer.get_annotated_glycans_hierarchical ( pdb_filename )

    from lxml import etree
    xml_tree = etree.fromstring ( xml )

    if not privateer.found_in_database ( first ) :
        first_residues  = xml_tree.xpath("glycan[contains(@root, '" + first + "')]")

        print "RESIDUE \t".rjust(18) + "SUGAR  \t".rjust(18) + "PHI  \t".rjust(18) + "PSI  ".rjust(18)

        for residue in first_residues :
            for linked in residue.xpath("./sugar[contains(@id, '" + second + "')]") :
                print residue.get('root').rjust(18) + "\t" + linked.get('id').rjust(18) + "\t" +\
                      linked.find('link').find('phi').text.rjust(18) + "\t" +\
                      linked.find('link').find('psi').text.rjust(18)

    else :
        first_residues  = xml_tree.xpath(".//sugar[contains(@id, '" + first + "')]")

        print "SUGAR 1  \t".rjust(18) + "SUGAR 2  \t".rjust(18) + "PHI  \t".rjust(18) + "PSI  ".rjust(18)

        for residue in first_residues :
            linked_sugars = residue.xpath("./sugar[contains(@id, '" + second + "')]")
            for linked in linked_sugars :
                print residue.get('id').rjust(18) + "\t" + linked.get('id').rjust(18) + "\t" +\
                      linked.find('link').find('phi').text.rjust(18) + "\t" +\
                      linked.find('link').find('psi').text.rjust(18)


def minimise_from_smiles ( smiles_string = "", n_conformers=50 ) :
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles( smiles_string )

    mol_with_h = AllChem.AddHs(mol)
    conformer_ids = AllChem.EmbedMultipleConfs ( mol_with_h, numConfs=n_conformers )

    best_conformer_id = 0
    lowest_energy = 5000

    for conformer_id in conformer_ids :
        force_field = AllChem.UFFGetMoleculeForceField ( mol_with_h, confId=conformer_id )
        success = force_field.Minimize()
        if success :
            if force_field.CalcEnergy() < lowest_energy :
                best_conformer_id = conformer_id
                lowest_energy = force_field.CalcEnergy()
                print "Energy: %f kcal/mol" % lowest_energy

    from rdkit.Chem import PDBWriter
    writer = PDBWriter ( "privateer-minimised.pdb" )
    writer.write ( mol = mol_with_h, confId = best_conformer_id )
    writer.close()

    with open("privateer-minimised.pdb", 'r+') as f:
        content = f.read()
        f.seek(0,0)
        f.write ( 'REMARK 300 THIS MOLECULE HAS BEEN PRODUCED BY PRIVATEER USING RDKIT\n')
        f.write ( 'REMARK 300 SMILES STRING: ' + smiles_string + '\n')
        f.write ( 'REMARK 300 NUMBER OF CONFORMERS EXPLORED: %i\n' % n_conformers )
        f.write ( 'REMARK 300 ENERGY AFTER MINIMISATION: %f KCAL/MOL\n' % lowest_energy )
        f.write ( content )


    print "lowest energy: %f" % lowest_energy

%}
