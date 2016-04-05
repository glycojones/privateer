
/* privateer-frontend_scripting.i */

%module privateer
%include "std_string.i"
%include "std_vector.i"

%{

#include "clipper-glyco_data.h"
#include "privateer-lib.h"
bool clipper::data::found_in_database ( std::string name );

%}

namespace clipper 
{
    namespace data
    {
        bool found_in_database ( std::string name );
    }
}

namespace privateer
{
    namespace scripting
    {
        std::string get_annotated_glycans ( std::string pdb_filename );
    }
}

%pythoncode
%{

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
            from rdkit.Chem import PDBWriter
            writer = PDBWriter ( "minimised_" + str(conformer_id) + ".pdb" )
            writer.write ( mol = mol_with_h, confId = conformer_id )

    from rdkit.Chem import PDBWriter
    writer = PDBWriter ( "lowest.pdb" )
    writer.write ( mol = mol_with_h, confId = best_conformer_id )
    
    print "lowest energy: %f" % lowest_energy
    
%}