from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PDBWriter

def minimise_from_smiles ( smiles_string = "", n_conformers=50 ) :

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
                print ("Energy: %f kcal/mol" % lowest_energy)

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

    print ("Lowest energy: %f" % lowest_energy)
