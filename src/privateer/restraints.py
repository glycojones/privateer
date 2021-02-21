from privateer.libprivateer import check_monlib_access
from privateer.libprivateer import CarbohydrateDictionary
from privateer.libprivateer import CarbohydrateLibrary
import os

def minimise_from_smiles ( smiles_string = "", n_conformers=50, filename="privateer-minimised.pdb" ) :
    # encapsulate imports here so they don't show up when importing the module
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import PDBWriter

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

    writer = PDBWriter ( filename )
    writer.write ( mol = mol_with_h, confId = best_conformer_id )
    writer.close()

    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0,0)
        f.write ( 'REMARK 300 THIS MOLECULE HAS BEEN PRODUCED BY PRIVATEER USING RDKIT\n')
        f.write ( 'REMARK 300 SMILES STRING: ' + smiles_string + '\n')
        f.write ( 'REMARK 300 NUMBER OF CONFORMERS EXPLORED: %i\n' % n_conformers )
        f.write ( 'REMARK 300 ENERGY AFTER MINIMISATION: %f KCAL/MOL\n' % lowest_energy )
        f.write ( content )

    print ("Lowest energy: %f" % lowest_energy)


def get_bond_params_from_monlib ( three_letter_code="", first_atom="", second_atom="" ):
    if three_letter_code != "" and first_atom != "" and second_atom != "" :
        path_to_monlib = check_monlib_access ( )
        if len(path_to_monlib) > 5 :
            dictionary = CarbohydrateDictionary()
            dictionary.read_from_monlib(three_letter_code)
            bond_params = dictionary.get_bond(first_atom, second_atom)
            return bond_params["length"], bond_params["esd"]
        else :
            return 0.0, 0.0


def get_bond_params_from_file ( filename="", first_atom="", second_atom="" ):
    if filename != "" and first_atom != "" and second_atom != "" :
        if os.path.isfile(filename) :
            dictionary = CarbohydrateDictionary()
            dictionary.read_from_file(filename)
            bond_params = dictionary.get_bond(first_atom, second_atom)
            return bond_params["length"], bond_params["esd"]
        else :
            return 0.0, 0.0


def unimodalise_from_file ( filename = "", also_1c4 = False ) :
    dictionary = CarbohydrateDictionary(filename)
    dictionary.restrain_rings_unimodal_from_conformer()
    if also_1c4 : dictionary.add_inverted_torsions()
    pathname, extension = os.path.splitext(filename)
    dictionary.write_to_file(pathname + "_unimodal.cif")


def unimodalise_from_monlib ( code = "", also_1c4 = False ) :
    dictionary = CarbohydrateDictionary()
    dictionary.read_from_monlib (code)
    dictionary.restrain_rings_unimodal_from_conformer()
    if also_1c4 : dictionary.add_inverted_torsions()
    dictionary.write_to_file(code + "_unimodal.cif")
