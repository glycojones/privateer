import os
import shutil
import sys
import privateer.restraints

def test_dictionaries ( ):
    test_output = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_output')
    if not os.path.exists ( test_output ) : os.makedirs ( test_output )

    print ("Testing access to the CCP4 monomer library")
    path_to_monlib = ""
    path_to_monlib = privateer.restraints.check_monlib_access ( )

    assert ( len(path_to_monlib) > 5 )

    dictionary = privateer.restraints.CarbohydrateDictionary(path_to_monlib + "/n/NAG.cif")
    assert(dictionary.get_chemcomp_id() == "NAG")

    dictionary.read_from_monlib("MAN")
    assert(dictionary.get_chemcomp_id() == "MAN")

    dictionary.read_from_monlib("XYP")
    assert(dictionary.get_chemcomp_id() == "XYP")

    dictionary.read_from_monlib("15L")
    assert(dictionary.get_chemcomp_id() == "15L")


def test_libraries ( ):
    dictionary = privateer.restraints.CarbohydrateDictionary()
    dictionary.read_from_monlib("GLC")
    library = privateer.restraints.CarbohydrateLibrary()
    assert (library.number_of_entries() == 0)
    library.add_dictionary (dictionary)
    assert (library.number_of_entries() == 1)
    dictionary.read_from_monlib("BGC")
    library.add_dictionary (dictionary)
    assert (library.number_of_entries() == 2)


def test_restraints ( ):
    test_output = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_output')
    if not os.path.exists ( test_output ) : os.makedirs ( test_output )
    dictionary = privateer.restraints.CarbohydrateDictionary()
    dictionary.read_from_monlib("GLC")
    dictionary.restrain_rings_unimodal()
    dictionary.print_torsion_restraints()
    dictionary.write_to_file(test_output + "/GLC_unimodal.cif")
    assert os.path.exists(test_output + "/GLC_unimodal.cif")


def test_chemistry ( ):
    bond_length, bond_esd = privateer.restraints.get_bond_params("GLC", "C1", "C2")
    assert (bond_length == 1.524)
    # to do: test esd's
