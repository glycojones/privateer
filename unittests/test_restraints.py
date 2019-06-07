import os
import shutil
import sys
import privateer.restraints

def test_svg_graphics ( ):

    '''
    Test restraint and conformer generation
    '''
    test_output = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_output')
    if not os.path.exists ( test_output ) : os.makedirs ( test_output )

    print ("Testing access to the CCP4 monomer library")
    path_to_monlib = ""
    path_to_monlib = privateer.restraints.check_monlib_access ( )

    assert ( len(path_to_monlib) > 5 )
