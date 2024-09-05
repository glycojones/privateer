import os
import shutil
import sys
from privateer import privateer_core as pvt

def test_svg_graphics ( ):

    '''
    Test SVG graphics in SNFG
    '''
    test_output = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'test_output')
    if not os.path.exists ( test_output ) : os.makedirs ( test_output )

    # GRAPHICS AND LIBRARY CALLS
    print ("Testing SVG graphics output demo (Essentials colour scheme)")

    pvt.svg_graphics_demo ( True, False )

    assert os.path.exists ( "privateer-glycoplot_demo.svg" )

    current_dir = cwd = os.getcwd()
    os.rename ( os.path.join(current_dir, "privateer-glycoplot_demo.svg"), os.path.join(test_output, "privateer-glycoplot_demo_original.svg") )

    print ("Testing SVG graphics output demo (Privateer colour scheme)")

    pvt.svg_graphics_demo ( False, False )

    assert os.path.exists ( "privateer-glycoplot_demo.svg" )
    os.rename ( os.path.join(current_dir, "privateer-glycoplot_demo.svg"), os.path.join(test_output, "privateer-glycoplot_demo_new.svg") )


    print ("Testing SVG graphics output demo (Essentials colour scheme, dark background)")

    pvt.svg_graphics_demo ( True, True )

    assert os.path.exists ( "privateer-glycoplot_demo.svg" )
    os.rename ( os.path.join(current_dir, "privateer-glycoplot_demo.svg"), os.path.join(test_output, "privateer-glycoplot_demo_original_dark.svg") )

    print ("Testing SVG graphics output demo (Privateer colour scheme, dark background)")

    pvt.svg_graphics_demo ( False, True )

    assert os.path.exists ( "privateer-glycoplot_demo.svg" )
    os.rename ( os.path.join(current_dir, "privateer-glycoplot_demo.svg"), os.path.join(test_output, "privateer-glycoplot_demo_new_dark.svg") )
