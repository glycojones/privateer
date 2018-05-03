import unittest
import os
import shutil
import sys
import privateer
import test_data
from datetime import datetime

class Test(unittest.TestCase):
    '''
    Example test script.  For all python code please write a
    corresponding unit test in a similar format to this example.  The test
    script should be named test_<module_name>.py and placed in same source
    directory as the module it is testing.  For further information on unit
    tests please see:

    https://docs.python.org/2/library/unittest.html
    '''

    def setUp(self):
        '''
        Always run at start of test, e.g. for creating directory to store
        temporary test data producing during unit test
        '''
        self.test_data_path = os.path.dirname(test_data.__file__)
        assert os.path.exists(self.test_data_path)
        self.test_output = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            'test_output')
        if not os.path.exists(self.test_output):
            os.makedirs ( self.test_output )
        

    def tearDown(self):
        '''
        Always run at end of test, e.g. to remove temporary data
        '''
        #shutil.rmtree(self.test_output)

    def test_annotated_output (self, verbose=False):
        
        '''
        Test Privateer
        '''
        
        pdb_input = os.path.join(self.test_data_path, "5fjj-high_mannose.pdb")
        
        assert os.path.exists(pdb_input)
        
        
        
        
        # TEST: SEQUENTIAL ANNOTATION
        
        print "Testing sequential annotation    (heaviest glycosylation in PDB)"
        tick = datetime.now()
        xml = privateer.get_annotated_glycans ( pdb_input, True, "fungal" )
        tock = datetime.now()
        
        diff = tock - tick
        print " -> executed in %f seconds" % diff.total_seconds()
        
        # this stuff will be used interactively, so let's clock it under 1 second
        assert ( diff.total_seconds() < 1.0 )
        
        from lxml import etree
        xml_tree = etree.fromstring ( xml )

        tick = datetime.now()
        with open ( "test_output/annotated_glycans_sequential.xml" , "w" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree, pretty_print=True ))
        tock = datetime.now()
        
        diff = tock - tick
        print " -> results parsed in %f seconds" % diff.total_seconds()
        
        # again, let's keep this within acceptable limits
        assert ( diff.total_seconds() < 0.05 )
        
        
        assert ( len(xml_tree.xpath("glycan/sugar[contains(@id, 'NAG')]")) == 69 )
        assert ( len(xml_tree.xpath("glycan/sugar[contains(@id, 'BMA')]")) == 22 )
        assert ( len(xml_tree.xpath("glycan/sugar[contains(@id, 'MAN')]")) == 66 )

        assert ( os.path.exists ("test_output/annotated_glycans_sequential.xml") )
        
        
        
        
        # TEST: HEAVIEST HIGH MANNOSE GLYCOSYLATION IN THE PDB
        
        print "Testing hierarchical annotation  (heaviest glycosylation in PDB)"
        
        tick = datetime.now()
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "fungal" )
        tock = datetime.now()
        
        diff = tock - tick
        print " -> executed in %f seconds" % diff.total_seconds()
        
        xml_tree = etree.fromstring ( xml )

        with open ( "test_output/test-high_mannose_ao.xml" , "w" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree, pretty_print=True ))
        
        assert ( os.path.exists ("test_output/test-high_mannose_ao.xml") )
        
        assert ( len(xml_tree.xpath("//sugar[contains(@id, 'NAG')]")) == 69 )
        assert ( len(xml_tree.xpath("//sugar[contains(@id, 'BMA')]")) == 22 )
        assert ( len(xml_tree.xpath("//sugar[contains(@id, 'MAN')]")) == 66 )


        print "Testing conformational code"
        
        assert ( xml_tree.xpath("//sugar[@id = '/B/1401(NAG)']")[0].find('anomer').text == 'beta' )
        assert ( xml_tree.xpath("//sugar[@id = '/B/1401(NAG)']")[0].find('conformation').text == '4c1' )
        assert ( xml_tree.xpath("//sugar[@id = '/B/1401(NAG)']")[0].find('hand').text == 'D' )
        assert ( xml_tree.xpath("//sugar[@id = '/B/1401(NAG)']")[0].find('cremer-pople_Q').text == '0.535968' )
        assert ( xml_tree.xpath("//sugar[@id = '/B/1401(NAG)']")[0].find('cremer-pople_Phi').text == '350.524' )
        assert ( xml_tree.xpath("//sugar[@id = '/B/1401(NAG)']")[0].find('cremer-pople_Theta').text == '50.2343' )

        print "Testing detection of stacked residues"
        
        assert ( xml_tree.xpath("//sugar[@id = '/A/1401(NAG)']")[0].find('stacked_against').find('residue').get('id') == '/A/431(TRP)' )
        assert ( xml_tree.xpath("//sugar[@id = '/B/1401(NAG)']")[0].find('stacked_against').find('residue').get('id') == '/B/431(TRP)' )
        assert ( xml_tree.xpath("//sugar[@id = '/C/1401(NAG)']")[0].find('stacked_against').find('residue').get('id') == '/C/431(TRP)' )
        assert ( xml_tree.xpath("//sugar[@id = '/D/1401(NAG)']")[0].find('stacked_against').find('residue').get('id') == '/D/431(TRP)' )
        
        
        
        
        # TEST: CLEAREST (MAYBE) HIGH MANNOSE GLYCOSYLATION IN THE PDB
        
        print "Testing hierarchical annotation  (high-mannose glycans with double conformations)"
        
        pdb_input = os.path.join(self.test_data_path, "5fji-high_mannose.pdb")
        
        tick = datetime.now()
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "fungal"  )
        tock = datetime.now()
        
        diff = tock - tick
        print " -> executed in %f seconds" % diff.total_seconds()
        
        xml_tree = etree.fromstring ( xml )

        with open ( "test_output/test-high_mannose_af.xml" , "w" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree, pretty_print=True ))
        
        assert ( os.path.exists ("test_output/test-high_mannose_af.xml") )



        
        # TEST: PLANT GLYCANS
        
        print "Testing hierarchical annotation  (plant glycans)"
        pdb_input = os.path.join(self.test_data_path, "5aog-plant_glycans.pdb")
        
        assert os.path.exists(pdb_input)
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "plant" )
        # write the file before parsing it in case there are problems
        
        with open ( "test_output/test-plant_glycans.xml" , "w" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree, pretty_print=True ))

        assert ( os.path.exists ("test_output/test-plant_glycans.xml") )


        xml_tree = etree.fromstring ( xml )


        

        # TEST: MAMMALIAN GLYCANS

        print "Testing hierarchical annotation  (mammalian glycans)"
        pdb_input = os.path.join(self.test_data_path, "5ajm-mammalian_glycans.pdb")
        
        assert os.path.exists(pdb_input)
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "mammalian" )

        from lxml import etree
        xml_tree = etree.fromstring ( xml )

        with open ( "test_output/test-mammalian_glycans.xml" , "w" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree, pretty_print=True ))
        
        assert ( os.path.exists ("test_output/test-mammalian_glycans.xml") )


        
        
        # TEST: HUMAN ANTIBODIES

        print "Testing hierarchical annotation (human antibodies, wrong links)"
        pdb_input = os.path.join(self.test_data_path, "3sgk-nglycans_antibodies.pdb")
        
        assert os.path.exists(pdb_input)
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "human" )

        from lxml import etree
        xml_tree = etree.fromstring ( xml )

        with open ( "test_output/test-human_antibodies.xml" , "w" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree, pretty_print=True ))
        
        assert ( os.path.exists ("test_output/test-human_antibodies.xml") )




        # TEST: HUMAN ANTIBODIES (SIALYLATED ANTENNAE)

        print "Testing hierarchical annotation (human antibodies with sialylated antennae)"
        pdb_input = os.path.join(self.test_data_path, "4byh-antibodies_sialylated_fc.pdb")
        
        assert os.path.exists(pdb_input)
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "human" )

        from lxml import etree
        xml_tree = etree.fromstring ( xml )

        with open ( "test_output/test-human_sialylated_antibodies.xml" , "w" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree, pretty_print=True ))

        assert ( os.path.exists ("test_output/test-human_sialylated_antibodies.xml" ) )



        
        # TEST: HIGH RESOLUTION CELLULASE DATA

        print "Testing hierarchical annotation  (high resolution cellulase)"
        pdb_input = os.path.join(self.test_data_path, "1kwf-ligand_cellulose_boat.pdb")
        
        assert os.path.exists(pdb_input)
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input )

        from lxml import etree
        xml_tree = etree.fromstring ( xml )

        #with open ( "annotated_glycans_hierarchical.xml" , "w" ) as xml_file :
        #    xml_file.write ( etree.tostring ( xml_tree, pretty_print=True ))
        
        #assert ( len(xml_tree.xpath("//sugar[contains(@id, 'BGC')]")) == 8 )

        #glucose_boat = xml_tree.xpath("//sugar[@id='BGC' and ]")




        # TEST: O-GLYCOSYLATION (X-RAY)

        print "Testing hierarchical annotation  (o-glycans)"
        pdb_input = os.path.join ( self.test_data_path, "4a5t-o_glycans.pdb" )
        
        assert os.path.exists ( pdb_input )
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "human" )

        from lxml import etree
        xml_tree = etree.fromstring ( xml )

        with open ( "test_output/test-o_glycans.xml" , "w" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree, pretty_print=True ))
        
        assert ( os.path.exists ("test_output/test-o_glycans.xml") )




        # TEST: NMR (HIGH-MANNOSE)
        
        print "Testing support for NMR models (high mannose n-glycan)"
        pdb_input = os.path.join(self.test_data_path, "1gya-nmr_n-glycan.pdb")
        
        assert os.path.exists(pdb_input)
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "human" )

        from lxml import etree
        xml_tree = etree.fromstring ( xml )

        with open ( "test_output/test-nmr_high_mannose.xml" , "w" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree, pretty_print=True ))
        
        assert ( os.path.exists ("test_output/test-nmr_high_mannose.xml") )
        


        
        # TEST: NMR (O-GLYCOSYLATION)
        
        print "Testing support for NMR models (o-glycosylation)"
        pdb_input = os.path.join(self.test_data_path, "2lhx-nmr_o-glycans.pdb")
        
        assert os.path.exists(pdb_input)
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "undefined" )

        from lxml import etree
        xml_tree = etree.fromstring ( xml )

        with open ( "test_output/test-nmr_o_glycans.xml" , "w" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree, pretty_print=True ))
        
        assert ( os.path.exists ("test_output/test-nmr_o_glycans.xml") )



        # GRAPHICS AND LIBRARY CALLS

        print "Testing SVG graphics output demo (Essentials colour scheme)"

        privateer.svg_graphics_demo ( True, False )

        assert os.path.exists ( "privateer-glycoplot_demo.svg" )
        os.rename ( "privateer-glycoplot_demo.svg", "test_output/privateer-glycoplot_demo_original.svg" )

        print "Testing SVG graphics output demo (Privateer colour scheme)"

        privateer.svg_graphics_demo ( False, False )

        assert os.path.exists ( "privateer-glycoplot_demo.svg" )
        os.rename ( "privateer-glycoplot_demo.svg", "test_output/privateer-glycoplot_demo_new.svg" )

        print "Testing SVG graphics output demo (Essentials colour scheme, dark background)"

        privateer.svg_graphics_demo ( True, True )

        assert os.path.exists ( "privateer-glycoplot_demo.svg" )
        os.rename ( "privateer-glycoplot_demo.svg", "test_output/privateer-glycoplot_demo_original_dark.svg" )

        print "Testing SVG graphics output demo (Privateer colour scheme, dark background)"

        privateer.svg_graphics_demo ( False, True )

        assert os.path.exists ( "privateer-glycoplot_demo.svg" )
        os.rename ( "privateer-glycoplot_demo.svg", "test_output/privateer-glycoplot_demo_new_dark.svg" )


        print "Testing nomenclature"

        assert ( privateer.carbname_of ( "GLC" ) == "Glc" )
        assert ( privateer.carbname_of ( "BGC" ) == "Glc" )
        assert ( privateer.carbname_of ( "SIA" ) == "Neu5Ac" )
        assert ( privateer.carbname_of ( "ALA" ) == "Unknown" )



if __name__ == '__main__':
    unittest.main()
