import unittest
import os
import shutil
import sys
import privateer
import test_data


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
        shutil.rmtree(self.test_output)

    def test_annotated_output (self, verbose=False):
        '''
        Test tools_coordinate_kicks
        '''
        pdb_input = os.path.join(self.test_data_path, '5fjj.pdb')
        
        assert os.path.exists(pdb_input)
        
        
        print "Testing sequential annotation..."
        
        xml = privateer.get_annotated_glycans ( pdb_input )

        from lxml import etree
        xml_tree = etree.fromstring ( xml )

        with open ( "annotated_glycans.xml" , "w" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree, pretty_print=True ))
        
        assert ( len(xml_tree.xpath("glycan/sugar[contains(@id, 'NAG')]")) == 69 )
        assert ( len(xml_tree.xpath("glycan/sugar[contains(@id, 'BMA')]")) == 22 )
        assert ( len(xml_tree.xpath("glycan/sugar[contains(@id, 'MAN')]")) == 66 )

        print "Testing hierarchical annotation..."
        
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input )

        from lxml import etree
        xml_tree = etree.fromstring ( xml )

        with open ( "annotated_glycans_hierarchical.xml" , "w" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree, pretty_print=True ))
        
        assert ( len(xml_tree.xpath("glycan/sugar[contains(@id, 'NAG')]")) == 69 )
        assert ( len(xml_tree.xpath("glycan/sugar[contains(@id, 'BMA')]")) == 22 )
        assert ( len(xml_tree.xpath("glycan/sugar[contains(@id, 'MAN')]")) == 66 )
        


        #print xml

        

if __name__ == '__main__':
    unittest.main()
