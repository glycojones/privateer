import unittest
import os
import shutil
import sys
import privateer
import test_data
from xml.etree import ElementTree as etree
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


    def test_nomenclature (self, verbose=True):
        '''
        Test Privateer's nomenclature translations
        '''

        print ("Testing nomenclature")

        assert ( privateer.found_in_database ("GLC") == True )
        assert ( privateer.found_in_database ("ALA") == False )

        assert ( privateer.carbname_of ( "GLC" ) == "Glc" )
        assert ( privateer.carbname_of ( "BGC" ) == "Glc" )
        assert ( privateer.carbname_of ( "SIA" ) == "Neu5Ac" )
        assert ( privateer.carbname_of ( "ALA" ) == "Unknown" )


    def test_sequentially_annotated_output (self, verbose=False):

        '''
        Test sequentially annotated XML
        '''

        pdb_input = os.path.join(self.test_data_path, "5fjj-high_mannose.pdb")
        assert os.path.exists(pdb_input)

        print ("Testing sequential annotation    (heaviest glycosylation in PDB)")
        tick = datetime.now()
        xml = privateer.get_annotated_glycans ( pdb_input, True, "fungal" )
        tock = datetime.now()

        diff = tock - tick
        print ( " -> executed in %f seconds" % diff.total_seconds() )

        # this stuff will be used interactively, so let's clock it under 1 second
        assert ( diff.total_seconds() < 1.0 )

        xml_tree = etree.fromstring ( xml )

        tick = datetime.now()
        with open ( os.path.join ( self.test_output, "annotated_glycans_sequential.xml" ) , "wb" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree ))
        tock = datetime.now()

        diff = tock - tick
        print (" -> results parsed in %f seconds" % diff.total_seconds() )

        # again, let's keep this within acceptable limits
        assert ( diff.total_seconds() < 0.08 )

        assert ( len(xml_tree.findall("glycan/sugar[@id='/B/1401(NAG)']")) == 1 )

        assert ( os.path.exists ( os.path.join ( self.test_output, "annotated_glycans_sequential.xml")) )


    def test_hierarchically_annotated_output (self, verbose=False):

        '''
        Test sequentially annotated XML
        '''

        pdb_input = os.path.join(self.test_data_path, "5fjj-high_mannose.pdb")
        assert os.path.exists(pdb_input)

        # TEST: HEAVIEST HIGH MANNOSE GLYCOSYLATION IN THE PDB

        print ("Testing hierarchical annotation  (heaviest glycosylation in PDB)")

        tick = datetime.now()
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "fungal" )
        tock = datetime.now()

        diff = tock - tick
        print (" -> executed in %f seconds" % diff.total_seconds())

        xml_tree = etree.fromstring ( xml )

        with open ( os.path.join ( self.test_output, "test-high_mannose_ao.xml" ) , "wb" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree ))

        assert ( os.path.exists (os.path.join ( self.test_output, "test-high_mannose_ao.xml" )) )

        assert ( len(xml_tree.findall("glycan/sugar[@id='/B/1401(NAG)']")) == 1 )


        print ("Testing conformational code")

        assert ( xml_tree.findall("glycan/sugar[@id = '/B/1401(NAG)']")[0].find('anomer').text == 'beta' )
        assert ( xml_tree.findall("glycan/sugar[@id = '/B/1401(NAG)']")[0].find('conformation').text == '4c1' )
        assert ( xml_tree.findall("glycan/sugar[@id = '/B/1401(NAG)']")[0].find('hand').text == 'D' )
        assert ( xml_tree.findall("glycan/sugar[@id = '/B/1401(NAG)']")[0].find('cremer-pople_Q').text == '0.535968' )
        assert ( xml_tree.findall("glycan/sugar[@id = '/B/1401(NAG)']")[0].find('cremer-pople_Phi').text == '350.524' )
        assert ( xml_tree.findall("glycan/sugar[@id = '/B/1401(NAG)']")[0].find('cremer-pople_Theta').text == '50.2343' )

        print ("Testing detection of stacked residues")

        assert ( xml_tree.findall("glycan/sugar[@id = '/A/1401(NAG)']")[0].find('stacked_against').find('residue').get('id') == '/A/431(TRP)' )
        assert ( xml_tree.findall("glycan/sugar[@id = '/B/1401(NAG)']")[0].find('stacked_against').find('residue').get('id') == '/B/431(TRP)' )
        assert ( xml_tree.findall("glycan/sugar[@id = '/C/1401(NAG)']")[0].find('stacked_against').find('residue').get('id') == '/C/431(TRP)' )
        assert ( xml_tree.findall("glycan/sugar[@id = '/D/1401(NAG)']")[0].find('stacked_against').find('residue').get('id') == '/D/431(TRP)' )


    def test_high_mannose_glycans (self, verbose=False):

        '''
        Test high-mannose glycans
        '''

        # TEST: CLEAREST (MAYBE) HIGH MANNOSE GLYCOSYLATION IN THE PDB
        print ("Testing hierarchical annotation  (high-mannose glycans with double conformations)")

        pdb_input = os.path.join(self.test_data_path, "5fji-high_mannose.pdb")
        assert os.path.exists(pdb_input)

        tick = datetime.now()
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "fungal"  )
        tock = datetime.now()

        diff = tock - tick
        print (" -> executed in %f seconds" % diff.total_seconds())

        xml_tree = etree.fromstring ( xml )

        with open ( os.path.join ( self.test_output, "test-high_mannose_af.xml" ) , "wb" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree ))

        assert ( os.path.exists (os.path.join(self.test_output, "test-high_mannose_af.xml")) )


    def test_plant_glycans (self, verbose=False):

        '''
        Test plant glycans
        '''

        # TEST: PLANT GLYCANS
        print ("Testing hierarchical annotation  (plant glycans)")
        pdb_input = os.path.join(self.test_data_path, "5aog-plant_glycans.pdb")
        assert os.path.exists(pdb_input)

        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "plant" )
        xml_tree = etree.fromstring ( xml )
        # write the file before parsing it in case there are problems

        with open ( os.path.join ( self.test_output, "test-plant_glycans.xml") , "wb" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree ))

        assert ( os.path.exists ( os.path.join(self.test_output, "test-plant_glycans.xml")) )

        xml_tree = etree.fromstring ( xml )


    def test_mammalian_glycans (self, verbose=False):

        '''
        Test mammalian glycans
        '''

        # TEST: MAMMALIAN GLYCANS
        print ("Testing hierarchical annotation  (mammalian glycans)")
        pdb_input = os.path.join(self.test_data_path, "5ajm-mammalian_glycans.pdb")

        assert os.path.exists(pdb_input)
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "mammalian" )

        xml_tree = etree.fromstring ( xml )

        with open ( os.path.join ( self.test_output, "test-mammalian_glycans.xml") , "wb" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree ))

        assert ( os.path.exists (os.path.join(self.test_output,"test-mammalian_glycans.xml")) )


    def test_human_glycans_simple (self, verbose=False):

        '''
        Test human glycans (antibodies)
        '''

        # TEST: HUMAN ANTIBODIES
        print ("Testing hierarchical annotation (human antibodies, wrong links)")
        pdb_input = os.path.join(self.test_data_path, "3sgk-nglycans_antibodies.pdb")

        assert os.path.exists(pdb_input)
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "human" )

        xml_tree = etree.fromstring ( xml )

        with open ( os.path.join ( self.test_output, "test-human_antibodies.xml") , "wb" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree ))

        assert ( os.path.exists (os.path.join(self.test_output,"test-human_antibodies.xml")) )


    def test_human_glycans_sialylated (self, verbose=False):

        '''
        Test human glycans (sialylated antibodies)
        '''

        # TEST: HUMAN ANTIBODIES (SIALYLATED ANTENNAE)
        print ("Testing hierarchical annotation (human antibodies with sialylated antennae)")
        pdb_input = os.path.join(self.test_data_path, "4byh-antibodies_sialylated_fc.pdb")

        assert os.path.exists(pdb_input)
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "human" )

        xml_tree = etree.fromstring ( xml )

        with open ( os.path.join ( self.test_output, "test-human_sialylated_antibodies.xml") , "wb" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree ))

        assert ( os.path.exists (os.path.join(self.test_output, "test-human_sialylated_antibodies.xml" )) )


    def test_cellulose_ligand (self, verbose=False):

        '''
        Test cellulose ligand
        '''

        # TEST: HIGH RESOLUTION CELLULASE DATA
        print ("Testing hierarchical annotation  (high resolution cellulase)")
        pdb_input = os.path.join(self.test_data_path, "1kwf-ligand_cellulose_boat.pdb")

        assert os.path.exists(pdb_input)
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input )

        xml_tree = etree.fromstring ( xml )



    def test_o_glycans (self, verbose=False):

        '''
        Test o-glycans
        '''

        # TEST: O-GLYCOSYLATION (X-RAY)

        print ("Testing hierarchical annotation  (o-glycans)")
        pdb_input = os.path.join ( self.test_data_path, "4a5t-o_glycans.pdb" )

        assert os.path.exists ( pdb_input )
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "human" )

        xml_tree = etree.fromstring ( xml )

        with open ( os.path.join ( self.test_output, "test-o_glycans.xml") , "wb" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree ))

        assert ( os.path.exists (os.path.join(self.test_output, "test-o_glycans.xml")) )



    def test_nmr_glycans_simple (self, verbose=False):

        '''
        Test glycans as determined by nmr
        '''

        # TEST: NMR (HIGH-MANNOSE)
        print ("Testing support for NMR models (high mannose n-glycan)")
        pdb_input = os.path.join(self.test_data_path, "1gya-nmr_n-glycan.pdb")

        assert os.path.exists(pdb_input)
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "human" )

        xml_tree = etree.fromstring ( xml )

        with open ( os.path.join ( self.test_output, "test-nmr_high_mannose.xml") , "wb" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree ))

        assert ( os.path.exists (os.path.join(self.test_output, "test-nmr_high_mannose.xml")) )



    def test_nmr_o_glycans (self, verbose=False):

        '''
        Test o-glycans as determined by nmr
        '''

        # TEST: NMR (O-GLYCOSYLATION)
        print ("Testing support for NMR models (o-glycosylation)")
        pdb_input = os.path.join(self.test_data_path, "2lhx-nmr_o-glycans.pdb")

        assert os.path.exists(pdb_input)
        xml = privateer.get_annotated_glycans_hierarchical ( pdb_input, True, "undefined" )

        xml_tree = etree.fromstring ( xml )

        with open ( os.path.join ( self.test_output, "test-nmr_o_glycans.xml") , "wb" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree ))

        assert ( os.path.exists (os.path.join(self.test_output, "test-nmr_o_glycans.xml")) )


if __name__ == '__main__':
    unittest.main()
