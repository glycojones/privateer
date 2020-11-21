import unittest
import os
import shutil
import sys
import privateer
import privateer.analysis
import test_data
import requests
import json
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

        assert ( privateer.analysis.found_in_database ("GLC") == True )
        assert ( privateer.analysis.found_in_database ("ALA") == False )

        assert ( privateer.analysis.carbname_of ( "GLC" ) == "Glc" )
        assert ( privateer.analysis.carbname_of ( "BGC" ) == "Glc" )
        assert ( privateer.analysis.carbname_of ( "SIA" ) == "Neu5Ac" )
        assert ( privateer.analysis.carbname_of ( "ALA" ) == "Unknown" )


    def test_sequentially_annotated_output (self, verbose=False):

        '''
        Test sequentially annotated XML
        '''

        pdb_input = os.path.join(self.test_data_path, "5fjj-high_mannose.pdb")
        assert os.path.exists(pdb_input)

        print ("Testing sequential annotation    (heaviest glycosylation in PDB)")
        tick = datetime.now()
        xml = privateer.analysis.get_annotated_glycans ( pdb_input, True, "fungal" )
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
        xml = privateer.analysis.get_annotated_glycans_hierarchical ( pdb_input, True, "fungal" )
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
        xml = privateer.analysis.get_annotated_glycans_hierarchical ( pdb_input, True, "fungal"  )
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

        xml = privateer.analysis.get_annotated_glycans_hierarchical ( pdb_input, True, "plant" )
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
        xml = privateer.analysis.get_annotated_glycans_hierarchical ( pdb_input, True, "mammalian" )

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
        xml = privateer.analysis.get_annotated_glycans_hierarchical ( pdb_input, True, "human" )

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
        xml = privateer.analysis.get_annotated_glycans_hierarchical ( pdb_input, True, "human" )

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
        xml = privateer.analysis.get_annotated_glycans_hierarchical ( pdb_input )

        xml_tree = etree.fromstring ( xml )



    def test_o_glycans (self, verbose=False):

        '''
        Test o-glycans
        '''

        # TEST: O-GLYCOSYLATION (X-RAY)

        print ("Testing hierarchical annotation  (o-glycans)")
        pdb_input = os.path.join ( self.test_data_path, "4a5t-o_glycans.pdb" )

        assert os.path.exists ( pdb_input )
        xml = privateer.analysis.get_annotated_glycans_hierarchical ( pdb_input, True, "human" )

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
        xml = privateer.analysis.get_annotated_glycans_hierarchical ( pdb_input, True, "human" )

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
        xml = privateer.analysis.get_annotated_glycans_hierarchical ( pdb_input, True, "undefined" )

        xml_tree = etree.fromstring ( xml )

        with open ( os.path.join ( self.test_output, "test-nmr_o_glycans.xml") , "wb" ) as xml_file :
            xml_file.write ( etree.tostring ( xml_tree ))

        assert ( os.path.exists (os.path.join(self.test_output, "test-nmr_o_glycans.xml")) )

        

    def test_wurcs_online (self, verbose=False):

        '''
        Test correctness of WURCS outputs of all test glycans
        '''

        print ("Testing WURCS output string correctness for a single model")
        tock = datetime.now()
        totalGlycansInDataset = 0
        model_input = os.path.join(self.test_data_path, "3sgk-nglycans_antibodies.pdb")
        assert os.path.exists(self.test_data_path)


        totalWURCS = privateer.analysis.print_wurcs(model_input)
        temporaryString = totalWURCS.split('\n', 1)[0]
        if temporaryString[0:21] == 'Total Glycans Found: ':
            totalGlycansInModel = int(temporaryString[21:])
        else:
            totalGlycansInModel = 0


        temporaryString = totalWURCS.splitlines()


        confirmNumGlycansInModel=0
        for line in temporaryString:
            if line[:6] == 'WURCS=':
                confirmNumGlycansInModel+=1
                totalGlycansInDataset+=1
                queryLink = 'https://api.glycosmos.org/glytoucan/sparql/wurcs2gtcids?wurcs=' + line
                serverResponse = requests.get(queryLink).json()
                for item in serverResponse:
                    assert (item['id']) and item['WURCS'] == line



        assert(totalGlycansInModel == confirmNumGlycansInModel)
        tick = datetime.now()
        diff = tick - tock

        print("Test duration: %f seconds" % diff.total_seconds())

    def test_wurcs_offline (self, verbose=False):

        '''
        Test correctness of WURCS outputs of all test glycans
        '''

        print ("Testing WURCS output string correctness for all test cases")

        def GetJSON(path):
            if path.endswith(".json"):
                with open(path) as json_file:
                    data = json.load(json_file)
                    return data

        def Find(list, key, value):
            for i, dic in enumerate(list):
                if dic[key] == value:
                    return i
            return "Not Found"


        tock = datetime.now()
        numModels = 0
        totalGlycansInDataset = 0
        totalCarbohydrateLigandsInDataset = 0
        glycosmosDir = os.path.join(self.test_data_path, "glycosmos")
        jsonFile = os.path.join(glycosmosDir, "glycosmos_data_2020-02-20.json")
        glycosmosData = GetJSON(jsonFile)

        for model_input in os.listdir(self.test_data_path):
            if model_input.endswith(".pdb") or model_input.endswith(".cif") or model_input.endswith(".mmcif"):


                file_input = os.path.join(self.test_data_path, model_input)
                assert os.path.exists(self.test_data_path)


                numModels+=1

                totalWURCS = privateer.print_wurcs(file_input)
                # print(totalWURCS)
                temporaryString = totalWURCS.split('\n', 1)[0]
                if temporaryString[0:21] == 'Total Glycans Found: ':
                    totalGlycansInModel = int(temporaryString[21:])
                else:
                    totalGlycansInModel = 0
                    totalCarbohydrateLigandsInDataset+=1


                temporaryListOfStrings = totalWURCS.splitlines()
                temporaryListOfStrings = temporaryListOfStrings[1:]


                confirmNumGlycansInModel=0
                for i in range(totalGlycansInModel):
                    privateerWURCS = temporaryListOfStrings[1]
                    indexMatch = Find(glycosmosData, "Sequence", privateerWURCS)
                    if(indexMatch != "Not Found"):
                        confirmNumGlycansInModel+=1
                        totalGlycansInDataset+=1
                        glycosmosWURCS = glycosmosData[indexMatch]["Sequence"]
                    else:
                        glycosmosWURCS = "Not Found"
                    assert(privateerWURCS == glycosmosWURCS)


                    temporaryListOfStrings = temporaryListOfStrings[2:]
                assert(confirmNumGlycansInModel == totalGlycansInModel)
        tick = datetime.now()
        diff = tick - tock

        print("Test duration: %f seconds" % diff.total_seconds())
        print("Total number of glycoprotein models tested: " + str(numModels))
        print("Total number of glycans in the dataset: " + str(totalGlycansInDataset))
        print("Total number of carbohydrates as ligands in the dataset: " + str(totalCarbohydrateLigandsInDataset))
        print("Total number of proteins containing glycosylation: " + str(numModels - totalCarbohydrateLigandsInDataset) + " out of " + str(numModels) + " models.")


if __name__ == '__main__':
    unittest.main()



if __name__ == '__main__':
    unittest.main()
