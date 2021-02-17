import xml.etree.ElementTree as etree
from privateer.privateer_core import get_annotated_glycans,get_annotated_glycans_hierarchical

def print_glycosidic_torsions ( pdb_filename = "", first="ASN", second="NAG" ) :

    xml = privateer.get_annotated_glycans_hierarchical ( pdb_filename )

    xml_tree = etree.fromstring ( xml )

    if not privateer.found_in_database ( first ) :
        first_residues  = xml_tree.xpath("glycan[contains(@root, '" + first + "')]")

        print (str("RESIDUE \t").rjust(18) + str("SUGAR  \t").rjust(18) + str("PHI  \t").rjust(18) + str("PSI  ").rjust(18))

        for residue in first_residues :
            for linked in residue.xpath("./sugar[contains(@id, '" + second + "')]") :
                print (residue.get('root').rjust(18) + "\t" + linked.get('id').rjust(18) + "\t" +\
                       linked.find('link').find('phi').text.rjust(18) + "\t" +\
                       linked.find('link').find('psi').text.rjust(18))

    else :
        first_residues  = xml_tree.xpath(".//sugar[contains(@id, '" + first + "')]")

        print (str("SUGAR 1  \t").rjust(18) + str("SUGAR 2  \t").rjust(18) + str("PHI  \t").rjust(18) + str("PSI  ").rjust(18))

        for residue in first_residues :
            linked_sugars = residue.xpath("./sugar[contains(@id, '" + second + "')]")
            for linked in linked_sugars :
                print (residue.get('id').rjust(18) + "\t" + linked.get('id').rjust(18) + "\t" +\
                       linked.find('link').find('phi').text.rjust(18) + "\t" +\
                       linked.find('link').find('psi').text.rjust(18) )
