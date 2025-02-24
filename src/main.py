import pandas as pd
import re
import os
import time
import tempfile

from . import privateer_core as pvt


def privateer_validation_wrapper(OutputFolderPath,m,i,display=False):
    timestr = time.strftime("%Y%m%d-%H%M%S")
    tempdirpath = tempfile.gettempdir()
    InputStructureFilePath = os.path.join(tempdirpath,f"temp_{timestr}.pdb")
    _write_pdb(m,InputStructureFilePath)
    dpath = os.path.dirname(os.path.abspath(__file__))
    zscorefilepath = os.path.join(dpath,"data","linkage_torsions","privateer_torsions_z_score_database.json")
    AllGlycans = pvt.validate(InputStructureFilePath,zscorefilepath)
    os.remove(InputStructureFilePath)
    if display:
        return AllGlycans
    else:
        for glycan in AllGlycans:
            svgstring = glycan["svg"]
            rootID = glycan["RootID"]
            rootID = rootID.replace("/","-")
            svgfile = open(os.path.join(OutputFolderPath,f"model-{i}_{rootID}.svg"),"w")
            svgfile.write(svgstring)
            svgfile.close()
        df = pd.DataFrame.from_dict(AllGlycans)
        df = df.drop(["svg"], axis=1)
        csv_out = os.path.join(OutputFolderPath,f"model-{i}_privateer-report.csv")
        df.to_csv(csv_out)  

_ATOM_FMT = ("ATOM  %5d %-4s%1s"                # serial, atom name, altloc
             "%-3s %1s%4s%1s   "                # res name, chain, seq, insert
             "%8.3f%8.3f%8.3f%6.2f%6.2f      "  # xyz, occupancy, bfactor
             "%4s%2s%2s")                       # segment, element, charge
def _write_pdb(m, filename):
    # Function for writing a pdb from chimeraX model
    # Taken from an example chimeraX bundle
    # Writing then deleting temp.pdb is kinda hacky -- look into proper way to do this converting model too a clipper::mmol object
    atoms = m.atoms
    coords = atoms.coords
    atom_names = atoms.names
    element_names = atoms.element_names
    residues = atoms.residues
    chain_ids = residues.chain_ids
    residue_names = residues.names
    residue_numbers = residues.numbers
    with open(filename, "w") as f:
        serial = 0
        for i in range(len(coords)):
            serial += 1
            atom_name = atom_names[i]
            alt_loc = ' '
            res_name = residue_names[i]
            chain_id = chain_ids[i]
            res_seq = residue_numbers[i]
            res_insert_code = ' '
            x, y, z = coords[i]
            occupancy = 1.0
            b_factor = 0.0
            segment = ' '
            element = element_names[i]
            if len(element) == 1:
                atom_name = ' ' + atom_name
            charge = ' '
            print(_ATOM_FMT % (serial, atom_name, alt_loc,
                               res_name, chain_id, res_seq, res_insert_code,
                               x, y, z, occupancy, b_factor,
                               segment, element, charge),
                               file=f)

def torsion_plot(sugar_1,atom_number_1,sugar_2,atom_number_2,phi,psi):
    import os
    dpath = os.path.dirname(os.path.abspath(__file__))
    torsiondatabasefilepath = os.path.join(dpath,"data","linkage_torsions","privateer_torsion_database.json")

    import json
    with open(torsiondatabasefilepath) as json_file:
        torsions = json.load(json_file)
    torsions = torsions["data"]
    phis = []
    psis = []
    for firsts in torsions:
        if firsts["first"] == sugar_1:
            for seconds in firsts["second"]:
                if seconds["sugar"] == sugar_2 and str(seconds["donor_position"]) == str(atom_number_2) and str(seconds["acceptor_position"]) == str(atom_number_1): #Do I need to swap these two positions?
                    for torsion_pairs in seconds["torsions"]:
                        phis.append(torsion_pairs["Phi"])
                        psis.append(torsion_pairs["Psi"])

    linkageString = f"{sugar_1}-{atom_number_1},{atom_number_2}-{sugar_2}"

    import matplotlib.pyplot as plt
    if sugar_1 == "ASN" and sugar_2 == "NAG":
        phimin = 0
        phimax = 360
        psimin = -180
        psimax = 180
    else:
        phimin = -180
        phimax = 180
        psimin = -180
        psimax = 180

    fig, axs = plt.subplots(1,1,figsize=(5,5))
    axs.hist2d(phis,psis,bins = 180,range=[[phimin,phimax],[psimin,psimax]], cmin = 1)
    axs.plot(phi,psi,"rx")
    axs.set_title(linkageString)
    axs.set_xlabel("$\phi$")
    axs.set_ylabel("$\psi$")
    fig.tight_layout()
    fig.show()




