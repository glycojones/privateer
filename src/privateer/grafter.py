import os
import shutil
import re
import sys
import argparse
import requests
from datetime import datetime
import warnings
import json
import numpy as np
import pandas as pd
import gemmi
import subprocess
from privateer import privateer_core as pvtcore
from privateer import privateer_modelling as pvtmodelling


def _run_refmac(mtz_in: str, pdb_in: str, mtz_out: str, pdb_out: str, other_out: str, restraint_file: str, ncycles: int): 
    # Taken and edited from ModelCraft 
    # Note: Need to make sure path to ccp4 is set up via: source /Applications/ccp4-8.0/bin/ccp4.setup-sh
    print("Running REFMAC with", mtz_in, pdb_in)
    _args = []
    _args += ["HKLIN", mtz_in]
    _args += ["XYZIN", pdb_in]
    
    _args += ["HKLOUT", mtz_out]
    _args += ["XYZOUT", pdb_out]
    _args += ["XMLOUT", f"{other_out}.xml"]
    #labin = "FP=F"
    #labin += " SIGFP=SIGF"
    #labin += " FREE=FREER"
    labin = "FP=FP"
    labin += " SIGFP=SIGFP"
    labin += " FREE=FREE"
    _stdin = []
    _stdin.append("LABIN " + labin)
    _stdin.append(f"NCYCLES {ncycles}")
    _stdin.append("WEIGHT AUTO")
    _stdin.append("MAKE HYDR NO")
    _stdin.append("MAKE NEWLIGAND NOEXIT")
    _stdin.append("PHOUT")
    _stdin.append("PNAME modelcraft")
    _stdin.append("DNAME modelcraft")
    if os.path.isfile(restraint_file):
        with open(restraint_file) as input_file: 
            for line in input_file:
                _stdin.append(line)
        _stdin.append("END")
    process = subprocess.Popen(
    args=["/Applications/ccp4-8.0/bin/refmac5"] + _args,
    #args=["/jarvis/programs/xtal/ccp4-8.0/bin/refmac5"] + _args,
    stdin=subprocess.PIPE if _stdin else None,
    # stdout=subprocess.PIPE,
    # stderr=subprocess.PIPE,
    encoding="utf8",
    env={**os.environ,},
    cwd=os.getcwd(),
    )
    if _stdin:
        stdin_str = '\n'.join(_stdin)
        process.communicate(input=stdin_str)
    return

def _run_refmacat(mtz_in: str, pdb_in: str, mtz_out: str, pdb_out: str, other_out: str, restraint_file: str, ncycles: int): 
    
    try:
        mtz = gemmi.read_mtz_file(mtz_in)
    except RuntimeError as e:
        print (f'{e} at running REFMACAT for {mtz_in}')
        return

    # print("Running CCP4 REFMAC with", mtz_in, pdb_in)
    _args = []
    _args += ["HKLIN", mtz_in]
    _args += ["XYZIN", pdb_in]
    
    _args += ["HKLOUT", mtz_out]
    _args += ["XYZOUT", pdb_out]
    _args += ["XMLOUT", f"{other_out}.xml"]
    
    if ('F' in mtz.column_labels()) and ('SIGF' in mtz.column_labels()):
        labin = "FP=F"
        labin += " SIGFP=SIGF"
    elif ('FP' in mtz.column_labels()) and ('SIGFP' in mtz.column_labels()):
        labin = "FP=FP"
        labin += " SIGFP=SIGFP"
    elif ('FMEAN' in mtz.column_labels()) and ('SIGFMEAN' in mtz.column_labels()):
        labin = "FP=FMEAN"
        labin += " SIGFP=SIGFMEAN"
    else:
        print(f'Refmac cannot find the amplitude at {mtz_in}')
        print(mtz.column_labels())
        return
    
    labin += " FREE=FreeR_flag"
    _stdin = []
    _stdin.append("LABIN " + labin)
    _stdin.append(f"NCYCLES {ncycles}")
    _stdin.append("WEIGHT AUTO")
    _stdin.append("PHOUT")
    if os.path.isfile(restraint_file):
        with open(restraint_file) as input_file: 
            for line in input_file:
                _stdin.append(line)
    _stdin.append("END")

    process = subprocess.Popen(
    args=["/jarvis/programs/xtal/ccp4-8.0/bin/refmacat"] + _args,
    #args=["/Applications/ccp4-8.0/bin/refmacat"] + _args,
    stdin=subprocess.PIPE if _stdin else None,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    encoding="utf8",
    env={**os.environ,},
    cwd=os.getcwd(),
    )
    if _stdin:
        stdin_str = '\n'.join(_stdin)
        process.communicate(input=stdin_str)
    return



def _import_list_of_uniprotIDs_to_glycosylate(inputFilePath):
    output = []
    with open(inputFilePath) as file:
        lines = file.readlines()
    for line in lines:
        lineSplitCommaList = line.split(",")
        if lineSplitCommaList:
            for split_line in lineSplitCommaList:
                cleanLine = re.sub("\W+", "", split_line)
                output.append(cleanLine)
        else:
            cleanLine = re.sub("\W+", "", line)
            output.append(cleanLine)

    return output


def _parse_json_for_grafting_instructions(uniprotIDListPath):
    with open(uniprotIDListPath) as json_file:
        data = json.load(json_file)

    return data


def _download_and_prepare_alphafoldDB_model(uniprotID, downloadLocation):
    outputFileName = uniprotID + ".pdb"
    outputFilePath = os.path.join(downloadLocation, outputFileName)
    requestURL = f"https://alphafold.ebi.ac.uk/files/AF-{uniprotID}-F1-model_v4.pdb"
    query = requests.get(requestURL, allow_redirects=True)

    outputLines = []
    downloadedLines = query.iter_lines()
    for line in downloadedLines:
        decodedLine = line.decode("utf-8")
        # if decodedLine[:5] != "MODEL":
        outputLines.append(decodedLine)

    with open(outputFilePath, "w") as file:
        file.writelines("%s\n" % l for l in outputLines)

    print(
        f"Successfully downloaded model from AlphaFoldDB with UniProt ID: '{uniprotID}' to {outputFilePath}"
    )
    return outputFilePath


def _query_uniprot_for_glycosylation_locations(uniprotID):
    uniprotRequestURL = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprotID}"
    uniprotResponse = requests.get(uniprotRequestURL,
                                   headers={"Accept": "application/json"})
    if not uniprotResponse.ok:
        uniprotResponse.raise_for_status()
        sys.exit()
    uniprotResponseJSON = uniprotResponse.json()
    uniprotSequence = uniprotResponseJSON["sequence"]
    uniprotFeatures = uniprotResponseJSON["features"]
    uniprotGlycosylations = []
    for item in uniprotFeatures:
        if item["type"] == "CARBOHYD":
            uniprotGlycosylations.append(item)
    # FLAG: When running in CMannosylation mode might want to only keep some of these targets
    outputSequence = uniprotSequence["sequence"]
    outputSequenceLength = uniprotSequence["length"]
    output = {
        "sequenceLength": outputSequenceLength,
        "sequence": outputSequence,
        "glycosylations": uniprotGlycosylations,
    }

    return output


def _get_sequences_in_receiving_model(receiverpath):
    builder_sequence_only = pvtmodelling.Builder(receiverpath, True)
    receiver_sequence = builder_sequence_only.get_receiving_model_sequence_info(
    )

    return receiver_sequence


# privateer::pymodelling::Builder::graft_glycan_to_receiver(int mglycanindex, int receiver_chain_index, int received_residue_index)
def _get_information_about_input_files(receiverpath, donorpath, uniprotID):
    if receiverpath is None and donorpath is None and uniprotID is None:
        raise ValueError(
            "'-info' flag needs to be used in conjuction with '-uniprotID' and/or '-donor_path' and/or '-local_receiver_path'. At least one of those flags need to be provided."
        )
    model_chains = []
    bulletIndex = 1
    if receiverpath is not None:
        print(
            f"\n\n\n{bulletIndex}. Information related to receiver model(glycans will be grafted to this PDB file). Associated with '-local_receiver_path' flag."
        )
        print(f"\nPath: {receiverpath}")
        bulletIndex += 1
        sequences = _get_sequences_in_receiving_model(receiverpath)
        for item in sequences:
            chainIndex = item["index"]
            chainID = item["ChainID"]
            chainSequence = item["Sequence"]
            chainResidues = item["Residues"]
            dict_for_uniprot = {
                "index": chainIndex,
                "chainID": chainID,
                "sequence": chainSequence,
            }
            model_chains.append(dict_for_uniprot)
            print(
                f"\nreceiver_chain_index: {chainIndex} \t\t\t(manual_grafting.json key:value pair)"
            )
            print(f"Chain ID: {chainID}")
            print(f"Chain Sequence: \n{chainSequence}")
            for residue in chainResidues:
                residueIndex = residue["index"]
                residueType = residue["residueType"]
                residueCode = residue["residueCode"]
                residueSeqnum = residue["residueSeqnum"]
                print(
                    f"\treceiving_aa_index: {residueIndex} \t\t(manual_grafting.json key:value pair)"
                )
                print(
                    f"\t\tResidue PDB code: {chainID}/{residueType}-{residueCode}/{residueSeqnum}"
                )
    if donorpath is not None:
        print(
            f"\n\n\n{bulletIndex}. Information related to donor model(glycans will be taken from this PDB file and translocated to another PDB file). Associated with '-donor_path' flag."
        )
        print(f"\nPath: {donorpath}")
        bulletIndex += 1
        glycosylation = pvtcore.GlycosylationComposition_memsafe(donorpath)
        glycans = glycosylation.get_summary_of_detected_glycans()
        for item in glycans:
            glycanIndex = item["GlycanID"]
            glycanWURCS = item["WURCS"]
            glycanType = item["GlycosylationType"]
            glycanRootInfo = item["RootInfo"]
            GlycosidicLinkageTorsions = item["ProteinGlycanLinkageTorsion"]
            print(
                f"\nglycan_index: {glycanIndex} \t\t\t(manual_grafting.json key:value pair)"
            )
            print(f"\tGlycan WURCS: {glycanWURCS}")
            print(f"\tGlycosylation Type: {glycanType}")
            print(
                f"\tGlycan in the donor linked to: {glycanRootInfo['ProteinChainID']}/{glycanRootInfo['ProteinResidueID']}-{glycanRootInfo['ProteinResidueType']} and is modelled as Chain {glycanRootInfo['RootSugarChainID']}"
            )
            print(
                f"\tGlycosidic linkage torsions in the donor model - Phi = {GlycosidicLinkageTorsions['Phi']}, Psi = {GlycosidicLinkageTorsions['Psi']}"
            )
    if uniprotID is not None:
        print(
            f"\n\n\n{bulletIndex}. Information related to glycosylation site data held by UniProt. Associated with '-uniprotID' flag."
        )
        print(f"\nUniProt ID: {uniprotID}")
        bulletIndex += 1
        uniprotQuery = _query_uniprot_for_glycosylation_locations(uniprotID)
        uniprotSequence = uniprotQuery["sequence"]
        uniprotSequenceLength = uniprotQuery["sequenceLength"]
        glycosylations = uniprotQuery["glycosylations"]
        print(f"\nProtein Sequence: \n{uniprotSequence}\n")
        print(f"Protein Sequence Length: {uniprotSequenceLength}\n")
        if model_chains and receiverpath is not None:
            search_result = next(
                (item for item in model_chains
                 if item["sequence"] == uniprotSequence),
                False,
            )
            if search_result == False:
                print(
                    f"WARNING: Unable to find a retrieved {uniprotID} sequence in {receiverpath}.\n"
                )
            else:
                print(
                    f"Successfully managed to find a sequence match for {uniprotID} in receiver_chain_index: {search_result['index']} which is modelled as Chain {search_result['chainID']} in '{receiverpath}'!\n"
                )
        for glycosylation in glycosylations:
            description = glycosylation["description"]
            residueindex = glycosylation["begin"]
            print(f"\n\tUniProt description: {description}")
            print(f"\tGlycosylated at residue seqnum: {residueindex}")


def _get_NGlycosylation_targets_via_consensus_seq(sequences):
    NGlycosylationConsensus = "[N][^P][ST]|[N][A-Z][C]"
    output = []

    for item in sequences:
        currentChainIndex = item["index"]
        currentChainID = item["ChainID"]
        currentSequence = item["Sequence"]

        glycosylationTargets = []

        for match in re.finditer(NGlycosylationConsensus, currentSequence):
            if (currentSequence[match.start()] == item["Residues"][
                    match.start()]["residueCode"]):
                glycosylationTargets.append({
                    "start": match.start(),
                    "end": match.end(),
                    "match": match.group()
                })

        output.append({
            "Sequence": currentSequence,
            "chainIndex": currentChainIndex,
            "currentChainID": currentChainID,
            "glycosylationTargets": glycosylationTargets,
        })

    return output

def _get_CMannosylation_targets_via_consensus_seq(sequences):
    CMannosylationConsensus = "[W]" #This needs to change as it is currently finding too many targets
    output = []
    for item in sequences:
        currentChainIndex = item["index"]
        currentChainID = item["ChainID"]
        currentSequence = item["Sequence"]

        glycosylationTargets = []

        for match in re.finditer(CMannosylationConsensus, currentSequence):
            if (currentSequence[match.start()] == item["Residues"][
                    match.start()]["residueCode"]):
                glycosylationTargets.append({
                    "start": match.start(),
                    "end": match.end(),
                    "match": match.group()
                })

        output.append({
            "Sequence": currentSequence,
            "chainIndex": currentChainIndex,
            "currentChainID": currentChainID,
            "glycosylationTargets": glycosylationTargets,
        })
    return output

def _get_CMannosylation_targets_manual(sequences):
    CMannosylationConsensus = "[W]" #This needs to change as it is currently finding too many targets
    output = []
    for item in sequences:
        currentChainIndex = item["index"]
        currentChainID = item["ChainID"]
        currentSequence = item["Sequence"]

        glycosylationTargets = []

        for match in re.finditer(CMannosylationConsensus, currentSequence):
            if (int(item["Residues"][match.start()]["residueSeqnum"]) == int(195)) and (str(currentChainID) == "A"):
                        glycosylationTargets.append({
                            "start": match.start(),
                            "end": match.end(),
                            "match": match.group()
                        })

        output.append({
            "Sequence": currentSequence,
            "chainIndex": currentChainIndex,
            "currentChainID": currentChainID,
            "glycosylationTargets": glycosylationTargets,
        })
    return output

def _get_CMannosylation_targets_via_blob_search(pdbfile, mtzfile,sequences, avg_dens_threshold = None):
    avglength = 6.4118
    #threshold = 0.08
    st = gemmi.read_structure(pdbfile)
    mtz = gemmi.read_mtz_file(mtzfile)
    grid = mtz.transform_f_phi_to_map('DELFWT', 'PHDELWT', sample_rate=2.0)
    grid.normalize()
    start = 1000 # arbitrary number
    if avg_dens_threshold == None:
        resolution = st.resolution
        threshold = 0.41868*resolution - 0.17116
    else:
        threshold = avg_dens_threshold

    pointlist = []; residuelist = []; chainlist = []
    if not os.path.exists(mtzfile):
        raise ValueError("Selected to find mannosylation sited via blob search but no mtz file")
    for model in st:
        for chain in model:
            for residue in chain:
                # GET ESTIMATED CENTROID WITH AVERAGE TRANSLATION LENGTH ~ 6.411 Å
                if (residue.name == 'TRP'):
                    ce3,cd1 = None,None
                    for atom in residue:
                        if atom.name == 'CE3': 
                            ce3 = atom.pos
                        elif atom.name == 'CD1': 
                            cd1 = atom.pos
                    if ce3 != None and cd1 != None:
                        vCED = cd1-ce3; norm = vCED.length(); uvCED = vCED/norm
                        translatedCED = uvCED*avglength
                        newpoint = translatedCED + ce3
                        gr = grid.clone()
                        gr.set_points_around(newpoint, radius=3, value=start)
                        pointgroup = np.argwhere(gr.array == start) # return positions in symmetry ???

                        # TRACE BACK ONTO ORIGINAL GRID --> SUM DENSITY VALUES
                        values = []
                        for point in pointgroup:
                            value = grid.get_value(point[0],point[1],point[2])
                            values.append(value)
                        avgdense = np.mean(values).round(3)
                        if avgdense > threshold: 
                            residuelist.append(residue.seqid.num)
                            chainlist.append(chain.name)
    CMannosylationConsensus = "[W]"
    output = []
    for item in sequences:
        currentChainIndex = item["index"]
        currentChainID = item["ChainID"]
        currentSequence = item["Sequence"]

        glycosylationTargets = []

        for match in re.finditer(CMannosylationConsensus, currentSequence):
            if (currentSequence[match.start()] == item["Residues"][match.start()]["residueCode"]):
                for i in range(len(residuelist)):
                    blob_chainID = chainlist[i]
                    blob_resID = residuelist[i]
                    if (item["Residues"][match.start()]["residueSeqnum"] == blob_resID) and (currentChainID == blob_chainID):
                        glycosylationTargets.append({
                            "start": match.start(),
                            "end": match.end(),
                            "match": match.group()
                        })
        output.append({
            "Sequence": currentSequence,
            "chainIndex": currentChainIndex,
            "currentChainID": currentChainID,
            "glycosylationTargets": glycosylationTargets,
        })
    return output

def _get_CMannosylation_targets_via_water_search(pdbfile, sequences):
    st = gemmi.read_structure(pdbfile)
    st.standardize_crystal_frame() # some structures don't have standard orientation of crystal frames
    ns = gemmi.NeighborSearch(st[0], st.cell, 5).populate(include_h=False)
    residuelist = []; chainlist = []
    for chain in st[0]:
        for residue in chain:
            if residue.name == 'TRP':
                # print(chain.name,residue)
                for atom in residue:
                    if atom.name=="CD1":
                        marks = ns.find_neighbors(atom, min_dist=0.1, max_dist=3.5) # search waters in distance of 4 Angstrom
                        for mark in marks:
                            cra = mark.to_cra(st[0])
                            # print(cra)
                            if mark.image_idx != 0: 
                                continue # image_idx == 0 ~ identical to model within symmetry operation    
                            if ((cra.residue.het_flag == 'H') and (cra.residue.name == 'HOH')): # check waters
                                residuelist.append(residue.seqid.num)
                                chainlist.append(chain.name)
    CMannosylationConsensus = "[W]"
    output = []
    for item in sequences:
        currentChainIndex = item["index"]
        currentChainID = item["ChainID"]
        currentSequence = item["Sequence"]

        glycosylationTargets = []

        for match in re.finditer(CMannosylationConsensus, currentSequence):
            if (currentSequence[match.start()] == item["Residues"][match.start()]["residueCode"]):
                for i in range(len(residuelist)):
                    blob_chainID = chainlist[i]
                    blob_resID = residuelist[i]
                    if (item["Residues"][match.start()]["residueSeqnum"] == blob_resID) and (currentChainID == blob_chainID):
                        glycosylationTargets.append({
                            "start": match.start(),
                            "end": match.end(),
                            "match": match.group()
                        })
        output.append({
            "Sequence": currentSequence,
            "chainIndex": currentChainIndex,
            "currentChainID": currentChainID,
            "glycosylationTargets": glycosylationTargets,
        })
    return output

def _remove_waters_close_to_TRP(pdb_in, pdb_out): # return pdb. file
    st = gemmi.read_structure(pdb_in)
    st.standardize_crystal_frame() # some structures don't have standard orientation of crystal frames
    ns = gemmi.NeighborSearch(st[0], st.cell, 5).populate(include_h=False)
    print(f'Number of atoms in original model: {st[0].count_atom_sites()}')
    for chain in st[0]:
        for residue in chain:
            if residue.name == 'TRP':
                # print(chain.name,residue)
                for atom in residue:
                    if atom.name=="CD1":
                        marks = ns.find_neighbors(atom, min_dist=0.1, max_dist=3.5) # search waters in distance of 4 Angstrom
                        for mark in marks:
                            cra = mark.to_cra(st[0])
                            # print(cra)
                            if mark.image_idx != 0: 
                                continue # image_idx == 0 ~ identical to model within symmetry operation      
                            if ((cra.residue.het_flag == 'H') and (cra.residue.name == 'HOH')): # check waters
                                nearest = st.cell.find_nearest_pbc_image(atom.pos, cra.atom.pos, mark.image_idx)
                                # print(hetatom)
                                deletedresidue = cra.residue
                                del deletedresidue[0] # delete [O] from water residue --> waters have only [O]
                                print(f'Deleted waters: {cra} within {nearest.dist()} Å to {chain.name}/{residue}')    
    print(f'Number of atoms after deleting local waters: {st[0].count_atom_sites()}')
    # write out file
    st.write_pdb(pdb_out)


def _make_connection_between_protein_and_glycan(filepath):
    st = gemmi.read_structure(filepath)
    ns = gemmi.NeighborSearch(st[0], st.cell, 5).populate(include_h=False) 
    for model in st:
        for chain in model:
            for residue in chain:
                hetatom = gemmi.find_tabulated_residue(residue.name)
                # print(residue,hetatom)
                if str(hetatom.kind) != 'ResidueKind.AA': residue.het_flag = 'H' # change flag ATOM or HETATM
                if str(hetatom.kind) != 'ResidueKind.PYR': continue
                print(chain,residue)
                try: c1 = residue['C1'][0].pos
                except RuntimeError as e:
                    print(f'{e} at {chain,residue}')
                    continue
                marks = ns.find_atoms(c1, '\0', radius=5)
                # marks = ns.find_neighbors(C1, min_dist=0.1, max_dist=3)
                # print(marks)
                # exit()
                for mark in marks:
                    cra = mark.to_cra(st[0])
                    if ((cra.residue.name == 'TRP' and cra.atom.name == 'CD1')):
                        #or (cra.residue.name == 'ASN' and cra.atom.name == 'ND2')
                        #or (cra.residue.name in ['SER','THR'] and cra.atom.name in ['OG','OG1'])):
                        # print(cra)
                        dist = (c1).dist(cra.atom.pos)
                        if dist >= 2.0: continue
                        con = gemmi.Connection()
                        con.name = f'new_covale'
                        con.type = gemmi.ConnectionType.Covale
                        con.asu = gemmi.Asu.Same
                        try:
                            con.partner1 = gemmi.make_address(cra.chain, cra.residue, cra.residue.sole_atom(cra.atom.name))
                        except:
                            con.partner1 = gemmi.make_address(cra.chain, cra.residue, cra.residue[cra.atom.name][0])
                        try:
                            con.partner2 = gemmi.make_address(chain, residue, residue.sole_atom('C1'))
                        except:
                            con.partner2 = gemmi.make_address(chain, residue, residue['C1'][0])
                        con.reported_distance = (c1).dist(cra.atom.pos)
                        st.connections.append(con)
    st.write_pdb(filepath)
     

def _copy_metadata(inputpdb, outputpdb, graftedGlycans):
    # Use input pdb to get meta data, remarks, links etc
    struct_in = gemmi.read_structure(inputpdb)
    struct_out = gemmi.read_structure(outputpdb)
    num_chains_in = len(struct_in)
    num_chains_out = len(struct_out)
    # Copy the grafted glycans into the input structure
    for i in range(len(graftedGlycans)):
        glycan = graftedGlycans[i]
        gly_chain_id = str(glycan["receiving_protein_residue_chain_PDBID"])
        gly_seq_num = int(glycan["donor_glycan_root_PDBID"])
        graft_status = glycan["GraftStatus"]
        num_sugars = glycan["donor_glycan_num_sugars"]
        for m in range(len(struct_in)):
                model = struct_in[m]
                model_out = struct_out[m]
                for c in range(len(model)):
                    chain = model[c]
                    chain_out = model_out[c]
                    st_chain_id = str(chain.name)
                    for r in range(len(chain_out)):
                        residue = chain_out[r]
                        st_seq_num = int(residue.seqid.num)
                        if st_chain_id == gly_chain_id and st_seq_num == gly_seq_num and graft_status:
                            for s in range(num_sugars):
                                residue_out = chain_out[r+s]
                                chain.add_residue(residue_out)
    struct_in.write_pdb(outputpdb)

def _generate_restraints(grafted_pdb, outputpath, resolution):
    glycosylation = pvtcore.GlycosylationComposition(grafted_pdb)
    restraints = glycosylation.return_external_restraints(resolution)
    pdbcode = os.path.basename(grafted_pdb).partition(".")[0]
    output_restraints = outputpath+f'/privateer-restraints_{pdbcode}.txt'
    with open(output_restraints, "w") as restraint_file:
        restraint_file.write(restraints)
    return output_restraints


def _refine_grafted_glycans(grafted_pdb, mtzfile, outputpath, pdbout, mtzout, ncycles, resolution):
    restraints_file = _generate_restraints(grafted_pdb, outputpath, resolution)
    filename = os.path.basename(grafted_pdb).partition(".")[0]
    otherdir =outputpath + "/temp"
    if not os.path.isdir(otherdir):
        os.mkdir(otherdir)
    other  = os.path.join(otherdir, filename)
    _run_refmacat(mtzfile, grafted_pdb, mtzout, pdbout, other, restraints_file, ncycles)
    shutil.rmtree(otherdir)
    if os.path.isfile(pdbout):
        os.remove(grafted_pdb)
        os.remove(restraints_file)
    else:
        print(f"Error refining structure {grafted_pdb}")
    return pdbout, mtzout

def _remove_waters_and_recalc_map(input_pdb, mtzfile, outputpath, pdbout, mtzout):
    st = gemmi.read_structure(input_pdb)
    st.remove_waters()
    st.write_pdb(pdbout)
    filename = os.path.basename(input_pdb).partition(".")[0]
    otherdir =outputpath + "/temp"
    pdb_out = outputpath + "/temp_out.pdb"
    if not os.path.isdir(otherdir):
        os.mkdir(otherdir)
    other  = os.path.join(otherdir, filename)
    _run_refmacat(mtzfile, pdbout, mtzout, pdb_out, other, "Not a file", 0)
    if os.path.isfile(pdb_out):
        os.remove(pdb_out)
    else:
        print(f"Error refining structure {input_pdb}")
    return pdbout, mtzout

def _calc_rscc_grafted_glycans(refined_pdb, original_mtz, graftedGlycans):
    mtz = gemmi.read_mtz_file(original_mtz)
    if ('F' in mtz.column_labels()) and ('SIGF' in mtz.column_labels()):
        glycosylation = pvtcore.GlycosylationComposition(refined_pdb, original_mtz, "F,SIGF")
    elif ('FP' in mtz.column_labels()) and ('SIGFP' in mtz.column_labels()):
        glycosylation = pvtcore.GlycosylationComposition(refined_pdb, original_mtz, "FP,SIGFP")
    elif ('FMEAN' in mtz.column_labels()) and ('SIGFMEAN' in mtz.column_labels()):
        glycosylation = pvtcore.GlycosylationComposition(refined_pdb, original_mtz, "FMEAN,SIGFMEAN")
    num_glycans = glycosylation.get_number_of_glycan_chains_detected()  
    for i in range(len(graftedGlycans)):
        graftedglycan = graftedGlycans[i]
        graftedGlycans[i]["RSCC"] = 0
        for glycan_num in range(num_glycans):
            glycan = glycosylation.get_glycan(glycan_num)
            numsugars = glycan.get_total_number_of_sugars()
            for j in range(numsugars):
                sugar = glycan.get_monosaccharide(j)
                root_info = glycan.get_root_info()
                summary = sugar.get_sugar_summary()
                if summary["sugar_name_short"] == graftedglycan["donor_glycan_root_type"] and root_info["ProteinResidueID"] == graftedglycan["receiving_protein_residue_monomer_PDBID"] and root_info["ProteinChainID"] == graftedglycan["receiving_protein_residue_chain_PDBID"]:
                    graftedGlycans[i]["RSCC"] = summary["RSCC"]
    return graftedGlycans

def _remove_grafted_glycans(refined_pdb, original_mtz, graftedGlycans, outputpath, rscc_threshold = 0.5):
    st = gemmi.read_structure(refined_pdb)
    ms = []
    cs = []
    rs = []
    for i in range(len(graftedGlycans)):
        glycan = graftedGlycans[i]
        # Add check of clashes here too if there is no RSCC yet
        if glycan["RSCC"] < rscc_threshold:
            graftedGlycans[i]["GraftStatus"] = False
            for m, model in enumerate(st):
                for c, chain in enumerate(model):
                    for r, residue in enumerate(chain):
                        st_seq_num = int(residue.seqid.num)
                        gly_seq_num = int(glycan["receiving_protein_residue_monomer_PDBID"])
                        st_chain_id = str(chain.name)
                        gly_chain_id = str(glycan["receiving_protein_residue_chain_PDBID"])
                        if st_seq_num == gly_seq_num and st_chain_id == gly_chain_id:
                            print(f"Deleting glycan grafted at {chain.name}-{residue.seqid.num} from {refined_pdb} due to poor RSCC")
                            ms.append(m)
                            cs.append(c)
                            rs.append(r)
    if len(rs) > 0:
        l = sorted(zip(rs, cs, ms))
        rs, cs, ms = zip(*l)
        for m, c, r in zip(ms[::-1], cs[::-1], rs[::-1]):
            del st[m][c][r]
    st.write_pdb(refined_pdb)
    count = 0
    for glycan in graftedGlycans:
        if not glycan["GraftStatus"]:
            count +=1
    if count >= len(graftedGlycans):
        print("Deleting output PDB file as no grafts occurred.")
        os.remove(refined_pdb)
    else:
        # FLAG: Add waters in here before final refinement???
        filename = os.path.basename(refined_pdb).partition("_")[0]
        pdbout = os.path.join(outputpath, filename + "_grafted.pdb")
        mtzout = os.path.join(outputpath, filename + "_grafted.mtz")
        final_pdb, final_mtz = _refine_grafted_glycans(refined_pdb, original_mtz, outputpath, pdbout, mtzout, 0, 0.1)
    return graftedGlycans

def _glycosylate_receiving_model_using_consensus_seq(
    receiverpath,
    donorpath,
    outputpath,
    glycosylationTargets,
    enableUserMessages,
    trimGlycanIfClashesDetected,
    removeGlycanIfClashesDetected,
):
    builder = pvtmodelling.Builder(
        receiverpath,
        donorpath,
        4, # nThreads -- set to -1 to use all available cores
        trimGlycanIfClashesDetected,
        removeGlycanIfClashesDetected,
        True, # ANY_search_policy
        enableUserMessages,
        False, # debug_output = True or False
    )
    for item in glycosylationTargets:
        chainIndex = item["chainIndex"]
        targets = item["glycosylationTargets"]
        for target in targets:
            currentTargetIndex = target["start"]
            try:
                builder.graft_glycan_to_receiver(0, chainIndex, currentTargetIndex)
            except:
                print(f"Unknown exception while grafting to chain {chainIndex} target {currentTargetIndex}. Skipping graft.")

    graftedGlycanSummary = builder.get_summary_of_grafted_glycans()
    builder.export_grafted_model(outputpath)

    return graftedGlycanSummary


def _glycosylate_receiving_model_using_uniprot_info(
    receiverpath,
    donorpath,
    outputpath,
    targets,
    enableUserMessages,
    trimGlycanIfClashesDetected,
):
    builder = pvtmodelling.Builder(
        receiverpath,
        donorpath,
        4, #nThreads
        trimGlycanIfClashesDetected,
        True, #ANY_search_policy
        enableUserMessages,
        False, #debug_output
    )
    for currentTarget in targets:
        chainIndex = 0
        builder.graft_glycan_to_receiver(0, chainIndex, currentTarget)

    graftedGlycanSummary = builder.get_summary_of_grafted_glycans()
    builder.export_grafted_model(outputpath)

    return graftedGlycanSummary


def _print_grafted_glycans_summary(graftedGlycans):
    for idx, graft in enumerate(graftedGlycans):
        proteinChainID = graft["receiving_protein_residue_chain_PDBID"]
        proteinPDBID = graft["receiving_protein_residue_monomer_PDBID"]
        proteinResidueType = graft["receiving_protein_residue_monomer_type"]
        graftedGlycanChainID = graft["glycan_grafted_as_chainID"]
        graftStatus = graft["GraftStatus"]

        if len(graft["ClashingResidues"]):
            averageTotalAtomicDistance = graft["AvgTotalAtomicDistance"]
            numberOfClashingResidues = len(graft["ClashingResidues"])
            if graftStatus:
                print(
                    f"{idx+1}/{len(graftedGlycans)}: Grafted donor glycan as chain {graftedGlycanChainID} to {proteinChainID}/{proteinResidueType}-{proteinPDBID}. The graft has resulted in {numberOfClashingResidues} clashes with an average atomic distance of from detected clashing residues: {averageTotalAtomicDistance}."
                )
            else:
                print(
                    f"{idx+1}/{len(graftedGlycans)}: Did not graft donor glycan to {proteinChainID}/{proteinResidueType}-{proteinPDBID} as the graft resulted in {numberOfClashingResidues} clashes with an average atomic distance of from detected clashing residues: {averageTotalAtomicDistance}."
                )
        else:
            print(
                f"{idx+1}/{len(graftedGlycans)}: Grafted donor glycan as chain {graftedGlycanChainID} to {proteinChainID}/{proteinResidueType}-{proteinPDBID}. The graft did not produce any clashes."
            )


def _store_grafted_glycans_summary(graftedGlycans, index, totalCount):
    for graft in graftedGlycans:
        proteinChainID = graft["receiving_protein_residue_chain_PDBID"]
        proteinPDBID = graft["receiving_protein_residue_monomer_PDBID"]
        proteinResidueType = graft["receiving_protein_residue_monomer_type"]
        graftedGlycanChainID = graft["glycan_grafted_as_chainID"]
        graftStatus = graft["GraftStatus"]

        output = ""
        if len(graft["ClashingResidues"]):
            averageTotalAtomicDistance = graft["AvgTotalAtomicDistance"]
            numberOfClashingResidues = len(graft["ClashingResidues"])
            if graftStatus:
                output = f"{index+1}/{totalCount}: Grafted donor glycan as chain {graftedGlycanChainID} to {proteinChainID}/{proteinResidueType}-{proteinPDBID}. The graft has resulted in {numberOfClashingResidues} clashes with an average atomic distance of from detected clashing residues: {averageTotalAtomicDistance}."
            else:
                output = f"{index+1}/{totalCount}: Did not graft donor glycan to {proteinChainID}/{proteinResidueType}-{proteinPDBID} as the graft resulted in {numberOfClashingResidues} clashes with an average atomic distance of from detected clashing residues: {averageTotalAtomicDistance}."
        else:
            output = f"{index+1}/{totalCount}: Grafted donor glycan as chain {graftedGlycanChainID} to {proteinChainID}/{proteinResidueType}-{proteinPDBID}. The graft did not produce any clashes."

    return output


def _local_input_model_pipeline(receiverpath, donorpath, outputpath,
                                uniprotID, mode, mtzfile, cryoEM):

    sequences = _get_sequences_in_receiving_model(receiverpath)
    print(f"Local Receiver Model Sequence corresponding to file {receiverpath}")
    if uniprotID is not None:
        outputFileName = uniprotID + ".pdb"
    else:
        outputFileName = os.path.basename(receiverpath).partition('.')[0] + '_grafted.pdb'
    outputpath = os.path.join(outputpath, outputFileName)
    if uniprotID is not None:
        uniprotQuery = _query_uniprot_for_glycosylation_locations(uniprotID)
        uniprotSequence = uniprotQuery["sequence"]
        receiverModelSequence = sequences[0]["Sequence"]
        if receiverModelSequence != uniprotSequence:
            raise ValueError(
                "Receiving model sequence does not match the sequence retrieved from UniProt. Please graft glycans using consensus sequence glycan grafting method"
            )
        else:
            uniprotGlycosylations = uniprotQuery["glycosylations"]
            targets = []
            for item in uniprotGlycosylations:
                if (item["description"][0] == "N"
                        or item["description"][0] == "O"
                        or item["description"][0] == "S"
                        or item["description"][0] == "C"
                        or item["description"][0] == "P"):
                    targets.append(int(item["begin"]) - 1)
            graftedGlycans = _glycosylate_receiving_model_using_uniprot_info(
                receiverpath, donorpath, outputpath, targets, True, False)
            _print_grafted_glycans_summary(graftedGlycans)

    elif mtzfile is not None:
        if mode == 'CMannosylation':
            if cryoEM:
                targets = _get_CMannosylation_targets_via_consensus_seq(sequences)
                removeclashes = False
            else:
                targets_1 = _get_CMannosylation_targets_via_blob_search(receiverpath, mtzfile, sequences)
                targets_2 = _get_CMannosylation_targets_via_water_search(receiverpath, sequences)
                for target_1 in targets_1:
                    for target_2 in targets_2[:]:
                        if target_1["chainIndex"]==target_2["chainIndex"] and target_1["currentChainID"]==target_2["currentChainID"]:
                            targets_2.remove(target_2)
                targets = targets_1 + targets_2
                removeclashes = True
        else:
            raise ValueError("Mode of operation not yet supported. Only CMannosylation and has blob search functionality.")
        graftedGlycans = _glycosylate_receiving_model_using_consensus_seq(
            receiverpath, donorpath, outputpath, targets, True, False, removeclashes)
        _print_grafted_glycans_summary(graftedGlycans)
    
    else:
        if mode == 'CMannosylation':
            targets = _get_CMannosylation_targets_via_consensus_seq(sequences)
            removeclashes = True
        elif mode == 'NGlycosylation':
            targets = _get_NGlycosylation_targets_via_consensus_seq(sequences)
            removeclashes = False
        else:
            raise ValueError(
                "Mode of operation not yet supported. Only CMannosylation and NGlycosylation are currently supported"
            )
        graftedGlycans = _glycosylate_receiving_model_using_consensus_seq(
            receiverpath, donorpath, outputpath, targets, True, False, removeclashes)
        _print_grafted_glycans_summary(graftedGlycans)
    count = 0
    if removeclashes:
        for glycan in graftedGlycans:
            if not glycan["GraftStatus"]:
                count +=1
    if len(graftedGlycans) == 0:
        print("Deleting output PDB file as no grafts occurred.")
        if os.path.isfile(outputpath):
            os.remove(outputpath)
    elif count >= len(graftedGlycans):
        print("Deleting output PDB file as no grafts occurred.")
        if os.path.isfile(outputpath):
            os.remove(outputpath)
    else:
        if mode == 'CMannosylation':
            try:
                _copy_metadata(receiverpath, outputpath, graftedGlycans)
            except:
                print(f"Failed to copy metadata for {outputpath}")
            try:
                _make_connection_between_protein_and_glycan(outputpath)
            except:
                print(f"Failed to generate link between TRP-MAN for file {outputpath}")
        if mtzfile is not None:
            outputlocation = outputpath.rpartition("/")[0]
            filename = os.path.basename(outputpath).partition("_")[0]
            pdbout = os.path.join(outputlocation, filename + "_refined.pdb")
            mmcifout = os.path.join(outputlocation, filename + "_refined.mmcif")
            mtzout = os.path.join(outputlocation, filename + "_refined.mtz")
            st = gemmi.read_structure(outputpath)
            resolution = st.resolution
            refined_pdb, refined_mtz = _refine_grafted_glycans(outputpath, mtzfile, outputlocation, pdbout, mtzout, 20, resolution)
            if os.path.isfile(refined_pdb):
                os.remove(refined_mtz)
                os.remove(mmcifout)
                if mode == 'CMannosylation':
                    try:
                        pdbout = os.path.join(outputlocation, filename + "_removed_waters.pdb")
                        _remove_waters_close_to_TRP(refined_pdb, pdbout)
                        os.remove(refined_pdb)
                    except:
                        pdbout = refined_pdb
                        print(f"Failed to remove waters close to TRP in {pdbout}")
                graftedGlycans = _calc_rscc_grafted_glycans(pdbout, mtzfile, graftedGlycans)
                graftedGlycans = _remove_grafted_glycans(pdbout, mtzfile, graftedGlycans, outputlocation)
    return graftedGlycans


# P27918 is a good test for TRP, as TRP is a bit more unique and requires different approach. C-mannosylation currently no worky.
def _online_input_model_pipeline(uniprotID, donorpath, defaultInputModelPath,
                                 outputLocation):
    outputFileName = uniprotID + ".pdb"
    outputpath = os.path.join(outputLocation, outputFileName)
    receiverpath = _download_and_prepare_alphafoldDB_model(
        uniprotID, defaultInputModelPath)
    uniprotGlycosylationQuery = _query_uniprot_for_glycosylation_locations(
        uniprotID)
    uniprotGlycosylations = uniprotGlycosylationQuery["glycosylations"]
    targets = []
    for item in uniprotGlycosylations:
        if (item["description"][0] == "N" or item["description"][0] == "O"
                or item["description"][0] == "S"
                or item["description"][0] == "C"
                or item["description"][0] == "P"):
            targets.append(int(item["begin"]) - 1)
    graftedGlycans = _glycosylate_receiving_model_using_uniprot_info(
        receiverpath, donorpath, outputpath, targets, True, False)
    _print_grafted_glycans_summary(graftedGlycans)
    return graftedGlycans


def glycosylate_receiving_model_using_manual_instructions(
        receiverpath,
        donorpath,
        outputpath,
        glycanIndex,
        receiverChainIndex,
        receiverResidueIndex,
        enableUserMessages=True,
        trimGlycanIfClashesDetected=False,
        nThreads=4):
    builder = pvtmodelling.Builder(
        receiverpath,
        donorpath,
        4, #nThreads
        trimGlycanIfClashesDetected,
        False, #Don't remove glycans if clashes detected
        True, #ANY_search_Policy
        enableUserMessages,
        False, #debug_output
    )
    try:
        builder.graft_glycan_to_receiver(glycanIndex, receiverChainIndex, receiverResidueIndex)
    except:
        print(f"Unknown exception while grafting to chain {receiverChainIndex} target {receiverResidueIndex}. Skipping graft.")

    graftedGlycanSummary = builder.get_summary_of_grafted_glycans()
    builder.export_grafted_model(outputpath)

    return graftedGlycanSummary


def glycosylate_receiving_model_using_schema(schema, print_messages=True):
    initialInputPath = schema["receiver_path"]
    initialOutputSubsequentInputOutputPath = schema["output_path"]
    glycosylations = schema["glycosylations"]
    graftedGlycansSummary = []
    detailed_graftedGlycansSummary = []
    for count, item in enumerate(glycosylations):
        donorPath = item["donor_path"]
        glycanIndex = item["glycan_index"]
        receivingChainIndex = item["receiving_chain_index"]
        receivingAminoAcidIndex = item["receiving_aa_index"]
        if count == 0:
            currentGraftedGlycanSummary = (
                glycosylate_receiving_model_using_manual_instructions(
                    initialInputPath, donorPath,
                    initialOutputSubsequentInputOutputPath, glycanIndex,
                    receivingChainIndex, receivingAminoAcidIndex, True, False,
                    -1))
            detailed_graftedGlycansSummary.append(currentGraftedGlycanSummary)
            messageString = _store_grafted_glycans_summary(
                currentGraftedGlycanSummary, count, len(glycosylations))
            graftedGlycansSummary.append(messageString)
        else:
            currentGraftedGlycanSummary = (
                glycosylate_receiving_model_using_manual_instructions(
                    initialOutputSubsequentInputOutputPath, donorPath,
                    initialOutputSubsequentInputOutputPath, glycanIndex,
                    receivingChainIndex, receivingAminoAcidIndex, True, False,
                    -1))
            detailed_graftedGlycansSummary.append(currentGraftedGlycanSummary)
            messageString = _store_grafted_glycans_summary(
                currentGraftedGlycanSummary, count, len(glycosylations))
            graftedGlycansSummary.append(messageString)
    if print_messages:
        print("\n")
        for message in graftedGlycansSummary:
            print(message + "\n")

    return detailed_graftedGlycansSummary, graftedGlycansSummary


if __name__ == "__main__":
    dt_string = datetime.now().strftime("%d_%m_%Y-%H_%M_%S")
    if os.getenv("PRIVATEERDATA", None) is not None:
        ROOTENV = os.getenv("PRIVATEERDATA", None)
        envbased = True
    elif os.getenv("CLIBD", None) is not None:
        ROOTENV = os.getenv("CLIBD", None)
        ROOTENV = os.path.join(ROOTENV, "privateer_data")
        envbased = True

    else:
        envbased = False
        ROOTENV = os.getcwd()
        defaultDonorLocation = None
        defaultInputModelLocation = None

    if envbased:
        defaultDonorPath = os.path.join(ROOTENV, "glycan_donor_repertoire",
                                        "High_Mannose", "man5", "cluster1.pdb")
        defaultInputModelPath = os.path.join(ROOTENV,
                                             "glycan_donor_repertoire")

    if os.getenv("PRIVATEERRESULTS", None) is not None:
        RESULTSENV = os.getenv("PRIVATEERRESULTS", None)
        defaultOutputModelDirectory = os.path.join(RESULTSENV,"grafter_job" + "__" +
                                                   dt_string)
    else:
        defaultOutputModelDirectory = "grafter_job" + "__" + dt_string

    defaultuniprotIDsListPath = os.path.join(ROOTENV, "uniprotIDinputs.txt")
    defaultJSONgrafting = os.path.join(ROOTENV, "manual_grafting.json")
    defaultUniprotID = "P29016"
    defaultmode = "NGlycosylation"
    defaultSaveSummary = False

    printInfo = False

    parser = argparse.ArgumentParser(
        prog="grafter.py",
        usage=
        "%(prog)s [options]. Most convenient usage: python grafter.py -import_uniprotIDs_from_file uniprotIDinputs.txt",
        description=
        f"Graft Glycans to AlphaFoldDB models using Privateer Modelling module.",
        epilog=
        f"If -local_receiver_path or -uniprotID are not provided, the script will default to using UniProtID: {defaultUniprotID} as default input. Will download the PDB from AlphaFoldDB and N-glycosylate according to UniProt data.",
    )
    parser.add_argument(
        "-uniprotID",
        action="store",
        default=None,
        dest="user_uniprotID",
        help=
        "If used with -local_receiver_path, N-glycosylate according to UniProt targets. If used without -local_receiver_path, this variable is used in the download of AlphaFoldDB .pdb file and N-Glycosylation according to UniProt targets.",
    )
    parser.add_argument(
        "-local_receiver_path",
        action="store",
        default=None,
        dest="user_localReceiverPath",
        help=
        f"Path to locally saved AlpfaFoldDB model on the computer. If -uniprotID is not provided, will carry out N-glycosylation according to regex consensus sequence of '[N][^P][ST]|[N][A-Z][C]'. The argument overrides default behaviour of downloading AlpfaFoldDB model from the server. WARNING: Ensure that \"MODEL 0\" line is deleted in the local file, as otherwise Privateer's MMBD dependency will not be able to import the model!",
    )
    parser.add_argument(
        "-local_receiver_path_mtz",
        action="store",
        default=None,
        dest="user_localReceiverPath_mtz",
        help=
        f"Path to locally saved mtz file on the computer. If provided will carry out mannosylation according to blob search. Overridden by -uniprotID as alphafold models don't have mtz data.",
    )
    parser.add_argument(
        "-donor_path",
        action="store",
        default=None,
        dest="user_donorPath",
        help=
        f"Path to the glycan that is to be grafted throughout AlphaFoldDB model. If not specified, the script will default to using glycan located in '{defaultDonorPath}'",
    )
    parser.add_argument(
        "-download_path",
        action="store",
        default=None,
        dest="user_inputModelDirectory",
        help=
        f"Specify download directory where original AlpfaFoldDB models downloaded from the server should be saved. If unspecified, the script will default to '{defaultOutputModelDirectory}'",
    )
    parser.add_argument(
        "-output_path",
        action="store",
        default=None,
        dest="user_outputPath",
        help=
        f"Specify output directory where AlpfaFoldDB models with grafted glycans should be saved. If unspecified, the script will default to '{defaultOutputModelDirectory}'. If the argument is used alongside -local_receiver_path, then the name of PDB output file should be provided, for example 'P29016_output.pdb'",
    )
    parser.add_argument(
        "-import_uniprotIDs_from_file",
        action="store",
        default=None,
        dest="user_uniprotIDsList",
        help=
        f"Glycosylate multiple AlphaFoldDB models from a list of UniProtIDs. Example file is located in '{defaultuniprotIDsListPath}' By default will download files from the server and save them localy in specified or default directory locations.",
    )

    parser.add_argument(
        "-info",
        action="store_true",
        default=False,
        dest="user_infoFlag",
        help=
        f"Print out relevant information about donor PDB(where glycans are taken from) and receiver PDB(where glycans are grafted to). To be used in conjuction with '-local_receiver_path' and/or '-donor_path' and/or '-uniprotID' flags. Usage of this flag overrides grafting functionality, i.e. no grafting will be carried out.",
    )
    parser.add_argument(
        "-manual_grafting",
        action="store",
        default=None,
        dest="user_JSONgrafting",
        help=
        f"Import a JSON file to manually graft glycans with total control over glycosylation sites. Example file is located at '{defaultJSONgrafting}'",
    )

    parser.add_argument(
        "-operation_mode",
        action="store",
        default=None,
        dest="user_mode",
        help=
        f"Specify the mode of operation for the grafter as either NGlycosylation or CMannosylation. If unspecified, the script will default to '{defaultmode}'",
    )

    parser.add_argument(
        "-save_summary",
        action="store",
        default=None,
        dest="save_summary",
        help=
        f"Save glycan summary to a csv. If unspecified, the script will default to '{defaultSaveSummary}'",
    )

    parser.add_argument(
        "-cryoEM",
        action="store",
        default=False,
        dest="cryoEM",
        help=
        f"If this is a strucutre from experiment, is this a cryoEM structure?",
    )

    args = parser.parse_args()

    if args.user_uniprotID is not None:
        uniprotID = args.user_uniprotID
    else:
        uniprotID = defaultUniprotID
    if args.user_donorPath is not None:
        donorPath = args.user_donorPath
    else:
        donorPath = defaultDonorPath
    if args.user_outputPath is not None:
        outputPath = args.user_outputPath
        if (os.path.isdir(args.user_localReceiverPath)
                and args.user_localReceiverPath is not None):
            raise ValueError(
                "ERROR: The combination of provided arguments requires -output_path argument to be a file name, rather than directory!"
            )
    else:
        outputPath = defaultOutputModelDirectory

    if os.path.isdir(outputPath) == False:
        os.mkdir(outputPath)

    if args.user_inputModelDirectory is not None:
        inputModelDirectory = args.user_inputModelDirectory
    else:
        inputModelDirectory = defaultInputModelPath

    if args.user_uniprotIDsList is not None:
        uniprotIDListPath = args.user_uniprotIDsList

    if args.user_JSONgrafting is not None:
        JSONgraftingPath = args.user_JSONgrafting

    if args.user_infoFlag == True and not None:
        printInfo = True
    
    if args.user_mode is not None:
        mode = args.user_mode
    else:
        mode = defaultmode
    
    if args.save_summary is not None:
        SaveSummary = args.save_summary
    else:
        SaveSummary = defaultSaveSummary


    if (args.user_localReceiverPath is not None and args.user_uniprotID is None
            and printInfo == False):
        uniprotID = None
        graftedGlycans = _local_input_model_pipeline(args.user_localReceiverPath, donorPath,
                                    outputPath, uniprotID, mode, args.user_localReceiverPath_mtz, args.cryoEM)
    elif (args.user_localReceiverPath is not None
          and args.user_uniprotID is not None and printInfo == False):
        graftedGlycans = _local_input_model_pipeline(args.user_localReceiverPath, donorPath,
                                    outputPath, uniprotID, mode, args.user_localReceiverPath_mtz, args.cryoEM)
    elif args.user_uniprotIDsList is not None and printInfo == False:
        uniprotIDList = _import_list_of_uniprotIDs_to_glycosylate(
            uniprotIDListPath)
        for idx, uniprotID in enumerate(uniprotIDList):
            graftedGlycans = _online_input_model_pipeline(uniprotID, donorPath,
                                         inputModelDirectory, outputPath)
            print(
                f"\n{idx+1}/{len(uniprotIDList)}: Successfully finished processing AlphaFoldDB model with UniProt ID of {uniprotID}.\n"
            )
    elif args.user_JSONgrafting is not None and printInfo == False:
        JSONGraftInstructions = _parse_json_for_grafting_instructions(
            JSONgraftingPath)
        initialInputPath = JSONGraftInstructions["receiver_path"]
        initialOutputSubsequentInputOutputPath = JSONGraftInstructions[
            "output_path"]
        glycosylations = JSONGraftInstructions["glycosylations"]
        graftedGlycansSummary = []
        graftedGlycans = []
        for count, item in enumerate(glycosylations):
            donorPath = item["donor_path"]
            glycanIndex = item["glycan_index"]
            receivingChainIndex = item["receiving_chain_index"]
            receivingAminoAcidIndex = item["receiving_aa_index"]
            if count == 0:
                currentGraftedGlycanSummary = (
                    glycosylate_receiving_model_using_manual_instructions(
                        initialInputPath, donorPath,
                        initialOutputSubsequentInputOutputPath, glycanIndex,
                        receivingChainIndex, receivingAminoAcidIndex, True,
                        False, -1))
                messageString = _store_grafted_glycans_summary(
                    currentGraftedGlycanSummary, count, len(glycosylations))
                graftedGlycansSummary.append(messageString)
            else:
                currentGraftedGlycanSummary = (
                    glycosylate_receiving_model_using_manual_instructions(
                        initialOutputSubsequentInputOutputPath, donorPath,
                        initialOutputSubsequentInputOutputPath, glycanIndex,
                        receivingChainIndex, receivingAminoAcidIndex, True,
                        False, -1))
                messageString = _store_grafted_glycans_summary(
                    currentGraftedGlycanSummary, count, len(glycosylations))
                graftedGlycans.append(currentGraftedGlycanSummary)
                graftedGlycansSummary.append(messageString)
        print("\n")
        for message in graftedGlycansSummary:
            print(message + "\n")

    elif printInfo == True:
        warnings.warn(
            "-info flag was provided, overriding all arguments regarding grafting and printing info only. Please remove -info flag if you actually want to graft glycans."
        )
        _get_information_about_input_files(args.user_localReceiverPath,
                                           args.user_donorPath,
                                           args.user_uniprotID)
    else:
        if printInfo == False:
            graftedGlycans = _online_input_model_pipeline(uniprotID, donorPath,
                                         inputModelDirectory, outputPath)
        else:
            warnings.warn(
                "-info flag was provided, overriding all arguments regarding grafting and printing info only. Please remove -info flag if you actually want to graft glycans."
            )
    
    if SaveSummary:
        df = pd.DataFrame.from_dict(graftedGlycans)
        if uniprotID is not None:
            outputFileName = uniprotID + "_graft_summary.csv"
        elif args.user_localReceiverPath is not None:
            outputFileName = os.path.basename(args.user_localReceiverPath).rpartition(".")[0] + "_graft_summary.csv"
        else:
            outputFileName = "graft_summary.csv"
        saveCSVto = os.path.join(outputPath, outputFileName)
        df.to_csv(saveCSVto)

    