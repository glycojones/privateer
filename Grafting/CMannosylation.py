import os
import sys
import matplotlib.pyplot as plt
import json
import gemmi
import pandas as pd
import urllib.request
import re
import argparse
import numpy as np
sys.path.append("/y/people/lah583/privateer/src/privateer")
import grafter
from privateer import privateer_core as pvtcore

def file_paths(root_directory, filetype=None):
    """
    Function to find all the files in the database
    """
    filepathlist = []
    for root, dirs, files in os.walk(root_directory):
        for f in files:
            if filetype is not None:
                if f.rpartition(".")[2] == filetype:
                    filepathlist.append(os.path.join(root,f))
            else:
                filepathlist.append(os.path.join(root,f))
    return filepathlist



def find_mtz_path(mtzdir,receiverdir,pdbcode, redo = False):
    if redo:
        #output_path_redo = os.path.join(receiverdir,pdbcode[1]+pdbcode[2],pdbcode, f"{pdbcode}_final.mtz")
        output_path_redo = os.path.join(mtzdir,pdbcode[1]+pdbcode[2],pdbcode,f"{pdbcode}_final.mtz")
        if os.path.exists(output_path_redo):
            return output_path_redo
        else:
            return ""
    else:
        output_path_1 = os.path.join(mtzdir,f"{pdbcode}.mtz")
        output_path_2 = os.path.join(receiverdir, f"{pdbcode}.mtz")
        url = f"https://edmaps.rcsb.org/coefficients/{pdbcode}.mtz"
        if os.path.exists(output_path_1):
            return output_path_1
        elif os.path.exists(output_path_2):
            return output_path_2
        else:
            try:
                urllib.request.urlretrieve(url, output_path_2)
                return output_path_2
            except Exception as e :
                return ""
            
def sf_mtz_path(mtzdir, pdbcode):
    dir = os.path.join(mtzdir,f"r{pdbcode}sf")
    path1 = os.path.join(dir,f"r{pdbcode}sf_fs_free.mtz")
    path2 = os.path.join(dir,f"r{pdbcode}sf_fs.mtz")
    path3 = os.path.join(dir,f"r{pdbcode}sf.mtz")
    if os.path.isfile(path1):
        return path1
    elif os.path.isfile(path2):
        return path2
    elif os.path.isfile(path3):
        return path3
    else:
        #print(f"Error finding mtz file for {pdbcode}")
        return ""

def get_RSCC_database(databasedir, pdbcode, protein_chain_ID, protein_res_ID):
    subdir = pdbcode[1] + pdbcode[2]
    jsonfile = f"{databasedir}/{subdir}/{pdbcode}.json"
    count = 0
    with open(jsonfile, "r") as f:
        data = json.load(f)
        glycandata = data["glycans"]
        cglycans = glycandata["c-glycan"]
        for cglycan in cglycans:
            if str(cglycan["proteinChainId"]) == str(protein_chain_ID) and int(cglycan["proteinResidueSeqnum"]) == int(protein_res_ID):
                for sugar in cglycan["sugars"]:
                    count +=1
                    RSCC = sugar["rscc"]
    if count > 1:
        print("Warning, original RSCC value for {pdbcode} protein chain ID {protein_chain_ID} and res ID {protein_res_ID} may be incorrect...")
    return RSCC

def get_consensus(inputchain:gemmi.Chain, inputresidue:gemmi.Residue) -> str:
    concat = ''
    rng = np.array([-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6])
    inputresidueNum = inputresidue.seqid.num
    residuerng = rng + inputresidueNum
    for seqid in residuerng:
        try:
            neighbour = inputchain[f'{seqid}'].__getitem__(0)
        except IndexError: 
            concat += '-'
            continue
        if ((neighbour != None)
            and (neighbour.label_seq != None) 
            and (gemmi.find_tabulated_residue(neighbour.name).is_amino_acid())
            and ((neighbour.label_seq - inputresidue.label_seq) in rng)
            and ((neighbour.seqid.num - inputresidue.seqid.num) in rng)):
            concat += gemmi.find_tabulated_residue(neighbour.name).one_letter_code
        else: 
            concat += '-'
    return concat

def check_consensus_sequence(sequence:str) -> bool:
    output = False
    for index in range(len(sequence)):
        segment = sequence[index:index+4]
        if re.match('W.{2}[W|C]',segment) != None:
            output = True
    return output

def check_expression_system_with_cif(receiverpath:str, pdbcode:str) -> tuple[list,bool]:
    if os.getenv("PRIVATEERDATA", None) is not None:
        datadir = os.getenv("PRIVATEERDATA", None)
    else:
        datadir = ''
    sourcefile = os.path.join(datadir,'taxon_summary.json') # include metazoan and toxoplasma taxonomy ids
    with open(sourcefile,'r') as f: 
        taxsrc = json.load(f)
        taxonids = taxsrc['taxonids']
        sciNames = taxsrc['sciNames']

    # pdb = os.path.basename(path).split('.')[0]

    chainlist,output = [], False

    doc = gemmi.cif.read(receiverpath)
    block = doc.sole_block()

    hostentities = block.find_values('_entity_src_gen.entity_id') # check which proteins are expressed in expression_system
    taxhosts = block.find_values('_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id')
    hostnames = block.find_values('_entity_src_gen.pdbx_host_org_scientific_name')

    natentities = block.find_values('_entity_src_nat.entity_id') # check which proteins are extracted directly from original_source
    taxnats = block.find_values('_entity_src_nat.pdbx_ncbi_taxonomy_id')
    natnames = block.find_values('_entity_src_nat.pdbx_organism_scientific')

    # totentities = block.find_values('_entity_poly.entity_id')
    polymers = block.find_values('_entity_poly.pdbx_seq_one_letter_code_can')
    chainforentity = block.find_values('_entity_poly.pdbx_strand_id')

    for (hent,taxhost,hostname) in zip(hostentities,taxhosts,hostnames): # PDB 3ccz host E.coli without taxid!

        try: taxhost = taxhost.replace(' ', ''); taxhost = int(taxhost)
        except (ValueError,TypeError,SyntaxError,RuntimeError): pass

        if (taxhost in taxonids) or (hostname.lower() in sciNames):
            # print(hent)
            try: proteinsequence = polymers[int(hent)-1]
            except (IndexError,ValueError,TypeError,SyntaxError) as e: 
                print(f'{pdbcode} with error = {e}')
                continue
            
            output = check_consensus_sequence(sequence=proteinsequence)
            chainname = chainforentity[int(hent)-1]
            if output == False: continue
            if len(chainname) > 1: chainname = chainname.split(',')
            if chainname not in chainlist: chainlist += chainname
        # else: print(f'{pdb} with taxhost = {taxhost} and hostname = {hostname}')

    for (natent,taxnat,natname) in zip(natentities,taxnats,natnames):

        try: taxnat = taxnat.replace(' ', ''); taxnat = int(taxnat)
        except (ValueError,TypeError,SyntaxError,RuntimeError): pass

        if (taxnat in taxonids) or (natname.lower() in sciNames):

            try: proteinsequence = polymers[int(natent)-1]
            except (IndexError,ValueError,TypeError,SyntaxError) as e: 
                print(f'{pdbcode} with error = {e}')
                continue

            output = check_consensus_sequence(sequence=proteinsequence)
            chainname = chainforentity[int(natent)-1]
            if output == False: continue
            if len(chainname) > 1: chainname = chainname.split(',')
            if chainname not in chainlist: chainlist += chainname

    return chainlist#,output


def get_targets_via_blob_search_and_consensus_sequence(ciffile:str,mtzfile:str,requestedchains:list,sequences:list,threshold = None) -> dict: # JUST DO BLOB_SEARCH AT PROTEINCHAINS HAVING WXXW|C
    avglength = 6.4118
    pdbid = os.path.basename(ciffile).split('.')[0]
    st = gemmi.read_structure(ciffile)
    if threshold == None:
        resolution = st.resolution
        threshold = 0.41868*resolution - 0.17116
    if not os.path.exists(mtzfile): 
        print(f"No mtzfile for {pdbid}")
        return None

    pointlist = []
    chainlist = []; residuelist = []; sum_density = []; mean_density = []; gridpoints = []
    consensus = []
    target_chainlist = []   ; target_residuelist = []
    for model in st:
        for chain in model:
            if requestedchains is not None:
                if chain.name not in requestedchains: 
                    continue
            for residue in chain:
                # GET ESTIMATED CENTROID WITH AVERAGE TRANSLATION LENGTH ~ 6.411 Å
                if (residue.name == 'TRP'):
                    if residue.label_seq != None: 
                        pentaseq = get_consensus(inputchain=chain,inputresidue=residue)
                        if re.search('W.{2}W',pentaseq[3:10]) == None: 
                            continue # JUST DO BLOB_SEARCH AT W RESIDUES FOLLOWING WXXW|C
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
                        pointlist.append(newpoint)
                        consensus.append(pentaseq)
                        residuelist.append(residue.seqid.num)
                        chainlist.append(chain.name)         
    # READ RECALCULATED MAP   
    try:
        mtz = gemmi.read_mtz_file(mtzfile)
    except RuntimeError as e:
        print(e)
        print(f'RuntimeError in reading mtz at {pdbid}')
        return
    grid = mtz.transform_f_phi_to_map('DELFWT', 'PHDELWT', sample_rate=2.0)
    grid.normalize()
    # COPY THE GRID --> SAMPLE GRIDPOINTS IN RADIUS = 3 Å ( BY SETTING VALUES THOSE GRIDPOINTS ON CLONED GRID ~ 10)
    # gr = grid.clone()
    start = 1000 # arbitary number
    for i, newpoint in enumerate(pointlist):
        gr = grid.clone()
        gr.set_points_around(newpoint, radius=3, value=start)
        pointgroup = np.argwhere(gr.array == start) # return positions in symmetry ???
        # TRACE BACK ONTO ORIGINAL GRID --> SUM DENSITY VALUES
        atomlist = []; values = []
        for point in pointgroup:
            atom = grid.get_point(point[0],point[1],point[2])
            atomlist.append(atom)
            value = grid.get_value(point[0],point[1],point[2])
            values.append(value)
        if atomlist: 
            s_density = np.sum(values).round(3)
            avg_density = np.mean(values).round(3)
        else: 
            s_density = 0
            avg_density = 0
        if avg_density > threshold: 
            target_chainlist.append(chainlist[i])
            target_residuelist.append(residuelist[i])
    CMannosylationConsensus = "[W]"
    output = []
    for item in sequences:
        currentChainIndex = item["index"]
        currentChainID = item["ChainID"]
        currentSequence = item["Sequence"]
        glycosylationTargets = []
        append = False
        for match in re.finditer(CMannosylationConsensus, currentSequence):
            if (currentSequence[match.start()] == item["Residues"][match.start()]["residueCode"]):
                for i in range(len(target_residuelist)):
                    blob_chainID = target_chainlist[i]
                    blob_resID = target_residuelist[i]
                    if (str(item["Residues"][match.start()]["residueSeqnum"]) == str(blob_resID)) and (str(currentChainID) == str(blob_chainID)):
                        glycosylationTargets.append({
                            "start": match.start(),
                            "end": match.end(),
                            "match": match.group()
                        })
                        append = True
        if append:
            output.append({
                "Sequence": currentSequence,
                "chainIndex": currentChainIndex,
                "currentChainID": currentChainID,
                "glycosylationTargets": glycosylationTargets,
            })
    return output

def find_and_delete_glycans_to_replace_database(databasedir,pdbmirrordir,mtzdir,receiverdir,donordir,outputdir,redo):
    filepathlist = file_paths(databasedir)
    data_out = {}
    pdbcodes = []
    failedlist = "failed.txt"
    with open(failedlist, "w") as myfile:
        myfile.write(f"A list of pdb codes that failed to find and delete potentially problematic c-glycans\n")
    for i in range(len(filepathlist)):
        jsonfile = filepathlist[i]
        jsonfilename = os.path.basename(jsonfile)
        pdbcode = jsonfilename.rpartition('.')[0]
        print(f"Checking glycans in {pdbcode}")
        if redo:
            pdbfile = pdbmirrordir + f"/{pdbcode[1]}{pdbcode[2]}/{pdbcode}/{pdbcode}_final.pdb"
            mmciffile = pdbmirrordir+f"/{pdbcode[1]}{pdbcode[2]}/{pdbcode}/{pdbcode}_final.cif"
            
        else:
            pdbfile = pdbmirrordir + f"/pdb/pdb{pdbcode}.ent.gz"
            mmciffile = pdbmirrordir+f"/mmCIF/{pdbcode}.cif.gz"
        try:
            st = gemmi.read_structure(pdbfile)
            ns = gemmi.NeighborSearch(model=st[0],cell=st.cell,max_radius=5.0).populate(include_h=False)
        except:
            try:
                st = gemmi.read_structure(mmciffile)
                ns = gemmi.NeighborSearch(model=st[0],cell=st.cell,max_radius=5.0).populate(include_h=False)
            except:
                print(f"Error opening structure {pdbcode}")
                with open(failedlist, "a") as myfile:
                    myfile.write(f"{pdbcode}\n")
        save_structure = False
        glycosylations = []
        ms = []
        cs = []
        rs = []
        with open(jsonfile, 'r') as f:
            data_in = json.load(f)
            glycandata = data_in["glycans"]
            cglycans = glycandata["c-glycan"]
            ligands = glycandata["ligand"]
            for cglycan in cglycans:
                for sugar in cglycan["sugars"]:
                    if (sugar["diagnostic"] != "yes" and  sugar["sugarId"].partition("-")[0] == "MAN") or sugar["sugarId"].partition("-")[0] == "BMA":
                        save_structure = True
                        sugarResId = sugar["sugarId"].rpartition("-")[2]
                        sugarChainID = sugar["sugarId"].partition("-")[2].partition("-")[0]
                        print(f"In {pdbcode} found sugar to fix chain {sugarChainID} res {sugarResId}")
                        for m, model in enumerate(st):
                            for c, chain in enumerate(model):
                                for r, residue in enumerate(chain):
                                    if str(chain.name) == str(sugarChainID) and int(residue.seqid.num) == int(sugarResId) and (str(residue.name) == "MAN" or str(residue.name) == "BMA"):
                                        ms.append(m)
                                        cs.append(c)
                                        rs.append(r) 
                        glycosylation = {}
                        glycosylation["donor_path"] = donordir + "/Alpha-D-Mannose.pdb"
                        glycosylation["glycan_index"] = 0
                        glycosylation["receiving_chain_index"] = cglycan["proteinChainId"]
                        glycosylation["receiving_res_index"] = cglycan["proteinResidueSeqnum"]
                        glycosylations.append(glycosylation)
            for ligand in ligands:
                for sugar in ligand["sugars"]:
                    if sugar["diagnostic"] != "yes":
                        if sugar["sugarId"].partition("-")[0] == "MAN" or sugar["sugarId"].partition("-")[0] == "BMA":
                            sugResId = sugar["sugarId"].rpartition("-")[2]
                            sugChainID = sugar["sugarId"].partition("-")[2].rpartition("-")[0]
                            try:
                                sugarResidue = st[0][sugChainID][sugResId]
                                for atom in sugarResidue:
                                    if atom.name == "C1":
                                        marks = ns.find_neighbors(atom)
                                        for mark in marks:
                                            cra = mark.to_cra(st[0])
                                            if cra.residue.name == "TRP" and cra.atom.name == "CD1":
                                                protseqnum = cra.residue.seqid.num
                                                protchainID = cra.chain.name
                                                save_structure = True
                                                for m, model in enumerate(st):
                                                    for c, chain in enumerate(model):
                                                        for r, residue in enumerate(chain):
                                                            # If the corresponding protein residue is a TRP
                                                            if str(chain.name) == str(sugChainID) and int(residue.seqid.num) == int(sugResId):
                                                                ms.append(m)
                                                                cs.append(c)
                                                                rs.append(r) 
                                                glycosylation = {}
                                                glycosylation["donor_path"] = donordir + "/Alpha-D-Mannose.pdb"
                                                glycosylation["glycan_index"] = 0
                                                glycosylation["receiving_chain_index"] = protchainID
                                                glycosylation["receiving_res_index"] = protseqnum
                                                glycosylations.append(glycosylation)
                            except:
                                save_structure = False
                                print(f"Error finding and deleting glycans in {pdbcode}")
                                with open(failedlist, "a") as myfile:
                                    myfile.write(f"{pdbcode}\n")
                                continue
        if save_structure and len(ms) == 0:
            print(f"Error finding and deleting glycans in {pdbcode}")     
            with open(failedlist, "a") as myfile:
                                    myfile.write(f"{pdbcode}\n")               
        if save_structure and len(ms) > 0:
            l = sorted(zip(rs, cs, ms))
            rs, cs, ms = zip(*l)
            for m, c, r in zip(ms[::-1], cs[::-1], rs[::-1]):
                del st[m][c][r]
            receiver_path = receiverdir + f"/{pdbcode}.pdb"
            #mtz_path = find_mtz_path(mtzdir,receiverdir,pdbcode,redo)
            mtz_path = sf_mtz_path(mtzdir,pdbcode)
            output_path = outputdir + f"/{pdbcode}.pdb"
            st.write_pdb(receiver_path)
            data_out[str(pdbcode)]={"receiver_path": receiver_path, "mtz_path": mtz_path, "output_path": output_path ,"glycosylations": glycosylations}
            pdbcodes.append(pdbcode)
    cwd = os.getcwd()
    output_file = os.path.join(cwd, "c-mannosylation_sites.json")
    with open(output_file, "w") as output_file: 
        json.dump(data_out, output_file)
    return pdbcodes

def find_and_delete_glycans_to_replace_privateer(pdbmirrordir,mtzdir,receiverdir,donordir,outputdir,redo):
    failedlist = "failed.txt"
    with open(failedlist, "w") as myfile:
        myfile.write(f"A list of pdb codes that failed to find and delete potentially problematic c-glycans\n")
    if redo:
        pdbfiles = file_paths(pdbmirrordir, ".pdb")
        mmcifiles = file_paths(pdbmirrordir, ".cif")
    else:
        pdbfiles = file_paths(os.path.join(pdbmirrordir , "pdb"))
        mmcifiles = file_paths(os.path.join(pdbmirrordir , "mmCIF"))
    data_out = {}
    pdbcodes = []
    for i in range(len(mmcifiles)):
        mmcifile = mmcifiles[i]
        filename = os.path.basename(mmcifile)
        if redo:
            pdbcode = filename.partition("_")[0]
            if mmcifile != os.path.join(pdbmirrordir,f"{pdbcode[1]}{pdbcode[2]}",f"{pdbcode}_final.cif"):
                continue
        else:
            pdbcode = filename.partition(".")[0]
        mtzfile = sf_mtz_path(mtzdir, pdbcode) # find_mtz_path(mtzdir,receiverdir,pdbcode,redo)
        if len(mtzfile) < 1:
            continue
        print(f"Looking for c-glycans to fix in {pdbcode}")
        try:
            st = gemmi.read_structure(mmcifile)
            ns = gemmi.NeighborSearch(model=st[0],cell=st.cell,max_radius=5.0).populate(include_h=False)
        except:
            print(f"Error opening structure {pdbcode}")
            with open(failedlist, "a") as myfile:
                myfile.write(f"{pdbcode}\n")
            continue
        save_structure = False
        glycosylations = []
        ms = []
        cs = []
        rs = []
        #try:
        glycosylation = pvtcore.GlycosylationComposition(mmcifile, mtzfile, "FP,SIGFP", nthreads=4)
        # except:
        #     try:
        #         glycosylation = pvtcore.GlycosylationComposition_memsafe(mmcifile,nthreads=4)
        #     except:
        #         print(f"Error running privateer on structure {pdbcode}")
        #         with open(failedlist, "a") as myfile:
        #             myfile.write(f"{pdbcode}\n")
        #         continue
            
        glycans = glycosylation.get_summary_of_detected_glycans()
        num_glycans = glycosylation.get_number_of_glycan_chains_detected() 
        for glycan in glycans:
            glycan_type = glycan["GlycosylationType"]
            if glycan_type == "c-glycan":
                num_sugars = glycan.get_total_number_of_sugars()
                root_info = glycan.get_root_info()
                for j in range(num_sugars):
                    sugar = glycan.get_monosaccharide(j)
                    summary = sugar.get_sugar_summary()
                    name = summary["sugar_name_short"]
                    diagnostic = sugar.get_privateer_diagnostic()
                    if name == "MAN" and diagnostic != "yes" or name == "BMA":
                        save_structure = True
                        for m, model in enumerate(st):
                            for c, chain in enumerate(model):
                                for r, residue in enumerate(chain):
                                    if str(chain.name) == str(summary["sugar_pdb_chain"]) and int(residue.seqid.num) == int(summary["sugar_seqnum"]):
                                        ms.append(m)
                                        cs.append(c)
                                        rs.append(r) 
                        glycosylation = {}
                        glycosylation["donor_path"] = donordir + "/Alpha-D-Mannose.pdb"
                        glycosylation["glycan_index"] = 0
                        glycosylation["receiving_chain_index"] = root_info["ProteinChainID"]
                        glycosylation["receiving_res_index"] = root_info["ProteinResidueSeqnum"]
                        glycosylations.append(glycosylation)
            if glycan_type == "ligand":
                num_sugars = glycan.get_total_number_of_sugars()
                root_info = glycan.get_root_info()
                for j in range(num_sugars):
                    sugar = glycan.get_monosaccharide(j)
                    summary = sugar.get_sugar_summary()
                    name = summary["sugar_name_short"]
                    sugChainID = summary["sugar_pdb_chain"]
                    sugResID = summary["sugar_seqnum"]
                    diagnostic = sugar.get_privateer_diagnostic()
                    if diagnostic != "yes":
                        if name == "MAN" or name == "BMA":
                            try:
                                sugarResidue = st[0][sugChainID][sugResID]
                                for atom in sugarResidue:
                                    if atom.name == "C1":
                                        marks = ns.find_neighbors(atom)
                                        for mark in marks:
                                            cra = mark.to_cra(st[0])
                                            if cra.residue.name == "TRP" and cra.atom.name == "CD1":
                                                protseqnum = cra.residue.seqid.num
                                                protchainID = cra.chain.name
                                                save_structure = True
                                                for m, model in enumerate(st):
                                                    for c, chain in enumerate(model):
                                                        for r, residue in enumerate(chain):
                                                            # If the corresponding protein residue is a TRP
                                                            if str(chain.name) == str(sugChainID) and int(residue.seqid.num) == int(sugResID):
                                                                ms.append(m)
                                                                cs.append(c)
                                                                rs.append(r) 
                                glycosylation = {}
                                glycosylation["donor_path"] = donordir + "/Alpha-D-Mannose.pdb"
                                glycosylation["glycan_index"] = 0
                                glycosylation["receiving_chain_index"] = root_info["ProteinChainID"]
                                glycosylation["receiving_res_index"] = root_info["ProteinResidueSeqnum"]
                                glycosylations.append(glycosylation)
                            except:
                                save_structure = False
                                print(f"Error finding and deleting glycans in {pdbcode}")
                                with open(failedlist, "a") as myfile:
                                    myfile.write(f"{pdbcode}\n")
                                continue
        if save_structure and len(ms) == 0:
            print(f"Error finding and deleting glycans in {pdbcode}")     
            with open(failedlist, "a") as myfile:
                                    myfile.write(f"{pdbcode}\n")  
        if save_structure and len(ms) > 0:
            l = sorted(zip(rs, cs, ms))
            rs, cs, ms = zip(*l)
            for m, c, r in zip(ms[::-1], cs[::-1], rs[::-1]):
                del st[m][c][r]
            receiver_path = receiverdir + f"/{pdbcode}.pdb"
            mtz_path = find_mtz_path(mtzdir,receiverdir,pdbcode)
            output_path = outputdir + f"/{pdbcode}.pdb"
            st.write_pdb(receiver_path)
            data_out[str(pdbcode)]={"receiver_path": receiver_path, "mtz_path": mtz_path, "output_path": output_path ,"glycosylations": glycosylations}
            pdbcodes.append(pdbcode)
    cwd = os.getcwd()
    output_file = os.path.join(cwd, "c-mannosylation_sites.json")
    with open(output_file, "w") as output_file: 
        json.dump(data_out, output_file)
    return pdbcodes

def fix_Cglycans(databasedir,pdbmirrordir,mtzdir,receiverdir,donordir,outputdir,redo):
    if databasedir is not None:
        pdbcodes = find_and_delete_glycans_to_replace_database(databasedir,pdbmirrordir,mtzdir,receiverdir,donordir,outputdir,redo)
    else:
        pdbcodes = find_and_delete_glycans_to_replace_privateer(pdbmirrordir,mtzdir,receiverdir,donordir,outputdir,redo)
    df_in = pd.read_json("c-mannosylation_sites.json")
    AllGlycans = []
    for pdbcode in pdbcodes:
        print(f"Regrafting structure {pdbcode}...")
        schema = df_in[pdbcode]
        receiverpath = schema["receiver_path"]
        outputpath = schema["output_path"]
        sequences = grafter._get_sequences_in_receiving_model(receiverpath)
        CMannosylationConsensus = "[W]"
        targets = []
        for item in sequences:
            currentChainIndex = item["index"]
            currentChainID = item["ChainID"]
            currentSequence = item["Sequence"]
            glycosylationTargets = []
            for match in re.finditer(CMannosylationConsensus, currentSequence):
                if (currentSequence[match.start()] == item["Residues"][match.start()]["residueCode"]):
                    for glycan in schema["glycosylations"]:
                        target_chainID = glycan["receiving_chain_index"]
                        target_resID = glycan["receiving_res_index"]
                        donorpath = glycan["donor_path"]
                        if (item["Residues"][match.start()]["residueSeqnum"] == target_resID) and (currentChainID == target_chainID):
                            glycosylationTargets.append({
                                "start": match.start(),
                                "end": match.end(),
                                "match": match.group()
                            })
            targets.append({
                "Sequence": currentSequence,
                "chainIndex": currentChainIndex,
                "currentChainID": currentChainID,
                "glycosylationTargets": glycosylationTargets,
            })
        try:
            graftedGlycans = grafter._glycosylate_receiving_model_using_consensus_seq(receiverpath, donorpath, outputpath, targets, True, False, False)
        except:
            print(f"Error grafting glycans to pdb {pdbcode}. Skipping graft...")
            continue
        try:
            grafter._copy_metadata(receiverpath, outputpath, graftedGlycans)
        except:
            print(f"Failed to copy metadata for {pdbcode}")
        try:
            grafter._make_connection_between_protein_and_glycan(outputpath)
        except:
            print(f"Failed to generate link between TRP-MAN for file {pdbcode}")
        grafted_pdb = schema["output_path"]
        mtzfile = schema["mtz_path"]
        outputloc = outputdir
        pdbout = outputdir+f"/{pdbcode}_grafted.pdb"
        mtzout = outputdir+f"/{pdbcode}_grafted.mtz"
        print(f"Refining grafted strucutre...")
        st = gemmi.read_structure(grafted_pdb)
        resolution = st.resolution
        try:
            refined_pdb, refined_mtz = grafter._refine_grafted_glycans(grafted_pdb, mtzfile, outputloc, pdbout, mtzout, 20, resolution)
        except:
            print(f"Error refining grafted structure {pdbcode}")
        if os.path.isfile(refined_pdb):
            print(f"Calculating RSCC for the grafted glycans...")
            try:
                graftedGlycans = grafter._calc_rscc_grafted_glycans(pdbout, mtzfile, graftedGlycans)
            except:
                print(f"Error calculating RSCC for grafted glycans in {pdbcode}")
                continue
        for graft in graftedGlycans:
            protein_chain_ID = graft["receiving_protein_residue_chain_PDBID"]
            protein_res_ID = graft["receiving_protein_residue_monomer_PDBID"]
            graft["OriginalRSCC"] = get_RSCC_database(databasedir, pdbcode, protein_chain_ID, protein_res_ID)
            graft["pdbcode"] = pdbcode
            AllGlycans.append(graft)
        df_temp = pd.DataFrame(graftedGlycans)
        df_temp.to_csv(outputdir+f"/{pdbcode}_graft_summary.csv")
    df_out = pd.DataFrame.from_dict(AllGlycans)
    output_csv = outputdir + "/full_graft_summary.csv"
    df_out.to_csv(output_csv)
    return

def find_and_graft_Cglycans(receiverdir,mtzdir,donordir,outputdir,redo,graftedlist,savesummary):
    donorpath = os.path.join(donordir, "Alpha-D-Mannose.pdb")
    receivers = file_paths(receiverdir)
    AllGlycans = []
    for receiverpath in receivers:
        if graftedlist is not None:
            with open(graftedlist) as myfile:
                if receiverpath in myfile.read():
                    continue # If the pdb has already been grafted, do not graft again
        filename = os.path.basename(receiverpath)
        if redo:
            if filename.partition("_")[2] != "final.pdb":
                continue
            pdbcode = filename.partition("_")[0]
            ciffile = receiverpath.rpartition(".")[0] + ".cif"
            expsysfile = f"/vault/pdb_mirror/data/structures/all/mmCIF/{pdbcode}.cif.gz"
        else:
            pdbcode = filename.partition(".")[0]
            ciffile = receiverpath.rpartition(".")[0] + ".cif"
            expsysfile = ciffile
        if "pdb" in pdbcode:
            pdbcode = filename.partition("pdb")[2]
        if redo and receiverpath != os.path.join(receiverdir,f"{pdbcode[1]}{pdbcode[2]}",f"{pdbcode}",f"{pdbcode}_final.pdb"):
            continue
        with open(graftedlist, "a") as myfile:
            myfile.write(receiverpath)
        mtzpath = find_mtz_path(mtzdir,receiverdir,pdbcode,redo)
        outputpath = os.path.join(outputdir,f"{pdbcode}.pdb")
        try:
            requestedchains = check_expression_system_with_cif(expsysfile,pdbcode)
            if not requestedchains: 
                with open(graftedlist, "a") as myfile:
                        myfile.write("\tWrong expression system")
                        myfile.write("\n")
                continue 
        except:
            print("Error checking expression system. Proceeding with grafting regardless of expression system...")
            requestedchains = None
        sequences = grafter._get_sequences_in_receiving_model(receiverpath)
        targets = get_targets_via_blob_search_and_consensus_sequence(ciffile, mtzpath, requestedchains, sequences, 0.08)
        #targets_2 = grafter._get_CMannosylation_targets_via_water_search(receiverpath, sequences) #FLAG: remove water search???
        #for target_1 in targets_1:
        #    for target_2 in targets_2[:]:
        #        if target_1["chainIndex"]==target_2["chainIndex"] and target_1["currentChainID"]==target_2["currentChainID"]:
        #            targets_2.remove(target_2)
        #targets = targets_1 + targets_2

        removeclashes = True
        if targets is None:
            with open(graftedlist, "a") as myfile:
                    myfile.write("\tNo C-Mannosylation Targets found")
                    myfile.write("\n")
            continue
        elif len(targets) < 1:
            with open(graftedlist, "a") as myfile:
                    myfile.write("\tNo C-Mannosylation Targets found")
                    myfile.write("\n")
            continue
        try:
            graftedGlycans = grafter._glycosylate_receiving_model_using_consensus_seq(
                receiverpath, donorpath, outputpath, targets, True, False, removeclashes)
            grafter._copy_metadata(receiverpath, outputpath, graftedGlycans)
            grafter._make_connection_between_protein_and_glycan(outputpath)
            if graftedlist is not None:
                with open(graftedlist, "a") as myfile:
                    myfile.write("\tCompleted")
                    myfile.write("\n")
        except:
            print(f"Error grafting glycans to pdb {pdbcode}. Skipping graft...")
            if graftedlist is not None:
                with open(graftedlist, "a") as myfile:
                    myfile.write("\tFailed")
                    myfile.write("\n")
            continue
        if len(graftedGlycans) < 1:
            with open(graftedlist, "a") as myfile:
                    myfile.write("\tNo C-Mannosylation Targets found")
                    myfile.write("\n")
            continue
        #FLAG: Here want a step removing cases with too many clashes??? Or just stick to removing grafts with any clashes???
        print(f"Refining grafted strucutre...")
        st = gemmi.read_structure(outputpath)
        resolution = st.resolution
        refined_pdb, refined_mtz = grafter._refine_grafted_glycans(outputpath, mtzpath, outputdir, outputdir+f"/{pdbcode}_refined.pdb", outputdir+f"/{pdbcode}_refined.mtz", 20, resolution)
        if os.path.isfile(refined_pdb):
            os.remove(refined_mtz)
            os.remove(outputdir+f"/{pdbcode}_refined.mmcif")
            try:
                pdbout = os.path.join(outputdir, f"{pdbcode}_removed_waters.pdb")
                grafter._remove_waters_close_to_TRP(refined_pdb, pdbout)
                os.remove(refined_pdb)
            except:
                pdbout = refined_pdb
                print(f"Failed to remove waters close to TRP in {pdbout}")
            graftedGlycans = grafter._calc_rscc_grafted_glycans(pdbout, mtzpath, graftedGlycans)
            graftedGlycans = grafter._remove_grafted_glycans(pdbout, mtzpath, graftedGlycans, outputdir, rscc_threshold = 0.5)
        for graft in graftedGlycans:
            graft["pdbcode"] = pdbcode
            AllGlycans.append(graft)
        if savesummary:
            df_single = pd.DataFrame(graftedGlycans)
            df_single.to_csv(outputdir+f"/{pdbcode}_graft_summary.csv")
            df_temp = pd.DataFrame.from_dict(AllGlycans)
            temp_csv = outputdir + "/full_graft_summary_temp.csv"
            df_temp.to_csv(temp_csv)
    if savesummary:
        df_all = pd.DataFrame.from_dict(AllGlycans)
        output_csv = outputdir + "/full_graft_summary.csv"
        df_all.to_csv(output_csv)
        temp_csv = outputdir + "/full_graft_summary_temp.csv"
        if os.path.isfile(temp_csv):
            os.remove(temp_csv)
    return 

def graft_Cglycans_from_csv(csvfile,receiverdir,mtzdir,donordir,outputdir,redo,graftedlist,savesummary):
    df_in = pd.read_csv(csvfile)
    pdbcodes = df_in["pdbid"].unique()
    donorpath = os.path.join(donordir, "Alpha-D-Mannose.pdb")
    for pdbcode in pdbcodes:
        if graftedlist is not None:
            with open(graftedlist) as myfile:
                if pdbcode in myfile.read():
                    continue
            with open(graftedlist, "a") as myfile:
                myfile.write(pdbcode)
        if redo:
            receiverpath = os.path.join(receiverdir,f"{pdbcode[1]}{pdbcode[2]}",f"{pdbcode}",f"{pdbcode}_final.pdb")
        else:
            cifpath = os.path.join(receiverdir,f"{pdbcode}.cif.gz")
            pdbpath = os.path.join(receiverdir,f"pdb{pdbcode}.ent.gz")
            if os.path.isfile(cifpath):
                receiverpath = cifpath
            elif os.path.isfile(pdbpath):
                receiverpath = pdbpath
            else:
                print(f"Error opening receiver structure {pdbcode}")
                if graftedlist is not None:
                    with open(graftedlist, "a") as myfile:
                        myfile.write("\tError opening receiver structure")
                        myfile.write("\n")
                continue
        mtzpath = sf_mtz_path(mtzdir, pdbcode)
        if len(mtzpath)<1:
            if graftedlist is not None:
                with open(graftedlist, "a") as myfile:
                    myfile.write("\tError opening receiver mtzfile")
                    myfile.write("\n")
        outputpath = os.path.join(outputdir,f"{pdbcode}.pdb")
        try:
            sequences = grafter._get_sequences_in_receiving_model(receiverpath)
        except:
            print(f"Failed to get sequences for {pdbcode}. Skipping...")
            if graftedlist is not None:
                with open(graftedlist, "a") as myfile:
                    myfile.write("\tFailed pvtmodelling.Builder()")
                    myfile.write("\n")
            continue
        CMannosylationConsensus = "[W]"
        targets = []
        for index, row in df_in.iterrows():
            if row["pdbid"] ==  pdbcode and row["status"] == "yes":
                resolution = row["resolution"]
                for item in sequences:
                    currentChainIndex = item["index"]
                    currentChainID = item["ChainID"]
                    currentSequence = item["Sequence"]
                    glycosylationTargets = []
                    for match in re.finditer(CMannosylationConsensus, currentSequence):
                        if (currentSequence[match.start()] == item["Residues"][match.start()]["residueCode"]):
                            target_chainID = row["authchain"]
                            target_resID = row["authresidue"]
                            if (item["Residues"][match.start()]["residueSeqnum"] == target_resID) and (currentChainID == target_chainID):
                                glycosylationTargets.append({
                                    "start": match.start(),
                                    "end": match.end(),
                                    "match": match.group()
                                })
                    targets.append({
                        "Sequence": currentSequence,
                        "chainIndex": currentChainIndex,
                        "currentChainID": currentChainID,
                        "glycosylationTargets": glycosylationTargets,
                    })
        removeclashes = True
        if targets is None:
            if graftedlist is not None:
                with open(graftedlist, "a") as myfile:
                        myfile.write("\tNo C-Mannosylation Targets found")
                        myfile.write("\n")
            continue
        elif len(targets) < 1:
            if graftedlist is not None:
                with open(graftedlist, "a") as myfile:
                        myfile.write("\tNo C-Mannosylation Targets found")
                        myfile.write("\n")
            continue
        try:
            graftedGlycans = grafter._glycosylate_receiving_model_using_consensus_seq(
                receiverpath, donorpath, outputpath, targets, True, False, removeclashes)
            grafter._copy_metadata(receiverpath, outputpath, graftedGlycans)
            grafter._make_connection_between_protein_and_glycan(outputpath)
            if graftedlist is not None:
                with open(graftedlist, "a") as myfile:
                    myfile.write("\tGrafted")
        except:
            print(f"Error grafting glycans to pdb {pdbcode}. Skipping graft...")
            if graftedlist is not None:
                with open(graftedlist, "a") as myfile:
                    myfile.write("\tFailed")
                    myfile.write("\n")
            continue
        if len(graftedGlycans) < 1:
            if graftedlist is not None:
                with open(graftedlist, "a") as myfile:
                        myfile.write("\tNo C-Mannosylation Targets found")
                        myfile.write("\n")
            continue
        print(f"Refining grafted strucutre...")
        refined_pdb, refined_mtz = grafter._refine_grafted_glycans(outputpath, mtzpath, outputdir, outputdir+f"/{pdbcode}_refined.pdb", outputdir+f"/{pdbcode}_refined.mtz", 20, resolution)
        if os.path.isfile(refined_pdb):
            if graftedlist is not None:
                with open(graftedlist, "a") as myfile:
                    myfile.write("\tRefined")
            os.remove(refined_mtz)
            os.remove(outputdir+f"/{pdbcode}_refined.mmcif")
            try:
                pdbout = os.path.join(outputdir, f"{pdbcode}_removed_waters.pdb")
                grafter._remove_waters_close_to_TRP(refined_pdb, pdbout)
                os.remove(refined_pdb)
            except:
                pdbout = refined_pdb
                print(f"Failed to remove waters close to TRP in {pdbout}")
            try:
                graftedGlycans = grafter._calc_rscc_grafted_glycans(pdbout, mtzpath, graftedGlycans)
                graftedGlycans = grafter._remove_grafted_glycans(pdbout, mtzpath, graftedGlycans, outputdir, rscc_threshold = 0.5)
                if graftedlist is not None:
                    with open(graftedlist, "a") as myfile:
                        myfile.write("\n")
            except:
                print(f"Failed to calculate RSCC in {pdbout}")
                if graftedlist is not None:
                    with open(graftedlist, "a") as myfile:
                        myfile.write("\tFailed RSCC calc")
                        myfile.write("\n")
        else: 
            if graftedlist is not None:
                with open(graftedlist, "a") as myfile:
                    myfile.write("\tFailed refinement")
                    myfile.write("\n")
        for graft in graftedGlycans:
            graft["pdbcode"] = pdbcode
        if savesummary:
            df_single = pd.DataFrame(graftedGlycans)
            df_single.to_csv(outputdir+f"/{pdbcode}_graft_summary.csv")
            full_csv = outputdir + "/full_graft_summary.csv"
            if os.path.isfile(full_csv):
                df_single.to_csv(full_csv, mode="a", index=False, header=False)
            else:
                df_single.to_csv(full_csv, index=False)
    return 



if __name__ == "__main__":
    defaultdatabasedir = None #"/vault/privateer_database/pdb" # Location the privateer database is (if using database). If not using database, set to None
    defaultpdbmirrordir = "/vault/pdb_mirror/data/structures/all" # Location the unedited pdb/mmcif files are stored (assumes pdb files are in subdir "pdb" and mmcif in subdir "mmcif")
    defaultmtzdir = "/vault/pdb_mtz_files"
    defaultcsvsites = None
    cwd = os.getcwd()
    defaultreceiverdir = os.path.join(cwd,"receivers") # Location you want to save edited pdb files with incorrect c-glycans removed
    if os.getenv("PRIVATEERDATA", None) is not None:
        rootdir = os.getenv("PRIVATEERDATA", None)
        defaultdonordir = os.path.join(rootdir, "glycan_donor_repertoire")
    else:
        defaultdonordir = "/y/people/lah583/privateer/data/glycan_donor_repertoire"
    defaultoutputdir = os.path.join(cwd,"output")
    defaultredo = False
    parser = argparse.ArgumentParser(
        prog="CMannosylation.py",
        usage=
        "%(prog)s [options]. Basic usage: 'python CMannosylation.py -mode fix' or 'python CMannnosylation.py -mode find'.",
        description=
        f"Either fix problematic cglycans in structures or find and graft potential cglycans in structures.",
        epilog=
        f"If mode of operation is not provided the script will not run. All other parameters are optional and if not specified will use default parameters.",
    )
    parser.add_argument(
        "-mode",
        action="store",
        default=None,
        dest="mode",
        help=
        f"Mode of operation: 'fix' to find and fix problematic cglycans in structures, 'find' to find potential c-mannosylation sited and graft mannose there.",
    )
    parser.add_argument(
        "-databasedir",
        action="store",
        default=defaultdatabasedir,
        dest="databasedir",
        help=
        f"Path to the locally stored privateer database. If not set, defaults to {defaultdatabasedir}. If -mode is set to 'find' this parameter is ignored.",
    )
    parser.add_argument(
        "-inputstructuredir",
        action="store",
        default=defaultpdbmirrordir,
        dest="pdbmirrordir",
        help=
        f"Path to the locally stored structures with cglycans to fix. If not set, defaults to {defaultpdbmirrordir}. If -mode is set to 'find' this parameter is ignored.",
    )
    parser.add_argument(
        "-mtzdirdir",
        action="store",
        default=defaultmtzdir,
        dest="mtzdir",
        help=
        f"Path to the locally stored mtz filed. If not set, defaults to {defaultmtzdir}.",
    )
    parser.add_argument(
        "-graftsitescsv",
        action="store",
        default=defaultcsvsites,
        dest="graftsitescsv",
        help=
        f"Path to csvfile which contains results from blob search to graft. If not set, defaults to {defaultcsvsites}.",
    )
    parser.add_argument(
        "-receiverdir",
        action="store",
        default=defaultreceiverdir,
        dest="receiverdir",
        help=
        f"If mode is set to 'fix' this is the location structures are saved once problematic cglycans are removed before grafting. If mode is set to 'find' this is the location of the original input structures to find potential c-mannosylation and graft. If not set, defaults to {defaultreceiverdir}.",
    )
    parser.add_argument(
        "-donordir",
        action="store",
        default=defaultdonordir,
        dest="donordir",
        help=
        f"Path to the locally stored alpha-D-mannose which will be grafted onto the structures. If not set, defaults to {defaultdonordir}.",
    )
    parser.add_argument(
        "-outputdir",
        action="store",
        default=defaultoutputdir,
        dest="outputdir",
        help=
        f"Path to save final fixed/grafted structures. If not set, defaults to {defaultoutputdir}.",
    )
    parser.add_argument(
        "--redo",
        action="store_true",
        dest="redo",
        help=
        f"Boolean to say whether running on pdbredo (True) or not (False). If not set, defaults to {defaultredo}.",
    )
    args = parser.parse_args()
    if args.mode == 'fix':
        fix_Cglycans(args.databasedir,args.pdbmirrordir,args.mtzdir,args.receiverdir,args.donordir,args.outputdir,args.redo)
    elif args.mode == 'find':
        if args.graftsitescsv is not None:
            graft_Cglycans_from_csv(args.graftsitescsv,args.receiverdir,args.mtzdir,args.donordir,args.outputdir,args.redo,"grafted_pdbs.txt",True)
        else:
            find_and_graft_Cglycans(args.receiverdir,args.mtzdir,args.donordir,args.outputdir,args.redo,"grafted_pdbs.txt",True)
    else:
        print("Mode of operation not specified. Exiting...")



    



