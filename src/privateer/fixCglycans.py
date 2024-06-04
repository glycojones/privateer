import os
import sys
import matplotlib.pyplot as plt
import json
import gemmi
import pandas as pd
import urllib.request
import re
sys.path.append("/y/people/lah583/privateer/src/privateer")
import grafter
from privateer import privateer_core as pvtcore

def file_paths(root_directory):
    """
    Function to find all the files in the database
    """
    filepathlist = []
    for root, dirs, files in os.walk(root_directory):
        for f in files:
            filepathlist.append(os.path.join(root,f))
    return filepathlist

def save_json(pdbcodes, receiver_paths, mtz_paths, output_paths, glycosylationss):
    assert(len(pdbcodes) == len(receiver_paths))
    assert(len(pdbcodes) == len(mtz_paths))
    assert(len(pdbcodes) == len(output_paths))
    #assert(len(pdbcodes) == len(glycosylationss[0]))
    output_file = f"/y/people/lah583/privateer/results/Fixing_c-glycans/c-mannosylation_sites.json"
    data = {}
    for pdbcode, receiver_path, mtz_path, output_path, glycosylations in zip(pdbcodes, receiver_paths, mtz_paths, output_paths, glycosylationss):
        data[str(pdbcode)]={"receiver_path": receiver_path, "mtz_path": mtz_path, "output_path": output_path ,"glycosylations": glycosylations}
    with open(output_file, "w") as output_file: 
        json.dump(data, output_file)

def find_mtz_path(pdbcode):
    output_path_1 = f"/vault/pdb_mtz_files/{pdbcode}.mtz"
    output_path_2 = f"/y/people/lah583/privateer/results/Fixing_c-glycans/receivers/{pdbcode}.mtz"
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

def get_RSCC(databasedir, pdbcode, protein_chain_ID, protein_res_ID):
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

def find_and_delete_glycans_to_replace_database(databasedir,pdbmirrordir,receiverdir,donordir,outputdir):
    filepathlist = file_paths(databasedir)
    data_out = {}
    pdbcodes = []
    #EM_structures = ["6dlw", "7nyc", "7nyd", "8de6", "8g04", "8u18"]
    for i in range(len(filepathlist)):
        jsonfile = filepathlist[i]
        jsonfilename = os.path.basename(jsonfile)
        pdbcode = jsonfilename.rpartition('.')[0]
        pdbfile = pdbmirrordir + f"/pdb/pdb{pdbcode}.ent.gz"
        mmciffile = pdbmirrordir+f"/mmCIF/{pdbcode}.cif.gz"
        try:
            st = gemmi.read_structure(pdbfile)
        except:
            st = gemmi.read_structure(mmciffile)
        save_structure = False
        glycosylations = []
        ms = []
        cs = []
        rs = []
        with open(jsonfile, 'r') as f:
            data_in = json.load(f)
            glycandata = data_in["glycans"]
            cglycans = glycandata["c-glycan"]
            for cglycan in cglycans:
                for sugar in cglycan["sugars"]:
                    if sugar["diagnostic"] != "yes" and  sugar["sugarId"].partition("-")[0] == "MAN" or sugar["sugarId"].partition("-")[0] == "BMA":
                        save_structure = True
                        sugarResId = sugar["sugarId"].rpartition("-")[2]
                        for m, model in enumerate(st):
                            for c, chain in enumerate(model):
                                for r, residue in enumerate(chain):
                                    if str(chain.name) == str(cglycan["rootSugarChainId"]) and int(residue.seqid.num) == int(sugarResId):
                                        ms.append(m)
                                        cs.append(c)
                                        rs.append(r) 
                        glycosylation = {}
                        glycosylation["donor_path"] = donordir + "/Alpha-D-Mannose.pdb"
                        glycosylation["glycan_index"] = 0
                        glycosylation["receiving_chain_index"] = cglycan["proteinChainId"]
                        glycosylation["receiving_res_index"] = cglycan["proteinResidueSeqnum"]
                        glycosylations.append(glycosylation)
        if save_structure:
            print(pdbcode)
            l = sorted(zip(rs, cs, ms))
            rs, cs, ms = zip(*l)
            for m, c, r in zip(ms[::-1], cs[::-1], rs[::-1]):
                print(m,c,r)
                del st[m][c][r]
            receiver_path = receiverdir + f"/{pdbcode}.pdb"
            mtz_path = find_mtz_path(pdbcode)
            output_path = outputdir + f"/{pdbcode}.pdb"
            st.write_pdb(receiver_path)
            data_out[str(pdbcode)]={"receiver_path": receiver_path, "mtz_path": mtz_path, "output_path": output_path ,"glycosylations": glycosylations}
            pdbcodes.append(pdbcode)
    cwd = os.getcwd()
    output_file = os.path.join(cwd, "c-mannosylation_sites.json")
    with open(output_file, "w") as output_file: 
        json.dump(data_out, output_file)
    return pdbcodes

def fix_glycans_using_database(databasedir,pdbmirrordir,receiverdir,donordir,outputdir):
    pdbcodes = find_and_delete_glycans_to_replace_database(databasedir,pdbmirrordir,receiverdir,donordir,outputdir)
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
        refined_pdb, refined_mtz = grafter._refine_grafted_glycans(grafted_pdb, mtzfile, outputloc, pdbout, mtzout, 10)
        if os.path.isfile(refined_pdb):
            print(f"Calculating RSCC for the grafted glycans...")
            try:
                graftedGlycans = grafter._calc_rscc_grafted_glycans(pdbout, mtzfile, graftedGlycans)
            except:
                print(f"Error calculating RSCC for grafted glycans in {pdbcode}")
        for graft in graftedGlycans:
            protein_chain_ID = graft["receiving_protein_residue_chain_PDBID"]
            protein_res_ID = graft["receiving_protein_residue_monomer_PDBID"]
            graft["OriginalRSCC"] = get_RSCC(databasedir, pdbcode, protein_chain_ID, protein_res_ID)
            graft["pdbcode"] = pdbcodes
            AllGlycans.append(graft)
        df_temp = pd.DataFrame(graftedGlycans)
        df_temp.to_csv(outputdir+f"/{pdbcode}_graft_summary.csv")
    df_out = pd.DataFrame.from_dict(AllGlycans)
    output_csv = outputdir + "/full_graft_summary.csv"
    df_out.to_csv(output_csv)


def fix_glycans_using_privateer(pdbmirrordir,receiverdir,donordir,outputdir):
    pdbfiles = file_paths(pdbmirrordir + "pdb")
    mmcifiles = file_paths(pdbmirrordir + "mmcif")
    for i in range(len(mmcifiles)):
        mmcifile = mmcifiles[i]
        filename = os.path.basename(mmcifile)
        pdbcode = filename.partition(".")[0]
        mtzfile = find_mtz_path(pdbcode)
        try:
            glycosylation = pvtcore.GlycosylationComposition(mmcifile, mtzfile, "FP,SIGFP")
        except:
            glycosylation = pvtcore.GlycosylationComposition_memsafe(mmcifile)
            mtzfile = None
        glycans = glycosylation.get_summary_of_detected_glycans()
        num_glycans = glycosylation.get_number_of_glycan_chains_detected() 
        for i in range(num_glycans):
            glycan_type = glycosylation.get_glyosylation_type(i)
            if glycan_type == "c-glycan":
                glycan = glycosylation.get_glycan(i)
                num_sugars = glycan.get_total_number_of_sugars()
                for j in range(num_sugars):
                    sugar = glycan.get_monosaccharide(j)
                    summary = sugar.get_sugar_summary()
                    name = summary["sugar_name_short"]
                    diagnostic = sugar.get_privateer_diagnostic()
                    if name == "MAN" and diagnostic != "yes" or name == "BMA":
                        





    




if __name__ == "__main__":
    databasedir = "/vault/privateer_database/pdb" # Location the privateer database is (if using database). If not using database, set to None
    pdbmirrordir = "/vault/pdb_mirror/data/structures/all" # Location the unedited pdb/mmcif files are stored (assumes pdb files are in subdir "pdb" and mmcif in subdir "mmcif")
    cwd = os.getcwd()
    receiverdir = os.path.join(cwd,"receivers") # Location you want to save edited pdb files with incorrect c-glycans removed
    if os.getenv("PRIVATEERDATA", None) is not None:
        rootdir = os.getenv("PRIVATEERDATA", None)
        donordir = os.path.join(rootdir, "glycan_donor_repertoire")
    else:
        donordir = "/y/people/lah583/privateer/data/glycan_donor_repertoire"
    outputdir = os.path.join(cwd,"output")
    if databasedir is not None:
        fix_glycans_using_database(databasedir,pdbmirrordir,receiverdir,donordir,outputdir)
    else:
        fix_glycans_using_privateer(pdbmirrordir,receiverdir,donordir,outputdir)

    



