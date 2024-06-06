import os
import sys
import matplotlib.pyplot as plt
import json
import gemmi
import pandas as pd
import urllib.request
import re
import argparse
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

def find_mtz_path(mtzdir,receiverdir,pdbcode):
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

def find_and_delete_glycans_to_replace_database(databasedir,pdbmirrordir,mtzdir,receiverdir,donordir,outputdir):
    filepathlist = file_paths(databasedir)
    data_out = {}
    pdbcodes = []
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

def find_and_delete_glycans_to_replace_privateer(pdbmirrordir,mtzdir,receiverdir,donordir,outputdir):
    pdbfiles = file_paths(pdbmirrordir + "pdb")
    mmcifiles = file_paths(pdbmirrordir + "mmcif")
    data_out = {}
    pdbcodes = []
    for i in range(len(mmcifiles)):
        mmcifile = mmcifiles[i]
        filename = os.path.basename(mmcifile)
        pdbcode = filename.partition(".")[0]
        mtzfile = find_mtz_path(mtzdir,receiverdir,pdbcode)
        st = gemmi.read_structure(mmcifile)
        save_structure = False
        glycosylations = []
        ms = []
        cs = []
        rs = []
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
        if save_structure:
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

def fix_Cglycans(databasedir,pdbmirrordir,mtzdir,receiverdir,donordir,outputdir):
    if databasedir is not None:
        pdbcodes = find_and_delete_glycans_to_replace_database(databasedir,pdbmirrordir,mtzdir,receiverdir,donordir,outputdir)
    else:
        pdbcodes = find_and_delete_glycans_to_replace_privateer(pdbmirrordir,mtzdir,receiverdir,donordir,outputdir)
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
        refined_pdb, refined_mtz = grafter._refine_grafted_glycans(grafted_pdb, mtzfile, outputloc, pdbout, mtzout, 20)
        if os.path.isfile(refined_pdb):
            print(f"Calculating RSCC for the grafted glycans...")
            try:
                graftedGlycans = grafter._calc_rscc_grafted_glycans(pdbout, mtzfile, graftedGlycans)
            except:
                print(f"Error calculating RSCC for grafted glycans in {pdbcode}")
        for graft in graftedGlycans:
            protein_chain_ID = graft["receiving_protein_residue_chain_PDBID"]
            protein_res_ID = graft["receiving_protein_residue_monomer_PDBID"]
            graft["OriginalRSCC"] = get_RSCC_database(databasedir, pdbcode, protein_chain_ID, protein_res_ID)
            graft["pdbcode"] = pdbcodes
            AllGlycans.append(graft)
        df_temp = pd.DataFrame(graftedGlycans)
        df_temp.to_csv(outputdir+f"/{pdbcode}_graft_summary.csv")
    df_out = pd.DataFrame.from_dict(AllGlycans)
    output_csv = outputdir + "/full_graft_summary.csv"
    df_out.to_csv(output_csv)
    return

def find_and_graft_Cglycans(receiverdir,mtzdir,donordir,outputdir):
    donorpath = os.path.join(donordir, "Alpha-D-Mannose.pdb")
    receivers = file_paths(receiverdir)
    AllGlycans = []
    for receiverpath in receivers:
        filename = os.path.basename(receiverpath)
        pdbcode = filename.partition(".")[0]
        if "pdb" in pdbcode:
            pdbcode = filename.partition("pdb")[2]
        mtzpath =find_mtz_path(mtzdir,receiverdir,pdbcode)
        outputpath = os.path.join(outputdir,f"{pdbcode}.pdb")
        sequences = grafter._get_sequences_in_receiving_model(receiverpath)
        #FLAG: Other criteria here to lower number of false positives. Sequence? Things from metadata about the structure?
        targets_1 = grafter._get_CMannosylation_targets_via_blob_search(receiverpath, mtzpath, sequences, avg_dens_threshold = 0.08) #FLAG: threshold in this function needs to change.
        targets_2 = grafter._get_CMannosylation_targets_via_water_search(receiverpath, sequences) #FLAG: remove water search???
        for target_1 in targets_1:
            for target_2 in targets_2[:]:
                if target_1["chainIndex"]==target_2["chainIndex"] and target_1["currentChainID"]==target_2["currentChainID"]:
                    targets_2.remove(target_2)
        targets = targets_1 + targets_2
        removeclashes = False
        try:
            graftedGlycans = grafter._glycosylate_receiving_model_using_consensus_seq(
                receiverpath, donorpath, outputpath, targets, True, False, removeclashes)
        except:
            print(f"Error grafting glycans to pdb {pdbcode}. Skipping graft...")
        #FLAG: Here want a step removing cases with too many clashes??? Or just stick to removing grafts with any clashes???
        try:
            grafter._copy_metadata(receiverpath, outputpath, graftedGlycans)
        except:
            print(f"Failed to copy metadata for {pdbcode}")
        try:
            grafter._make_connection_between_protein_and_glycan(outputpath)
        except:
            print(f"Failed to generate link between TRP-MAN for file {pdbcode}")
        print(f"Refining grafted strucutre...")
        refined_pdb, refined_mtz = grafter._refine_grafted_glycans(outputpath, mtzpath, outputdir, outputdir+f"/{pdbcode}_refined.pdb", outputdir+f"/{pdbcode}_refined.mtz", 20)
        if os.path.isfile(refined_pdb):
            os.remove(refined_mtz)
            os.remove(outputdir+f"{pdbcode}_refined.mmcif")
            try:
                pdbout = os.path.join(outputdir, f"{pdbcode}_removed_waters.pdb")
                grafter._remove_waters_close_to_TRP(refined_pdb, pdbout)
                os.remove(refined_pdb)
            except:
                pdbout = refined_pdb
                print(f"Failed to remove waters close to TRP in {pdbout}")
            graftedGlycans = grafter._calc_rscc_grafted_glycans(pdbout, mtzpath, graftedGlycans)
            graftedGlycans = grafter._remove_grafted_glycans(pdbout, mtzpath, graftedGlycans, outputdir, rscc_threshold = 0.8)
        for graft in graftedGlycans:
            graft["pdbcode"] = pdbcode
            AllGlycans.append(graft)
        df_temp = pd.DataFrame(graftedGlycans)
        df_temp.to_csv(outputdir+f"/{pdbcode}_graft_summary.csv")
    df_out = pd.DataFrame.from_dict(AllGlycans)
    output_csv = outputdir + "/full_graft_summary.csv"
    df_out.to_csv(output_csv)
    return 




if __name__ == "__main__":
    defaultdatabasedir = "/vault/privateer_database/pdb" # Location the privateer database is (if using database). If not using database, set to None
    defaultpdbmirrordir = "/vault/pdb_mirror/data/structures/all" # Location the unedited pdb/mmcif files are stored (assumes pdb files are in subdir "pdb" and mmcif in subdir "mmcif")
    defaultmtzdir = "/vault/pdb_mtz_files"
    cwd = os.getcwd()
    defaultreceiverdir = os.path.join(cwd,"receivers") # Location you want to save edited pdb files with incorrect c-glycans removed
    if os.getenv("PRIVATEERDATA", None) is not None:
        rootdir = os.getenv("PRIVATEERDATA", None)
        defaultdonordir = os.path.join(rootdir, "glycan_donor_repertoire")
    else:
        defaultdonordir = "/y/people/lah583/privateer/data/glycan_donor_repertoire"
    defaultoutputdir = os.path.join(cwd,"output")
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
    args = parser.parse_args()
    if args.mode == 'fix':
        fix_Cglycans(args.databasedir,args.pdbmirrordir,args.mtzdir,args.receiverdir,args.donordir,args.outputdir)
    elif args.mode == 'find':
        find_and_graft_Cglycans(args.receiverdir,args.mtzdir,args.donordir,args.outputdir)
    else:
        print("Mode of operation not specified. Exiting...")



    



