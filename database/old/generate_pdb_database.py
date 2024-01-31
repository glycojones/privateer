from privateer import privateer_core as pvt
from datetime import datetime
import json
import gzip
import shutil
import os
import urllib.request
from typing import Optional
from multiprocessing import Pool
import multiprocessing
from tqdm import tqdm 


def find_pdb_path(pdb_code: str) -> str:

    gzipped_pdb = f"/vault/pdb_mirror/data/structures/all/pdb/pdb{pdb_code}.ent.gz"
    unzipped_pdb = f"/vault/tmp_extracted_pdbs/pdb{pdb_code}.ent"

    if os.path.exists(unzipped_pdb):
        return unzipped_pdb

    with gzip.open(gzipped_pdb, 'rb') as f_in:
        with open(unzipped_pdb, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return unzipped_pdb

def delete_tmp_pdb(pdb_path: str) -> None:
    os.remove(pdb_path)

def find_mtz_path(pdb_code: str) -> Optional[str]:
    output_path = f"/vault/pdb_mtz_files/{pdb_code.lower()}.mtz"
    url = f"https://edmaps.rcsb.org/coefficients/{pdb_code.lower()}.mtz"
    
    try:
        urllib.request.urlretrieve(url, output_path)
    except Exception as e :
        return None
    return output_path

def worker(pdb_code: str) -> None: 

    pdb_path = find_pdb_path(pdb_code)
    mtz_path = find_mtz_path(pdb_code)

    if not pdb_path: 
        return 
    
    try:
        if mtz_path:
            glycosylation = pvt.GlycosylationComposition(pdb_path, mtz_path, "FP,SIGFP")
        else:
            glycosylation = pvt.GlycosylationComposition(pdb_path)
    except:
        return

    middlefix = pdb_code.lower()[1:3]
    base_dir = f"/vault/privateer_database/pdb/{middlefix}"

    output = f"{base_dir}/{pdb_code.lower()}.json"
    
    # if os.path.exists(output):
    #     return (True, "Already in DB")

    data = {}

    metadata = {
        "dateAdded": datetime.today().strftime('%Y-%m-%d'), 
        "experimentalDataAvailable": "yes" if mtz_path else "no"
    }

    glycans = {
        "n-glycan": [],
        "o-glycan": [],
        "s-glycan": [],
        "c-glycan": [],
        "ligand": []
    }

    number_of_glycans = glycosylation.get_number_of_glycan_chains_detected()
    # print("No glycans", number_of_glycans)
    for glycanNo in range(number_of_glycans):
        glycan = glycosylation.get_glycan(glycanNo)

        numsugars = glycan.get_total_number_of_sugars()
        glycosylation_type = glycan.get_glycosylation_type()
        wurcs = glycan.get_wurcs_notation()
        root_info = glycan.get_root_info()

        snfg_data = glycan.get_SNFG_strings()
        snfg = snfg_data["SNFG"]
        sugars= []

        for j in range(numsugars):
            sugar = glycan.get_monosaccharide(j)

            summary = sugar.get_sugar_summary()

            sugar_id = f'{summary["sugar_name_short"]}-{summary["sugar_pdb_chain"]}-{summary["sugar_seqnum"]}'
            Q = summary["Q"]
            Phi = summary["Phi"]
            Theta = summary["Theta"]
            RSCC = summary["RSCC"]
            detected_type = sugar.get_denomination()
            conformation = sugar.get_conformation_name()
            mFo = summary["mFo"]
            BFactor = sugar.get_bfactor()
            diagnostic = sugar.get_privateer_diagnostic()

            sugar_data = {
                "sugarId": sugar_id, 
                "q": Q,
                "phi": Phi, 
                "theta": Theta, 
                "rscc": RSCC,
                "detectedType": detected_type,
                "conformation": conformation, 
                "bFactor": BFactor,
                "mFo": mFo,
                "diagnostic": diagnostic
            }
            sugars.append(sugar_data)
         
        root_info_format = {
            "proteinResidueType": root_info["ProteinResidueType"],
            "proteinResidueId": root_info["ProteinResidueID"],
            "proteinResidueSeqnum": root_info["ProteinResidueSeqnum"],
            "proteinChainId": root_info["ProteinChainID"],
            "rootSugarChainId": root_info["RootSugarChainID"]
        }

        glycan_data = {
            **root_info_format, 
            "numberOfSugars": numsugars,
            "wurcs": wurcs, 
            "snfg": snfg, 
            "sugars": sugars,
        }

        if glycosylation_type in glycans:
            glycans[glycosylation_type].append(glycan_data)

    data["metadata"] = metadata
    data["glycans"] = glycans

    if not any(
        [glycans["n-glycan"], glycans["ligand"], glycans["o-glycan"], glycans["c-glycan"], glycans["s-glycan"]]
    ):
        return

    if not os.path.isdir(base_dir):
        os.makedirs(base_dir, exist_ok=True)


    with open(output, 'w') as fout:
        # print("writing ", output)
        fout.write(json.dumps(data))

    delete_tmp_pdb(pdb_path)

    del data

    return


def captured_worker(pdb_code):
    try: 
        result = worker(pdb_code)
        # if not result[0]:
        #     with lock:
        #         if result[1] == "No Glycans":
        #             with open("no_glycan_db.txt", "a") as failed_file: 
        #                 failed_file.write(f"{pdb_code}\n")
        #         elif result[1] == "Can't find PDB": 
        #             with open("no_pdb.txt", "a") as failed_file: 
        #                 failed_file.write(f"{pdb_code}\n")
            
    except Exception as e:
        # print("Exception", e)
        # with lock:
        #     with open("failed_db.txt", "a") as failed_file: 
        #         failed_file.write(f"{pdb_code},{e}\n")
        return


def get_pdb_list(pdb_vault_dir: str): 
    return [x[3:-7] for x in os.listdir(pdb_vault_dir)]


def main(): 
    pdb_vault_dir = "/vault/pdb_mirror/data/structures/all/pdb/"
    # pdb_list = get_pdb_list(pdb_vault_dir)[13:14]
    pdb_list = get_pdb_list(pdb_vault_dir)

    # print(pdb_list)
    # captured_worker("5fjj")
    # quit()
    # lock = Lock()

    # with open("failed_db.txt", "w") as f_:
    #     f_.write("PDB,Reason\n")

    # with open("no_glycan_db.txt", "w") as f_:
    #     f_.write("PDB\n")

    # with open("no_pdb.txt", "w") as f_:
    #     f_.write("PDB\n")

    with Pool(100, maxtasksperchild=1) as pool_:
        x = list(tqdm(pool_.map(worker, pdb_list), total=len(pdb_list)))
    # try:
    #     pool = multiprocessing.pool.ThreadPool(4)
    #     x = pool.map(worker, pdb_list)
    # finally: # To make sure processes are closed in the end, even if errors happen
    #     pool.close()
    #     pool.join()

if __name__ == "__main__":
    main()

    # pdb_vault_dir = "/vault/pdb_mirror/data/structures/all/pdb/"
    # pdb_list = get_pdb_list(pdb_vault_dir)[13:14]
    # pdb_list = get_pdb_list(pdb_vault_dir)[13:50]

    # print(pdb_list)
    # # captured_worker("5fjj")
    # # quit()
    # # lock = Lock()

    # # with open("failed_db.txt", "w") as f_:
    # #     f_.write("PDB,Reason\n")

    # # with open("no_glycan_db.txt", "w") as f_:
    # #     f_.write("PDB\n")

    # # with open("no_pdb.txt", "w") as f_:
    # #     f_.write("PDB\n")

    # # with Pool(1) as pool_:
    # #     x = list(tqdm(pool_.map(worker, pdb_list), total=len(pdb_list)))
    # try:
    #     pool = Pool(4) # on 8 processors
    #     x = pool.map(worker, pdb_list)
    # finally: # To make sure processes are closed in the end, even if errors happen
    #     pool.close()
    #     pool.join()