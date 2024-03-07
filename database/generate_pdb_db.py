
import gzip
import shutil
import os
import urllib.request
from typing import Optional
from multiprocessing import Pool
import multiprocessing
from tqdm import tqdm 
import subprocess
import requests
import resource

def find_pdb_path(pdb_code: str) -> str:

    gzipped_pdb = f"/vault/pdb_mirror/data/structures/all/mmCIF/{pdb_code}.cif.gz"
    unzipped_pdb = f"/vault/tmp_extracted_mmcif/{pdb_code}.mmcif"

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
    
    if os.path.exists(output_path):
        return output_path
    
    try:
        urllib.request.urlretrieve(url, output_path)
    except Exception as e :
        return ""
    return output_path
    
def find_information(pdb_code: str) -> None: 
    url = f" https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_code}"
    response = requests.get(url)
    json = response.json()
    if pdb_code not in json: 
        return None
    
    for entry in json[pdb_code]:
        if "Electron Microscopy" in entry["experimental_method"]:
            related = entry["related_structures"]
            for rel in related:
                if rel["resource"] == "EMDB":
                    emdb_accession = rel["accession"]
                    
                    res_response = requests.get(f"https://www.ebi.ac.uk/emdb/api/entry/{emdb_accession[4:]}")
                    res_json = res_response.json()
                    for method in res_json["structure_determination_list"]["structure_determination"]:
                        for im in method["image_processing"]:
                            res = im["final_reconstruction"]["resolution"]["valueOf_"]

                            return ("EM", emdb_accession, res)
                        
def find_map_path(em_code: str) -> Optional[str]:
    output_path = f"/vault/pdb_map_files/{em_code.lower()}.map.gz"
    output_path_ungz = f"/vault/pdb_map_files/{em_code.lower()}.map"
    if os.path.exists(output_path_ungz):
        return output_path_ungz

    url = f"https://files.rcsb.org/pub/emdb/structures/{em_code}/map/{em_code.lower().replace('-', '_')}.map.gz"
    try:
        urllib.request.urlretrieve(url, output_path)
        with gzip.open(output_path, 'rb') as f_in:
            with open(output_path_ungz, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    except Exception as e :
        print(e)
        print(url)
        return ""
    return output_path_ungz

def get_pdb_list(pdb_vault_dir: str): 
    return [x[:4] for x in os.listdir(pdb_vault_dir)]


def limit_virtual_memory():
    # The tuple below is of the form (soft limit, hard limit). Limit only
    # the soft part so that the limit can be increased later (setting also
    # the hard limit would prevent that).
    # When the limit cannot be changed, setrlimit() raises ValueError.
    MAX_VIRTUAL_MEMORY = 10 * 1024 * 1024 * 1024 # 1Gb

    resource.setrlimit(resource.RLIMIT_AS, (MAX_VIRTUAL_MEMORY, resource.RLIM_INFINITY))

def new_worker(pdb_code):
    pdbpath = find_pdb_path(pdb_code)
    mtzpath = find_mtz_path(pdb_code)
    mappath = None
    mapres = None
    if mtzpath == "":
        info = find_information(pdb_code)
        print(info)
        if info:
            if info[0] == "EM":
                mappath = find_map_path(info[1])
                mapres = float(info[2])

    middlefix = pdb_code.lower()[1:3]
    base_dir = f"/vault/privateer_database/pdb/{middlefix}"

    output = f"{base_dir}/{pdb_code.lower()}.json"
   
    
    try: 
        if mappath and mapres:
            subprocess.run(["python", "database/calculate.py", "-pdbpath", pdbpath, "-mappath", mappath, "-mapres", mapres, "-basedir", base_dir, "-output", output]
            ,timeout=180,
            stdout=subprocess.DEVNULL,
            preexec_fn=limit_virtual_memory
        )
        else:
            subprocess.run(["python", "database/calculate.py", "-pdbpath", pdbpath, "-mtzpath", mtzpath, "-basedir", base_dir, "-output", output]
                ,timeout=180,
                stdout=subprocess.DEVNULL,
                preexec_fn=limit_virtual_memory
            )
    except subprocess.TimeoutExpired as e: 
        with open(f"database/timeouts/{pdb_code}.txt", "w") as x_:
            x_.write("")
        return

def worker(pdb_code):
    pdbpath = find_pdb_path(pdb_code)
    mtzpath = find_mtz_path(pdb_code)

    middlefix = pdb_code.lower()[1:3]
    base_dir = f"/vault/privateer_database/pdb/{middlefix}"

    output = f"{base_dir}/{pdb_code.lower()}.json"
    try: 
        cmd = ["python", "/y/people/jsd523/dev/privateer/database/calculate.py", "-pdbpath", pdbpath, "-mtzpath", mtzpath, "-basedir", base_dir, "-output", output]
        subprocess.run(cmd
            ,timeout=180,
            stdout=subprocess.DEVNULL,
            preexec_fn=limit_virtual_memory
        )
    except subprocess.TimeoutExpired as e: 
        with open(f"database/timeouts/{pdb_code}.txt", "w") as x_:
            x_.write("")
        return

def main(): 
    pdb_vault_dir = "/vault/pdb_mirror/data/structures/all/mmCIF/"
    pdb_list = get_pdb_list(pdb_vault_dir)

    # worker("6plh")

    with Pool(128) as pool_:
        x = list(tqdm(pool_.imap_unordered(worker, pdb_list), total=len(pdb_list)))
 

if __name__ == "__main__":
    main()
    # find_information("8wtd")