
import gzip
import shutil
import os
import urllib.request
from typing import Optional
from multiprocessing import Pool
import multiprocessing
from tqdm import tqdm 
import subprocess
import resource

def find_pdb_path(pdb_code: str) -> str:

    # gzipped_pdb = f"/vault/pdb_mirror/data/structures/all/pdb/pdb{pdb_code}.ent.gz"
    unzipped_pdb = f"/vault/pdb-redo/{pdb_code[1:3]}/{pdb_code}/{pdb_code}_final.pdb"

    if os.path.exists(unzipped_pdb):
        return unzipped_pdb
    else:
        return ""

    # with gzip.open(gzipped_pdb, 'rb') as f_in:
    #     with open(unzipped_pdb, 'wb') as f_out:
    #         shutil.copyfileobj(f_in, f_out)

    # return unzipped_pdb

def delete_tmp_pdb(pdb_path: str) -> None:
    os.remove(pdb_path)

def find_mtz_path(pdb_code: str) -> Optional[str]:
    output_path = f"/vault/pdb_mtz_files/{pdb_code.lower()}.mtz"
    url = f"https://edmaps.rcsb.org/coefficients/{pdb_code.lower()}.mtz"
    
    try:
        urllib.request.urlretrieve(url, output_path)
    except Exception as e :
        return ""
    return output_path

def get_pdb_list(pdb_vault_dir: str): 

    pdbs = []
    for x in os.listdir(pdb_vault_dir):
        if len(x) == 2: 
            for pdb in os.listdir(os.path.join(pdb_vault_dir, x)):
                pdbs.append(pdb)
    return pdbs


def limit_virtual_memory():
    # The tuple below is of the form (soft limit, hard limit). Limit only
    # the soft part so that the limit can be increased later (setting also
    # the hard limit would prevent that).
    # When the limit cannot be changed, setrlimit() raises ValueError.
    MAX_VIRTUAL_MEMORY = 10 * 1024 * 1024 * 1024 # 1Gb

    resource.setrlimit(resource.RLIMIT_AS, (MAX_VIRTUAL_MEMORY, resource.RLIM_INFINITY))

def worker(pdb_code):
    pdbpath = find_pdb_path(pdb_code)
    mtzpath = find_mtz_path(pdb_code)


    middlefix = pdb_code.lower()[1:3]
    base_dir = f"/vault/privateer_database/pdbredo/{middlefix}"

    output = f"{base_dir}/{pdb_code.lower()}.json"
    try: 
        subprocess.run(["python", "database/calculate.py", "-pdbpath", pdbpath, "-mtzpath", mtzpath, "-basedir", base_dir, "-output", output]
            ,timeout=180,
            stdout=subprocess.DEVNULL,
            preexec_fn=limit_virtual_memory
        )
    except subprocess.TimeoutExpired as e: 
        with open(f"database/timeouts/{pdb_code}.txt", "w") as x_:
            x_.write("")
        return

def main(): 
    pdb_vault_dir = "/vault/pdb-redo"
    pdb_list = get_pdb_list(pdb_vault_dir)


    with Pool(100) as pool_:
        x = list(tqdm(pool_.imap_unordered(worker, pdb_list), total=len(pdb_list)))
 

if __name__ == "__main__":
    main()
