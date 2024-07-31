import gemmi
import os
import numpy as np
import pandas as pd
import json
import glob,re
import requests
import time
from datetime import timedelta
import subprocess
import gzip,shutil
import sys

# run code on 21-Jul-2024 

def blob_search(cifpath:str,requestedchains:list, keepfile:bool = False) -> tuple:

    avglength = 6.4118
    pdbid = os.path.basename(cifpath).split('.')[0]
    # print(pdbid)
    sfcif = f'/vault/pdb_mirror/data/structures/all/structure_factors/r{pdbid}sf.ent.gz'
    st = gemmi.read_structure(cifpath)
    st.remove_ligands_and_waters() # remove ligands and waters

    pointlist = []
    chainlist = []; residuelist = []; sum_density = []
    consensus = []; status = []; mean_density = []; gridpoints = []

    if not os.path.exists(sfcif):
        return chainlist,residuelist,sum_density,consensus,status,mean_density,gridpoints
    

    for model in st:
        for chain in model:
            if chain.name not in requestedchains: continue

            for residue in chain:
                # GET ESTIMATED CENTROID WITH AVERAGE TRANSLATION LENGTH ~ 6.411 Å
                if (residue.name == 'TRP'):
                    
                    if residue.label_seq != None: pentaseq = get_consensus(inputchain=chain,inputresidue=residue)
                    else: pentaseq = residue.name

                    if re.search('W.{2}W',pentaseq[3:10]) == None: continue

                    ce3,cd1 = None,None
                    for atom in residue:
                        if atom.name == 'CE3': ce3 = atom.pos
                        elif atom.name == 'CD1': cd1 = atom.pos
                    if ce3 != None and cd1 != None:
                        vCED = cd1-ce3; norm = vCED.length(); uvCED = vCED/norm
                        translatedCED = uvCED*avglength
                        newpoint = translatedCED + ce3
                        pointlist.append(newpoint)
                        consensus.append(pentaseq)
                        residuelist.append(residue.seqid.num)
                        chainlist.append(chain.name)


    resolution = st.resolution
    threshold = resolution*(0.3877) - 0.2417
    if not pointlist: 
        # print('Exist here')
        return chainlist,residuelist,sum_density,consensus,status,mean_density,gridpoints

    # recalculate map
    cmtz,stfolder = convert_sf_to_mtz(sfcif=sfcif)
    if None in [cmtz,stfolder]: return [],[],[],[],[],[],[]
    # editedpdb = stfolder + f'/{pdbid}_edited.pdb'
    # st.write_pdb(editedpdb)
    if 'cif' in os.path.basename(cifpath):
       cif_block = gemmi.cif.read(cifpath)[0]
       groups = gemmi.MmcifOutputGroups(True)
       st.update_mmcif_block(block=cif_block,groups=groups)
       editedcif = stfolder + f'/{pdbid}_edited.cif'
       cif_block.write_file(editedcif) 
    rfmtz = stfolder + f'/{pdbid}_refine.mtz'
    otherout = stfolder + f'/{pdbid}'
    pdbout = stfolder + f'/{pdbid}_refine.pdb'
    run_refmac(mtz_in=cmtz,pdb_in=editedcif,mtz_out=rfmtz,pdb_out=pdbout,other_out=otherout)
    
    try:
        mtz = gemmi.read_mtz_file(rfmtz)
    except RuntimeError as e:
        print(e)
        print(f'RuntimeError in reading mtz at {pdbid}')
        return chainlist,residuelist,sum_density,consensus,status,mean_density,gridpoints
    
    grid = mtz.transform_f_phi_to_map('DELFWT', 'PHDELWT', sample_rate=2.0)
    grid.normalize()

    # COPY THE GRID --> SAMPLE GRIDPOINTS IN RADIUS = 3 Å ( BY SETTING VALUES THOSE GRIDPOINTS ON CLONED GRID ~ 10000)
    start = 1000 # arbitary number
    for newpoint in pointlist:
        gr = grid.clone()
        gr.set_points_around(newpoint, radius=3, value=start)
        pointgroup = np.argwhere(gr.array == start) # return positions in symmetry ???

        # TRACE BACK ONTO ORIGINAL GRID --> SUM DENSITY VALUES
        atomlist = []; values = []
        for point in pointgroup:
            atom = grid.get_point(point[0],point[1],point[2])
            atom = grid.get_position(point[0],point[1],point[2])
            atomlist.append(atom)
            value = grid.get_value(point[0],point[1],point[2])
            values.append(value)
            
        s_density = np.sum(values).round(3)
        avg_density = np.mean(values).round(3)
        
        if avg_density > threshold: status.append('yes')
        else: status.append('no')
        sum_density.append(s_density)
        mean_density.append(avg_density)
        gridpoints.append(float(len(pointgroup)))

    # print(chainlist,residuelist,sum_density,consensus,status,mean_density,gridpoints)
    if all(x=='no' for x in status): return [],[],[],[],[],[],[]

    if keepfile == False: # keep only updated .cif, converted .mtz, REFMAC .pdb, REFMAC .mtz
        for file in glob.glob(stfolder + '/*'):
            if file not in [cmtz,rfmtz]: os.remove(file)

    return chainlist,residuelist,sum_density,consensus,status,mean_density,gridpoints

def run_cif2mtz(hklin:str,hklout:str):

    # print("Running CCP4 CIF2MTZ with", hklin)
    _args = []
    _args += ["hklin", hklin]
    _args += ["hklout", hklout]

    _stdin = []
    _stdin.append("END")

    process = subprocess.Popen(
    args=["/jarvis/programs/xtal/ccp4-8.0/bin/cif2mtz"] + _args,
    stdin=subprocess.PIPE if _stdin else None,
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
    encoding="utf8",
    env={**os.environ,},
    cwd=os.getcwd(),
    )
    if _stdin:
        stdin_str = '\n'.join(_stdin)
        process.communicate(input=stdin_str)

def convert_sf_to_mtz(sfcif:str):

    filename = os.path.basename(sfcif).split('.')[0]
    subfolder = os.path.join('/y/people/tpp508/results/structures/from_sf',f'{filename}')
    if not os.path.exists(subfolder): os.makedirs(subfolder)

    jsfcif = os.path.join(subfolder,f'r{filename}sf.ent')
    with gzip.open(sfcif, 'rb') as f_in:
        with open(jsfcif, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    output = subfolder + f'/{filename}.mtz'
    run_cif2mtz(hklin=jsfcif,hklout=output)

    try:
        mtz = gemmi.read_mtz_file(output) 
    except RuntimeError as e:
        try:
            doc = gemmi.cif.read(jsfcif)
            rblocks = gemmi.as_refln_blocks(doc).__getitem__(0)
            cif2mtz = gemmi.CifToMtz()
            mtz = cif2mtz.convert_block_to_mtz(rblocks)
            freeR = mtz.column_with_label('FreeR_flag')
            freeR.label = 'FREE'
            mtz.write_to_file(output)
        except RuntimeError as e:
            print('Cannot convert cif to mtz with cif2mtz and Gemmi')
            print(e)
            return None,None
    labels = mtz.column_labels()
    # print(labels)

    if (('FP' in labels) or ('F' in labels)):
        freeR_infile = output
    else:
        args = []
        ctr_infile = output; outfile1 = subfolder + f'/{filename}_fs.mtz'
        _args = ["-hklin",ctr_infile,"-hklout",outfile1]

        
        if 'I' in mtz.column_labels() and 'SIGI' in mtz.column_labels():
            args += ["-colin", "/*/*/[I,SIGI]"]
        elif 'F(+)' in mtz.column_labels() and 'SIGF(+)' in mtz.column_labels():
            args += ["-colano", "/*/*/[F(+),SIGF(+),F(-),SIGF(-)]", "-amplitudes"]
        elif 'IMEAN' in mtz.column_labels() and 'SIGIMEAN' in mtz.column_labels():
            args += ["-colin", "/*/*/[IMEAN,SIGIMEAN]"]
        elif 'I(+)' in mtz.column_labels() and 'SIGI(+)' in mtz.column_labels():
            args += ["-colano", "/*/*/[I(+),SIGI(+),I(-),SIGI(-)]"]
        else:
            print(f'No intensity or amplitudes found at {sfcif}')
            print(mtz.column_labels())
            return None,None
        
        if 'FREE' in labels: args += ["-freein", "/*/*/[FREE]"] # keep FREE columns
        
        _args += args
        # print('Running CCP4 CTRUNCATE') # in case converted .mtz only has I, SIGI
        try:
            command=["/jarvis/programs/xtal/ccp4-8.0/bin/ctruncate"] + _args 
            subprocess.run(command,capture_output=True)
            freeR_infile = outfile1
            mtz = gemmi.read_mtz_file(outfile1)
        except (RuntimeError,subprocess.CalledProcessError) as e:
            print(f'{e} at CTRUNCATE {outfile1}')
            print(mtz.column_labels())
            return None,None

    # labels = mtz.column_labels()
    # for idx in range(len(labels))[::-1]:
    #     if labels[idx] in ['FC','PHIC']: mtz.remove_column(idx)
    # mtz.write_to_file(output) # convert structure factor from cif to mtz

    if ('FreeR_flag' not in mtz.column_labels()): # if new (converted/ctruncate) mtz doesn't have 'FREERFLAG', run FREERFLAG
        outfile2 = subfolder + f'/{filename}_fs_free.mtz'
        if not os.path.exists(freeR_infile): return None,None
        run_freeR_flag(mtzin=freeR_infile,mtzout=outfile2)
        finalmtz = outfile2
    else: finalmtz = freeR_infile

    return finalmtz,subfolder

def run_freeR_flag(mtzin: str, mtzout: str):

    # print("Running CCP4 FREERFLAG with", mtzin)
    _args = []
    _args += ["hklin", mtzin]
    _args += ["hklout", mtzout]

    _stdin = []
    _stdin.append("END")

    process = subprocess.Popen(
    args=["/jarvis/programs/xtal/ccp4-8.0/bin/freerflag"] + _args,
    stdin=subprocess.PIPE if _stdin else None,
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
    encoding="utf8",
    env={**os.environ,},
    cwd=os.getcwd(),
    )
    if _stdin:
        stdin_str = '\n'.join(_stdin)
        process.communicate(input=stdin_str)


# Taken and edited from ModelCraft # written by Jordan Dialpuri 2024
def run_refmac(mtz_in: str, pdb_in: str, mtz_out: str, pdb_out: str, other_out: str): 
    
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
    _stdin.append("NCYCLES 0")
    _stdin.append("WEIGHT AUTO")
    _stdin.append("PHOUT")
    _stdin.append("END")

    process = subprocess.Popen(
    args=["/jarvis/programs/xtal/ccp4-8.0/bin/refmacat"] + _args,
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


def get_method_and_resolution(path:str) -> tuple[float,str]:

    st = gemmi.read_structure(path)

    resolution = round(st.resolution,2)
    for key, value in st.info.items():
        if key == '_exptl.method': method = value
    
    return resolution,method


def get_consensus(inputchain:gemmi.Chain, inputresidue:gemmi.Residue) -> str:
    
    concat = ''
    rng = np.array([-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6])

    inputresidueNum = inputresidue.seqid.num
    residuerng = rng + inputresidueNum
    
    for seqid in residuerng:
        try:
            neighbour = inputchain[f'{seqid}'].__getitem__(0)
        except IndexError: concat += '-'; continue
        
        if ((neighbour != None)
            and (neighbour.label_seq != None) 
            and (gemmi.find_tabulated_residue(neighbour.name).is_amino_acid())
            and ((neighbour.label_seq - inputresidue.label_seq) in rng)
            and ((neighbour.seqid.num - inputresidue.seqid.num) in rng)):
            concat += gemmi.find_tabulated_residue(neighbour.name).one_letter_code
        else: concat += '-'
    
    return concat


def check_consensus_sequence(sequence:str) -> bool:

    output = False

    for index in range(len(sequence)):
        segment = sequence[index:index+4]
        if re.match('W.{2}[W|C]',segment) != None:
            output = True
    
    return output


def check_expression_system_with_cif(path:str) -> tuple[list,bool]:

    sourcefile = '/y/people/tpp508/results/C_sites/taxon_summary.json' # include metazoan and toxoplasma taxonomy ids
    with open(sourcefile,'r') as f: 
        taxsrc = json.load(f)
        taxonids = taxsrc['taxonids']
        sciNames = taxsrc['sciNames']

    chainlist,output = [], False

    doc = gemmi.cif.read(path)
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
        except (ValueError,TypeError,SyntaxError,RuntimeError): 
            # print(f'{pdb} with taxhost = {taxhost} and hostname = {hostname}')
            pass
            
        if (taxhost in taxonids) or (hostname.lower() in sciNames):
            # print(hent)
            try: proteinsequence = polymers[int(hent)-1]
            except (IndexError,ValueError,TypeError,SyntaxError) as e: 
                # print(f'{pdb} with error = {e}')
                continue
            
            output = check_consensus_sequence(sequence=proteinsequence)
            chainname = chainforentity[int(hent)-1]
            if output == False: continue
            if len(chainname) > 1: chainname = chainname.split(',')
            if chainname not in chainlist: chainlist += chainname
            

    for (natent,taxnat,natname) in zip(natentities,taxnats,natnames): # PDB 2r2b host Mus musculus without taxid!

        try: taxnat = taxnat.replace(' ', ''); taxnat = int(taxnat)
        except (ValueError,TypeError,SyntaxError,RuntimeError): 
            # print(f'{pdb} with taxnat = {taxnat} and natname = {natname}')
            pass
            
        if (taxnat in taxonids) or (natname.lower() in sciNames):

            try: proteinsequence = polymers[int(natent)-1]
            except (IndexError,ValueError,TypeError,SyntaxError) as e: 
                # print(f'{pdb} with error = {e}')
                continue

            output = check_consensus_sequence(sequence=proteinsequence)
            chainname = chainforentity[int(natent)-1]
            if output == False: continue
            if len(chainname) > 1: chainname = chainname.split(',')
            if chainname not in chainlist: chainlist += chainname

    return chainlist#,output

start = time.time()

folder = '/vault/pdb_mirror/data/structures/all/mmCIF'

print('search_grid_pdb_new_ver.py')
for path in glob.glob(folder + '/*'):
    pdb = os.path.basename(path).split('.')[0]
    
    output_file = f"/y/people/tpp508/results/C_sites/pdb_json_files/{pdb}.json"
    if os.path.exists(output_file): continue
    
    # print('Found file', pdb)
    # if pdb != '6rv6': continue
    resolution,method = get_method_and_resolution(path=path)
    if method != 'X-RAY DIFFRACTION': continue 

    requestedchains = check_expression_system_with_cif(path=path)
    if not requestedchains: continue

    # print(method,requestedchains)
    authchains,authresidues,density_sum,sequons,status,density_mean,gridnumpoints = blob_search(cifpath=path,requestedchains=requestedchains)
    if not all([authchains,authresidues,density_sum,sequons,status,density_mean,gridnumpoints]): 
        folder = f'/y/people/tpp508/results/structures/from_sf/r{pdb}sf'
        if os.path.exists(folder):
            shutil.rmtree(f'/y/people/tpp508/results/structures/from_sf/r{pdb}sf')
        continue
    
    resolutions = [resolution]*len(authresidues); techniques = [method]*len(authresidues)

    newdict = {'pdbid':pdb,'authchains':authchains,'authresidues':authresidues,'sumdensities':density_sum,
               'consensus':sequons,'resolution':resolutions,'methods':techniques, 'status':status,'mean_density':density_mean,'numpoint':gridnumpoints}
    with open(output_file, "w") as json_file: 
        json.dump(newdict, json_file)

collated_data = {}
for path in os.scandir("/y/people/tpp508/results/C_sites/pdb_json_files/"):
    with open(path.path, "r") as json_file: 
        data = json.load(json_file)
    for k, v in data.items():
        collated_data.setdefault(k, [])
        collated_data[k].append(v)

with open('/y/people/tpp508/results/C_sites/csvfile/blob_search_whole_pdb_new_ver.json','w') as f:
    json.dump(collated_data,f) 

end = time.time()
elapsed = end-start
print(str(timedelta(seconds=elapsed)))
print('END')