import gemmi
import os
import numpy as np
import pandas as pd
import json
import re
import requests
import time
from datetime import timedelta
from pathlib import Path
from privateer import privateer_core as pvt
import sys
import subprocess
import gzip
import shutil
from multiprocessing.pool import ThreadPool
from multiprocessing import cpu_count
import argparse
import sys
from pathlib import Path
# sys.path.append('/Users/thaophuongpham/Dev/privateer')

def download_emmap(args:tuple) -> tuple[bool,str,float]:
    
    t0 = time.time()
    EMMAP = args[0]; savepath = args[1]
    # basename1 = os.path.basename(savepath)
    basename2 = os.path.basename(EMMAP)

    emmap = requests.get(EMMAP, allow_redirects=True)
    if emmap.ok:
        if not os.path.exists(savepath):
            with open (savepath,'wb') as s: s.write(emmap.content)
        return True,basename2,time.time() - t0
    else:
        return False,basename2,time.time() - t0


# FOLLOWED THE INSTRUCTION FROM KONRAD HAFEN (2023) 
# WEBSITE: 'https://towardsdatascience.com/use-python-to-download-multiple-files-or-urls-in-parallel-1759da9d6535' 
def download_em_halfmaps_parallel(pdbid:str,emdbid:str) -> tuple[str,str,str,str]:

    """DOWNLOAD HALF MAP 1, HALF MAP 2, MASKED MAP, AND PRIMMARY MAP FROM EMDB"""

    emdbid1 = (emdbid.replace('-','_')).lower()

    halfmap1 = f'https://ftp.ebi.ac.uk/pub/databases/emdb/structures/{emdbid}/other/{emdbid1}_half_map_1.map.gz'
    halfmap2 = f'https://ftp.ebi.ac.uk/pub/databases/emdb/structures/{emdbid}/other/{emdbid1}_half_map_2.map.gz'
    maskedmap1 = f'https://ftp.ebi.ac.uk/pub/databases/emdb/structures/{emdbid}/masks/{emdbid1}_msk_1.map'
    primmarymap = f'https://ftp.ebi.ac.uk/pub/databases/emdb/structures/{emdbid}/map/{emdbid1}.map.gz'

    urls = []

    outfolder = os.getcwd()
    output1 = os.path.join(outfolder,f'{pdbid}_half_map_1.map.gz')
    output2 = os.path.join(outfolder,f'{pdbid}_half_map_2.map.gz')
    output3 = os.path.join(outfolder,f'{pdbid}_msk_1.map')
    output4 = os.path.join(outfolder,f'{pdbid}_primmap.map.gz')

    fns = []

    extract_dir1 = os.path.join(outfolder,f'{pdbid}_half_map_1.map')
    extract_dir2 = os.path.join(outfolder,f'{pdbid}_half_map_2.map')
    extract_dir4 = os.path.join(outfolder,f'{pdbid}_primmap.map')

    if not os.path.exists(extract_dir1): fns.append(output1); urls.append(halfmap1)
    if not os.path.exists(extract_dir2): fns.append(output2); urls.append(halfmap2)
    if not os.path.exists(output3): fns.append(output3); urls.append(maskedmap1)
    if not os.path.exists(extract_dir4): fns.append(output4); urls.append(primmarymap)

    inputs = zip(urls,fns)

    # DOWNLOAD IN PARALLEL
    start = time.time()
    print(f'Start downloading EM maps at: {start}')
    cpus = cpu_count() 
    results = ThreadPool(cpus - 1).imap_unordered(download_emmap, inputs) 
    for result in results: 
        print('\ndownloaded_status', result[0], 'file', result[1] , 'time (s):', result[2])
    print(f'Map downloading finsishes in total {time.time()-start}')

    # UNPACK EM.MAP.ZIP
    print('Start unzipping the maps')
    if os.path.exists(output1):
        extract_dir1 = os.path.join(outfolder,f'{pdbid}_half_map_1.map')
        with gzip.open(output1, 'rb') as f_in:
            with open(extract_dir1, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(output1)
    elif not os.path.exists(extract_dir1): extract_dir1 = None
    
    if os.path.exists(output2):
        # extract_dir2 = os.path.join(outfolder,f'{pdbid}_half_map_2.map')
        with gzip.open(output2, 'rb') as f_in:
            with open(extract_dir2, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(output2)
    elif not os.path.exists(extract_dir2): extract_dir2 = None
    
    if not os.path.exists(output3): output3 = None

    if os.path.exists(output4):
        extract_dir4 = os.path.join(outfolder,f'{pdbid}_primmap.map')
        with gzip.open(output4, 'rb') as f_in:
            with open(extract_dir4, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(output4)
    elif not os.path.exists(extract_dir4): extract_dir4 = None

    print('Finish unzipping the maps')

    return extract_dir1,extract_dir2,output3,extract_dir4

def cube_surround_checkpoint(checkpoint) -> list:

    a,b,c = [],[],[]
    X,Y,Z = [],[],[]; newpositions = []; arr = []
    for i in range(-3,4):
        x = checkpoint.x; y = checkpoint.y; z= checkpoint.z
        a.append(x+i); b.append(y+i); c.append(z+i)

    for suba in a:
        for subb in b:
            for subc in c:
                newpos = gemmi.Position(suba,subb,subc)
                newpositions.append(newpos); arr.append([suba,subb,subc])
                X.append(suba); Y.append(subb); Z.append(subc)
    
    return newpositions

def blob_search_cryo(cifpath:str,emmap:str) -> list:

    """Find potential C-mannosylation sites"""
    
    avglength = 6.4118

    residuelist = []; chainlist = []
    consensus = []; status = []
    eclashes = []; pclashes = []
    resiaddress = []
    finallist = []

    st = gemmi.read_structure(cifpath)
    ns = gemmi.NeighborSearch(st[0], st.cell, 5).populate(include_h=False)

    pointlist = []; trplist = []

    chainidx = 0
    for model in st:
        for chain in model:
            residx = 0
            for residue in chain:

                # GET ESTIMATED CENTROID WITH AVERAGE TRANSLATION LENGTH ~ 6.411 Å
                if (residue.name == 'TRP'):
                    # print(chain,chainidx,residue,residx)
                    
                    pentaseq = get_consensus(inputchain=chain,inputresidue=residue)
                    
                    if re.search('W.{2}[W|C]',pentaseq[3:10]) != None:

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
                            trplist.append(residue)
                            # residx = residue.label_seq -1
                            address = {'chainidx':chainidx,'residx':residx}
                            resiaddress.append(address)
                residx += 1
            chainidx += 1
    
    if pointlist:

        # READ RECALCULATED MAP      
        threshold = 308.656

        map = gemmi.read_ccp4_map(emmap,setup=True)
        grid = map.grid
        grid.normalize()

        for (newpoint,trp,trpaddress) in zip(pointlist,trplist,resiaddress):

            print(f'HERE: {newpoint,trp}')
            estimatedvalues = np.nan
            if newpoint == None: 
                status.append(np.nan)
                eclashes.append(np.nan)
                pclashes.append(np.nan)
                continue

            values = []; clashes = 0; removedpos = []

            # MAKE A CUBE OF GRIDPOINTS, SURROUNDING CENTRAL POINT
            poslist = cube_surround_checkpoint(checkpoint=newpoint)
            for position in poslist:

                # cannot classify which densities from amino acids --> check if cube's points span proteinresidues
                marks = ns.find_atoms(position, '\0', radius=0.7)
                countclashes = 0
                for mark in marks:
                    # print(mark)
                    if mark.image_idx != 0: continue
                    cra = mark.to_cra(st[0])
                    # print(cra)
                    if (cra.residue.name in ['MAN','BMA']): continue
                    else:                
                        countclashes += 1 
                if countclashes > 0: clashes += 1 # clashes : scaling 1 position ~ 1 clash

                value = grid.interpolate_value(position)
                values.append(value)
            
            for pos in removedpos: poslist.remove(pos)  
            estimatedvalues = clashes/len(poslist)

            sumdense = np.sum(values).round(3)
            print(f'Sum of density = {sumdense} and Clash percent = {estimatedvalues}')
            if ((sumdense >= threshold) and (estimatedvalues <= 0.05)):
                print('Pass') 
                finallist.append(trpaddress)
            else: print('Not pass')

    return finallist

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

def check_expression_system_with_cif(cifpath:str) -> tuple[list,bool]:
    
    """Check expression system of the structures 
        (Only if the input structures are the experimental structures)"""
    
    PRIVATEERDIR = os.getenv("PRIVATEER_USEFUL_DATA", None)
    if PRIVATEERDIR is not None:
        sourcefile = os.path.join(PRIVATEERDIR,'taxon_summary.json')# include metazoan and toxoplasma taxonomy ids
    else:
        raise 'Cannot get Privateer environment --> check Path'
    
    with open(sourcefile,'r') as f: 
        taxsrc = json.load(f)
        taxonids = taxsrc['taxonids']
        sciNames = taxsrc['sciNames']

    pdb = os.path.basename(cifpath).split('.')[0]

    chainlist,output = [], False

    try:
        doc = gemmi.cif.read(cifpath)
        block = doc.sole_block()
    except (FileExistsError,RuntimeError) as e: 
        print (e)
        return chainlist

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
        
        if hostname in ['Drosophila melanogaster'] or taxhost in [7227]: continue
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
            print(f'{pdb} with taxnat = {taxnat} and natname = {natname}')
            pass
        
        if natname in ['Drosophila melanogaster'] or taxnat in [7227]: continue
        if (taxnat in taxonids) or (natname.lower() in sciNames):

            try: proteinsequence = polymers[int(natent)-1]
            except (IndexError,ValueError,TypeError,SyntaxError) as e: 
                print(f'{pdb} with error = {e}')
                continue

            output = check_consensus_sequence(sequence=proteinsequence)
            chainname = chainforentity[int(natent)-1]
            if output == False: continue
            if len(chainname) > 1: chainname = chainname.split(',')
            if chainname not in chainlist: chainlist += chainname

    return chainlist#,output

def get_method_and_resolution(cifpath:str) -> tuple[float,str]:

    """Check method and resolution of the input structure"""
    
    doc = gemmi.cif.read(cifpath)
    block = doc.sole_block()
    method,resolution = None,None
    method = block.find_value('_exptl.method')
    resolution = block.find_value('_em_3d_reconstruction.resolution')
    if resolution == None:
       resolution = block.find_value('_refine.ls_d_res_high')
    
    try: resolution = float(resolution)
    except (RuntimeError,SyntaxError,ValueError,TypeError):
        resolution = None
    
    return method,resolution

def get_emdbid (cifpath:str) -> tuple[str,str]:

    """Get EMDB id of Cryo-EM structures"""

    doc = gemmi.cif.read(cifpath)
    block = doc.sole_block()
    pdbid = block.name.lower()

    emid = None
    databases = block.find_values('_database_2.database_id')
    idcodes =  block.find_values('_database_2.database_code')
    for (database,idcode) in zip(databases,idcodes):
        if database == 'EMDB': emid = idcode
        
    if emid == None:
        dbnames = block.find_values('_pdbx_database_related.db_name')
        tags = block.find_values('_pdbx_database_related.content_type')
        dbids = block.find_values('_pdbx_database_related.db_id')

        for (dbname,tag,dbid) in zip(dbnames,tags,dbids):
            if ((dbname == 'EMDB') and (tag == "'associated EM volume'")):
                emid = dbid
    
    return pdbid,emid

def make_jsonfile_for_grafting(pdbfile:str, graftedmodel:str, residuelist:list):
    
    PRIVATEERDIR = os.getenv("PRIVATEER_USEFUL_DATA", None)
    if PRIVATEERDIR is not None:
        donorsugar = PRIVATEERDIR + '/glycan_donor_repertoire/Alpha-D-Mannose.pdb'
    else: 
        raise 'Cannot get environment of Privateer --> check Path'
    
    glycosylations = []
    for subdict in residuelist:
        glycosylation = {
            "donor_path": donorsugar,
            "glycan_index": 0,
            "receiving_chain_index": subdict['chainidx'],
            "receiving_aa_index": subdict['residx']
        }
        glycosylations.append(glycosylation)
    
    exportdict = {
        "receiver_path": pdbfile,
        "output_path": graftedmodel,
        "glycosylations": glycosylations
    }

    manualgrafting = 'manual_grafting.json'
    with open(manualgrafting,'w') as f: json.dump(exportdict,f)

    return manualgrafting

def try_grafting_interface(pdbfile:str, outputgraft:str, residuelist:list):

    """Running python grafter.py"""

    PRIVATEERSRC = os.getenv("PRIVATEERSRC", None)
    if PRIVATEERSRC is None:
        raise 'Cannot get Privateer environment --> check Path'
    
    manualgrafting = make_jsonfile_for_grafting(pdbfile=pdbfile,graftedmodel=outputgraft,residuelist=residuelist)

    name = os.path.basename(pdbfile).split('.')[0]
    print(f'Start grafting {name}')
    command = f"""python {PRIVATEERSRC}/grafter.py -manual_grafting {manualgrafting}"""# -save_summary True"""# >/dev/null 2>&1"""
    os.system(command)
    print(f'Finish grafting {name}')


def make_connnection_C_Man(path:str): # return .pdb file
    
    """Return PDB file/
    This function need re-considering to make N-linkages and O-linkages"""
    
    if os.path.isfile(path): pass
    else: raise ('The input not existing')

    st = gemmi.read_structure(path)
    st.assign_serial_numbers(numbered_ter=True)

    ns = gemmi.NeighborSearch(st[0], st.cell, 5).populate(include_h=False) 
    
    linkidx = 1
    for model in st:
        for chain in model:
            for residue in chain:

                hetatom = gemmi.find_tabulated_residue(residue.name)

                if str(hetatom.kind) != 'ResidueKind.AA': residue.het_flag = 'H' # change flag ATOM or HETATM
                if str(hetatom.kind) != 'ResidueKind.PYR': continue
                try: c1 = residue['C1'][0].pos
                except RuntimeError as e:
                    print(f'{e} at {chain,residue}')

                if residue.name != 'MAN': continue

                for atom in residue:
                    if atom.name == 'C1':

                        marks = ns.find_atoms(atom.pos, '\0', radius=5)
                        for mark in marks:
                            cra = mark.to_cra(st[0])
                            if ((cra.residue.name == 'TRP' and cra.atom.name == 'CD1')
                                # or (cra.residue.name == 'ASN' and cra.atom.name == 'ND2')
                                # or (cra.residue.name in ['SER','THR'] and cra.atom.name in ['OG','OG1'])
                                ):
                                dist = (c1).dist(cra.atom.pos)
                                if dist >= 2.0: continue

                    # write LINK record            
                                con = gemmi.Connection()
                                con.name = f'newcovale{linkidx}'
                                con.type = gemmi.ConnectionType.Covale
                                con.asu = gemmi.Asu.Same
                                print(f'Cra {cra.chain, cra.residue, cra.residue.sole_atom(cra.atom.name)}')
                                con.partner1 = gemmi.make_address(cra.chain, cra.residue, cra.residue.sole_atom(cra.atom.name))
                                print(f"atom {chain, residue, residue.sole_atom('C1')}")
                                con.partner2 = gemmi.make_address(chain, residue, residue.sole_atom('C1'))
                                con.reported_distance = (c1).dist(cra.atom.pos)
                                st.connections.append(con)
                    
                    # write CONECT records
                                st.add_conect(cra.atom.serial,atom.serial,order=1)
                                st.add_conect(atom.serial,residue['C2'][0].serial,order=1)
                                st.add_conect(atom.serial,residue['O5'][0].serial,order=1)
                                linkidx += 1

                    elif atom.name == 'C2':
                        st.add_conect(atom.serial,residue['C3'][0].serial,order=1) 
                        st.add_conect(atom.serial,residue['O2'][0].serial,order=1)     

                    elif atom.name == 'C3':
                        st.add_conect(atom.serial,residue['C4'][0].serial,order=1) 
                        st.add_conect(atom.serial,residue['O3'][0].serial,order=1) 
                    
                    elif atom.name == 'C4':
                        st.add_conect(atom.serial,residue['C5'][0].serial,order=1) 
                        st.add_conect(atom.serial,residue['O4'][0].serial,order=1) 
                    
                    elif atom.name == 'C5':
                        st.add_conect(atom.serial,residue['C6'][0].serial,order=1) 
                        st.add_conect(atom.serial,residue['O5'][0].serial,order=1) 
                    
                    elif atom.name == 'C6':
                        st.add_conect(atom.serial,residue['O6'][0].serial,order=1)                      

    # name = Path(path).stem
    st.write_pdb(path,gemmi.PdbWriteOptions(preserve_serial=True,conect_records=True))

    print('Finish making .pdb files with LINK and CONECT records')

def make_external_restraints(pdbfile:str,outfile:str,resolution:float): # return .txt file

    """Return man 1C4 restraints in .txt file """
    glycosylation = pvt.GlycosylationComposition_memsafe(pdbfile) 
    restraints = glycosylation.return_external_restraints(resolution)
    with open(outfile,'w') as f: f.write(restraints)

def run_servalcat_refine_spa(pdbfile:str,resolution:str,outputname:str,restraints:str,
                                halfmap1:str=None,halfmap2:str=None,
                                maskmap:str=None, primmap:str=None): # return .pdb/.cif directories
    
    """Run refinement with Servalcat"""
    print('Start running Servalcat')
    pdbid = os.path.basename(pdbfile).split('.')[0]

    _args = []
    _args += ['--ncycle', '20']
    
    _args += ['--model', pdbfile]

    if halfmap1 != None and halfmap2 != None:
        _args += ['--halfmaps', halfmap1,halfmap2]
    elif primmap != None:
        _args += ['--map', primmap]
    
    if maskmap != None: # cryo-EM structures might not have masked map
        _args += ['--mask_for_fofc',maskmap]
        _args += ['--trim_fofc_mtz']

    _args += ['--resolution', f'{resolution}']
    # _args += ['--pixel_size', '0.9815']
    _args += ['--output_prefix', outputname]
    _args += ['--hydrogen', 'no']
    _args += ['--hout']
    _args += ['--keyword_file', restraints]

    args = ['/Applications/ccp4-8.0/bin/servalcat', 'refine_spa'] + _args
    
    outdirectory = os.getcwd()
    
    # Don't want to print Servalcat run in Terminal? -> change capture_output = False
    try: subprocess.run(args=args,check=True,capture_output=True)
    except subprocess.CalledProcessError as e:
        print(e.stderr.decode("utf-8"))
        if 'Error: Model is out of mask.' in str(e.stderr.decode("utf-8")):
            print(f'Running no check mask with model')
            args += ['--no_check_mask_with_model']
            subprocess.run(args=args,check=True,capture_output=True)
        else:
            sys.exit(1)
    
    print(f'Finish running Servalcat for PDB ID: {pdbid}')

    pdbfile = os.path.join(outdirectory,f'{outputname}.pdb')
    ciffile = os.path.join(outdirectory,f'{outputname}.mmcif')
    # finalmtz = os.path.join(outdirectory,f'diffmap.mtz')

    return pdbfile,ciffile

def calculate_em_rscc(emmap:str,pdbfile:str,resolution:float,outfile:str): # return .csv file with nested list_dict

    """Calculate RSCC and make table summary"""
    print('Run RSCC')
    glycosylation = pvt.GlycosylationComposition(pdbfile,emmap,float(resolution))

    number_of_glycans = glycosylation.get_number_of_glycan_chains_detected()
    sugarsite,sugarchain,rscc,sugarname,sugarconf,glycantype = [],[],[],[],[],[]
    proteinchain, proteinresidue,proteinname = [],[],[]
    status = []

    for glycanNo in range(number_of_glycans):
        glycan = glycosylation.get_glycan(glycanNo)
        chain = glycan.get_root_info()
        gtype = glycan.get_glycosylation_type()

        numsugars = glycan.get_total_number_of_sugars()
            
        for j in range(numsugars):
            sugar = glycan.get_monosaccharide(j)
            sugardict = sugar.get_sugar_summary()
            conf = sugar.get_conformation_name()

            if chain['ProteinResidueType'] in ['TRP','ASN']:

                sta = 'yes'
                if gtype == 'c-glycan' and conf != '1c4': sta = 'no'
                elif gtype == 'n-glycan' and conf != '4c1': sta = 'no'
                elif gtype == 'o-glycan': sta = 'check'
                
                sugarsite.append(sugardict['sugar_seqnum'])
                sugarchain.append(sugardict['sugar_pdb_chain'])
                rscc.append(sugardict['RSCC'].__round__(3))
                sugarname.append(sugardict['sugar_name_short'])

                proteinchain.append(chain['ProteinChainID'])
                proteinresidue.append(chain['ProteinResidueSeqnum'])
                proteinname.append(chain['ProteinResidueType'])

                glycantype.append(gtype)
                sugarconf.append(conf)
                status.append(sta)

                sugardict = {'glycantype':gtype} | sugardict
                sugardict['sugarconf'] = conf
                sugardict['sigma'] = (12/resolution)-1; sugardict['resolution'] = resolution
                sugardict['status'] = sta
            # outlist.append(sugardict)

    # with open(outfile,'w') as f:json.dump(outlist,f)

    # outfile = outfile.replace('.json', '.csv')
    pdbid = os.path.basename(pdbfile).split('.')[0]
    sigma = (12/resolution)-1
    pdbids = [pdbid]*len(sugarsite)
    sig = [sigma]*len(sugarsite)
    data = {'pdbid':pdbids,'glycantype':glycantype,'proteinname':proteinname,'proteinchain':proteinchain,'proteinresidue':proteinresidue,'sugarname':sugarname,
            'sugarchain':sugarchain,'sugarsite':sugarsite,'sugarconformation':sugarconf,'rscc':rscc,'sigma':sig,'resolution':resolution,'status':status}
    df = pd.DataFrame(data)
    df.to_csv(outfile)

def process_grafted_file(inputfile:str) -> str:

    """Reindex grafted Mannose residues"""

    st = gemmi.read_structure(inputfile)
    for chain in st[0]:
        chain.name = chain.name[0]
        resilist = []
        for residue in chain:
            oldidx = 0
            if residue.name == 'MAN':
                # print(residue['C1'].__len__())
                # print(residue['C1'][0])
                resitag = residue.seqid.num
                resiname = residue.name
                checklen = residue['C1'].__len__()
                if checklen > 1:
                    batches = []
                    atomname = []
                    num = 0
                    for atom in residue:
                        if len(atomname) < 11:
                            atomname.append(atom)
                        elif len(atomname) == 11:
                            batches.append(atomname)
                            atomname = []
                            atomname.append(atom)
                    
                    for batch in batches:
                        # print(batch)
                        newresi = gemmi.Residue()
                        newresi.name = resiname
                        newid = str(resitag) + str(num)
                        # print(newid)
                        newresi.seqid.num = int(newid)
                        for subatom in batch:
                            newresi.add_atom(subatom)
                        chain.add_residue(newresi)
                        num +=1
                    
                    resilist.append(oldidx)
        
        if resilist:
            for cidx in resilist[::-1]:
                del chain[cidx]

    st.write_pdb(inputfile)


if __name__ == "__main__":
    
    start = time.time()

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-input", 
                        action="store",
                        dest = "InputPath",
                        required= True,
                        help = "Input of structure expected to be mmCIF file")
    
    parser.add_argument("-output", 
                        action="store",
                        default= None,
                        dest = "OutputPath",
                        help = "Directory to save all output files")
    
    args = parser.parse_args()

    if (args.InputPath == None) or (not os.path.exists(args.InputPath)):
        print('Expected a existing file')
        sys.exit(0)
    else:
        inputfile = args.InputPath
        if Path(inputfile).suffix != '.cif':
            print('Expected mmCIF as inputfile')
            sys.exit(0)

    filename = Path(inputfile).stem

    if args.OutputPath == None:
        print('The result will be saved to the current directory with the name of the program')
        outputpath = os.getenv("PRIVATEERRESULTS",None)
        if outputpath is None: raise 'Cannot get Privateer environment --> check path'
        outputpath = os.path.join(outputpath,f'{filename}_cmannosylation')
    else:
        outputpath = args.OutputPath
    
    if not os.path.exists(outputpath): os.makedirs(outputpath)
    os.chdir(outputpath)

    # Check if Cryo-EM structures
    method,resolution = get_method_and_resolution(cifpath=inputfile)
    if (method != "'ELECTRON MICROSCOPY'"):
        print('Not Cryo-EM structure')
        sys.exit(0)
    if resolution == None: 
        print(f'No resolution in file')
        sys.exit(0)
    
    # Check expression system
    requestedchains = check_expression_system_with_cif(cifpath=inputfile)
    if not requestedchains:
        print('No proteinchain in the structures satisfied the metazoan expression system')
        sys.exit(0)

    # Check if EM structures have EMmap
    pdbid,emdbid = get_emdbid(cifpath=inputfile)
    if emdbid == None:
        print('No EMDB code of EMMAP?')
        sys.exit(0)        

    # download all emmaps of the structures
    hf1,hf2,maskedmap,primmap = download_em_halfmaps_parallel(pdbid=pdbid,emdbid=emdbid)
    if not os.path.exists(primmap): 
        print("No primmary/full map to run downstream process")
        sys.exit(0)

    # check which are potential C-mannosylated tryptophan
    ctrp = blob_search_cryo(cifpath=inputfile,emmap=primmap)
    print(ctrp)

    # graft C-mannose
    graftedmodel = os.path.join(outputpath,f'{pdbid}_grafted.pdb')
    if not os.path.exists(graftedmodel):
        try_grafting_interface(pdbfile=inputfile,outputgraft=graftedmodel,residuelist=ctrp)
        process_grafted_file(inputfile=graftedmodel)
        make_connnection_C_Man(path=graftedmodel)
    # make external restraint for Servalcat refinement
    restraintfile = f'{pdbid}_external_restraint.txt'
    make_external_restraints(pdbfile=graftedmodel,
                            outfile=restraintfile,
                            resolution=resolution)
    
    # run Servalcat refinement
    outputname = f'{filename}_refined'
    refinepdb,refinecif = run_servalcat_refine_spa(pdbfile=graftedmodel,restraints=restraintfile,outputname=outputname,
                            resolution=resolution,halfmap1=hf1,halfmap2=hf2,maskmap=maskedmap,primmap=primmap)
    
    # calculate RSCC
    rsccfile = f'{pdbid}_rscc_data.csv'
    calculate_em_rscc(emmap=primmap,pdbfile=refinepdb,resolution=resolution,outfile=rsccfile)
    print('Recommend: sugars with RSCC < 0.2 should be discarded')
    end = time.time()
    elapsed = end-start
    print(str(timedelta(seconds=elapsed)))
    print('END')