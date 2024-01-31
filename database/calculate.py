from privateer import privateer_core as pvt
from datetime import datetime
import os
import json
import argparse 

def main(args): 
    
    if not os.path.exists(args.pdb): 
        raise FileNotFoundError(f"Can't find a PDB path")
    
    try:
        
        if args.mtz:
            if os.path.exists(args.mtz):
                glycosylation = pvt.GlycosylationComposition(args.pdb, args.mtz, "FP,SIGFP")
        elif args.map and args.mapres:        
            if os.path.exists(args.map):
                glycosylation = pvt.GlycosylationComposition(args.pdb, path_to_mrc_file=args.map, resolution=float(args.mapres), nThreads=4, debug_output=False)
        else:
            glycosylation = pvt.GlycosylationComposition(args.pdb)

    except Exception as pe:
        raise RuntimeError(f"Error in Privater call {pe}")

    
    # if os.path.exists(args.output):
    #     return 

    data = {}

    metadata = {
        "dateAdded": datetime.today().strftime('%Y-%m-%d'), 
        "experimentalDataAvailable": "yes" if args.mtz else "no"
    }

    glycans = {
        "n-glycan": [],
        "o-glycan": [],
        "s-glycan": [],
        "c-glycan": [],
        "ligand": []
    }

    number_of_glycans = glycosylation.get_number_of_glycan_chains_detected()

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

    if not os.path.isdir(args.basedir):
        os.makedirs(args.basedir, exist_ok=True)
        
    with open(args.output, 'w') as fout:
        fout.write(json.dumps(data))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="PrivateerCalculate",
        usage="Calculate Validation Report using Privateer"
    )
    parser.add_argument(
        "-pdb",
        "-pdbpath"
    )
    parser.add_argument(
        "-mtz",
        "-mtzpath"
    )
    parser.add_argument("-map",
                        "-mappath")
    parser.add_argument("-mapres")
    parser.add_argument(
        "-basedir"
    )
    parser.add_argument(
        "-output"
    )

    args = parser.parse_args()

    try:
        main(args)
    except Exception as e:
        raise RuntimeError(f"Error in main call {e}")