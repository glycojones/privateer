import os
import sys
import argparse
from datetime import datetime
import pandas as pd
from privateer import privateer_core as pvtcore
from privateer import privateer_modelling as pvtmodelling
sys.path.append("/y/people/lah583/privateer/src/privateer")
import grafter

if __name__ == "__main__":
    dt_string = datetime.now().strftime("%d_%m_%Y-%H_%M_%S")
    if os.getenv("PRIVATEERDATA", None) is not None:
        DATAENV = os.getenv("PRIVATEERDATA", None)
        defaultDonorPath = os.path.join(DATAENV, "glycan_donor_repertoire",
                                        "Alpha-D-Mannose.pdb")
        defaultInputModelPath = os.path.join(DATAENV,
                                        "glycan_donor_repertoire")
    else:
        defaultDonorLocation = None
        defaultInputModelLocation = None

    if os.getenv("PRIVATEERRESULTS", None) is not None:
        RESULTSENV = os.getenv("PRIVATEERRESULTS", None)
        defaultOutputModelDirectory = os.path.join(RESULTSENV,"grafter_job" + "__" +
                                               dt_string)
    else:
        defaultOutputModelDirectory = "grafter_job" + "__" + dt_string

    defaultmode = "CMannosylation"
    defaultSaveSummary = True
    uniprotID = None

    parser = argparse.ArgumentParser(
        prog="mannose_grafting_text.py",
        usage=
        "%(prog)s [options]. Most convenient usage: python mannose_grafting_test.py -input_pdbs_dir PATH_TO_PDB_FILES_TO_GRAFT_MANNOSE_ONTO",
        description=
        f"Graft Mannoses to TRP in pdb files where doing so does not produce any clashes using Privateer Modelling module.",
        epilog=
        f"If -output_path is not provided, the script will default to saving output files in: {defaultOutputModelDirectory}.",
    )

    parser.add_argument(
        "-input_pdbs_dir",
        action="store",
        default=None,
        dest="input_pdbs_dir",
        help=
        f"Path to input receiver files that will be grafted. Should be a directory containing at least one pdb file. This is a required argument",
    )
    parser.add_argument(
        "-input_mtzs_dir",
        action="store",
        default=None,
        dest="input_mtzs_dir",
        help=
        f"Path to input density maps for blob search=. Should be a directory containing at least one mtz file. This is a required argument",
    )
    parser.add_argument(
        "-donor_path",
        action="store",
        default=None,
        dest="user_donorPath",
        help=
        f"Path to the glycan that is to be grafted throughout receiver model pdb. If not specified, the script will default to using glycan located in '{defaultDonorPath}'",
    )

    parser.add_argument(
        "-output_path",
        action="store",
        default=None,
        dest="user_outputPath",
        help=
        f"Specify output directory where models with grafted glycans should be saved along with other output files. If unspecified, the script will default to '{defaultOutputModelDirectory}'. If the argument is used alongside -local_receiver_path, then the name of PDB output file should be provided, for example 'P29016_output.pdb'",
    )

    parser.add_argument(
        "-save_summary",
        action="store",
        default=None,
        dest="save_summary",
        help=
        f"Save glycan summary to a csv. If unspecified, the script will default to '{defaultSaveSummary}'",
    )

    args = parser.parse_args()
    if args.input_pdbs_dir is not None:
        receiverPath = args.input_pdbs_dir
        if not os.path.isdir(receiverPath):
            raise ValueError("ERROR: The provided directory of input pdb files to graft to does not exist!")
    else:
        raise ValueError("ERROR: An input directory of pdb files to graft to must be provided!")
    
    if args.input_mtzs_dir is not None:
        mtzPath = args.input_mtzs_dir
        if not os.path.isdir(mtzPath):
            raise ValueError("ERROR: The provided directory of input pdb files to graft to does not exist!")
    else:
        mtzPath = None

    
    if args.user_donorPath is not None:
        donorPath = args.user_donorPath
    else:
        donorPath = defaultDonorPath
    
    if args.user_outputPath is not None:
        outputPath = args.user_outputPath
    else:
        outputPath = defaultOutputModelDirectory
    if not os.path.isdir(outputPath):
        os.mkdir(outputPath)
    
    if args.save_summary is not None:
        SaveSummary = args.save_summary
    else:
        SaveSummary = defaultSaveSummary

    pdbFilePathList = []
    for root, dirs, files in os.walk(receiverPath):
        for f in files:
            if f.partition(".")[2] == "pdb":
                pdbFilePathList.append(os.path.join(root,f))
            elif f.partition(".")[2] == "ent.gz":
                pdbFilePathList.append(os.path.join(root,f))
    
    if mtzPath is not None:
        mtzFilePathList = []
        for root, dirs, files in os.walk(mtzPath):
            for f in files:
                if f.partition(".")[2] == "mtz":
                    mtzFilePathList.append(os.path.join(root,f))

    AllGlycans = []
    offset = 0
    for receiverFile in pdbFilePathList:
        if mtzPath is not None:
            for mtzFile in mtzFilePathList:
                if os.path.basename(mtzFile).partition(".")[0] == os.path.basename(receiverFile).partition(".")[0].partition("pdb")[2]:
                    print(f"Grafting {receiverFile}")
                    try:
                        graftedGlycans = grafter._local_input_model_pipeline(receiverFile, donorPath,
                                                outputPath, None, "CMannosylation", mtzFile)
                    except:
                        print(f"Graft failed for {receiverFile}. Continuuing to next file...")
                        continue
                    outputFileName = os.path.basename(receiverFile).partition('.')[0] + '_grafted.pdb'
                    outFileName = os.path.join(outputPath, outputFileName)
                    if os.path.isfile(outFileName):
                        try:
                            glycosylation = pvtcore.GlycosylationComposition(outFileName, mtzFile, "FP,SIGFP")
                            glycans = []
                            num_glycans = glycosylation.get_number_of_glycan_chains_detected()  
                            for i in range(len(graftedGlycans)):
                                graftedglycan = graftedGlycans[i]
                                for glycan_num in range(num_glycans):
                                    glycan = glycosylation.get_glycan(glycan_num)
                                    numsugars = glycan.get_total_number_of_sugars()
                                    for j in range(numsugars):
                                        sugar = glycan.get_monosaccharide(j)
                                        root_info = glycan.get_root_info()
                                        summary = sugar.get_sugar_summary()
                                        if summary["sugar_name_short"] == "MAN" and root_info["ProteinResidueID"] == graftedglycan["receiving_protein_residue_monomer_PDBID"] and root_info["ProteinChainID"] == graftedglycan["receiving_protein_residue_chain_PDBID"]:
                                            graftedGlycans[i]["RSCC"] = summary["RSCC"]
                        except:
                            for i in range(len(graftedGlycans)):
                                graftedGlycans[i]["RSCC"] = 0
                        for i in range(len(graftedGlycans)):
                            graftedGlycans[i]["pdbCode"] = os.path.basename(receiverFile).partition('.')[0]
                        df = pd.DataFrame.from_dict(graftedGlycans)
                        cols = df.columns.tolist()
                        cols = cols[-1:] + cols[:-1]
                        df = df[cols]
                        if len(df.index) > 0:
                            if SaveSummary:
                                outputFileName = os.path.basename(receiverFile).partition(".")[0] + "_graft_summary.csv"
                                saveCSVto = os.path.join(outputPath, outputFileName)
                                try:
                                    df.to_csv(saveCSVto)
                                except:
                                    continue
                            df.reset_index(drop=True)
                            df.index.name = "index"
                            df["index"] = range(offset, len(df.index)+offset)
                            df = df.set_index("index", append=True)
                            AllGlycans.append(df)
                            offset += len(df.index)

    full_df = pd.concat(AllGlycans)
    outputFileName = "full_graft_summary.csv"
    saveCSVto = os.path.join(outputPath, outputFileName)
    full_df.to_csv(saveCSVto)
    indx_to_drop = []
    for index, row in full_df.iterrows():
        if not row['GraftStatus']:
            indx_to_drop.append(index)
    compact_df = full_df.drop(indx_to_drop)
    outputFileName = "compact_graft_summary.csv"
    saveCSVto = os.path.join(outputPath, outputFileName)
    compact_df.to_csv(saveCSVto)
    
