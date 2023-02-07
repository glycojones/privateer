import os
import re
import sys
import argparse
import requests
import datetime
import warnings
import json
from privateer import privateer_core as pvtcore
from privateer import privateer_modelling as pvtmodelling


def _import_list_of_uniprotIDs_to_glycosylate(inputFilePath):
    output = []
    with open(inputFilePath) as file:
        lines = file.readlines()
    for line in lines:
        lineSplitCommaList = line.split(",")
        if lineSplitCommaList:
            for split_line in lineSplitCommaList:
                cleanLine = re.sub("\W+", "", split_line)
                output.append(cleanLine)
        else:
            cleanLine = re.sub("\W+", "", line)
            output.append(cleanLine)

    return output


def _parse_json_for_grafting_instructions(uniprotIDListPath):
    with open(uniprotIDListPath) as json_file:
        data = json.load(json_file)

    return data


def _download_and_prepare_alphafoldDB_model(uniprotID, downloadLocation):
    outputFileName = uniprotID + ".pdb"
    outputFilePath = os.path.join(downloadLocation, outputFileName)
    requestURL = f"https://alphafold.ebi.ac.uk/files/AF-{uniprotID}-F1-model_v1.pdb"
    query = requests.get(requestURL, allow_redirects=True)

    outputLines = []
    downloadedLines = query.iter_lines()
    for line in downloadedLines:
        decodedLine = line.decode("utf-8")
        if decodedLine[:5] != "MODEL":
            outputLines.append(decodedLine)

    with open(outputFilePath, "w") as file:
        file.writelines("%s\n" % l for l in outputLines)

    print(
        f"Successfully downloaded model from AlphaFoldDB with UniProt ID: '{uniprotID}' to {outputFilePath}"
    )
    return outputFilePath


def _query_uniprot_for_glycosylation_locations(uniprotID):
    uniprotRequestURL = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprotID}"
    uniprotResponse = requests.get(uniprotRequestURL,
                                   headers={"Accept": "application/json"})
    if not uniprotResponse.ok:
        uniprotResponse.raise_for_status()
        sys.exit()
    uniprotResponseJSON = uniprotResponse.json()
    uniprotSequence = uniprotResponseJSON["sequence"]
    uniprotFeatures = uniprotResponseJSON["features"]
    uniprotGlycosylations = []
    for item in uniprotFeatures:
        if item["type"] == "CARBOHYD":
            uniprotGlycosylations.append(item)

    outputSequence = uniprotSequence["sequence"]
    outputSequenceLength = uniprotSequence["length"]
    output = {
        "sequenceLength": outputSequenceLength,
        "sequence": outputSequence,
        "glycosylations": uniprotGlycosylations,
    }

    return output


def _get_sequences_in_receiving_model(receiverpath):
    builder_sequence_only = pvtmodelling.Builder(receiverpath, True)
    receiver_sequence = builder_sequence_only.get_receiving_model_sequence_info(
    )

    return receiver_sequence


# privateer::pymodelling::Builder::graft_glycan_to_receiver(int mglycanindex, int receiver_chain_index, int received_residue_index)
def _get_information_about_input_files(receiverpath, donorpath, uniprotID):
    if receiverpath is None and donorpath is None and uniprotID is None:
        raise ValueError(
            "'-info' flag needs to be used in conjuction with '-uniprotID' and/or '-donor_path' and/or '-local_receiver_path'. At least one of those flags need to be provided."
        )
    model_chains = []
    bulletIndex = 1
    if receiverpath is not None:
        print(
            f"\n\n\n{bulletIndex}. Information related to receiver model(glycans will be grafted to this PDB file). Associated with '-local_receiver_path' flag."
        )
        print(f"\nPath: {receiverpath}")
        bulletIndex += 1
        sequences = _get_sequences_in_receiving_model(receiverpath)
        for item in sequences:
            chainIndex = item["index"]
            chainID = item["ChainID"]
            chainSequence = item["Sequence"]
            chainResidues = item["Residues"]
            dict_for_uniprot = {
                "index": chainIndex,
                "chainID": chainID,
                "sequence": chainSequence,
            }
            model_chains.append(dict_for_uniprot)
            print(
                f"\nreceiver_chain_index: {chainIndex} \t\t\t(manual_grafting.json key:value pair)"
            )
            print(f"Chain ID: {chainID}")
            print(f"Chain Sequence: \n{chainSequence}")
            for residue in chainResidues:
                residueIndex = residue["index"]
                residueType = residue["residueType"]
                residueCode = residue["residueCode"]
                residueSeqnum = residue["residueSeqnum"]
                print(
                    f"\treceiving_aa_index: {residueIndex} \t\t(manual_grafting.json key:value pair)"
                )
                print(
                    f"\t\tResidue PDB code: {chainID}/{residueType}-{residueCode}/{residueSeqnum}"
                )
    if donorpath is not None:
        print(
            f"\n\n\n{bulletIndex}. Information related to donor model(glycans will be taken from this PDB file and translocated to another PDB file). Associated with '-donor_path' flag."
        )
        print(f"\nPath: {donorpath}")
        bulletIndex += 1
        glycosylation = pvtcore.GlycosylationComposition_memsafe(donorpath)
        glycans = glycosylation.get_summary_of_detected_glycans()
        for item in glycans:
            glycanIndex = item["GlycanID"]
            glycanWURCS = item["WURCS"]
            glycanType = item["GlycosylationType"]
            glycanRootInfo = item["RootInfo"]
            GlycosidicLinkageTorsions = item["ProteinGlycanLinkageTorsion"]
            print(
                f"\nglycan_index: {glycanIndex} \t\t\t(manual_grafting.json key:value pair)"
            )
            print(f"\tGlycan WURCS: {glycanWURCS}")
            print(f"\tGlycosylation Type: {glycanType}")
            print(
                f"\tGlycan in the donor linked to: {glycanRootInfo['ProteinChainID']}/{glycanRootInfo['ProteinResidueID']}-{glycanRootInfo['ProteinResidueType']} and is modelled as Chain {glycanRootInfo['RootSugarChainID']}"
            )
            print(
                f"\tGlycosidic linkage torsions in the donor model - Phi = {GlycosidicLinkageTorsions['Phi']}, Psi = {GlycosidicLinkageTorsions['Psi']}"
            )
    if uniprotID is not None:
        print(
            f"\n\n\n{bulletIndex}. Information related to glycosylation site data held by UniProt. Associated with '-uniprotID' flag."
        )
        print(f"\nUniProt ID: {uniprotID}")
        bulletIndex += 1
        uniprotQuery = _query_uniprot_for_glycosylation_locations(uniprotID)
        uniprotSequence = uniprotQuery["sequence"]
        uniprotSequenceLength = uniprotQuery["sequenceLength"]
        glycosylations = uniprotQuery["glycosylations"]
        print(f"\nProtein Sequence: \n{uniprotSequence}\n")
        print(f"Protein Sequence Length: {uniprotSequenceLength}\n")
        if model_chains and receiverpath is not None:
            search_result = next(
                (item for item in model_chains
                 if item["sequence"] == uniprotSequence),
                False,
            )
            if search_result == False:
                print(
                    f"WARNING: Unable to find a retrieved {uniprotID} sequence in {receiverpath}.\n"
                )
            else:
                print(
                    f"Successfully managed to find a sequence match for {uniprotID} in receiver_chain_index: {search_result['index']} which is modelled as Chain {search_result['chainID']} in '{receiverpath}'!\n"
                )
        for glycosylation in glycosylations:
            description = glycosylation["description"]
            residueindex = glycosylation["begin"]
            print(f"\n\tUniProt description: {description}")
            print(f"\tGlycosylated at residue seqnum: {residueindex}")


def _get_NGlycosylation_targets_via_consensus_seq(sequences):
    NGlycosylationConsensus = "[N][^P][ST]|[N][A-Z][C]"
    output = []

    for item in sequences:
        currentChainIndex = item["index"]
        currentChainID = item["ChainID"]
        currentSequence = item["Sequence"]

        glycosylationTargets = []

        for match in re.finditer(NGlycosylationConsensus, currentSequence):
            if (currentSequence[match.start()] == item["Residues"][
                    match.start()]["residueCode"]):
                glycosylationTargets.append({
                    "start": match.start(),
                    "end": match.end(),
                    "match": match.group()
                })

        output.append({
            "Sequence": currentSequence,
            "chainIndex": currentChainIndex,
            "currentChainID": currentChainID,
            "glycosylationTargets": glycosylationTargets,
        })

    return output


def _glycosylate_receiving_model_using_manual_instructions(
    receiverpath,
    donorpath,
    outputpath,
    glycanIndex,
    receiverChainIndex,
    receiverResidueIndex,
    enableUserMessages,
    trimGlycanIfClashesDetected,
):
    builder = pvtmodelling.Builder(
        receiverpath,
        donorpath,
        -1,
        trimGlycanIfClashesDetected,
        True,
        enableUserMessages,
        False,
    )

    builder.graft_glycan_to_receiver(glycanIndex, receiverChainIndex,
                                     receiverResidueIndex)

    graftedGlycanSummary = builder.get_summary_of_grafted_glycans()
    builder.export_grafted_model(outputpath)

    return graftedGlycanSummary


def _glycosylate_receiving_model_using_consensus_seq(
    receiverpath,
    donorpath,
    outputpath,
    glycosylationTargets,
    enableUserMessages,
    trimGlycanIfClashesDetected,
):
    builder = pvtmodelling.Builder(
        receiverpath,
        donorpath,
        -1,
        trimGlycanIfClashesDetected,
        True,
        enableUserMessages,
        False,
    )
    for item in glycosylationTargets:
        chainIndex = item["chainIndex"]
        targets = item["glycosylationTargets"]
        for target in targets:
            currentTargetIndex = target["start"]
            builder.graft_glycan_to_receiver(0, chainIndex, currentTargetIndex)

    graftedGlycanSummary = builder.get_summary_of_grafted_glycans()
    builder.export_grafted_model(outputpath)

    return graftedGlycanSummary


def _glycosylate_receiving_model_using_uniprot_info(
    receiverpath,
    donorpath,
    outputpath,
    targets,
    enableUserMessages,
    trimGlycanIfClashesDetected,
):
    builder = pvtmodelling.Builder(
        receiverpath,
        donorpath,
        -1,
        trimGlycanIfClashesDetected,
        True,
        enableUserMessages,
        False,
    )
    for currentTarget in targets:
        chainIndex = 0
        builder.graft_glycan_to_receiver(0, chainIndex, currentTarget)

    graftedGlycanSummary = builder.get_summary_of_grafted_glycans()
    builder.export_grafted_model(outputpath)

    return graftedGlycanSummary


def _print_grafted_glycans_summary(graftedGlycans):
    for idx, graft in enumerate(graftedGlycans):
        proteinChainID = graft["receiving_protein_residue_chain_PDBID"]
        proteinPDBID = graft["receiving_protein_residue_monomer_PDBID"]
        proteinResidueType = graft["receiving_protein_residue_monomer_type"]
        graftedGlycanChainID = graft["glycan_grafted_as_chainID"]

        if len(graft["ClashingResidues"]):
            averageTotalAtomicDistance = graft["AvgTotalAtomicDistance"]
            numberOfClashingResidues = len(graft["ClashingResidues"])
            print(
                f"{idx+1}/{len(graftedGlycans)}: Grafted donor glycan as chain {graftedGlycanChainID} to {proteinChainID}/{proteinResidueType}-{proteinPDBID}. The graft has resulted in {numberOfClashingResidues} clashes with an average atomic distance of from detected clashing residues: {averageTotalAtomicDistance}."
            )
        else:
            print(
                f"{idx+1}/{len(graftedGlycans)}: Grafted donor glycan as chain {graftedGlycanChainID} to {proteinChainID}/{proteinResidueType}-{proteinPDBID}. The graft did not produce any clashes."
            )


def _store_grafted_glycans_summary(graftedGlycans, index, totalCount):
    for graft in graftedGlycans:
        proteinChainID = graft["receiving_protein_residue_chain_PDBID"]
        proteinPDBID = graft["receiving_protein_residue_monomer_PDBID"]
        proteinResidueType = graft["receiving_protein_residue_monomer_type"]
        graftedGlycanChainID = graft["glycan_grafted_as_chainID"]

        output = ""
        if len(graft["ClashingResidues"]):
            averageTotalAtomicDistance = graft["AvgTotalAtomicDistance"]
            numberOfClashingResidues = len(graft["ClashingResidues"])
            output = f"{index+1}/{totalCount}: Grafted donor glycan as chain {graftedGlycanChainID} to {proteinChainID}/{proteinResidueType}-{proteinPDBID}. The graft has resulted in {numberOfClashingResidues} clashes with an average atomic distance of from detected clashing residues: {averageTotalAtomicDistance}."

        else:
            output = f"{index+1}/{totalCount}: Grafted donor glycan as chain {graftedGlycanChainID} to {proteinChainID}/{proteinResidueType}-{proteinPDBID}. The graft did not produce any clashes."

    return output


def _local_input_model_pipeline(receiverpath, donorpath, outputpath,
                                uniprotID):
    sequences = _get_sequences_in_receiving_model(receiverpath)
    if uniprotID is not None:
        uniprotQuery = _query_uniprot_for_glycosylation_locations(uniprotID)
        uniprotSequence = uniprotQuery["sequence"]
        receiverModelSequence = sequences[0]["Sequence"]
        if receiverModelSequence != uniprotSequence:
            raise ValueError(
                "Receiving model sequence does not match the sequence retrieved from UniProt. Please graft glycans using consensus sequence glycan grafting method"
            )
        else:
            uniprotGlycosylations = uniprotQuery["glycosylations"]
            targets = []
            for item in uniprotGlycosylations:
                if (item["description"][0] == "N"
                        or item["description"][0] == "O"
                        or item["description"][0] == "S"
                        or item["description"][0] == "C"
                        or item["description"][0] == "P"):
                    targets.append(int(item["begin"]) - 1)
            graftedGlycans = _glycosylate_receiving_model_using_uniprot_info(
                receiverpath, donorpath, outputpath, targets, True, False)
            _print_grafted_glycans_summary(graftedGlycans)

    else:
        targets = _get_NGlycosylation_targets_via_consensus_seq(sequences)
        graftedGlycans = _glycosylate_receiving_model_using_consensus_seq(
            receiverpath, donorpath, outputpath, targets, True, False)
        _print_grafted_glycans_summary(graftedGlycans)


# P27918 is a good test for TRP, as TRP is a bit more unique and requires different approach. C-mannosylation currently no worky.
def _online_input_model_pipeline(uniprotID, donorpath, defaultInputModelPath,
                                 outputLocation):
    outputFileName = uniprotID + ".pdb"
    outputpath = os.path.join(outputLocation, outputFileName)
    receiverpath = _download_and_prepare_alphafoldDB_model(
        uniprotID, defaultInputModelPath)
    uniprotGlycosylationQuery = _query_uniprot_for_glycosylation_locations(
        uniprotID)
    uniprotGlycosylations = uniprotGlycosylationQuery["glycosylations"]
    targets = []
    for item in uniprotGlycosylations:
        if (item["description"][0] == "N" or item["description"][0] == "O"
                or item["description"][0] == "S"
                or item["description"][0] == "C"
                or item["description"][0] == "P"):
            targets.append(int(item["begin"]) - 1)
    graftedGlycans = _glycosylate_receiving_model_using_uniprot_info(
        receiverpath, donorpath, outputpath, targets, True, False)
    _print_grafted_glycans_summary(graftedGlycans)


if __name__ == "__main__":
    dt_string = datetime.now().strftime("%d_%m_%Y-%H_%M_%S")
    if os.getenv("PRIVATEERDATA", None) is not None:
        ROOTENV = os.getenv("PRIVATEERDATA", None)
    elif os.getenv("CLIBD", None) is not None:
        ROOTENV = os.getenv("CLIBD", None)
        ROOTENV = os.path.join(ROOTENV, "privateer_data")

    else:
        envbased = False
        ROOTENV = os.getcwd()
        defaultDonorLocation = None
        defaultInputModelLocation = None

    if envbased:
        defaultDonorPath = os.path.join(ROOTENV, "glycan_donor_repertoire",
                                        "High_Mannose", "man5", "cluster1.pdb")
        defaultInputModelPath = os.path.join(ROOTENV,
                                             "glycan_donor_repertoire",
                                             "P29016.pdb")

    if os.getenv("PRIVATEERRESULTS", None) is not None:
        RESULTSENV = os.getenv("PRIVATEERRESULTS", None)
        defaultOutputModelDirectory = os.path.join("grafter_job" + "__" +
                                                   dt_string)
    else:
        defaultOutputModelDirectory = "grafter_job" + "__" + dt_string

    defaultuniprotIDsListPath = os.path.join(ROOTENV, "uniprotIDinputs.txt")
    defaultJSONgrafting = os.path.join(ROOTENV, "manual_grafting.json")
    defaultUniprotID = "P29016"

    printInfo = False

    parser = argparse.ArgumentParser(
        prog="grafter.py",
        usage=
        "%(prog)s [options]. Most convenient usage: python grafter.py -import_uniprotIDs_from_file uniprotIDinputs.txt",
        description=
        f"Graft Glycans to AlphaFoldDB models using Privateer Modelling module.",
        epilog=
        f"If -local_receiver_path or -uniprotID are not provided, the script will default to using UniProtID: {defaultUniprotID} as default input. Will download the PDB from AlphaFoldDB and N-glycosylate according to UniProt data.",
    )
    parser.add_argument(
        "-uniprotID",
        action="store",
        default=None,
        dest="user_uniprotID",
        help=
        "If used with -local_receiver_path, N-glycosylate according to UniProt targets. If used without -local_receiver_path, this variable is used in the download of AlphaFoldDB .pdb file and N-Glycosylation according to UniProt targets.",
    )
    parser.add_argument(
        "-local_receiver_path",
        action="store",
        default=None,
        dest="user_localReceiverPath",
        help=
        f"Path to locally saved AlpfaFoldDB model on the computer. If -uniprotID is not provided, will carry out N-glycosylation according to regex consensus sequence of '[N][^P][ST]|[N][A-Z][C]'. The argument overrides default behaviour of downloading AlpfaFoldDB model from the server. WARNING: Ensure that \"MODEL 0\" line is deleted in the local file, as otherwise Privateer's MMBD dependency will not be able to import the model!",
    )
    parser.add_argument(
        "-donor_path",
        action="store",
        default=None,
        dest="user_donorPath",
        help=
        f"Path to the glycan that is to be grafted throughout AlphaFoldDB model. If not specified, the script will default to using glycan located in '{defaultDonorPath}'",
    )
    parser.add_argument(
        "-download_path",
        action="store",
        default=None,
        dest="user_inputModelDirectory",
        help=
        f"Specify download directory where original AlpfaFoldDB models downloaded from the server should be saved. If unspecified, the script will default to '{defaultOutputModelDirectory}'",
    )
    parser.add_argument(
        "-output_path",
        action="store",
        default=None,
        dest="user_outputPath",
        help=
        f"Specify output directory where AlpfaFoldDB models with grafted glycans should be saved. If unspecified, the script will default to '{defaultOutputModelDirectory}'. If the argument is used alongside -local_receiver_path, then the name of PDB output file should be provided, for example 'P29016_output.pdb'",
    )
    parser.add_argument(
        "-import_uniprotIDs_from_file",
        action="store",
        default=None,
        dest="user_uniprotIDsList",
        help=
        f"Glycosylate multiple AlphaFoldDB models from a list of UniProtIDs. Example file is located in '{defaultuniprotIDsListPath}' By default will download files from the server and save them localy in specified or default directory locations.",
    )

    parser.add_argument(
        "-info",
        action="store_true",
        default=False,
        dest="user_infoFlag",
        help=
        f"Print out relevant information about donor PDB(where glycans are taken from) and receiver PDB(where glycans are grafted to). To be used in conjuction with '-local_receiver_path' and/or '-donor_path' and/or '-uniprotID' flags. Usage of this flag overrides grafting functionality, i.e. no grafting will be carried out.",
    )
    parser.add_argument(
        "-manual_grafting",
        action="store",
        default=None,
        dest="user_JSONgrafting",
        help=
        f"Import a JSON file to manually graft glycans with total control over glycosylation sites. Example file is located at '{defaultJSONgrafting}'",
    )

    args = parser.parse_args()

    if args.user_uniprotID is not None:
        uniprotID = args.user_uniprotID
    else:
        uniprotID = defaultUniprotID
    if args.user_donorPath is not None:
        donorPath = args.user_donorPath
    else:
        donorPath = defaultDonorPath
    if args.user_outputPath is not None:
        outputPath = args.user_outputPath
        if (os.path.isdir(args.user_localReceiverPath)
                and args.user_localReceiverPath is not None):
            raise ValueError(
                "ERROR: The combination of provided arguments requires -output_path argument to be a file name, rather than directory!"
            )
    else:
        outputPath = defaultOutputModelDirectory
    if args.user_inputModelDirectory is not None:
        inputModelDirectory = args.user_inputModelDirectory
    else:
        inputModelDirectory = defaultInputModelPath

    if args.user_uniprotIDsList is not None:
        uniprotIDListPath = args.user_uniprotIDsList

    if args.user_JSONgrafting is not None:
        JSONgraftingPath = args.user_JSONgrafting

    if args.user_infoFlag == True and not None:
        printInfo = True

    if (args.user_localReceiverPath is not None and args.user_uniprotID is None
            and printInfo == False):
        uniprotID = None
        _local_input_model_pipeline(args.user_localReceiverPath, donorPath,
                                    outputPath, uniprotID)
    elif (args.user_localReceiverPath is not None
          and args.user_uniprotID is not None and printInfo == False):
        _local_input_model_pipeline(args.user_localReceiverPath, donorPath,
                                    outputPath, uniprotID)
    elif args.user_uniprotIDsList is not None and printInfo == False:
        uniprotIDList = _import_list_of_uniprotIDs_to_glycosylate(
            uniprotIDListPath)
        for idx, uniprotID in enumerate(uniprotIDList):
            _online_input_model_pipeline(uniprotID, donorPath,
                                         inputModelDirectory, outputPath)
            print(
                f"\n{idx+1}/{len(uniprotIDList)}: Successfully finished processing AlphaFoldDB model with UniProt ID of {uniprotID}.\n"
            )
    elif args.user_JSONgrafting is not None and printInfo == False:
        JSONGraftInstructions = _parse_json_for_grafting_instructions(
            JSONgraftingPath)
        initialInputPath = JSONGraftInstructions["receiver_path"]
        initialOutputSubsequentInputOutputPath = JSONGraftInstructions[
            "output_path"]
        glycosylations = JSONGraftInstructions["glycosylations"]
        graftedGlycansSummary = []
        for count, item in enumerate(glycosylations):
            donorPath = item["donor_path"]
            glycanIndex = item["glycan_index"]
            receivingChainIndex = item["receiving_chain_index"]
            receivingAminoAcidIndex = item["receiving_aa_index"]
            if count == 0:
                currentGraftedGlycanSummary = (
                    _glycosylate_receiving_model_using_manual_instructions(
                        initialInputPath,
                        donorPath,
                        initialOutputSubsequentInputOutputPath,
                        glycanIndex,
                        receivingChainIndex,
                        receivingAminoAcidIndex,
                        True,
                        False,
                    ))
                messageString = _store_grafted_glycans_summary(
                    currentGraftedGlycanSummary, count, len(glycosylations))
                graftedGlycansSummary.append(messageString)
            else:
                currentGraftedGlycanSummary = (
                    _glycosylate_receiving_model_using_manual_instructions(
                        initialOutputSubsequentInputOutputPath,
                        donorPath,
                        initialOutputSubsequentInputOutputPath,
                        glycanIndex,
                        receivingChainIndex,
                        receivingAminoAcidIndex,
                        True,
                        False,
                    ))
                messageString = _store_grafted_glycans_summary(
                    currentGraftedGlycanSummary, count, len(glycosylations))
                graftedGlycansSummary.append(messageString)
        print("\n")
        for message in graftedGlycansSummary:
            print(message + "\n")

    elif printInfo == True:
        warnings.warn(
            "-info flag was provided, overriding all arguments regarding grafting and printing info only. Please remove -info flag if you actually want to graft glycans."
        )
        _get_information_about_input_files(args.user_localReceiverPath,
                                           args.user_donorPath,
                                           args.user_uniprotID)
    else:
        if printInfo == False:
            _online_input_model_pipeline(uniprotID, donorPath,
                                         inputModelDirectory, outputPath)
        else:
            warnings.warn(
                "-info flag was provided, overriding all arguments regarding grafting and printing info only. Please remove -info flag if you actually want to graft glycans."
            )