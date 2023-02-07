import os
import json
from datetime import date
from re import split

amino_acids = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
]


def split_initial_list_into_two_categories(importedList):
    amino_acid_containing_list = []
    sugar_sugar_only_list = []

    for item in importedList:
        firstResidue = item[:3]
        if firstResidue in amino_acids:
            amino_acid_containing_list.append(item)
        else:
            sugar_sugar_only_list.append(item)

    return sorted(amino_acid_containing_list), sorted(sugar_sugar_only_list)


def get_unique_keys_for_output_dictionary(contains_amino_acid_list,
                                          sugar_sugar_only_list):
    output = []

    for item in contains_amino_acid_list:
        currentAminoAcidCode = item[:3]
        if currentAminoAcidCode not in output:
            output.append(currentAminoAcidCode)
        else:
            continue

    for item in sugar_sugar_only_list:
        currentSugarCode = item[:3]
        if currentSugarCode not in output:
            output.append(currentSugarCode)
        else:
            continue

    return output


def populate_torsions(jsonTemplate, contains_amino_acid_list,
                      sugar_sugar_only_list, unprocessedDataPath):
    output = {
        "date_last_updated": date.today().strftime("%m/%d/%Y"),
        "database_name": "torsion_database",
        "data": [],
    }

    for item in jsonTemplate:
        currentKey = item["first"]
        currentValues = item["second"]
        if currentKey in amino_acids:
            currentKeyDict = {
                "first": currentKey,
                "type": "protein-sugar",
                "second": [],
            }
            for second_sugar in currentValues:
                currentFileName = f"{currentKey}-{second_sugar}_reduced.json"
                if currentFileName in contains_amino_acid_list:
                    currentFullPath = os.path.join(unprocessedDataPath,
                                                   currentFileName)
                    with open(currentFullPath) as jsonFile:
                        imported_torsions_for_current_pair = json.load(
                            jsonFile)
                    torsions_for_current_pair = []
                    for torsion in imported_torsions_for_current_pair:
                        torsions_for_current_pair.append({
                            "Phi":
                            float(torsion["Phi"]),
                            "Psi":
                            float(torsion["Psi"])
                        })
                    current_second_sugar_dict = {
                        "sugar": second_sugar,
                        "torsions": torsions_for_current_pair,
                    }
                    currentKeyDict["second"].append(current_second_sugar_dict)
                else:
                    raise ValueError(
                        f"Something went wrong with the script when attempting to populate torsions - {currentFileName} not found in {contains_amino_acid_list}"
                    )
            output["data"].append(currentKeyDict)
        else:
            currentKeyDict = {
                "first": currentKey,
                "type": "sugar-sugar",
                "second": []
            }
            for second_sugar in currentValues:
                currentFileName = f"{currentKey}-{second_sugar}_reduced.json"
                if currentFileName in sugar_sugar_only_list:
                    currentFullPath = os.path.join(unprocessedDataPath,
                                                   currentFileName)
                    with open(currentFullPath) as jsonFile:
                        imported_torsions_for_current_pair = json.load(
                            jsonFile)
                    torsions_for_current_pair = []
                    for torsion in imported_torsions_for_current_pair:
                        torsions_for_current_pair.append({
                            "Phi":
                            float(torsion["Phi"]),
                            "Psi":
                            float(torsion["Psi"])
                        })
                    current_second_sugar_dict = {
                        "sugar": second_sugar,
                        "torsions": torsions_for_current_pair,
                    }
                    currentKeyDict["second"].append(current_second_sugar_dict)
                else:
                    raise ValueError(
                        f"Something went wrong with the script when attempting to populate torsions - {currentFileName} not found in {sugar_sugar_only_list}"
                    )
            output["data"].append(currentKeyDict)

    return output


if __name__ == "__main__":
    if os.getenv("PRIVATEERDATA", None) is not None:
        ROOTPATH = os.getenv("PRIVATEERDATA", None)
    else:
        ROOTPATH = os.getenv("CLIBD", None)
        if ROOTPATH is None:
            raise EnvironmentError(
                "Unable to retrieve 'PRIVATEERDATA' nor 'CLIBD' environment variable. Please try sourcing the ccp4.envsetup-sh file again"
            )
        ROOTPATH = os.path.join(ROOTPATH, "privateer_data")

    unprocessedDataPath = os.path.join(ROOTPATH,
                                       "linkage_torsions/unprocessed_files")
    processedDataPath = os.path.join(ROOTPATH, "linkage_torsions")

    if os.path.isdir(unprocessedDataPath):
        files = os.listdir(unprocessedDataPath)
        contains_amino_acid, sugar_sugar_only = split_initial_list_into_two_categories(
            files)
        dict_keys = get_unique_keys_for_output_dictionary(
            contains_amino_acid, sugar_sugar_only)
        jsonTemplate = []
        for key in dict_keys:
            pairs_for_current_key = []
            for entry in contains_amino_acid:
                currentCode = entry[:3]
                if key == currentCode:
                    pairs_for_current_key.append(entry[4:7])
                else:
                    continue

            for entry in sugar_sugar_only:
                currentCode = entry[:3]
                if key == currentCode:
                    pairs_for_current_key.append(entry[4:7])
                else:
                    continue

            currentDict = {"first": key, "second": pairs_for_current_key}
            jsonTemplate.append(currentDict)

        exportJSON = populate_torsions(jsonTemplate, contains_amino_acid,
                                       sugar_sugar_only, unprocessedDataPath)
        print(jsonTemplate)
        with open(
                os.path.join(processedDataPath, "torsion_database.json"),
                "w",
                encoding="utf-8",
        ) as export_json_file:
            json.dump(exportJSON,
                      export_json_file,
                      indent=4,
                      ensure_ascii=False)
        print(
            f'Torsions database has been successfully updated, and saved locally as: {os.path.join(processedDataPath, "torsion_database.json")}'
        )
    else:
        raise ValueError(
            f"Expected {unprocessedDataPath} to be a directory, for some reason the script is not detecting it as a directory."
        )
