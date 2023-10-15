from datetime import date
from email.mime import base
import json
import csv
import os
from pkgutil import get_data
from typing import Dict, List


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

def get_datafile_list(ROOT_PATH) -> List[str]: 

    relative_path = 'linkage_torsions/unprocessed_files'
    base_path = os.path.join(ROOT_PATH, relative_path)

    file_list = []

    for file in os.scandir(base_path):
        file_list.append(file.path)

    return file_list

def get_torsion_dict(file_path) -> List[Dict[str,str]]:

    return_list = []

    with open(file_path, 'r') as data_file: 
        data = csv.reader(data_file)

        for index, entry in enumerate(data):
            if index == 0: continue

            try:
                temp_dict = { 
                    "Phi": float(entry[4]), 
                    "Psi": float(entry[5])
                }
            except ValueError: 
                continue

            return_list.append(temp_dict)
    
    return return_list

def main(): 
    if os.getenv("PRIVATEERDATA", None) is not None:
        ROOTPATH = os.getenv("PRIVATEERDATA", None)
    else:
        ROOTPATH = os.getenv("CLIBD", None)
        if ROOTPATH is None:
            raise EnvironmentError(
                "Unable to retrieve 'PRIVATEERDATA' nor 'CLIBD' environment variable. Please try sourcing the ccp4.envsetup-sh file again"
            )
        ROOTPATH = os.path.join(ROOTPATH, "privateer_data")
    
    file_list = get_datafile_list(ROOTPATH)

    temp_data_json = []

    for file_path in file_list: 
        
        file_name_with_extension = file_path.split('/')[-1]

        file_name = file_name_with_extension.split('.')[0]
        
        first_sugar = file_name.split('-')[0]
        second_sugar = file_name.split('-')[2]

        linkage = file_name.split('-')[1]
        first_linkage_number = linkage.split(',')[0]
        second_linkage_number = linkage.split(',')[1]
        
        torsion_list = get_torsion_dict(file_path)

        temp_data_json_item = { 
            "linkage_name": file_name,
            "first": first_sugar, 
            "second": second_sugar, 
            "donor_position": first_linkage_number, 
            "acceptor_position": second_linkage_number, 
            "torsions": torsion_list,
            "type": "protein-sugar" if first_sugar in amino_acids or second_sugar in amino_acids else "sugar-sugar"
        }
        
        if temp_data_json_item["first"] == "ASN" and temp_data_json_item["second"] == "NAG" and temp_data_json_item["donor_position"] == "1" and temp_data_json_item["acceptor_position"] == "1":
            temp_data_json_item["acceptor_position"] = "2"

        temp_data_json.append(temp_data_json_item)
    
    added_files = []

    all_torsion_data = []

    for item in temp_data_json: 
        if item['linkage_name'] not in added_files: 
            output = {
                "first": item['first'], 
                "type": item['type'],
                "second": []
            }
            for link in temp_data_json:
                if item['first'] == link['first']:
                    output['second'].append({
                        "sugar": link['second'], 
                        "donor_position": link['donor_position'],
                        "acceptor_position": link['acceptor_position'],
                        "torsions": link['torsions']
                    })  

                    added_files.append(link['linkage_name'])

            all_torsion_data.append(output)

    output = { 
        "date_last_updated": date.today().strftime("%m/%d/%Y"),
        "database_name": "torsion_database",
        "data": all_torsion_data
    }
    output_json_path = os.path.join(ROOTPATH, "linkage_torsions", "privateer_torsion_database.json")
    with open(output_json_path, 'w', encoding="utf-8",) as output_file: 
        json.dump(output, output_file, indent=3, ensure_ascii=False)

if __name__ == "__main__":
    main()