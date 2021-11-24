import requests
import shutil
import os
import csv
import re
import json
from datetime import date


def return_response_from_glyconnect_api(glytoucanID):
    headers = {
        "accept": "application/json",
        "Content-Type": "application/json",
    }
    data = f'{{ "glytoucan_id": "{glytoucanID}"}}'
    response = requests.post(
        "https://glyconnect.expasy.org/api/structures/search/glytoucan",
        headers=headers,
        data=data,
    )
    return response


def return_glyconnect_id_from_gtc_external_id(glytoucanID):
    jsonObject = return_json(
        f"https://sparqlist.glycosmos.org/sparqlist/api/gtc_external_id?accNum={glytoucanID}"
    )
    if "GlyConnect" in jsonObject:
        listResponse = jsonObject["GlyConnect"]["list"]
        for item in listResponse:
            if "id" in item:
                glyconnectID = item["id"]
                return glyconnectID
            else:
                return "Error: return_glyconnect_id_from_gtc_external_id() found GlyConnect, but did not find id."
    else:
        return "NotFound"


def parse_json_file(path):
    with open(path) as json_file:
        data = json.load(json_file)

    return data


def is_json(myjson):
    try:
        json_object = json.loads(myjson)
    except ValueError as e:
        return False
    return True


def CreateFolder(path):
    if not os.path.exists(path):
        print(f"\n...{path} does not exist. Creating directory...\n")
        os.makedirs(path)
    else:
        print(
            f"\n...{path} already exists. Removing old directory to create a new directory...\n"
        )
        shutil.rmtree(path)  # Removes all the subdirectories!
        os.makedirs(path)


def return_json(URL):
    response = requests.get(URL)

    try:
        response.raise_for_status()
    except requests.exceptions.HTTPError as e:
        # Whoops it wasn't a 200
        return "Error: " + str(e)

    # Must have been a 200 status code
    json_obj = response.json()
    return json_obj


def return_data_from_glyconnect(outputFolder):
    def return_glyconnect_id_from_list(list):
        glytoucan_pattern = "[G]\d{5}[A-Z][A-Z]"
        for item in list:
            regex_result = re.search(glytoucan_pattern, item)
            if regex_result is not None:
                return regex_result.group()

        return None

    glyconnect_address = (
        "https://glyconnect.expasy.org/browser/export?type=structure&query=id:"
    )
    _output_file = os.path.join(outputFolder, "glyconnect_query.json")
    print(f"Downloading GlyConnect data and making backup in {_output_file}")

    exceptions = {8849: "G49108TO", 8850: "G70323CJ"}

    _ranges = [{"start": 0, "end": 3687}, {"start": 8847, "end": 11050}]
    # _ranges = [{"start": 0, "end": 25}, {"start": 8845, "end": 8897}]
    _structure_entries_begin = 7
    output = {}

    for currentRange in _ranges:
        for idx in range(currentRange["start"], currentRange["end"]):
            if idx not in exceptions:
                URL = f"{glyconnect_address}{idx}"
                print(URL)
                download = requests.get(URL)
                decoded_content = download.content.decode("utf-8")

                split_csv = csv.reader(decoded_content.splitlines(), delimiter=";")
                downloaded_list = list(split_csv)
                condensed_list = downloaded_list[_structure_entries_begin:]
                for row in condensed_list:
                    glytoucan_id = return_glyconnect_id_from_list(row)
                    if glytoucan_id is not None:
                        if glytoucan_id not in output:
                            output.update({glytoucan_id: str(idx)})
                        else:
                            if output[glytoucan_id] != str(idx):
                                if type(output[glytoucan_id]) is str:
                                    old_value = output[glytoucan_id]
                                    output[glytoucan_id] = []
                                    output[glytoucan_id].append(str(old_value))
                                    output[glytoucan_id].append(str(idx))
                                else:
                                    if str(idx) not in output[glytoucan_id]:
                                        output[glytoucan_id].append(str(idx))
                            else:
                                continue
            else:
                output.update({exceptions[idx]: str(idx)})

    with open(_output_file, "w") as file:
        file.write(json.dumps(output))

    print(f"Finished obtaining GlyConnect Data")

    return output


def find_glytoucan_id_in_glyconnect_data(glyconnect_data, glytoucan_id):
    if glytoucan_id in glyconnect_data:
        return glyconnect_data[glytoucan_id]
    else:
        return "NotFound"


if __name__ == "__main__":
    fileName = "privateer_glycomics_database.json"

    PRIVATEERDATAPATH = os.getenv("PRIVATEERDATA")
    if PRIVATEERDATAPATH is not None:
        outputFolder = os.path.join(PRIVATEERDATAPATH, "glycomics")
        fullPath = os.path.join(outputFolder, fileName)
    else:
        raise EnvironmentError(
            "Unable to retrieve 'PRIVATEERDATA' environment variable. Please try sourcing the ccp4.envsetup-sh file again"
        )

    import_backup = False

    if import_backup == False:
        glyconnect_data = return_data_from_glyconnect(outputFolder)
    else:
        backup_file = os.path.join(outputFolder, "glyconnect_query.json")
        with open(backup_file) as f:
            glyconnect_data = json.load(f)

    print(glyconnect_data)

    print(f"Downloading GlyTouCan data as: {fileName} in {fullPath}")

    jsonObject = return_json(
        "https://api.glycosmos.org/glytoucan/sparql/glytoucan-data"
    )

    print(f"Finished downloading GlyTouCan data.")

    jsonOutput = {
        "date_last_updated": date.today().strftime("%m/%d/%Y"),
        "database_name": "glycomics_database",
        "data": [],
    }
    array_of_entries = []
    for count, line in enumerate(jsonObject):
        print(
            f'Currently processing: {line["AccessionNumber"]}\nProgress: {count} out of {len(jsonObject)}\t Progress - {int((count/len(jsonObject)*100))}%'
        )
        glytoucanID = line["AccessionNumber"]

        line["glyconnect"] = find_glytoucan_id_in_glyconnect_data(
            glyconnect_data, glytoucanID
        )

        array_of_entries.append(line)

    if os.path.exists(fullPath):
        os.remove(fullPath)

    jsonOutput["data"] = array_of_entries

    with open(fullPath, "w", encoding="utf-8") as file:
        json.dump(jsonOutput, file, ensure_ascii=False, indent=4)

    print(
        f"Finished downloading and appending GlyConnect data. \nAbsolute path of the output: {fullPath}"
    )
    print("The script has finished running, bye bye.")
