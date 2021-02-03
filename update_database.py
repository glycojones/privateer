import requests
import shutil
import os
import json
from datetime import datetime

def return_response_from_glyconnect_api(glytoucanID):
    headers = {
        'accept': 'application/json',
        'Content-Type': 'application/json',
    }
    data = f'{{ "glytoucan_id": "{glytoucanID}"}}'
    response = requests.post('https://glyconnect.expasy.org/api/structures/search/glytoucan', headers=headers, data=data)
    return response
   

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
        print(f"\n...{path} already exists. Removing old directory to create a new directory...\n")
        shutil.rmtree(path)           # Removes all the subdirectories!
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

glyconnectHTTP500Exceptions = { "G82348BZ": {"id": 3247,
                                            "core": "Core 0",
                                            "type": "Pauci-Mannose",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3247 for full details" },
                                "G57321FI": {"id": 2305,
                                            "core": "Core 0",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/2305 for full details" },
                                "G25559RS": {"id": 3277,
                                            "core": "No-core",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3277 and https://glyconnect.expasy.org/browser/structures/2470 for full details" },
                                "G17413FR": {"id": 3314,
                                            "core": "Complex",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3314 and https://glyconnect.expasy.org/browser/structures/2459 for full details" },
                                "G96153DN": {"id": 1954,
                                            "core": "High-Mannose",
                                            "type": "N-Linked",
                                            "iupac": "Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)[Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1954 and https://glyconnect.expasy.org/browser/structures/2702 for full details" },                                            
                                "G57930PW": {"id": 3617,
                                            "core": "High-Mannose",
                                            "type": "N-Linked",
                                            "iupac": "Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)[Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3617 and https://glyconnect.expasy.org/browser/structures/2332 for full details" },
                                "G10514PF": {"id": 1348,
                                            "core": "Core 2",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1348 and https://glyconnect.expasy.org/browser/structures/674 for full details" },
                                "G83161QT": {"id": 654,
                                            "core": "High-Mannose",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/654 and https://glyconnect.expasy.org/browser/structures/2327 for full details" },
                                "G53402KW": {"id": 915,
                                            "core": "No-core",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/915 and https://glyconnect.expasy.org/browser/structures/2440 for full details" },
                                "G51851RG": {"id": 1475,
                                            "core": "No-core",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1475 and https://glyconnect.expasy.org/browser/structures/1759 for full details" },
                                "G56749GV": {"id": 3302,
                                            "core": "Complex",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3302 and https://glyconnect.expasy.org/browser/structures/311 for full details" },                                            
                                "G18948TG": {"id": 1602,
                                            "core": "High-Mannose",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1602 and https://glyconnect.expasy.org/browser/structures/2203 for full details" },
                                "G52880ZN": {"id": 1655,
                                            "core": "Complex",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1655 and https://glyconnect.expasy.org/browser/structures/211 for full details" },
                                "G77865AH": {"id": 2625,
                                            "core": "No-core",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/2625 and https://glyconnect.expasy.org/browser/structures/1188 for full details" },
                                "G00031MO": {"id": 1947,
                                            "core": "Core 1",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1947 for full details"},
                                "G68668TB": {"id": 3516,
                                            "core": "High-Mannose",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3516 and https://glyconnect.expasy.org/browser/structures/3373 for full details"}}

date = datetime.now()
date = date.strftime('%Y-%m-%d')
fileName = "database.json"

rootDir = os.getcwd()
outputFolder = os.path.join(rootDir, "src/privateer")
fullPath = os.path.join(outputFolder, fileName)

print(f"Downloading GlyTouCan data as: {fileName} in {fullPath}")

jsonObject = return_json("https://api.glycosmos.org/glytoucan/sparql/glytoucan-data")

print(f"Finished downloading GlyTouCan data.")

numFailedToMatch = 0

jsonObjectNew = []
for count, line in enumerate(jsonObject):
    print(f'Currently processing: {line["AccessionNumber"]}\nProgress: {count} out of {len(jsonObject)}\t Progress - {int((count/len(jsonObject)*100))}%')
    glytoucanID = line['AccessionNumber']
    responseGlyConnect = return_response_from_glyconnect_api(glytoucanID)
    print("Current HTTP response: " + str(responseGlyConnect.status_code))
    if responseGlyConnect.status_code == 200:
        line['glyconnect'] = responseGlyConnect.json()
    elif responseGlyConnect.status_code == 404: 
        line['glyconnect'] = "NotFound"
    elif responseGlyConnect.status_code == 500:
        if glytoucanID in glyconnectHTTP500Exceptions:
            line['glyconnect'] = glyconnectHTTP500Exceptions[glytoucanID]
        else:
            line['glyconnect'] = "Unable to match GlyTouCan ID in Glyconnect database as HTTP 500 error was returned, nor it is described in glyconnectHTTP500Exceptions dict. Please report this to hb1115@york.ac.uk"
            numFailedToMatch += 1
    else: 
        try:
            responseGlyConnect.raise_for_status()
        except requests.exceptions.HTTPError as e:
            line['glyconnect'] = "Error: " + str(e)
            numFailedToMatch += 1
    jsonObjectNew.append(line)

if os.path.exists(fullPath):
    os.remove(fullPath)

with open(fullPath, 'w', encoding='utf-8') as file:
    json.dump(jsonObjectNew, file, ensure_ascii=False, indent=4)

print(f"Finished downloading and appending GlyConnect data. \nAbsolute path of the output: {fullPath}")
print(f"Failed to match {numFailedToMatch} GlyConnect IDs in GlyConnect database out of {len(jsonObject)}. Please CTRL + F for \"Unable\" in output file and modify the exception list if possible.")
print("The script has finished running, bye bye.")
