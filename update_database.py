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
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3516 and https://glyconnect.expasy.org/browser/structures/3373 for full details"},
                                "G96881BQ": {"id": 1898,
                                            "core": "No-core",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1898 for full details"},
                                "G47613FH": {"id": 3625,
                                            "core": "Hybrid",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3625 for full details"},
                                "G47318KU": {"id": 2543,
                                            "core": "Undefined core",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/2543 for full details"},
                                "G22140GZ": {"id": 3471,
                                            "core": "Complex",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3471 for full details"},            
                                "G27122CE": {"id": 3567,
                                            "core": "Hybrid",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3567 for full details"},                                            
                                "G43157UW": {"id": 1147,
                                            "core": "Complex",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1147 for full details"},                                            
                                "G65575TV": {"id": 3431,
                                            "core": "Hybrid",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3431 for full details"},                                            
                                "G08146BT": {"id": 1519,
                                            "core": "Complex",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1519 for full details"},                                            
                                "G58338QG": {"id": 2491,
                                            "core": "Hybrid",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/2491 for full details"},                                            
                                "G42493LZ": {"id": 939,
                                            "core": "Hybrid",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/939 for full details"},                                            
                                "G06830JN": {"id": 372,
                                            "core": "Complex",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/372 for full details"},                                            
                                "G74587YW": {"id": 1153,
                                            "core": "Complex",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1153 and https://glyconnect.expasy.org/browser/structures/2652 for full details"},                                            
                                "G55220VL": {"id": 1443,
                                            "core": "High-Mannose",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1443 for full details"},                                            
                                "G38979OV": {"id": 2056,
                                            "core": "No-core",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/2056 for full details"},
                                "G57818FI": {"id": 2642,
                                            "core": "No-core",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/2642 for full details"},            
                                "G49108TO": {"id": 8849,
                                            "core": "Undefined core",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/8849 for full details"},                                            
                                "G27608TI": {"id": 446,
                                            "core": "No-core",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/446 for full details"},                                            
                                "G50236GJ": {"id": 3394,
                                            "core": "No-core",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3394 for full details"},                                            
                                "G17199HN": {"id": 1648,
                                            "core": "Hybrid",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1648 for full details"},                                            
                                "G91297KX": {"id": 2358,
                                            "core": "Complex",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/2358 for full details"},                                            
                                "G17409FQ": {"id": 2497,
                                            "core": "Complex",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/2497 for full details"},                                            
                                "G77079RV": {"id": 1446,
                                            "core": "Lactose",
                                            "type": "Free",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1446 and https://glyconnect.expasy.org/browser/structures/1701 for full details"},                                            
                                "G70323CJ": {"id": 3570,
                                            "core": "No-core",
                                            "type": "O-Linked/C-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3570 for full details"},                                            
                                "G88947NR": {"id": 1413,
                                            "core": "High-Mannose",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1413 for full details"},                                            
                                "G88751GL": {"id": 1707,
                                            "core": "Hybrid",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1707 for full details"},            
                                "G24693UW": {"id": 129,
                                            "core": "No-core",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/129 for full details"},                                            
                                "G42358LZ": {"id": 3353,
                                            "core": "Complex",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3353 for full details"},                                            
                                "G39431MA": {"id": 3513,
                                            "core": "Hybrid",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3513 for full details"},                                            
                                "G33780DA": {"id": 432,
                                            "core": "Complex",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/432 for full details"},                                            
                                "G22573RC": {"id": 2007,
                                            "core": "No-core",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/2007 for full details"},                                            
                                "G00281HB": {"id": 2477,
                                            "core": "Hybrid",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/2477 for full details"},                                            
                                "G56903ZB": {"id": 2643,
                                            "core": "Complex",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/2643 for full details"},                                            
                                "G00395TQ": {"id": 1879,
                                            "core": "No-core",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1879 for full details"},                                            
                                "G10651WD": {"id": 223,
                                            "core": "Core 2",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/223 for full details"},                                            
                                "G15021LG": {"id": 534,
                                            "core": "No-core",
                                            "type": "N-Linked/O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/534 for full details"},                                            
                                "G22768VO": {"id": 11,
                                            "core": "Pauci-Mannose",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/11 for full details"},                                           
                                "G35609TJ": {"id": 627,
                                            "core": "Core 1",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/627 for full details"},                                            
                                "G60230HH": {"id": 2472,
                                            "core": "High-Mannose",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/2472 for full details"},                                            
                                "G64581RP": {"id": 3168,
                                            "core": "Truncated if N-Linked",
                                            "type": "N-Linked/O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3168 for full details"},                                            
                                "G70649KP": {"id": 927,
                                            "core": "Core 2",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/927 for full details"},                                            
                                "G84467IZ": {"id": 1378,
                                            "core": "Complex",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1378 for full details"},                                            
                                "G48809BS": {"id": 2820,
                                            "core": "Core 2",
                                            "type": "O-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/2820 for full details"},                                            
                                "G12511VU": {"id": 398,
                                            "core": "Hybrid",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/398 for full details"},                                            
                                "G10471FG": {"id": 3485,
                                            "core": "High-Mannose",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3485 for full details"},                                                                                        
                                "G15683CN": {"id": 20,
                                            "core": "Hybrid",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/20 for full details"},                                            
                                "G13165FV": {"id": 788,
                                            "core": "Complex",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/788 for full details"},                                                                                        
                                "G42440OJ": {"id": 1180,
                                            "core": "Complex",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/1180 for full details"},                                            
                                "G59784AY": {"id": 3309,
                                            "core": "Complex",
                                            "type": "N-Linked",
                                            "comment": "Incomplete information in .json file, please visit: https://glyconnect.expasy.org/browser/structures/3309 for full details"}
                                            }

date = datetime.now()
date = date.strftime('%Y-%m-%d')
fileName = "privateer_database.json"

rootDir = os.getcwd()
outputFolder = os.path.join(rootDir, "src/privateer")
fullPath = os.path.join(outputFolder, fileName)

print(f"Downloading GlyTouCan data as: {fileName} in {fullPath}")

jsonObject = return_json("https://api.glycosmos.org/glytoucan/sparql/glytoucan-data")

print(f"Finished downloading GlyTouCan data.")

numFailedToMatch = 0

jsonOutput = { "privateer_database": None}
array_of_entries = []
for count, line in enumerate(jsonObject):
    print(f'Currently processing: {line["AccessionNumber"]}\nProgress: {count} out of {len(jsonObject)}\t Progress - {int((count/len(jsonObject)*100))}%')
    glytoucanID = line['AccessionNumber']
    responseGlyConnect = return_response_from_glyconnect_api(glytoucanID)
    print("Current HTTP response: " + str(responseGlyConnect.status_code))
    if responseGlyConnect.status_code == 200:
        tempjson = responseGlyConnect.json()
        line['glyconnect'] = str(tempjson["id"])
    elif responseGlyConnect.status_code == 404: 
        line['glyconnect'] = "NotFound"
    elif responseGlyConnect.status_code == 500:
        if glytoucanID in glyconnectHTTP500Exceptions:
            line['glyconnect'] = str(glyconnectHTTP500Exceptions[glytoucanID]["id"])
        else:
            line['glyconnect'] = "Unable to match GlyTouCan ID in Glyconnect database as HTTP 500 error was returned, nor it is described in glyconnectHTTP500Exceptions dict. Please report this to hb1115@york.ac.uk"
            numFailedToMatch += 1
    else: 
        try:
            responseGlyConnect.raise_for_status()
        except requests.exceptions.HTTPError as e:
            line['glyconnect'] = "Error: " + str(e)
            numFailedToMatch += 1
    array_of_entries.append(line)

if os.path.exists(fullPath):
    os.remove(fullPath)

jsonOutput["privateer_database"] = {"entries": array_of_entries}

with open(fullPath, 'w', encoding='utf-8') as file:
    json.dump(jsonOutput, file, ensure_ascii=False, indent=4)

print(f"Finished downloading and appending GlyConnect data. \nAbsolute path of the output: {fullPath}")
print(f"Failed to match {numFailedToMatch} GlyConnect IDs in GlyConnect database out of {len(jsonObject)}. Please CTRL + F for \"Unable\" in output file and modify the exception list if possible.")
print("The script has finished running, bye bye.")
