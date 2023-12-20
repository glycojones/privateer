import json

def main(): 
    data_file = "glycomics/privateer_glycomics_database.json"

    with open(data_file, "r") as data_file: 

        data = json.load(data_file)
    
    output = {}
    d = data["data"]

    for x in d:
        output[x["Sequence"]] = {"GlyToucan": x["AccessionNumber"], "GlyConnect": x["glyconnect"]}

    output_path = "glycomics/privateer_glycomics_database_slim.json"
    with open(output_path, "w") as output_file: 
        json.dump(output, output_file)


main()