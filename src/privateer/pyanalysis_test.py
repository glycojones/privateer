from privateer import libprivateer as pvt


classone = pvt.GlycosylationStructure("/home/harold/Dev/privateer_python/tests/test_data/5fjj.pdb")
print(classone.get_expression_system_used())
print(classone.get_path_of_model_file_used())
print(classone.get_number_of_glycan_chains_detected())

listOfDetectedGlycans = classone.get_summary_of_detected_glycans()

for entry in listOfDetectedGlycans:
    for key, value in entry.items():
        print('{}: {}'.format(key, value))
    print("_______________________")

