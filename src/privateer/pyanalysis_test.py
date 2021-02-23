from privateer import libprivateer as pvt


glycosylation = pvt.GlycosylationComposition("/home/harold/Dev/privateer_python/tests/test_data/5fjj.pdb")
print(glycosylation.get_expression_system_used())
print(glycosylation.get_path_of_model_file_used())
print(glycosylation.get_number_of_glycan_chains_detected())

listOfDetectedGlycans = glycosylation.get_summary_of_detected_glycans()

for entry in listOfDetectedGlycans:
    for key, value in entry.items():
        print('{}: {}'.format(key, value))
    print("_______________________")

glycan = glycosylation.get_glycan(31)

print(glycan.get_glycan_summary())
print(glycan.get_total_of_glycosidic_bonds())
print(glycan.get_unique_monosaccharides())

sugar = glycan.get_monosaccharide(2)
print(sugar.get_sugar_id())
print(sugar.get_glycan_id())
print(sugar.get_type())
print(sugar.get_anomer())
