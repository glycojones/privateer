from privateer import privateer_core as pvt

print("____________________GLYCOSYLATION LEVEL_______________")
glycosylation = pvt.GlycosylationComposition("/home/harold/Dev/privateer_python/tests/test_data/5fjj.pdb")
# glycosylation = pvt.GlycosylationComposition("/home/harold/Dev/privateer_python/tests/test_data/2h6o_carbremediation.pdb")
# print("Expression system defined: " + glycosylation.get_expression_system_used())
# print("Path of the input PDB: " + glycosylation.get_path_of_model_file_used())
# print("Number of Glycan chains detected in the model: " + str(glycosylation.get_number_of_glycan_chains_detected()))

listOfDetectedGlycans = glycosylation.get_summary_of_detected_glycans()

print("\n")
for entry in listOfDetectedGlycans:
    for key, value in entry.items():
        print('{}: {}'.format(key, value))
    print("_______________________")

glycan = glycosylation.get_glycan(13)

print("\n")
print("____________________GLYCAN CHAIN LEVEL_______________")
print("Summary of the glycan: " + str(glycan.get_glycan_summary()))
print("Number of glycosidic bonds in the glycan: " + str(glycan.get_total_of_glycosidic_bonds()))
print("Unique monosaccharide codes in the glycan: " + str(glycan.get_unique_monosaccharide_codes()))

print("\n")
print("____________________INDIVIDUAL SUGAR LEVEL_______________")
sugar = glycan.get_monosaccharide(0)
print("get_sugar_id() = " + str(sugar.get_sugar_id()))
print("get_glycan_id() = " + str(sugar.get_glycan_id()))
print("get_sugar_seqnum() = " + str(sugar.get_seqnum()))
print("get_sugar_pdb_id() = " + str(sugar.get_sugar_pdb_id()))
print("get_sugar_pdb_chain() = " + sugar.get_sugar_pdb_chain())
print("get_conformation_name() = " + sugar.get_conformation_name())
print("get_conformation_name_iupac() = " + sugar.get_conformation_name_iupac())
print("get_puckering_amplitude() = " + str(sugar.get_puckering_amplitude()))
print("get_anomer() = " + sugar.get_anomer())
print("get_handedness() = " + sugar.get_handedness())
print("get_denomination() = " + sugar.get_denomination())
print("get_ring_cardinality() = " + str(sugar.get_ring_cardinality()))
print("get_cremer_pople_params() = " + str(sugar.get_cremer_pople_params()))
print("is_sane() = " + str(sugar.is_sane()))
print("get_name_full() = " + sugar.get_name_full())
print("get_name_short() = " + sugar.get_name_short())
print("get_type() = " + sugar.get_type())
print("get_ring_angles() = " + str(sugar.get_ring_angles()))
print("get_ring_bonds() = " + str(sugar.get_ring_bonds()))
print("get_ring_torsion() = " + str(sugar.get_ring_torsion()))
print("get_ring_bond_rmsd() = " + str(sugar.get_ring_bond_rmsd()))
print("get_ring_angle_rmsd() = " + str(sugar.get_ring_angle_rmsd()))
print("get_bfactor() = " + str(sugar.get_bfactor()))
print("is_supported() = " + str(sugar.is_supported()))
print("ok_with_ring() = " + str(sugar.ok_with_ring()))
print("ok_with_bonds_rmsd() = " + str(sugar.ok_with_bonds_rmsd()))
print("ok_with_angles_rmsd() = " + str(sugar.ok_with_angles_rmsd()))
print("ok_with_anomer() = " + str(sugar.ok_with_anomer()))
print("ok_with_chirality() = " + str(sugar.ok_with_chirality()))
print("ok_with_conformation() = " + str(sugar.ok_with_conformation()))
print("ok_with_puckering() = " + str(sugar.ok_with_puckering()))
print("get_glycosylation_type() = " + sugar.get_glycosylation_type())

sugars = glycan.get_all_monosaccharides()
print("Sugars: " + str(sugars))
print("len(sugars) = " + str(len(sugars)))

sugarSummary = sugars[2].get_sugar_summary()

print("\n")

for key, value in sugarSummary.items():
    print('{}: {}'.format(key, value))
print("_______________________")

# XRayData = pvt.XRayData("/Users/haroldas/Dev/privateer_python/tests/test_data/5fjj.mtz", "/Users/haroldas/Dev/privateer_python/tests/test_data/5fjj.pdb", "NONE", 2.5, -1)
XRayData = pvt.XRayData("/home/harold/Dev/privateer_python/tests/test_data/5fjj.mtz", "/home/harold/Dev/privateer_python/tests/test_data/5fjj.pdb", "NONE", 2.5, -1)

experimentaldata = XRayData.get_sugar_summary_with_experimental_data()
print("\n")
for entry in experimentaldata:
    for key, value in entry.items():
        print('{}: {}'.format(key, value))
    print("_______________________")

XRayData.print_cpp_console_output_summary()


glycosylation.update_with_experimental_data(XRayData)
status = glycosylation.check_if_updated_with_experimental_data()
print(f'Glycosylation object updated with experimental data: {status}')
glycanUpdated = glycosylation.get_glycan(32)
sugarUpdated = glycanUpdated.get_monosaccharide(2)
sugarUpdatedTwo = glycanUpdated.get_monosaccharide(4)
print("get_sugar_rscc() = " + str(sugarUpdated.get_sugar_rscc()))
print("get_sugar_accum() = " + str(sugarUpdated.get_sugar_accum()))
print("get_sugar_occupancy_check() = " + str(sugarUpdated.get_sugar_occupancy_check()))
print("get_sugar_rscc()Updated = " + str(sugarUpdatedTwo.get_sugar_rscc()))
print("get_sugar_accum()Updated = " + str(sugarUpdatedTwo.get_sugar_accum()))
print("get_sugar_occupancy_check() Updated = " + str(sugarUpdatedTwo.get_sugar_occupancy_check()))
print(sugarUpdated.get_sugar_summary())
print(sugarUpdatedTwo.get_sugar_summary())

ligands = glycosylation.get_ligands()
print("Ligands.len() = " + str(len(ligands)))

glycosylationMasterConstructor = pvt.GlycosylationComposition("/home/harold/Dev/privateer_python/tests/test_data/5fjj.pdb", "/home/harold/Dev/privateer_python/tests/test_data/5fjj.mtz", "NONE", -1, 2.5, "undefined")
discoteque = glycosylationMasterConstructor.get_summary_of_detected_glycans()

print("\n")
for entry in discoteque:
    for key, value in entry.items():
        print('{}: {}'.format(key, value))
    print("_______________________")





# glycosylation = pvt.GlycosylationComposition("/home/harold/Dev/privateer_python/tests/test_data/2h6o_carbremediation.pdb")

# listOfDetectedGlycans = glycosylation.get_summary_of_detected_glycans()

# print("\n")
# for entry in listOfDetectedGlycans:
#     for key, value in entry.items():
#         print('{}: {}'.format(key, value))
#     print("_______________________")

# glycan = glycosylation.get_glycan(5)

# json = pvt.OfflineDatabase()

# glycoproteomicsoutput = glycan.query_offline_database(json, True, True)

# print(glycoproteomicsoutput)

# permutations = glycoproteomicsoutput['permutations']
# for entry in permutations:
#     for key, value in entry.items():
#         print('{}: {}'.format(key, value))
#     print("_______________________")

# print("Now getting from class without having to repeat compute intensive functon\n\n\n\n")

# glycoproteomicsoutputrepeated = glycan.query_offline_database(json, True, True)

# permutations = glycoproteomicsoutputrepeated['permutations']
# for entry in permutations:
#     for key, value in entry.items():
#         print('{}: {}'.format(key, value))
#     print("_______________________")

# print("Summary of the glycan: " + str(glycan.get_glycan_summary()))

# SNFGstrings = glycan.get_SNFG_strings()
# # print(SNFGstrings)

# ligands = glycosylation.get_ligands()

# for count, item in enumerate(ligands):
#     print('{}: {}'.format(count, item.get_sugar_summary()))
#     print("_______________________")
