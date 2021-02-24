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

glycan = glycosylation.get_glycan(0)

print(glycan.get_glycan_summary())
print(glycan.get_total_of_glycosidic_bonds())
print(glycan.get_unique_monosaccharides())

sugar = glycan.get_monosaccharide(1)
print("get_sugar_id() = " + str(sugar.get_sugar_id()))
print("get_glycan_id() = " + str(sugar.get_glycan_id()))
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
print("get_glycosylation_context() = " + sugar.get_glycosylation_context())