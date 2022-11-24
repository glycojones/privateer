from privateer import privateer_core as pvtcore

interactions = pvtcore.GlycosylationInteractions(
    path_to_model_file="test_data/3pic.pdb",
    path_to_output_file="test_data/3pic_hydrogenated.pdb",
    enableHBonds=True )

interactions_summary=interactions.get_all_detected_interactions()

print(interactions_summary)

for d in interactions_summary:
    print(d['CH_Pi'])
