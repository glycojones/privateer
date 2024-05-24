from privateer import privateer_core as pvtcore
from pprint import pprint
interactions = pvtcore.GlycosylationInteractions(
    path_to_model_file="test_data/5fji.pdb",
    path_to_output_file="test_data/5fji_hydrogenated.pdb",
    enableHBonds=True )

interactions_summary=interactions.get_all_detected_interactions()

#pprint(interactions_summary)

for d in interactions_summary:
    print(d['CH_Pi'])
