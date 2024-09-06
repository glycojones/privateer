from privateer import privateer_core as pvtcore
import json
import pytest

# Tests TODO
# .def("get_all_detected_interactions",  &pa::GlycosylationInteractions::get_all_detected_interactions)
# .def("get_all_detected_hbonds",  &pa::GlycosylationInteractions::get_all_detected_hbonds)
# .def("get_all_detected_chpibonds",  &pa::GlycosylationInteractions::get_all_detected_chpibonds)
# .def("get_all_interactions_for_specific_glycan",  &pa::GlycosylationInteractions::get_all_interactions_for_specific_glycan)
# .def("get_hbonds_for_specific_glycan", &pa::GlycosylationInteractions::get_hbonds_for_specific_glycan)
# .def("get_chpibonds_for_specific_glycan", &pa::GlycosylationInteractions::get_chpibonds_for_specific_glycan)
#       
model_input = "tests/test_data/5fji.pdb"
model_output = "tests/test_data/5fji_hydrogenated.pdb"

interactions = pvtcore.GlycosylationInteractions(
    path_to_model_file=model_input, 
    path_to_output_file=model_output,
    enableHBonds=True )

@pytest.mark.fast
def test_get_path_of_model_file_used(): 
    result = interactions.get_path_of_model_file_used()
    assert result == model_input

@pytest.mark.fast
def test_get_protein_sequence_for_entire_model(): 
    result = interactions.get_protein_sequence_information_for_entire_model()
    with open("tests/expected_outputs/GlycosylationInteractions/get_protein_sequence_for_entire_model.json", "r") as in_file:
        expected = json.load(in_file)

    assert result == expected

@pytest.mark.fast
def test_get_protein_sequence_information_for_single_chain():
    result = interactions.get_protein_sequence_information_for_single_chain(0)
    with open("tests/expected_outputs/GlycosylationInteractions/get_protein_sequence_information_for_single_chain.json", "r") as in_file:
        expected = json.load(in_file)

    assert result == expected

@pytest.mark.slow
def test_get_neighborhood_for_specific_glycan():
    result = interactions.get_neighborhood_for_specific_glycan(0)
    with open("tests/expected_outputs/GlycosylationInteractions/get_neighborhood_for_specific_glycan.json", "r") as in_file:
        expected = json.load(in_file)

    assert result == expected

@pytest.mark.slow
def test_get_all_glycan_neighborhoods():
    result = interactions.get_all_glycan_neighborhoods()
    with open("tests/expected_outputs/GlycosylationInteractions/get_all_glycan_neighborhoods.json", "r") as in_file:
        expected = json.load(in_file)

    assert result == expected

def try_func(): 
    result = interactions.get_all_glycan_neighborhoods()
    with open("tests/expected_outputs/GlycosylationInteractions/get_all_glycan_neighborhoods.json", "w") as in_file:
        json.dump(result, in_file)

    # assert result == expected
if __name__ == "__main__":
    try_func()