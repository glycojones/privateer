from blobsearch import convert_sf_to_mtz

sf_dir = "/vault/pdb_mirror/data/structures/all/structure_factors/"
path_one = sf_dir + "r6plhsf.ent.gz"
path_two = sf_dir + "r6s0bsf.ent.gz"
convert_sf_to_mtz(path_one)
convert_sf_to_mtz(path_two)
