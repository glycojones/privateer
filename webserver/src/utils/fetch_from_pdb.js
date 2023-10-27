
export async function fetch_pdb(PDBCode) { 
    if (PDBCode == null) {return}
    console.log("Fetching PDB ", PDBCode)
    let pdb_url = `https://files.rcsb.org/download/${PDBCode.toUpperCase()}.pdb`

    let file = fetch(pdb_url)
    .then(response => {
        if (!response.ok) {
            throw new Error('Network error');
          }
        return response.text()
    })
    .then(file => {
        return Promise.resolve(file)
    })
    .catch((e) => {
        throw new Error("PDB Not Found")
    })
    return file
}

export async function fetch_map(PDBCode) { 
    if (PDBCode == null) {return}
    console.log("Fetching MTZ ", PDBCode)
    let mtz_url = `https://www.ebi.ac.uk/pdbe/entry-files/${PDBCode.toLowerCase()}.ccp4`

    let file = fetch(mtz_url)
    .then(response => {
        if (!response.ok) {
            throw new Error('Network error');
          }
        return response.arrayBuffer()
    })
    .then(file => {
        return Promise.resolve(file)
    })
    .catch((e) => {
        throw new Error("Map Not Found")
    })
    return file
}
