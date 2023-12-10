
export async function fetch_pdb(PDBCode: string) {
    if (PDBCode == null) { return }
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

export async function fetch_map(PDBCode: string) {
    if (PDBCode == null) { return }
    console.log("Fetching MTZ ", PDBCode)
    let mtz_url = `https://www.ebi.ac.uk/pdbe/entry-files/${PDBCode.toLowerCase()}.ccp4`

    try {
        let controller = new AbortController();
        let timeId = setTimeout(() => {
            controller.abort()
        }, 60000)

        let file = fetch(mtz_url, { signal: controller.signal })
            .then(response => {
                clearTimeout(timeId)

                if (!response.ok) {
                    throw new Error('Network error');
                }
                return response.arrayBuffer()
            })
            .then(file => {
                return Promise.resolve(file)
            })
            .catch((e) => {
                return Promise.reject("Map not found")
                throw new Error("Map Not Found")
                
            })
        return file

    }
    catch {
        return
    }

}
