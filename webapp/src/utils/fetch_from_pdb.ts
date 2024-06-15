import { detect } from 'detect-browser';

async function fetchPDBFile(PDBCode: string): Promise<string | void> {
    console.warn(
        'The CIF file for this PDB could not be found, trying for the PDB'
    );
    const pdbURL = `https://files.rcsb.org/download/${PDBCode.toUpperCase()}.pdb`;

    const file = fetch(pdbURL)
        .then(async (response) => {
            if (!response.ok) {
                throw new Error('Network error');
            }
            return await response.text();
        })
        .then(async (file) => {
            return file;
        })
        .catch(async (error) => {
            return await Promise.reject(error);
        });
    return await file;
}

export async function fetchPDB(PDBCode: string): Promise<string | void> {
    if (PDBCode == null) {
        return;
    }
    // first try fetching the cif
    console.log('Fetching PDB ', PDBCode);
    const pdbURL = `https://files.rcsb.org/download/${PDBCode.toUpperCase()}.cif`;

    // FIXME
    const browser = detect(); // FireFox doesn't work with CIF files, get the PDB.
    if (browser.name === 'firefox') {
        try {
            return await fetchPDBFile(PDBCode);
        } catch (e) {
            return await Promise.reject(e);
        }
    }

    const file = fetch(pdbURL)
        .then(async (response) => {
            if (!response.ok) {
                throw new Error('Network error');
            }
            return await response.text();
        })
        .then(async (file) => {
            return file;
        })
        .catch(async () => {
            // if we can't find the cif, try the PDB, if that fails, then report the failure.
            try {
                return await fetchPDBFile(PDBCode);
            } catch (e) {
                return await Promise.reject(e);
            }
        });
    return await file;
}

export async function fetchMap(PDBCode: string): Promise<ArrayBuffer | void> {
    if (PDBCode == null) {
        await Promise.resolve();
        return;
    }
    console.log('Fetching MTZ ', PDBCode);
    const mtzURL = `https://www.ebi.ac.uk/pdbe/entry-files/${PDBCode.toLowerCase()}.ccp4`;

    try {
        const controller = new AbortController();
        const timeId = setTimeout(() => {
            controller.abort();
        }, 60000);

        const file = fetch(mtzURL, { signal: controller.signal })
            .then(async (response) => {
                clearTimeout(timeId);

                if (!response.ok) {
                    throw new Error('Network error');
                }
                return await response.arrayBuffer();
            })
            .then(async (file) => {
                return file;
            })
            .catch(async () => {
                // throw 'Map not found'
                throw new Error('Map Not Found');
            });
        return await file;
    } catch (e) {
        return await Promise.reject(e);
    }
}
