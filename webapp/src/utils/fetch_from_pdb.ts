export async function fetchPDB (PDBCode: string): Promise<string | void> {
  if (PDBCode == null) {
    return;
  }
  console.log('Fetching PDB ', PDBCode);
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
    .catch(() => {
      throw new Error('PDB Not Found');
    });
  return await file;
}

export async function fetchMap (PDBCode: string): Promise<ArrayBuffer | void> {
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
  } catch {}
}
