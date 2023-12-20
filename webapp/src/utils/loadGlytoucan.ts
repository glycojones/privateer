import { type TableDataEntry } from '../interfaces/types';
import pako from 'pako';

async function getGlytoucanId(wurcs: string): Promise<any> {
    const url =
        'https://api.glycosmos.org/sparqlist/wurcs2gtcids?wurcs=' +
        encodeURIComponent(wurcs);

    return await fetch(url, {
        method: 'GET', // default, so we can ignore
    })
        .then(async (response) => await response.json())
        .catch((error) => {
            console.log(error);
        });
}

export default async function loadGlytoucan(
    tableData: TableDataEntry[]
): Promise<TableDataEntry[]> {
    const promises: Array<Promise<any>> = [];

    for (let i = 0; i < tableData.length; i++) {
        const item = tableData[i];
        promises.push(getGlytoucanId(item.wurcs));
    }

    const data = await Promise.all(promises);

    data.forEach((data, index) => {
        tableData[index].glytoucan_id = data[0].id;
    });
    return tableData;
}

export async function loadGlytoucanFromFile(
    tableData: TableDataEntry[]
): Promise<TableDataEntry[]> {
    const response = await fetch(
        'privateer_glycomics_database_slim.json.gzip',
        {
            headers: new Headers({ 'content-type': 'application/gzip' }),
            mode: 'no-cors',
        }
    );

    if (!response.ok) {
        throw new Error(`Failed to fetch the file. Status: ${response.status}`);
    }
    const gzippedData = await response.arrayBuffer();
    const output = pako.inflate(gzippedData, { to: 'string' });
    const glycomicsData = JSON.parse(output);

    tableData.forEach((data, index) => {
        if (!glycomicsData.hasOwnProperty(data.wurcs)) { 
            return
        }
        const glycomicsResult = glycomicsData[data.wurcs];

        // Neaten up NotFound -> Not Found
        if (glycomicsResult.GlyConnect === 'NotFound') {
            glycomicsResult.GlyConnect = 'Not Found';
        }

        tableData[index].glytoucan_id = glycomicsResult.GlyToucan;
        tableData[index].glyconnect_id = glycomicsResult.GlyConnect;
    });
    return tableData;
}
