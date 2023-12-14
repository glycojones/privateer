import { type TableDataEntry } from '../interfaces/types';

async function getGlytoucanId (wurcs: string): Promise<any> {
  const url =
    'https://api.glycosmos.org/sparqlist/wurcs2gtcids?wurcs=' +
    encodeURIComponent(wurcs);

  return await fetch(url, {
    method: 'GET' // default, so we can ignore
  })
    .then(async (response) => await response.json())
    .catch((error) => {
      console.log(error);
    });
}

export default async function loadGlytoucan (
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
