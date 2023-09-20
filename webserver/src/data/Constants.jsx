export const GITHUB_REPO = "https://github.com/glycojones/privateer"
export const GENERAL_CITATION = ""

export const COLUMNS = [
    {
        Header: 'Chain',
        accessor: 'chain',
    },
    {
        Header: 'ID',
        accessor: 'id',
    },
    {
        Header: 'SNFG',
        accessor: 'svg',
        Cell: tableProps => {
            return <img src={`data:image/svg+xml;base64,${btoa(tableProps.row.original.svg)}`} alt="" width="300"
                        height="300"/>
        }
    },
    {
        Header: 'GlytoucanID',
        accessor: 'glytoucan_id',
    },
    // {
    //     Header: "Explore",
    //     Cell: tableProps => {
    //         return <button onClick={tableProps.row.onClick}>Visualise</button>
    //     }
    // }
    // {
    //     Header: 'GlyConnect ID',
    //     accessor: 'glyconnect_id',
    // },
    // {
    //     Header: 'WURCS',
    //     accessor: 'wurcs',
    // },

];

export const linkage_db = { 
    "BMA-1,6-MAN" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/BMA-1,6-MAN.json',
    "ASN-1,2-NAG" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/ASN-1,1-NAG.json',
    "MAN-1,6-MAN" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/MAN-1,6-MAN.json',
    "NAG-1,3-FUC" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/NAG-1,3-FUC.json',
    "NAG-1,4-NAG" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/NAG-1,4-NAG.json',
    "NAG-1,6-FUC" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/NAG-1,6-FUC.json',
    "NAG-1,4-BMA" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/NAG-1,4-BMA.json',
    "BMA-1,3-MAN" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/BMA-1,3-MAN.json',
    "NAG-1,4-GAL" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/NAG-1,4-GAL.json',
    "MAN-1,3-MAN" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/MAN-1,3-MAN.json',
    "MAN-1,2-MAN" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/MAN-1,2-MAN.json',
    "MAN-1,2-NAG" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/MAN-1,2-NAG.json'
  }

  export const bin_db = {
    "ASN-1,2-NAG": {start: 0, end: 360,size: 4},
    "BMA-1,6-MAN" : {start: -180, end: 180, size: 4},
    "BMA-1,6-MAN" : {start: -180, end: 180, size: 4},
    "MAN-1,6-MAN" : {start: -180, end: 180, size: 4},
    "NAG-1,3-FUC" : {start: -180, end: 180, size: 4},
    "NAG-1,4-NAG" : {start: -180, end: 180, size: 4},
    "NAG-1,6-FUC" : {start: -180, end: 180, size: 4},
    "NAG-1,4-BMA" : {start: -180, end: 180, size: 4},
    "BMA-1,3-MAN" : {start: -180, end: 180, size: 4},
    "NAG-1,4-GAL" : {start: -180, end: 180, size: 4},
    "MAN-1,3-MAN" : {start: -180, end: 180, size: 4},
    "MAN-1,2-MAN" : {start: -180, end: 180, size: 4},
    "MAN-1,2-NAG" : {start: -180, end: 180, size: 4},
  }
