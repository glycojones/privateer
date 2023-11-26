export const GITHUB_REPO = "https://github.com/glycojones/privateer"
export const GENERAL_CITATION = "https://doi.org/10.1038/nsmb.3115"
export const TORSIONS_CITATION = "http://dx.doi.org/10.1107/S2059798323003510"
export const GLYCOMICS_CITATION = "https://doi.org/10.3762/bjoc.16.204"

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

export const DatabaseColumns = [
    {
        Header: 'Chain',
        accessor: 'chain',
    },
    {
        Header: 'SNFG',
        accessor: 'SNFG',
        Cell: tableProps => {
            return (
                
                <div className="mt-4 py-4 flex justify-end" id='svgContainer' dangerouslySetInnerHTML={{
                    __html: tableProps.row.original.SNFG
                }} />
                    )
        }
    },
    {
        Header: "WURCS", 
        accessor: 'WURCS'
    }
]

export const SugarListColumns = [
    {
        Header: "Sugar ID", 
        accessor: "Sugar ID"
    }, 
    {
        Header: "Q", 
        accessor: "Q"
    }, 
    {
        Header: "Phi", 
        accessor: "Phi",
    }, 
    {
        Header: "Theta", 
        accessor: "Theta"
    }, 
    {
        Header: "RSCC", 
        accessor: "RSCC"
    }, 
    {
        Header: "B Factor", 
        accessor: "BFactor"
    }, 
    {
        Header: "Detected Type", 
        accessor: "Detected Type"
    }, 
    {
        Header: "mFo", 
        accessor: "mFo"
    }, 
    {
        Header: "Type", 
        accessor: "type"
    }
]

export const linkage_db = { 
    "GAL-1,4-NAG" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/GAL-1,4-NAG.json',
    "NAG-1,2-MAN" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/NAG-1,2-MAN.json',
    "MAN-1,6-MAN" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/MAN-1,6-MAN.json',
    "FUC-1,3-NAG" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/FUC-1,3-NAG.json',
    "MAN-1,3-BMA" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/MAN-1,3-BMA.json',
    "NAG-1,4-NAG" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/NAG-1,4-NAG.json',
    "ASN-1,2-NAG" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/ASN-1,2-NAG.json',
    "MAN-1,6-BMA" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/MAN-1,6-BMA.json',
    "BMA-1,4-NAG" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/BMA-1,4-NAG.json',
    "FUC-1,6-NAG" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/FUC-1,6-NAG.json',
    "MAN-1,3-MAN" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/MAN-1,3-MAN.json',
    "MAN-1,2-MAN" : 'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/MAN-1,2-MAN.json',
  }

export const bin_db = {
    "ASN-1,2-NAG" : {start: 0, end: 360,size: 4},
    "NAG-1,2-MAN" : {start: -180, end: 180, size: 4},
    "MAN-1,6-MAN" : {start: -180, end: 180, size: 4},
    "FUC-1,3-NAG" : {start: -180, end: 180, size: 4},
    "MAN-1,3-BMA" : {start: -180, end: 180, size: 4},
    "NAG-1,4-NAG" : {start: -180, end: 180, size: 4},
    "GAL-1,4-NAG" : {start: -180, end: 180, size: 4},
    "MAN-1,6-BMA" : {start: -180, end: 180, size: 4},
    "BMA-1,4-NAG" : {start: -180, end: 180, size: 4},
    "FUC-1,6-NAG" : {start: -180, end: 180, size: 4},
    "MAN-1,3-MAN" : {start: -180, end: 180, size: 4},
    "MAN-1,2-MAN" : {start: -180, end: 180, size: 4},
  }
