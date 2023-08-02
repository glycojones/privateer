export const GITHUB_REPO = "https://github.com/glycojones/privateer"
export const GENERAL_CITATION = ""

export const COLUMNS = [
    {
        Header: 'Chain ID',
        accessor: 'chain',
    },
    {
        Header: 'SNFG',
        accessor: 'svg',
        Cell: tableProps => {
            return <img src={`data:image/svg+xml;base64,${btoa(tableProps.row.original.svg)}`} alt="" width="300" height="300" />
        } 
    },
    {
        Header: 'GlytoucanID',
        accessor: 'glytoucan_id',
    },
    // {
    //     Header: 'GlyConnect ID',
    //     accessor: 'glyconnect_id',
    // },
    // {
    //     Header: 'WURCS',
    //     accessor: 'wurcs',
    // },
    
];