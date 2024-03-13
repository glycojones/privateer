import React from 'react';

export const GITHUB_REPO = 'https://github.com/glycojones/privateer';
export const GENERAL_CITATION = 'https://doi.org/10.1038/nsmb.3115';
export const TORSIONS_CITATION = 'http://dx.doi.org/10.1107/S2059798323003510';
export const GLYCOMICS_CITATION = 'https://doi.org/10.3762/bjoc.16.204';

export const GlycanDetailLayout = {
    xl: [
        { i: 'info', x: 0, y: 0, w: 1, h: 0.5, isResizable: false },
        { i: 'snfg', x: 1, y: 0, w: 1, h: 1, isResizable: false },
        { i: 'moorhen', x: 2, y: 0, w: 1, h: 2, isResizable: false },
        { i: 'torsions', x: 3, y: 0, w: 1, h: 1.75, isResizable: false },
    ],
    lg: [
        { i: 'info', x: 0, y: 0, w: 1, h: 1, isResizable: false },
        { i: 'snfg', x: 1, y: 2, w: 1, h: 1, isResizable: false },
        { i: 'moorhen', x: 0, y: 2, w: 1, h: 2, isResizable: false },
        { i: 'torsions', x: 1, y: 2, w: 1, h: 2, isResizable: false },
    ],
    md: [
        { i: 'info', x: 0, y: 0, w: 1, h: 1, isResizable: false },
        { i: 'snfg', x: 1, y: 0, w: 1, h: 1, isResizable: false },
        { i: 'moorhen', x: 0, y: 2, w: 1, h: 2, isResizable: false },
        { i: 'torsions', x: 1, y: 1, w: 1, h: 2, isResizable: false },
    ],
    sm: [
        { i: 'info', x: 0, y: 0, w: 1, h: 0.75, isResizable: false },
        { i: 'snfg', x: 0, y: 2, w: 1, h: 1, isResizable: false },
        { i: 'moorhen', x: 0, y: 4, w: 1, h: 2, isResizable: false },
        { i: 'torsions', x: 0, y: 6, w: 1, h: 1.75, isResizable: false },
    ],
};

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
        Cell: (tableProps) => {
            return (
                <div
                    className="mt-4 py-4 flex justify-end w-24 sm:w-64 md:w-96"
                    id="svgContainer"
                    dangerouslySetInnerHTML={{
                        __html: tableProps.row.original.svg,
                    }}
                />
            );
            // return <img src={`data:image/svg+xml;base64,${btoa(tableProps.row.original.svg)}`} alt="" width="300"
            //             height="300"/>
        },
    },
    {
        Header: 'Glytoucan ID',
        accessor: 'glytoucan_id',
    },
    {
        Header: 'GlyConnect ID',
        accessor: 'glyconnect_id',
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
        accessor: 'snfg',
        Cell: (tableProps) => {
            // console.log("SNFG", tableProps.row.original.SNFG)

            return (
                <div
                    className="mt-4 py-4 flex justify-end"
                    id="svgContainer"
                    dangerouslySetInnerHTML={{
                        __html: tableProps.row.original.SNFG,
                    }}
                />
            );
        },
    },
    {
        Header: 'WURCS',
        accessor: 'wurcs',
        Cell: (tableProps) => {
            return (
                <button
                    className="justify-center"
                    title="Copy WURCS to clipboard"
                    onClick={() => {
                        navigator.clipboard
                            .writeText(tableProps.row.original.WURCS as string)
                            .then(
                                () => {},
                                () => {}
                            );
                    }}
                >
                    {/* <svg className="h-6 w-6 mx-auto my-auto" xmlns="http://www.w3.org/2000/svg" height="16" width="16" viewBox="0 0 448 512">
                    <path d="M208 0H332.1c12.7 0 24.9 5.1 33.9 14.1l67.9 67.9c9 9 14.1 21.2 14.1 33.9V336c0 26.5-21.5 48-48 48H208c-26.5 0-48-21.5-48-48V48c0-26.5 21.5-48 48-48zM48 128h80v64H64V448H256V416h64v48c0 26.5-21.5 48-48 48H48c-26.5 0-48-21.5-48-48V176c0-26.5 21.5-48 48-48z" />
                </svg> */}
                    <svg
                        className="h-6 w-6 mx-auto my-auto"
                        xmlns="http://www.w3.org/2000/svg"
                        height="16"
                        width="14"
                        viewBox="0 0 448 512"
                    >
                        <path d="M384 336H192c-8.8 0-16-7.2-16-16V64c0-8.8 7.2-16 16-16l140.1 0L400 115.9V320c0 8.8-7.2 16-16 16zM192 384H384c35.3 0 64-28.7 64-64V115.9c0-12.7-5.1-24.9-14.1-33.9L366.1 14.1c-9-9-21.2-14.1-33.9-14.1H192c-35.3 0-64 28.7-64 64V320c0 35.3 28.7 64 64 64zM64 128c-35.3 0-64 28.7-64 64V448c0 35.3 28.7 64 64 64H256c35.3 0 64-28.7 64-64V416H272v32c0 8.8-7.2 16-16 16H64c-8.8 0-16-7.2-16-16V192c0-8.8 7.2-16 16-16H96V128H64z" />
                    </svg>
                </button>
            );
        },
    },
];


function getColour(diagnostic) {
    return diagnostic !== "yes" ? "text-red-600" : (diagnostic === "check" ? "text-yellow-600" : "")
}


function extracted(props, accessor) {
    const typeValue = props.row.original[accessor];
    const diagnostic = props.row.original.diagnostic;

    const colour = getColour(diagnostic)
    return (
        <span
            className={colour}>
            {typeValue}
        </span>
    );
}

export const SugarListColumns = [
    {
        Header: 'Sugar ID',
        accessor: 'sugarId',
        Cell: (props) => extracted(props, 'sugarId'),
    },
    {
        Header: 'Conformation',
        accessor: 'conformation',
        Cell: (props: {
            row: { original: { conformation: any; diagnostic: any } };
        }) => {
            const typeValue = props.row.original.conformation;
            const diagnostic = props.row.original.diagnostic;
            const colour = getColour(diagnostic)

            const regex = /([a-zA-Z]*?\d*)([a-zA-Z])(\d*)/;
            const formattedString = typeValue.replace(
                regex,
                (_, before, letter, after) => {
                    let string = '';
                    if (before) {
                        string += '<sup>' + before + '</sup>';
                    }
                    string += letter;
                    if (after) {
                        string += '<sub>' + after + '</sub>';
                    }
                    return string;
                }
            );
            return (
                <div
                    className={colour}
                    dangerouslySetInnerHTML={{ __html: formattedString }}
                ></div>
            );
        },
    },
    {
        Header: 'Q',
        accessor: 'q',
        Cell: (props) => extracted(props, 'q'),
    },
    {
        Header: 'Phi',
        accessor: 'phi',
        Cell: (props) => extracted(props, 'phi'),
    },
    {
        Header: 'Theta',
        accessor: 'theta',
        Cell: (props) => extracted(props, 'theta'),
    },

    {
        Header: 'RSCC',
        accessor: 'rscc',
        Cell: (props) => extracted(props, 'rscc'),
    },
    {
        Header: 'B Factor',
        accessor: 'bFactor',
        Cell: (props) => extracted(props, 'bFactor'),
    },
    {
        Header: 'Detected Type',
        accessor: 'detectedType',
        Cell: (props) => extracted(props, 'detectedType'),
    },
    // {
    //     Header: 'mFo',
    //     accessor: 'mFo',
    // },
    {
        Header: 'Type',
        accessor: 'type',
        Cell: (props) => extracted(props, 'type'),
    },
    {
        Header: 'Diagnostic',
        accessor: 'diagnostic',
        Cell: (props) => extracted(props, 'diagnostic'),
    },
];

export const linkageDB: Record<string, string> = {
    'GAL-1,4-NAG':
        'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/GAL-1,4-NAG.json',
    'NAG-1,2-MAN':
        'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/NAG-1,2-MAN.json',
    'MAN-1,6-MAN':
        'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/MAN-1,6-MAN.json',
    'FUC-1,3-NAG':
        'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/FUC-1,3-NAG.json',
    'MAN-1,3-BMA':
        'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/MAN-1,3-BMA.json',
    'NAG-1,4-NAG':
        'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/NAG-1,4-NAG.json',
    'ASN-1,2-NAG':
        'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/ASN-1,2-NAG.json',
    'MAN-1,6-BMA':
        'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/MAN-1,6-BMA.json',
    'BMA-1,4-NAG':
        'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/BMA-1,4-NAG.json',
    'FUC-1,6-NAG':
        'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/FUC-1,6-NAG.json',
    'MAN-1,3-MAN':
        'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/MAN-1,3-MAN.json',
    'MAN-1,2-MAN':
        'https://raw.githubusercontent.com/Dialpuri/N-glycanTorsionDatabase/main/data_json/MAN-1,2-MAN.json',
};

export const binDB = {
    'ASN-1,2-NAG': { start: 0, end: 360, size: 4 },
    'NAG-1,2-MAN': { start: -180, end: 180, size: 4 },
    'MAN-1,6-MAN': { start: -180, end: 180, size: 4 },
    'FUC-1,3-NAG': { start: -180, end: 180, size: 4 },
    'MAN-1,3-BMA': { start: -180, end: 180, size: 4 },
    'NAG-1,4-NAG': { start: -180, end: 180, size: 4 },
    'GAL-1,4-NAG': { start: -180, end: 180, size: 4 },
    'MAN-1,6-BMA': { start: -180, end: 180, size: 4 },
    'BMA-1,4-NAG': { start: -180, end: 180, size: 4 },
    'FUC-1,6-NAG': { start: -180, end: 180, size: 4 },
    'MAN-1,3-MAN': { start: -180, end: 180, size: 4 },
    'MAN-1,2-MAN': { start: -180, end: 180, size: 4 },
};

export const sugarLinkageMap: Record<string, string[]> = {
    'N-glycans': [
        'Any',
        'NAG-1,2-ASN',
        'NAG-1,4-NAG',
        'BMA-1,4-NAG',
        'MAN-1,3-BMA',
        'MAN-1,6-BMA',
        'FUC-1,6-NAG',
        'MAN-1,3-MAN',
        'MAN-1,6-MAN',
        'MAN-1,2-MAN',
        'MAN-1,4-NAG',
        'NAG-1,2-MAN',
        'FUC-1,3-NAG',
        'BMA-1,3-BMA',
        'GAL-1,4-NAG',
        'FUL-1,6-NAG',
        'BMA-1,6-BMA',
        'NDG-1,4-NAG',
        'NAG-1,1-ASN',
        'XYP-1,2-BMA',
        'NAG-1,3-NAG',
        'MAN-1,4-BMA',
        'BMA-1,6-MAN',
        'MAN-1,4-MAN',
        'BMA-1,3-NAG',
        'FUL-1,3-NAG',
        'BMA-1,4-NDG',
        'BMA-1,4-MAN',
        'BMA-1,3-MAN',
        'NAG-1,6-NAG',
        'BMA-1,4-BMA',
        'NAG-1,4-MAN',
        'MAN-1,3-NAG',
        'MAN-1,4-NDG',
        'SIA-2,6-GAL',
        'MAN-1,6-NAG',
        'BMA-1,2-MAN',
        'NDG-1,2-ASN',
        'NAG-1,2-BMA',
        'MAN-1,2-BMA',
        'NAG-1,4-BMA',
        'FCA-1,3-NAG',
        'NAG-1,2-ARG',
        'GLC-1,3-MAN',
        'NAG-1,2-LYS',
        'NAG-1,4-NDG',
        'NDG-1,2-MAN',
        'NAG-1,6-MAN',
        'MAN-2,6-BMA',
        'NAG-3,4-NAG',
        'GAL-1,4-NDG',
        'FUC-1,4-NAG',
        'FCA-1,6-NAG',
        'NAG-5,4-NAG',
        'BMA-1,6-NAG',
        'BMA-1,2-BMA',
        'MAN-1,2-ASN',
        'BGC-1,2-ASN',
        'XYS-1,2-BMA',
        'SIA-2,3-GAL',
        'NAG-1,3-MAN',
        'NAG-4,4-NAG',
        'GLA-1,4-NAG',
        'GAL-1,4-MAN',
        'NDG-1,4-NDG',
        'NDG-1,2-BMA',
        'HSQ-1,2-ASN',
        'BMA-1,2-ASN',
        'MAN-2,3-BMA',
        'NAG-1,G-ASN',
        'MAN-1,3-XXR',
        'KDO-2,1-ASN',
        'RM4-1,4-XYP',
        'GLC-1,4-GLC',
        'FUL-1,6-NDG',
        'NAG-1,3-BMA',
        'M6D-1,2-MAN',
        'MAN-3,4-NAG',
        'NDG-1,1-ASN',
        'GLA-1,2-FUC',
        'MAN-3,3-MAN',
        'NAG-8,6-NAG',
        'B6D-1,2-ASN',
        'TGY-2,6-NAG',
        'NAG-1,B-ASN',
        'GAL-1,4-FUC',
        'XYP-1,4-FUC',
        'GLC-1,4-NAG',
        'XYP-1,4-BGC',
        'FUC-5,6-NAG',
        'NAG-1,1-MAN',
        'XXR-1,3-FUC',
        'FUC-1,3-BGC',
        'XYS-1,2-MAN',
        'NAG-1,Z-LYS',
        'NGK-1,4-NAG',
        'NAG-7,8-NAG',
        'FUC-1,2-GAL',
        'FRU-2,2-LYS',
        'NAG-4,1-NAG',
        'BMA-4,4-NAG',
        '7CV-1,2-RM4',
        'XYL-1,2-BMA',
        'BGC-1,4-NAG',
        'NAG-1,4-ARG',
        'BMA-1,3-SHD',
        'GMH-1,5-KDO',
        'XYP-1,2-MAN',
        'GL0-1,4-NGA',
        'MAN-3,3-BMA',
        'NGZ-1,3-B6D',
        'NAG-A,A-ASN',
        'GLA-1,6-BMA',
        'MAN-1,3-XYS',
        'GUP-1,4-NAG',
        'MAN-1,3-GLC',
        'NAG-1,4-GAL',
        'BMA-1,6-SHD',
        'NAA-1,3-NAG',
        'MAN-1,1-ASN',
        'BMA-5,4-NAG',
        'G6D-1,6-BGC',
        'FUC-1,6-BMA',
        'NDG-1,4-MAN',
        'Z9N-1,6-NAG',
        'GLC-1,2-ASN',
        'NAA-1,2-ASN',
        'MAN-1,1-ARG',
        'GUP-1,3-BMA',
        'BMA-1,6-GLC',
        'BMA-1,3-NDG',
        'MAN-3,6-BMA',
        'MAN-1,3-ARG',
        'GCS-1,1-ARG',
        'MAN-1,4-ARG',
        'LXZ-1,2-ASN',
        'MAN-2,3-NAG',
        'JHM-1,4-IDS',
        'NGZ-1,4-LXB',
        'NAG-6,6-NAG',
        'SHD-1,4-NAG',
        'GLC-1,3-GLC',
        'BGC-1,3-A2G',
        'GXL-1,6-LXB',
        'RIB-1,1-ARG',
        'NAG-2,4-NAG',
        'NDG-1,6-BMA',
        'BGC-1,6-BGC',
        'BMA-5,6-MAN',
        'NAG-3,3-NAG',
        'MAN-1,3-GAL',
        'MAN-1,6-GUP',
        'SGN-1,4-IDS',
        'GLC-1,6-GL0',
        'JHM-1,1-LYS',
        'BGC-1,2-BGC',
        'FUC-1,6-NDG',
        'GLC-1,4-NDG',
        'XYS-1,6-BMA',
        'NAG-1,1-ARG',
        'BGC-1,3-LXB',
        'IDS-1,4-JHM',
        'IDS-4,1-SGN',
        'XYZ-1,6-MAN',
        'GAL-1,3-GAL',
        'RIP-1,6-NAG',
        'BGC-1,3-LXZ',
        'MAN-1,3-LGU',
        'BGC-1,3-BGC',
        'A2G-1,4-A2G',
        'XYS-1,3-BMA',
        'NGA-1,4-LXZ',
        'NAG-1,2-GUP',
        'GLC-4,1-GLC',
        'XYP-1,2-ASN',
        'LGU-1,4-NAG',
        'SIA-2,6-GLA',
        'MAN-4,6-MAN',
        'IDS-1,4-SGN',
        'MAN-1,4-BGC',
        'BGC-1,1-LYS',
        'NAG-2,3-NAG',
        'MAN-1,3-GUP',
        'BDF-2,2-LYS',
        'A2G-1,3-B6D',
        'GAL-1,2-ASN',
        'GAL-1,3-MAN',
        'KDO-2,6-GCS',
        'LGU-1,6-LGU',
        'GUP-1,2-BMA',
        'NAG-4,3-MAN',
        'MAN-4,2-MAN',
        'GAL-1,6-NAG',
        'NAG-7,7-NAG',
        'FUC-4,3-NAG',
        'GLA-1,3-GAL',
        'SGN-4,1-IDS',
        'GLC-1,4-ASN',
        'MAN-1,Z-LYS',
        'GLC-1,6-LXZ',
        'NAG-1,4-HSQ',
        'BMA-6,3-BMA',
        'BMA-2,4-NAG',
        'SGN-1,4-ARG',
        'MAN-1,3-NDG',
        'MAN-1,3-ASN',
        'GL0-1,4-NGZ',
        'BGC-1,2-BMA',
        'BGC-1,3-GL0',
        'IDS-1,4-ARG',
        'BGS-1,1-LYS',
        'NDG-8,3-NDG',
        'LXB-1,2-ASN',
        'BGC-1,4-BGC',
        'GUP-1,2-MAN',
        'BGC-1,3-BMA',
        'KDO-2,4-KDO',
        'BMA-1,6-ARG',
    ],
    'C-glycans': [
        'Any',
        'MAN-1,1-TRP',
        'BMA-1,1-TRP',
        'GAL-1,1-TRP',
        'NAG-1,4-TRP',
        'BGC-1,1-TRP',
        'NAG-1,1-TRP',
        'NAG-1,2-TRP',
        'NAG-4,1-NAG',
    ],
    'O-glycans': [
        'Any',
        'MAN-1,G-SER',
        'MAN-1,1-THR',
        'A2G-1,1-THR',
        'BGC-1,G-SER',
        'FUC-1,1-THR',
        'NAG-1,G-SER',
        'FUC-1,G-SER',
        'NAG-1,1-THR',
        'A2G-1,G-SER',
        'GAL-1,3-A2G',
        'BGC-1,3-FUC',
        'NGA-1,1-THR',
        'GLC-1,G-SER',
        'RAM-1,4-MAN',
        'G2F-1,2-GLU',
        'NAG-1,2-GLU',
        'BGC-1,4-BGC',
        'MAN-1,2-MAN',
        'G2F-1,1-GLU',
        'GCU-1,2-MAN',
        'NAG-1,2-ASP',
        'XYS-1,3-BGC',
        'NGA-1,G-SER',
        'XYP-1,4-GCU',
        'BGC-1,4-G2F',
        'BGC-1,2-ASP',
        'NDG-1,1-THR',
        'RIP-1,G-SER',
        'BMA-1,4-NAG',
        'GLC-1,4-GLC',
        'NAG-1,2-THR',
        'SIA-2,6-A2G',
        'B9D-1,1-ASP',
        'NAG-1,2-SER',
        'XYS-1,3-XYS',
        'MXY-1,4-XYP',
        'NAG-1,2-TYR',
        'BMA-1,G-SER',
        'GAL-1,3-NGA',
        'NAG-1,4-NAG',
        'SIA-2,3-GAL',
        'MAN-A,A-SER',
        'DFX-1,2-GLU',
        'BMA-1,1-THR',
        '2DG-1,2-GLU',
        'BGC-1,3-BGC',
        'XYP-1,4-DFX',
        'BGC-1,H-TYR',
        'GLC-1,2-GLU',
        'XYP-1,4-X2F',
        'BGC-1,4-GLC',
        'GLC-1,1-GLU',
        'GLC-1,1-ASP',
        'GLC-1,4-BGC',
        'SIA-2,6-NDG',
        'X2F-1,2-GLU',
        'MAN-A,A-THR',
        'AC1-1,4-GLC',
        'NBG-1,1-GLU',
        'XYS-1,6-BGC',
        'GLA-1,3-DT6',
        'MAN-A,G-SER',
        'FUL-1,G-SER',
        'GLA-1,3-MXZ',
        'MXZ-1,4-XYP',
        'FUL-1,1-THR',
        'AC1-1,4-ASO',
        'NAG-1,4-ASP',
        'BGC-1,4-MXZ',
        'MAN-1,A-SER',
        'NAG-1,6-A2G',
        '7JZ-1,2-ASP',
        'NAG-1,3-FUC',
        '2FG-1,2-GLU',
        'XYS-1,G-SER',
        'ASO-1,2-ASP',
        'BGC-1,4-MXY',
        'GLC-1,2-ASP',
        'BGC-1,2-B9D',
        'MAN-A,1-THR',
        'GLC-4,1-GLC',
        'RIP-1,1-THR',
        'GAL-1,3-NDG',
        'EPG-1,1-GLU',
        'G4D-1,3-MXY',
        'GLA-1,3-NAG',
        'X2F-1,1-GLU',
        'XYS-1,1-THR',
        'BGC-1,4-NBG',
        'BGC-1,3-G2F',
        'DT6-1,G-SER',
        'G2F-1,1-ASP',
        'GLC-1,H-TYR',
        'XYS-1,2-GLU',
        'XYS-1,4-GCU',
        'NAG-1,H-TYR',
        'XYF-1,5-ASP',
        'NGA-1,3-NGA',
        'MFU-1,4-XYP',
        'GAL-1,1-THR',
        'XYP-1,3-BXF',
        'NAG-1,4-G2F',
        'GLC-1,4-THR',
        'XYP-1,3-XYP',
        'GAL-1,H-TYR',
        'GAF-1,2-GLU',
        'BDP-1,1-SER',
        'G2F-1,2-ASP',
        'GAL-1,4-GLC',
        'MAN-1,6-BMA',
        'XYP-2,2-XYS',
        'GLC-1,6-BGC',
        'GAL-1,3-NAG',
        'GLC-1,4-ASP',
        '8B9-1,1-THR',
        'EBG-1,1-GLU',
        'GLC-1,4-SHG',
        'NGA-1,H-TYR',
        'NAG-1,2-MAN',
        'SHG-4,1-BGC',
        'FUC-1,3-NAG',
        'SHG-1,1-ASP',
        'BMA-1,4-GLU',
        'A2G-1,1-GLU',
        'MAN-1,1-GLU',
        'BGC-1,1-SER',
        'BMA-2,1-MAN',
        'BXF-1,1-GLU',
        'XYS-1,1-GLU',
        'RAM-1,1-THR',
        'RIB-1,2-GLU',
        'BGC-1,3-NBG',
        'XYS-1,4-XYS',
        '3HD-1,1-THR',
        'NAG-1,4-MUB',
        'GLA-1,3-B6D',
        'XYP-2,1-GCV',
        'GLA-1,G-SER',
        'XYP-1,3-BGC',
        'AHR-1,1-GLU',
        'GCV-1,2-SER',
        '289-1,G-SER',
        'MAN-1,4-MAN',
        'AMP-1,9-ASP',
        'B8D-1,4-BGC',
        'SIA-2,6-8B9',
        'NAG-1,4-GLU',
        'BGC-4,1-BGC',
        'NAG-1,6-NGA',
        'GCU-1,1-SER',
        'GXL-1,4-NDG',
        'C3X-1,1-GLU',
        'SIA-2,1-THR',
        'XYL-4,1-XYP',
        'BGC-1,3-GLC',
        'GAL-1,4-NAG',
        'GXL-1,2-GAL',
        'MAN-1,3-GLU',
        'MAN-1,2-ASP',
        'XYP-1,4-XYS',
        'C5X-1,1-GLU',
        'ARB-1,2-MAN',
        'NGA-1,1-SER',
        'SHG-1,G-SER',
        'G4D-1,4-GLC',
        'BGC-1,4-GLU',
        'NDG-1,6-A2G',
        'ARA-1,1-THR',
        'GLC-1,2-FRU',
        'BGC-1,1-THR',
        'NAG-1,3-A2G',
        'GAL-1,4-BGC',
        'DFX-1,1-GLU',
        'B6D-1,G-SER',
        'NDG-1,H-TYR',
        'FRU-2,2-GLU',
        'XYP-4,1-XYP',
        'GLC-1,1-THR',
        'GAL-1,1-GLU',
        'MAN-1,2-BMA',
        'MUB-1,2-GLU',
    ],
    'S-glycans': [
        'Any',
        'NAG-1,G-CYS',
        'KDO-2,G-CYS',
        'BGC-1,G-CYS',
        'MAN-1,G-CYS',
        'A2G-1,G-CYS',
        'MAN-1,2-MAN',
    ],
    Ligands: [
        'Any',
        'NAG-1,2-ASN',
        'NAG-1,4-NAG',
        'NAG-1,4-BMA',
        'NAG-1,4-NAG',
        'NAG-1,4-NAG',
    ],
};
