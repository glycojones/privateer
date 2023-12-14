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
        { i: 'info', x: 0, y: 0, w: 1, h: 0.5, isResizable: false },
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
        accessor: 'WURCS',
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

export const SugarListColumns = [
    {
        Header: 'Sugar ID',
        accessor: 'Sugar ID',
    },
    {
        Header: 'Q',
        accessor: 'Q',
    },
    {
        Header: 'Phi',
        accessor: 'Phi',
    },
    {
        Header: 'Theta',
        accessor: 'Theta',
    },
    {
        Header: 'RSCC',
        accessor: 'RSCC',
    },
    {
        Header: 'B Factor',
        accessor: 'BFactor',
    },
    {
        Header: 'Detected Type',
        accessor: 'Detected Type',
    },
    {
        Header: 'mFo',
        accessor: 'mFo',
    },
    {
        Header: 'Type',
        accessor: 'type',
    },
    {
        Header: 'Diagnostic',
        accessor: 'Diagnostic',
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
