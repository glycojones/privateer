import React, {
    type Dispatch,
    type SetStateAction,
    useEffect,
    useState,
} from 'react';
import DatabaseSearchTable from '../DatabaseSearchTable/DatabaseSearchTable.tsx';
import Loading from '../../shared/Loading/Loading.tsx';

function ViewAllEntriesButton(props: {
    text: string;
    label: any;
    onClickMethod: any;
    secondLabel: string | undefined;
    onSecondClick: any;
}) {
    return (
        <div
            id="submit"
            className="flex flex-col m-12 px-12 pt-8 w-64
            items-center text-center justify-between h-64 border-2
            transition-all border-gray-300 rounded-lg
            bg-gray-50 flex-grow-0 border-dashed"
        >
            <svg
                xmlns="http://www.w3.org/2000/svg"
                className="w-6 h-6 mb-4 text-gray-500 dark:text-gray-400"
                viewBox="0 0 576 512"
            >
                <path d="M64 32C64 14.3 49.7 0 32 0S0 14.3 0 32v96V384c0 35.3 28.7 64 64 64H256V384H64V160H256V96H64V32zM288 192c0 17.7 14.3 32 32 32H544c17.7 0 32-14.3 32-32V64c0-17.7-14.3-32-32-32H445.3c-8.5 0-16.6-3.4-22.6-9.4L409.4 9.4c-6-6-14.1-9.4-22.6-9.4H320c-17.7 0-32 14.3-32 32V192zm0 288c0 17.7 14.3 32 32 32H544c17.7 0 32-14.3 32-32V352c0-17.7-14.3-32-32-32H445.3c-8.5 0-16.6-3.4-22.6-9.4l-13.3-13.3c-6-6-14.1-9.4-22.6-9.4H320c-17.7 0-32 14.3-32 32V480z" />
            </svg>
            <p className="mb-2 text-md text-gray-500 dark:text-gray-400">
                <span className="font-semibold">{props.text}</span>
            </p>
            <div className="flex w-full">
                <button
                    type="button"
                    id="fetch"
                    className="bg-gray-50 border font-bold mb-8 border-gray-300 text-gray-900 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block w-full p-2.5"
                    onClick={() => {
                        props.onClickMethod((value: any) => !value);
                    }}
                >
                    {props.label}
                </button>
                {props.secondLabel !== undefined ? (
                    <button
                        type="button"
                        id="fetch"
                        className="bg-gray-50 border ml-2 font-bold mb-8 border-gray-300 text-gray-900 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block w-full p-2.5"
                        onClick={() => {
                            props.onSecondClick((value: any) => !value);
                        }}
                    >
                        {props.secondLabel}
                    </button>
                ) : (
                    <></>
                )}
            </div>
        </div>
    );
}

function TypeFilterBox(props: { selected: string; onClickMethod: any }) {
    const sugars = ['N-glycans', 'O-glycans', 'S-glycans', 'C-glycans'];
    return (
        <div
            id="submit"
            className="flex flex-col m-12 px-2 pt-6 w-64
            items-center text-center justify-between h-64 border-2
            transition-all border-gray-300 rounded-lg
            bg-gray-50 flex-grow-0 border-dashed"
        >
            <svg
                xmlns="http://www.w3.org/2000/svg"
                viewBox="0 0 448 512"
                className="w-6 h-6 mb-2 text-gray-500 dark:text-gray-400"
            >
                <path d="M192 64v64c0 17.7 14.3 32 32 32h64c17.7 0 32-14.3 32-32V64c0-17.7-14.3-32-32-32H224c-17.7 0-32 14.3-32 32zM82.7 207c-15.3 8.8-20.5 28.4-11.7 43.7l32 55.4c8.8 15.3 28.4 20.5 43.7 11.7l55.4-32c15.3-8.8 20.5-28.4 11.7-43.7l-32-55.4c-8.8-15.3-28.4-20.5-43.7-11.7L82.7 207zM288 192c-17.7 0-32 14.3-32 32v64c0 17.7 14.3 32 32 32h64c17.7 0 32-14.3 32-32V224c0-17.7-14.3-32-32-32H288zm64 160c-17.7 0-32 14.3-32 32v64c0 17.7 14.3 32 32 32h64c17.7 0 32-14.3 32-32V384c0-17.7-14.3-32-32-32H352zM160 384v64c0 17.7 14.3 32 32 32h64c17.7 0 32-14.3 32-32V384c0-17.7-14.3-32-32-32H192c-17.7 0-32 14.3-32 32zM32 352c-17.7 0-32 14.3-32 32v64c0 17.7 14.3 32 32 32H96c17.7 0 32-14.3 32-32V384c0-17.7-14.3-32-32-32H32z" />
            </svg>

            <p className="mb-1 text-md text-gray-500 dark:text-gray-400">
                <span className="font-semibold">Filter Sugars</span>
            </p>

            <div className="w-full flex flex-wrap my-auto overflow-x-hidden justify-center align-items-center ">
                {sugars.map((item, index) => {
                    return (
                        <button
                            type="button"
                            key={index}
                            id="fetch"
                            className={
                                props.selected === item
                                    ? 'bg-gray-50 mx-2 my-2 border w-24 font-bold border-gray-300 text-gray-900 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block p-2'
                                    : 'bg-gray-50 mx-2 my-2 border w-24 border-gray-300 text-gray-900 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block p-2'
                            }
                            onClick={() => {
                                props.onClickMethod(item);
                            }}
                        >
                            {item}
                        </button>
                    );
                })}
            </div>
        </div>
    );
}

function LinkageFilterBox(props: {
    selected: string;
    sugars: string[];
    onClickMethod: Dispatch<SetStateAction<unknown>>;
}) {
    return (
        <div
            id="submit"
            className="flex flex-col m-12 px-2 pt-6 w-64
            items-center text-center justify-between h-64 border-2
            transition-all border-gray-300 rounded-lg
            bg-gray-50 flex-grow-0 border-dashed"
        >
            <svg
                xmlns="http://www.w3.org/2000/svg"
                viewBox="0 0 640 512"
                className="w-6 h-6 mb-2 text-gray-500 dark:text-gray-400"
            >
                <path d="M579.8 267.7c56.5-56.5 56.5-148 0-204.5c-50-50-128.8-56.5-186.3-15.4l-1.6 1.1c-14.4 10.3-17.7 30.3-7.4 44.6s30.3 17.7 44.6 7.4l1.6-1.1c32.1-22.9 76-19.3 103.8 8.6c31.5 31.5 31.5 82.5 0 114L422.3 334.8c-31.5 31.5-82.5 31.5-114 0c-27.9-27.9-31.5-71.8-8.6-103.8l1.1-1.6c10.3-14.4 6.9-34.4-7.4-44.6s-34.4-6.9-44.6 7.4l-1.1 1.6C206.5 251.2 213 330 263 380c56.5 56.5 148 56.5 204.5 0L579.8 267.7zM60.2 244.3c-56.5 56.5-56.5 148 0 204.5c50 50 128.8 56.5 186.3 15.4l1.6-1.1c14.4-10.3 17.7-30.3 7.4-44.6s-30.3-17.7-44.6-7.4l-1.6 1.1c-32.1 22.9-76 19.3-103.8-8.6C74 372 74 321 105.5 289.5L217.7 177.2c31.5-31.5 82.5-31.5 114 0c27.9 27.9 31.5 71.8 8.6 103.9l-1.1 1.6c-10.3 14.4-6.9 34.4 7.4 44.6s34.4 6.9 44.6-7.4l1.1-1.6C433.5 260.8 427 182 377 132c-56.5-56.5-148-56.5-204.5 0L60.2 244.3z" />
            </svg>

            <p className="mb-2 text-md text-gray-500 dark:text-gray-400">
                <span className="font-semibold">Filter Linkages</span>
            </p>

            <div className="w-full flex flex-wrap overflow-scroll my-auto overflow-x-hidden justify-center align-items-center ">
                {props.sugars.map((item, index) => {
                    return (
                        <button
                            key={index}
                            type="button"
                            id="fetch"
                            className={
                                props.selected === item
                                    ? 'bg-gray-50 mx-2 my-2 border w-32 font-bold border-gray-300 text-gray-900 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block p-2'
                                    : 'bg-gray-50 mx-2 my-2 border w-32 border-gray-300 text-gray-900 text-sm rounded-lg focus:ring-blue-500 focus:border-blue-500 block p-2'
                            }
                            onClick={() => {
                                props.onClickMethod(item);
                            }}
                        >
                            {item}
                        </button>
                    );
                })}
            </div>
        </div>
    );
}

function FilterZone(props: { setSearchBegin: any }) {
    const [search, setSearch] = useState<boolean>(false);
    const [linkage, setLinkage] = useState<string>('NAG-1,2-ASN');
    const [type, setType] = useState<string>('N-glycans');
    const [text, setText] = useState<string>('');
    const sugars = {
        'N-glycans': [
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
            'NAG-1,G-CYS',
            'KDO-2,G-CYS',
            'BGC-1,G-CYS',
            'MAN-1,G-CYS',
            'A2G-1,G-CYS',
            'MAN-1,2-MAN',
        ],
        Ligands: [
            'NAG-1,2-ASN',
            'NAG-1,4-NAG',
            'NAG-1,4-BMA',
            'NAG-1,4-NAG',
            'NAG-1,4-NAG',
        ],
    };

    useEffect(() => {
        let textString = 'Find';

        if (type === 'Any') {
            textString += ' any carbohydrate';
        } else {
            textString += ' ' + type;
        }
        textString += ' with ';
        if (linkage === 'Any') {
            textString += ' any linkage type';
        } else {
            textString += ' ' + linkage + ' linkages';
        }

        setText(textString);
    }, [linkage, type]);

    useEffect(() => {
        if (!search) {
            return;
        }

        let url =
            'https://raw.githubusercontent.com/Dialpuri/PrivateerDatabase/master/linkages/';
        if (type === 'N-glycans') {
            url += 'n-glycan/';
        }
        if (type === 'O-glycans') {
            url += 'o-glycan/';
        }
        if (type === 'S-glycans') {
            url += 's-glycan/';
        }
        if (type === 'C-glycans') {
            url += 'c-glycan/';
        }
        if (type === 'Ligands') {
            url += 'ligand/';
        }
        url += linkage.replace(',', '%2C');
        url += '.json';

        void fetch(url)
            .then(async (response) => await response.json())
            .then((json) => {
                const formattedData: Record<string, string | number> = json.map(
                    (item) => {
                        return {
                            pdb: item.pdb,
                            count: item.count,
                            resolution: item.resolution,
                            type: 'n-glycan',
                            link:
                                'https://privateer.york.ac.uk/database?pdb=' +
                                item.pdb,
                        };
                    }
                );
                setData(formattedData);
            })
            .catch(() => {});
    }, [search]);

    const [data, setData] = useState<Record<string, string | number> | null>(
        null
    );

    return !search ? (
        <>
            {ViewAllEntriesButton({
                label: 'Back',
                text,
                onClickMethod: props.setSearchBegin,
                secondLabel: 'Search',
                onSecondClick: setSearch,
            })}
            {TypeFilterBox({ onClickMethod: setType, selected: type })}
            {LinkageFilterBox({
                onClickMethod: setLinkage,
                sugars: sugars[type],
                selected: linkage,
            })}
        </>
    ) : data === null ? (
        <Loading loadingText={'Crunching Numbers'} />
    ) : (
        <div className="flex flex-col w-full text-center">
            <span className="font-semibold">Query: {text}</span>
            <button
                className="align-self-start w-32 mb-2"
                onClick={() => {
                    setSearch(false);
                }}
            >
                &#8592; Back to Search
            </button>
            <DatabaseSearchTable data={data} />
        </div>
    );
}

export default function DatabaseSearch(props: {
    setSearchBegin: any;
    searchBegin: boolean;
}) {
    return !props.searchBegin
        ? ViewAllEntriesButton({
              label: 'View',
              text: (
                  <>
                      View all entries in the <i>Privateer</i> Database
                  </>
              ),
              onClickMethod: props.setSearchBegin,
          })
        : FilterZone({ setSearchBegin: props.setSearchBegin });
}
