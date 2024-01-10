import React, { useEffect, useState } from 'react';
import { type DatabaseResultProps } from '../../interfaces/types.ts';
import CremerPopleGraph from '../DatabaseComponents/CremerPopleGraph.tsx';
import BFactorVsRSCC from '../DatabaseComponents/BFactorVsRSCC.tsx';
import SNFGList from '../DatabaseComponents/SNFGList.tsx';
import SugarList from '../DatabaseComponents/SugarList.tsx';

export default function DatabaseResult(props: DatabaseResultProps) {
    const [pdbShown, setPDBShown] = useState<boolean>(true);
    const [selectedData, setSelectedData] = useState<any>();
    const [isChecked, setIsChecked] = useState(false);

    const handleCheckboxChange = () => {
        setIsChecked(!isChecked);
        toggleDataSource();
    };
    useEffect(() => {
        if (props.pdbResults === null) return;
        if (props.pdbRedoResults === null) return;

        console.log('running use effect');
        setSelectedData(props.pdbResults);
    }, []);

    function toggleDataSource() {
        if (pdbShown) {
            setSelectedData(props.pdbRedoResults);
            setPDBShown(false);
        } else {
            setSelectedData(props.pdbResults);
            setPDBShown(true);
        }
    }

    return (
        <>
            {selectedData !== undefined ? (
                <div className="flex flex-col space-y-6">
                    <h2 className="text-center">
                        Validation Report - {props.PDBCode}
                        {props.pdbRedoResults !== '' ? toggleSwitch() : <></>}
                    </h2>
                    <div className="flex flex-wrap text-center">
                        <CremerPopleGraph data={selectedData} />
                        <BFactorVsRSCC data={selectedData} />
                    </div>
                    <SNFGList data={selectedData} />
                    <SugarList data={selectedData} />
                </div>
            ) : (
                <></>
            )}
        </>
    );

    function pdbRedoSVG() {
        return (
            <div className="w-12">
                <svg
                    version="1.1"
                    id="Layer_1"
                    xmlns="http://www.w3.org/2000/svg"
                    xmlnsXlink="http://www.w3.org/1999/xlink"
                    x="0px"
                    y="0px"
                    viewBox="0 0 160 120"
                    enableBackground="new 0 0 160 120"
                    xmlSpace="preserve"
                >
                    <g>
                        <g>
                            <path
                                fill="#3366CC"
                                d="M47.144,92.745c0,0.372-0.335,0.559-1.005,0.559h-4.731c-0.5,0-1.049-0.119-1.646-0.357
                    c-0.599-0.238-1.055-0.526-1.37-0.864l-3.352-3.458c-1.38-1.421-3.049-2.132-5.007-2.132h-4.811v5.59
                    c0,0.338-0.181,0.626-0.542,0.864s-0.799,0.357-1.311,0.357h-3.707c-0.512,0-0.952-0.119-1.32-0.357s-0.552-0.526-0.552-0.864
                    V75.897c0-0.338,0.177-0.628,0.532-0.871c0.354-0.242,0.788-0.363,1.301-0.363h16.264c3.01,0,5.478,0.532,7.403,1.599
                    c1.925,1.065,2.888,2.505,2.888,4.315c0,2.522-2.478,4.104-7.433,4.745c0.789,0.19,1.519,0.518,2.188,0.981
                    c0.671,0.464,1.439,1.138,2.307,2.021l3.588,3.757C47.038,92.324,47.144,92.546,47.144,92.745z M25.222,82.813h8.398
                    c1.354,0,2.546-0.178,3.578-0.533c1.031-0.355,1.547-0.923,1.547-1.703c0-0.779-0.516-1.348-1.547-1.703
                    c-1.032-0.354-2.225-0.532-3.578-0.532h-8.398V82.813z"
                            />
                            <path
                                fill="#3366CC"
                                d="M50.416,92.056V75.871c0-0.338,0.184-0.623,0.552-0.857s0.809-0.351,1.32-0.351h22.73
                    c0.513,0,0.949,0.121,1.311,0.363c0.361,0.243,0.543,0.533,0.543,0.871v1.222c0,0.338-0.182,0.627-0.543,0.865
                    s-0.798,0.357-1.311,0.357H57.848v3.782h14.785c0.513,0,0.95,0.119,1.312,0.357s0.542,0.526,0.542,0.865v1.222
                    c0,0.338-0.181,0.626-0.542,0.864s-0.799,0.357-1.312,0.357H57.848v3.835h17.132c0.499,0,0.926,0.117,1.281,0.351
                    c0.354,0.234,0.532,0.521,0.532,0.858v1.248c0,0.338-0.185,0.626-0.552,0.864c-0.368,0.238-0.809,0.357-1.321,0.357H52.288
                    c-0.512,0-0.952-0.119-1.32-0.357S50.416,92.411,50.416,92.056z"
                            />
                            <path
                                fill="#3366CC"
                                d="M80.5,92.107V75.924c0-0.355,0.184-0.654,0.552-0.897c0.368-0.242,0.809-0.363,1.32-0.363h11.021
                    c5.704,0,9.854,0.752,12.449,2.255c2.596,1.504,3.894,3.861,3.894,7.071c0,3.211-1.298,5.566-3.894,7.065
                    c-2.596,1.5-6.745,2.249-12.449,2.249H82.372c-0.512,0-0.952-0.115-1.32-0.345S80.5,92.445,80.5,92.107z M87.932,89.625h5.244
                    c3.535,0,5.94-0.425,7.215-1.274c1.275-0.849,1.912-2.303,1.912-4.361c0-2.058-0.637-3.514-1.912-4.367
                    c-1.274-0.854-3.68-1.28-7.215-1.28h-5.244V89.625z"
                            />
                            <path
                                fill="#3366CC"
                                d="M139.867,91.373c-2.688,1.556-6.631,2.334-11.828,2.334c-5.198,0-9.145-0.778-11.838-2.334
                    c-2.694-1.556-4.042-4.019-4.042-7.39s1.348-5.832,4.042-7.384c2.693-1.551,6.64-2.327,11.838-2.327
                    c5.197,0,9.141,0.776,11.828,2.327c2.688,1.552,4.031,4.013,4.031,7.384S142.555,89.817,139.867,91.373z M134.476,79.343
                    c-1.327-1.04-3.473-1.561-6.437-1.561s-5.112,0.521-6.446,1.561s-2.001,2.589-2.001,4.646c0,2.059,0.667,3.605,2.001,4.642
                    c1.334,1.035,3.482,1.553,6.446,1.553s5.109-0.518,6.437-1.553c1.327-1.036,1.991-2.583,1.991-4.642
                    C136.467,81.932,135.803,80.383,134.476,79.343z"
                            />
                        </g>
                    </g>
                    <g>
                        <path
                            fill="#FF9933"
                            d="M23.095,52.253v19.542H18.2V24.766h20.604c6.278,0,11.138,1.254,14.578,3.764
                c3.441,2.507,5.161,5.83,5.161,9.964c0,4.178-1.72,7.515-5.161,10.014c-3.44,2.499-8.3,3.746-14.578,3.746H23.095z M23.095,48.959
                h15.709c4.936,0,8.643-0.985,11.124-2.956s3.723-4.452,3.723-7.445c0-3.015-1.235-5.519-3.702-7.51
                c-2.468-1.993-6.184-2.988-11.145-2.988H23.095V48.959z"
                        />
                        <path
                            fill="#FF9933"
                            d="M60.034,71.795V24.766h17.396c7.319,0,13.363,1.916,18.136,5.749c4.771,3.833,7.154,8.753,7.154,14.761
                v6.041c0,6.029-2.384,10.95-7.154,14.762c-4.772,3.813-10.816,5.717-18.136,5.717H60.034z M64.928,28.06V68.5h12.503
                c5.921,0,10.794-1.622,14.618-4.873c3.825-3.249,5.737-7.358,5.737-12.329v-6.132c0-4.906-1.919-8.983-5.757-12.232
                c-3.839-3.249-8.706-4.874-14.599-4.874H64.928z"
                        />
                        <path
                            fill="#FF9933"
                            d="M104.607,71.795V24.766h17.478c6.223,0,11.069,1.029,14.537,3.087c3.468,2.057,5.202,5.153,5.202,9.29
                c0,2.414-0.87,4.514-2.61,6.304c-1.742,1.788-4.107,3.026-7.096,3.716c3.646,0.518,6.607,1.863,8.884,4.037
                c2.274,2.172,3.413,4.681,3.413,7.521c0,4.198-1.733,7.426-5.201,9.687c-3.47,2.259-8.109,3.389-13.921,3.389H104.607z
                 M109.5,45.857h13.778c4.33,0,7.683-0.747,10.054-2.245c2.371-1.496,3.558-3.698,3.558-6.605c0-2.972-1.263-5.206-3.784-6.702
                c-2.522-1.498-6.195-2.245-11.021-2.245H109.5V45.857z M109.5,49.151V68.5h15.793c4.414,0,7.882-0.852,10.403-2.56
                c2.522-1.705,3.784-4.092,3.784-7.162c0-2.812-1.146-5.118-3.434-6.922c-2.29-1.803-5.614-2.705-9.973-2.705H109.5z"
                        />
                    </g>
                </svg>
            </div>
        );
    }

    function toggleSwitch(): React.ReactNode {
        return (
            <div className="mt-4">
                <label className="themeSwitcherTwo shadow-card relative inline-flex cursor-pointer select-none items-center justify-center rounded-md bg-white p-1">
                    <input
                        type="checkbox"
                        className="sr-only"
                        checked={isChecked}
                        onChange={handleCheckboxChange}
                    />
                    <span
                        className={`flex items-center space-x-[6px] rounded py-2 px-[18px] text-sm font-medium ${
                            !isChecked
                                ? 'text-primary bg-[#f4f7ff]'
                                : 'text-body-color'
                        }`}
                    >
                        PDB
                    </span>
                    <span
                        className={`flex items-center space-x-[6px] rounded py-2 px-[18px] text-sm font-medium ${
                            isChecked
                                ? 'text-primary bg-[#f4f7ff]'
                                : 'text-body-color'
                        }`}
                    >
                        {/* {pdbRedoSVG()} */}
                        PDB-REDO
                    </span>
                </label>
            </div>
        );
    }
}
