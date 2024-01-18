import GlycansVsYear from './GlycansVsYear.tsx';
import React, { Suspense, useEffect, useState } from 'react';
import Loading from '../../shared/Loading/Loading.tsx';
import ErrorsVsYear from './ErrorsVsYear.tsx';
import ErrorsVsResolution from './ErrorsVsResolution.tsx';
export default function Graphs() {
    const [lastUpdated, setLastUpdated] = useState<string>();

    useEffect(() => {
        const url =
            'https://raw.githubusercontent.com/Dialpuri/PrivateerDatabase/master/stats/last_updated.json';

        fetch(url)
            .then(async (response) => await response.json())
            .then((data) => {
                setLastUpdated(data.date as string);
            })
            .catch((error) => {
                console.log('Error reading last updated date', error);
            });
    }, []);

    const [isChecked, setIsChecked] = useState(false);

    const handleCheckboxChange = () => {
        setIsChecked(!isChecked);
    };

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

    return (
        <>
            <Suspense
                fallback={<Loading loadingText={'Getting latest data'} />}
            >
                <div className="flex flex-col items-center justify-center">
                    <h2 className="w-full text-center sm:text-left pl-2 sm:pl-12">
                        Protein Data Bank Statistics
                    </h2>
                    <h4 className="w-full text-center sm:text-left pl-2 sm:pl-12">
                        Last Updated - {lastUpdated}
                    </h4>

                    <GlycansVsYear />
                    <div className="flex flex-col sm:flex-row flex-wrap sm:justify-between w-full items-center align-center">
                        <div>
                            <h2 className="w-full text-center sm:text-left pl-2 sm:pl-12">
                                Validation Statistics
                            </h2>
                            <h4 className="w-full text-center sm:text-left pl-2 sm:pl-12">
                                Last Updated - {lastUpdated}
                            </h4>
                        </div>
                        <div className="sm:px-12">{toggleSwitch()}</div>
                    </div>
                    <ErrorsVsYear database={isChecked ? 'pdbredo' : 'pdb'} />
                    <ErrorsVsResolution
                        database={isChecked ? 'pdbredo' : 'pdb'}
                    />
                </div>
            </Suspense>
        </>
    );
}
