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
                        PDB-REDO
                    </span>
                </label>
            </div>
        );
    }
}
