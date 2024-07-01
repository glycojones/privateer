import React, { type Dispatch, type SetStateAction } from 'react';

export default function toggleSwitch(props: {
    checkState: boolean;
    setCheckState: Dispatch<SetStateAction<boolean>>;
}) {
    const handleCheckboxChange = () => {
        props.setCheckState(!props.checkState);
    };

    return (
        <div className="mt-4">
            <label className="themeSwitcherTwo shadow-card relative inline-flex cursor-pointer select-none items-center justify-center rounded-md bg-white p-1">
                <input
                    type="checkbox"
                    className="sr-only"
                    checked={props.checkState}
                    onChange={handleCheckboxChange}
                />
                <span
                    className={`flex items-center space-x-[6px] rounded py-2 px-[18px] text-sm font-medium ${
                        !props.checkState
                            ? 'text-primary bg-[#f4f7ff]'
                            : 'text-body-color'
                    }`}
                >
                    PDB
                </span>
                <span
                    className={`flex items-center space-x-[6px] rounded py-2 px-[18px] text-sm font-medium ${
                        props.checkState
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