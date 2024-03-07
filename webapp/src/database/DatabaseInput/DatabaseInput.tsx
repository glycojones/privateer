import React, { useState } from 'react';
import DatabaseFetch from '../DatabaseFetch/DatabaseFetch';
import UploadButton from '../../shared/Upload/UploadButton.tsx';
import Submit from '../../shared/Submit/Submit.tsx';
import PDBFetch from '../../shared/PDBFetch/PDBFetch.tsx';
import DatabaseSearch from '../DatabaseSearch/DatabaseSearch.tsx';

export default function DatabaseInput(props: {
    PDBCode: string;
    setPDBCode: any;
    submitPressed: boolean;
}) {
    const [showSearchBegin, setSearchBegin] = useState(false);

    return (
        <>
            <div
                id="upload"
                className="flex flex-wrap align-middle items-center justify-center transition-none"
            >
                {
                    <DatabaseSearch
                        searchBegin={showSearchBegin}
                        setSearchBegin={setSearchBegin}
                    />
                }

                {!showSearchBegin ? (
                    <div className="mx-6 w-full lg:w-6 sm:w-full text-center">
                        OR
                    </div>
                ) : (
                    <></>
                )}
                {!showSearchBegin ? (
                    <DatabaseFetch
                        PDBCode={props.PDBCode}
                        setPDBCode={props.setPDBCode}
                        submitPressed={props.submitPressed}
                    />
                ) : (
                    <></>
                )}
            </div>
        </>
    );
}
