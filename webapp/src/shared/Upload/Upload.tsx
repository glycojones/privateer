import React, { type ReactElement, useEffect, useState } from 'react';
import UploadButton from './UploadButton.tsx';
import Submit from '../Submit/Submit.tsx';
import PDBFetch from '../PDBFetch/PDBFetch';
import {
    type HeaderProps,
    type UploadButtonProps,
} from '../../interfaces/types.ts';

export default function Upload(props: HeaderProps): ReactElement {
    const [showUploadAgain, setShowUploadAgain] = useState(true);
    const [showSubmit, setShowSubmit] = useState(false);
    const [allowSubmit, setAllowSubmit] = useState(false);
    const [showPDBFetch, setShowPDBFetch] = useState(true);

    useEffect(() => {
        setShowPDBFetch(true);
    }, [props.resetApp]);

    useEffect(() => {
        if (props.coordinateFile !== null && props.reflectionFile != null) {
            setShowSubmit(true);
            setShowUploadAgain(false);
            setAllowSubmit(true);
            setShowPDBFetch(false);
        }

        if (props.coordinateFile !== null && props.reflectionFile === null) {
            setShowSubmit(true);
            setShowUploadAgain(true);
            setAllowSubmit(true);
            setShowPDBFetch(false);
        }

        if (props.coordinateFile === null && props.reflectionFile !== null) {
            setShowSubmit(true);
            setShowUploadAgain(true);
            setAllowSubmit(false);
            setShowPDBFetch(false);
        }

        if (props.coordinateFile === null && props.reflectionFile === null) {
            setShowSubmit(false);
            setShowUploadAgain(true);
            setAllowSubmit(false);
            setShowPDBFetch(true);
        }
    }, [props.coordinateFile, props.reflectionFile]);

    const uploadButtonProps: UploadButtonProps = {
        setCoordinateFile: props.setCoordinateFile,
        setReflectionFile: props.setReflectionFile,
    };

    return (
        <div
            id="upload"
            className="flex flex-wrap align-middle items-center justify-center"
        >
            {showUploadAgain ? <UploadButton {...uploadButtonProps} /> : <></>}
            {showSubmit ? (
                <Submit
                    coordinateFile={props.coordinateFile}
                    reflectionFile={props.reflectionFile}
                    submitPressed={props.setSubmit}
                    setResetApp={props.setResetApp}
                    allowSubmit={allowSubmit}
                />
            ) : (
                <></>
            )}

            {showPDBFetch ? (
                <div className="mx-6 w-full lg:w-6 sm:w-full text-center">
                    OR
                </div>
            ) : (
                <></>
            )}
            {showPDBFetch ? (
                <PDBFetch
                    PDBCode={props.PDBCode}
                    setPDBCode={props.setPDBCode}
                    submitPressed={props.setSubmit}
                />
            ) : (
                <></>
            )}
        </div>
    );
}