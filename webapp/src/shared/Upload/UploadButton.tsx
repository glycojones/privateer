import React, { type ReactElement, useState } from 'react';
import { type UploadButtonProps } from '../../interfaces/types.ts';
export default function UploadButton(props: UploadButtonProps): ReactElement {
    const fileUpload: (e: any) => void = (e: any) => {
        if (e.target.files !== null) {
            for (let i = 0; i < e.target.files.length; i++) {
                const splitname = e.target.files[i].name.split('.');
                const ext = splitname[splitname.length - 1];

                if (ext === 'pdb' || ext === 'cif' || ext === 'mmcif') {
                    props.setCoordinateFile(e.target.files[i] as File);
                } else if (ext === 'mtz') {
                    props.setReflectionFile(e.target.files[i] as File);
                } else {
                    console.error('Unknown File Type.');
                }
            }
        }
    };

    const [dragActivate, setDragActive] = useState(false);

    const handleDragOver: (e: any) => void = (e: any) => {
        e.preventDefault();
        e.stopPropagation();
        setDragActive(true);
    };

    const handleDrop: (e: any) => void = (e: any) => {
        e.preventDefault();
        e.stopPropagation();
        if (e.dataTransfer.files[0] !== null) {
            props.setCoordinateFile(e.dataTransfer.files[0] as File);
        }
        setDragActive(false);
    };

    return (
        <div
            id="uploadButton"
            className="flex items-center justify-center m-12 w-64 hover:border-slate-500"
        >
            {!dragActivate ? (
                <label className="flex flex-col items-center justify-center w-full p-12 h-64 border-2 border-gray-300 border-dashed rounded-lg cursor-pointer hover:bg-hover border-gray-600">
                    <div
                        className="flex flex-col items-center justify-center pt-5 pb-6 text-center"
                        onDragOver={handleDragOver}
                        onDrop={handleDrop}
                        onDragEnter={() => {
                            console.log('DRAG ENTER');
                        }}
                        onDragExit={(e) => {
                            console.log('LEAVE', e);
                        }}
                    >
                        <svg
                            className="w-8 h-8 mb-4 text-gray-500 dark:text-gray-400"
                            aria-hidden="true"
                            xmlns="http://www.w3.org/2000/svg"
                            fill="none"
                            viewBox="0 0 20 16"
                        >
                            <path
                                stroke="currentColor"
                                strokeLinecap="round"
                                strokeLinejoin="round"
                                strokeWidth="2"
                                d="M13 13h3a3 3 0 0 0 0-6h-.025A5.56 5.56 0 0 0 16 6.5 5.5 5.5 0 0 0 5.207 5.021C5.137 5.017 5.071 5 5 5a4 4 0 0 0 0 8h2.167M10 15V6m0 0L8 8m2-2 2 2"
                            />
                        </svg>
                        <p className="mb-2 text-md text-gray-500 dark:text-gray-400">
                            <span className="font-semibold">Choose a file</span>
                        </p>
                        <p className="text-sm text-gray-500 dark:text-gray-400">
                            PDB, mmCIF or MTZ.
                            <br />
                            Files will never be sent externally.
                        </p>
                    </div>
                    <input
                        id="dropzone-file"
                        type="file"
                        className="hidden"
                        accept=".pdb,.mmcif,.cif,.mtz"
                        onChange={fileUpload}
                    />
                </label>
            ) : (
                <label className="flex flex-col items-center justify-center w-full p-12 h-64 border-2 border-gray-300  rounded-lg cursor-pointer hover:bg-hover border-gray-600">
                    <div
                        className="flex flex-col items-center justify-center pt-5 pb-6 text-center"
                        onDragOver={handleDragOver}
                        onDrop={handleDrop}
                        onDragEnter={() => {
                            console.log('DRAG ENTER');
                        }}
                        onDragExit={(e) => {
                            console.log('LEAVE', e);
                        }}
                    >
                        <svg
                            className="w-8 h-8 mb-4 text-gray-500 dark:text-gray-400"
                            aria-hidden="true"
                            xmlns="http://www.w3.org/2000/svg"
                            fill="none"
                            viewBox="0 0 20 16"
                        >
                            <path
                                stroke="currentColor"
                                strokeLinecap="round"
                                strokeLinejoin="round"
                                strokeWidth="2"
                                d="M13 13h3a3 3 0 0 0 0-6h-.025A5.56 5.56 0 0 0 16 6.5 5.5 5.5 0 0 0 5.207 5.021C5.137 5.017 5.071 5 5 5a4 4 0 0 0 0 8h2.167M10 15V6m0 0L8 8m2-2 2 2"
                            />
                        </svg>
                        <p className="mb-2 text-md text-gray-500 dark:text-gray-400">
                            <span className="font-semibold">Drop here!</span>
                        </p>
                    </div>
                </label>
            )}
        </div>
    );
}
