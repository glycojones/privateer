import React, { lazy, useEffect, useState } from 'react';
import { MainHeader } from '../../main/Header/MainHeader.tsx';
import { Information } from '../../shared/Information/Information.tsx';

// @ts-expect-error: Emscripten Generated JS does not conform to typescript conventions
import privateer_module from '../../wasm/privateer.js';
import { loadGlytoucanFromFile } from '../../utils/loadGlytoucan.ts';

import { fetchMap, fetchPDB } from '../../utils/fetch_from_pdb.ts';

import {
    type HeaderProps,
    type ResultsEntry,
    type TorsionEntry,
} from '../../interfaces/types';
import { detect } from 'detect-browser';

const Footer = lazy(async () => await import('../../layouts/Footer.tsx'));
const BorderElement = lazy(
    async () => await import('../../layouts/BorderElement.tsx')
);

export default function Home(): Element {
    const [coordinateFile, setCoordinateFile] = useState<File | null>(null);
    const [reflectionFile, setReflectionFile] = useState<File | null>(null);
    const [PDBCode, setPDBCode] = useState<string>('');
    const [fileContent, setFileContent] = useState<string | ArrayBuffer>('');
    const [mtzData, setMtzData] = useState<Uint8Array | null>(null);
    const [submit, setSubmit] = useState<boolean>(false);
    const [tableData, setTableData] = useState<ResultsEntry[] | null>(null);
    const [loadingText, setLoadingText] = useState<string>(
        'Validating Glycans...'
    );
    const [resetApp, setResetApp] = useState<boolean>(false);
    const [fallback, setFallBack] = useState<boolean>(false);
    const [failureText, setFailureText] = useState<string>('');

    function sanitizeID(id: string): string {
        const regex = /: *32/g;
        return id.replace(regex, '');
    }

    async function runPrivateer(
        Module: any,
        fileContent: string | ArrayBuffer | null,
        name: string | null
    ): Promise<void> {
        if (fileContent === null || name === null) {
            await Promise.resolve();
            return;
        }
        setFileContent(fileContent);

        const splitFileName = name.split('.');
        const fileExtension = splitFileName[splitFileName.length - 1];
        const newFileName = '/coordinates.' + fileExtension;
        Module.FS.writeFile(newFileName, fileContent);

        const disallowedExtensions = ["cif", "mmcif"]
        const browser = detect(); // FireFox doesn't work with CIF files, get the PDB.
        if (browser.name === 'firefox' && disallowedExtensions.includes(fileExtension)) {
            setFailureText(
                'CIF files are currently unsupported on FireFox. Please consider an alternate browser.'
            );
            setFallBack(true);
            return;
        }

        const results = Module.validate(newFileName, name);
        const data: ResultsEntry[] = [];
        const resultSize = results.size();
        for (let i = 0; i < resultSize; i++) {
            const entry = results.get(i);

            const collectedTorsions: TorsionEntry[] = [];
            const torsionsVec = entry.torsions;
            const torsionSize = torsionsVec.size();
            for (let j = 0; j < torsionSize; j++) {
                collectedTorsions.push(torsionsVec.get(j) as TorsionEntry);
            }
            torsionsVec.delete();

            const regex: RegExp = /: *32/g;
            const sanitisedSVG = entry.svg.replace(regex, '');
            const sanitisedID = sanitizeID(entry.id as string);
            // @ts-expect-error
            const entryJS: ResultsEntry = {
                torsions: collectedTorsions,
                svg: sanitisedSVG,
                id: sanitisedID,
                chain: entry.chain,
                wurcs: entry.wurcs,
                glyconnect_id: entry.glyconnect_id,
                glytoucan_id: entry.glytoucan_id,
                torsion_err: entry.torsion_err,
                anomer_err: entry.anomer_err,
                conformation_err: entry.conformation_err,
                puckering_err: entry.puckering_err,
                chirality_err: entry.chirality_err,
            };

            data.push(entryJS);
        }

        results.delete();

        if (data.length === 0) {
            setFailureText('There were no detected glycans in this model.');
            setFallBack(true);
        }

        // Get Glyconnect ID from WURCS
        setLoadingText('Querying Glytoucan...');
        await loadGlytoucanFromFile(data);

        setTableData(data);
    }

    useEffect(() => {
        async function handleLoad(): Promise<void> {
            if (PDBCode !== '') {
                setLoadingText(
                    `Fetching ${PDBCode.toUpperCase()} from the PDB`
                );

                try {
                    const mapResponse = await fetchMap(PDBCode);
                    const mapArray = new Uint8Array(mapResponse);
                    setMtzData(mapArray);
                } catch (err) {
                    console.log('No map found, continuing...');
                }

                fetchPDB(PDBCode)
                    .then((response: [string] | void) => {
                        // @ts-expect-error
                        const [fileContent, fileExtension]: [string, string] =
                            response;
                        setFileContent(fileContent);
                        setLoadingText('Validating Glycans...');

                        const fileName = PDBCode + fileExtension;
                        privateer_module().then(
                            async (Module: any) => {
                                await runPrivateer(
                                    Module,
                                    fileContent,
                                    fileName
                                );
                            },
                            () => {}
                        );
                    })
                    .catch(() => {
                        setFailureText('This PDB code could not be found');
                        setLoadingText(
                            'There were no detected glycans in this file.'
                        );
                        setFallBack(true);
                    });
            } else {
                privateer_module()
                    .then(
                        (
                            Module: Record<
                                string,
                                (
                                    arg0: string,
                                    arg1: string,
                                    arg2: Uint8Array,
                                    arg3: boolean,
                                    arg4: boolean,
                                    arg5: boolean
                                ) => void
                            >
                        ) => {
                            const coordinateReader = new FileReader();
                            const reflectionReader = new FileReader();

                            coordinateReader.onload = () => {
                                runPrivateer(
                                    Module,
                                    coordinateReader.result,
                                    coordinateFile.name
                                ).then(
                                    () => {},
                                    () => {}
                                );
                            };

                            if (coordinateFile !== null) {
                                coordinateReader.readAsText(coordinateFile);
                                console.log(coordinateFile);
                            }

                            reflectionReader.onload = async () => {
                                const mapData = new Uint8Array(
                                    reflectionReader.result
                                );
                                setMtzData(mapData);

                                Module.FS_createDataFile(
                                    '/',
                                    'input.mtz',
                                    mapData,
                                    true,
                                    true,
                                    true
                                );
                            };

                            if (reflectionFile !== null) {
                                reflectionReader.readAsArrayBuffer(
                                    reflectionFile
                                );
                            }
                        }
                    )
                    .catch((e: any) => {
                        console.log(e);
                    });
            }
        }
        handleLoad()
            .then(
                () => {},
                () => {}
            )
            .catch(() => {});
    }, [submit]);

    useEffect(() => {
        setReflectionFile(null);
        setCoordinateFile(null);
        setSubmit(false);
        setTableData(null);
        setFallBack(false);
        setResetApp(false);
        setPDBCode('');
    }, [resetApp]);

    const mainProps: HeaderProps = {
        resetApp,
        setResetApp,
        PDBCode,
        setPDBCode,
        coordinateFile,
        setCoordinateFile,
        reflectionFile,
        setReflectionFile,
        submit,
        setSubmit,
        tableData,
        loadingText,
        fileContent,
        fallback,
        mtzData,
        failureText,
    };

    return (
        <>
            <MainHeader {...mainProps} />
            <BorderElement
                topColor={'#D6D9E5'}
                bottomColor={'#F4F9FF'}
                reverse={false}
            ></BorderElement>
            <Information />
            <BorderElement
                topColor={'#F4F9FF'}
                bottomColor={'#D6D9E5'}
                reverse={true}
            ></BorderElement>
            <Footer></Footer>
        </>
    );
}
