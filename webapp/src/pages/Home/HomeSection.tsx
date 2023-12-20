import React, { lazy, useEffect, useState } from 'react';
import { Header } from '../../layouts/Header';
import { Information } from '../../shared/Information/Information.tsx';

// @ts-expect-error: Emscripten Generated JS does not conform to typescript conventions
import privateer_module from '../../wasm/privateer.js';
import loadGlytoucan, {
    loadGlytoucanFromFile,
} from '../../utils/loadGlytoucan.ts';

import { fetchMap, fetchPDB } from '../../utils/fetch_from_pdb.ts';

import { type HeaderProps, type TableDataEntry } from '../../interfaces/types';

const Footer = lazy(async () => await import('../../layouts/Footer.tsx'));
const BorderElement = lazy(
    async () => await import('../../layouts/BorderElement.tsx')
);

export default function HomeSection(): Element {
    const [coordinateFile, setCoordinateFile] = useState<File | null>(null);
    const [reflectionFile, setReflectionFile] = useState<File | null>(null);
    const [PDBCode, setPDBCode] = useState<string>('');
    const [fileContent, setFileContent] = useState<string | ArrayBuffer>('');
    const [mtzData, setMtzData] = useState<Uint8Array | null>(null);
    const [submit, setSubmit] = useState<boolean>(false);
    const [tableData, setTableData] = useState<TableDataEntry[] | null>(null);
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

        const x = Module.read_structure_to_table(fileContent, name);

        const tableData: TableDataEntry[] = [];
        for (let i = 0; i < x.size(); i++) {
            const tableEntry = x.get(i);

            tableEntry.id = sanitizeID(tableEntry.id as string);

            const collectedTorsions: any[] = [];
            for (let j = 0; j < tableEntry.torsions.size(); j++) {
                collectedTorsions.push(tableEntry.torsions.get(j));
            }
            tableEntry.torsions = collectedTorsions;
            const regex = /: *32/g;
            tableEntry.svg = tableEntry.svg.replace(regex, '');
            tableData.push(tableEntry as TableDataEntry);
        }

        if (x.size() === 0) {
            setFailureText('There were no detected glycans in this model.');
            setFallBack(true);
        }

        // Get Glyconnect ID from WURCS
        setLoadingText('Querying Glytoucan...');
        await loadGlytoucanFromFile(tableData);

        setTableData(tableData);
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
                    .then((response: ArrayBuffer) => {
                        setFileContent(response);
                        setLoadingText('Validating Glycans...');

                        privateer_module().then(
                            async (Module: any) => {
                                await runPrivateer(Module, response, PDBCode);
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
            <Header {...mainProps} />
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
