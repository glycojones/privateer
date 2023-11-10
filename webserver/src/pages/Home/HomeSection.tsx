import { lazy, useEffect, useState } from "react";

import { Header } from '../../layouts/Header';
import { Information } from '../../components/Information/Information.tsx';

// @ts-ignore
import privateer_module from "../../wasm/privateer.js"
import loadGlytoucan from "../../utils/loadGlytoucan.ts"

const Footer = lazy(() => import('../../layouts/Footer.tsx'));
const BorderElement = lazy(() => import('../../layouts/BorderElement.tsx'));

import { fetch_map, fetch_pdb } from "../../utils/fetch_from_pdb.ts"
 
import {TableDataEntry, HeaderProps} from "../../interfaces/types"


export default function HomeSection() {
    const [coordinateFile, setCoordinateFile] = useState<File | null>(null);
    const [reflectionFile, setReflectionFile] = useState<File| null>(null);
    const [PDBCode, setPDBCode] = useState<string>("")
    const [fileContent, setFileContent] = useState<string | ArrayBuffer>("")
    const [mtzData, setMtzData] = useState<Uint8Array | null>(null)
    const [submit, setSubmit] = useState<boolean>(false);
    const [tableData, setTableData] = useState<Array<TableDataEntry> | null>(null);
    const [loadingText, setLoadingText] = useState<string>("Validating Glycans...");
    const [resetApp, setResetApp] = useState<boolean>(false)
    const [fallback, setFallBack] = useState<boolean>(false)
    const [failureText, setFailureText] = useState<string>("")

    let sanitize_id = (id: string) => {
        const regex = /: *32/g;
        const new_id = id.replace(regex, "")
        return new_id
    }

    async function run_privateer(Module: any, fileContent: string | ArrayBuffer, name: string) {

        setFileContent(fileContent)

        let x = Module.read_structure_to_table(fileContent, name)

        let table_data = [];
        for (var i = 0; i < x.size(); i++) {
            let table_entry = x.get(i)

            table_entry.id = sanitize_id(table_entry.id)

            let collected_torsions = []
            for (var j = 0; j < table_entry.torsions.size(); j++) {
                collected_torsions.push(table_entry.torsions.get(j));
            }
            table_entry.torsions = collected_torsions
            const regex = /: *32/g;
            table_entry.svg = table_entry.svg.replace(regex, "")
            table_data.push(table_entry)

        }

        if (x.size() == 0) {
            setLoadingText("There were no detected glycans in this file.")
            setFallBack(true)
        }

        // Get Glyconnect ID from WURCS
        setLoadingText("Querying Glytoucan...")
        await loadGlytoucan(table_data)


        setTableData(table_data);
    }

    useEffect(() => {
        if (PDBCode != "") {
            setLoadingText(`Fetching ${PDBCode.toUpperCase()} from the PDB`)

            fetch_map(PDBCode).then((response: ArrayBuffer) => {
                let array = new Uint8Array(response)
                setMtzData(array)
            }).catch((e: any) => {
                setLoadingText("MTZ not found, continuing...")
            })

            fetch_pdb(PDBCode).then((response: ArrayBuffer) => {
                setFileContent(response)
                setLoadingText("Validating Glycans...")

                privateer_module().then((Module: any) => run_privateer(Module, response, PDBCode))

            }).catch((e: any) => {
                setFailureText("This PDB code could not be found")
                setLoadingText("There were no detected glycans in this file.")
                setFallBack(true)
            })

        } else {
            privateer_module().then((Module: { [x: string]: (arg0: string, arg1: string, arg2: Uint8Array, arg3: boolean, arg4: boolean, arg5: boolean) => void; }) => {


                var coordinateReader = new FileReader();
                var reflectionReader = new FileReader();

                coordinateReader.onload = () => { run_privateer(Module, coordinateReader.result, coordinateFile!.name) }

                if (coordinateFile) {
                    coordinateReader.readAsText(coordinateFile);
                }

                reflectionReader.onload = async () => {
                    let map_data = new Uint8Array(reflectionReader.result);
                    setMtzData(map_data)

                    Module['FS_createDataFile']('/', "input.mtz", map_data, true, true, true)
                }

                if (reflectionFile) {
                    reflectionReader.readAsArrayBuffer(reflectionFile)
                }

            }).catch((e: any) => console.log(e));
        }


    }, [submit])

    useEffect(() => {
        setReflectionFile(null)
        setCoordinateFile(null)
        setSubmit(false)
        setTableData(null)
        setFallBack(false)
        setResetApp(false)
        setPDBCode("")
    }, [resetApp])

    const main_props: HeaderProps = {
        resetApp: resetApp,
        setResetApp: setResetApp,
        PDBCode: PDBCode,
        setPDBCode: setPDBCode,
        coordinateFile: coordinateFile,
        setCoordinateFile: setCoordinateFile,
        reflectionFile: reflectionFile,
        setReflectionFile: setReflectionFile,
        submit: submit,
        setSubmit: setSubmit,
        tableData: tableData,
        loadingText: loadingText,
        fileContent: fileContent,
        fallback: fallback,
        mtzData: mtzData,
        failureText: failureText
    }

    return (
        <>
            <Header {...main_props} />
            <BorderElement topColor={"#D6D9E5"} bottomColor={"#F4F9FF"} reverse={false}></BorderElement>
            <Information />
            <BorderElement topColor={"#F4F9FF"} bottomColor={"#D6D9E5"} reverse={true}></BorderElement>
            <Footer></Footer>
        </>
    )
}