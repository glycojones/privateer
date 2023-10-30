import {lazy, useEffect, useState} from "react";

import {Header} from '../../layouts/Header';
import {Information} from '../../components/Information/Information';

import privateer_module from "../../wasm/privateer.js"
import loadGlytoucan from "../../utils/loadGlytoucan"

const Footer = lazy(() => import('../../layouts/Footer'));
const BorderElement = lazy(() => import('../../layouts/BorderElement'));

import {fetch_map, fetch_pdb} from "../../utils/fetch_from_pdb"

export default function HomeSection() {
    const [coordinateFile, setCoordinateFile] = useState(null);
    const [reflectionFile, setReflectionFile] = useState(null);
    const [PDBCode, setPDBCode] = useState("")
    const [fileContent, setFileContent] = useState(null)
    const [mtzData, setMtzData] = useState()
    const [submit, setSubmit] = useState(null);
    const [tableData, setTableData] = useState(null);
    const [loadingText, setLoadingText] = useState("Validating Glycans...");
    const [resetApp, setResetApp] = useState(false)
    const [fallback, setFallBack] = useState(false)
    const [failureText, setFailureText] = useState(false)

    let sanitize_id = (id) => { 
        const regex = /: *32/g;
        const new_id = id.replace(regex, "")
        return new_id
    }

    async function run_privateer(Module, fileContent, name) {
        setFileContent(fileContent)

        let x = Module.read_structure_to_table(fileContent, name)

        let table_data = [];
        for (var i = 0; i < x.size(); i++) {
            let table_entry = x.get(i)

            table_entry.id = sanitize_id(table_entry.id)

            let collected_torsions = []
            for(var j = 0; j < table_entry.torsions.size(); j++) { 
                collected_torsions.push(table_entry.torsions.get(j)); 
            }
            table_entry.torsions = collected_torsions
            const regex = /: *32/g;
            table_entry.svg = table_entry.svg.replace(regex, "")
            table_data.push(table_entry)

        }
        
        if (x.size() == 0 ) { 
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
        
            fetch_map(PDBCode).then((response) => { 
                let array = new Uint8Array(response)
                setMtzData(array)
            }).catch((e) => {
                setLoadingText("MTZ not found, continuing...")
            }) 

            fetch_pdb(PDBCode).then((response) => { 
                setFileContent(response)
                setLoadingText("Validating Glycans...")

                privateer_module().then(Module => run_privateer(Module, response, PDBCode))

            }).catch(e => {
                setFailureText("This PDB code could not be found")
                setLoadingText("There were no detected glycans in this file.")
                setFallBack(true)})
           
         } else {
            privateer_module().then((Module) => {


                var coordinateReader = new FileReader();
                var reflectionReader = new FileReader();
                
                coordinateReader.onload = () => {run_privateer(Module, coordinateReader.result, coordinateFile.name)}
    
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
    
            }).catch((e) => console.log(e));
         }

        
    }, [submit])

    useEffect(() => {
        setReflectionFile(null)
        setCoordinateFile(null)
        setSubmit(null)
        setTableData(null)
        setFallBack(false)
        setResetApp(false)
        setPDBCode("")
    }, [resetApp])

    const main_props = {
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
            <Header {...main_props}/>
            <BorderElement topColor={"#D6D9E5"} bottomColor={"#F4F9FF"}></BorderElement>
            <Information/>
            <BorderElement topColor={"#F4F9FF"} bottomColor={"#D6D9E5"} reverse={true}></BorderElement>
            <Footer></Footer>
        </>
    )
}