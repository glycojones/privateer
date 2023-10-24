import {lazy, useEffect, useState} from "react";

import {Header} from '../../layouts/Header';
import {Information} from '../../components/Information/Information';

import privateer_module from "../../wasm/privateer.js"
import loadGlytoucan from "../../utils/loadGlytoucan"
// import Footer from "../../layouts/Footer"
// import BorderElement from '../../layouts/BorderElement';

const Footer = lazy(() => import('../../layouts/Footer'));
const BorderElement = lazy(() => import('../../layouts/BorderElement'));
// const Header = lazy(() => import('../../layouts/Header'));
// const Information = lazy(() => import('../../components/Information/Information'));

async function fetch_pdb(PDBCode) { 
    if (PDBCode == null) {return}
    console.log("Fetching PDB ", PDBCode)
    let mtz_url = `https://edmaps.rcsb.org/coefficients/${PDBCode.toLowerCase()}.mtz`
    let pdb_url = `https://files.rcsb.org/download/${PDBCode.toUpperCase()}.pdb`

    let file = fetch(pdb_url)
    .then(response => {
        if (!response.ok) {
            throw new Error('Network error');
          }
        return response.text()
    })
    .then(file => {
        return Promise.resolve(file)
    })
    .catch((e) => {
        throw new Error("PDB Not Found")
    })
    return file
}

async function fetch_mtz(PDBCode) { 
    if (PDBCode == null) {return}
    console.log("Fetching MTZ ", PDBCode)
    let mtz_url = `https://edmaps.rcsb.org/coefficients/${PDBCode.toLowerCase()}.mtz`

    let file = fetch(mtz_url)
    .then(response => {
        if (!response.ok) {
            throw new Error('Network error');
          }
        return response.arrayBuffer()
    })
    .then(file => {
        return Promise.resolve(file)
    })
    .catch((e) => {
        throw new Error("PDB Not Found")
    })
    return file
}

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

    useEffect(() => {            
        if (PDBCode != "") {
            setLoadingText(`Fetching ${PDBCode.toUpperCase()} from the PDB`)
        
            fetch_mtz(PDBCode).then((response) => { 
                let array = new Uint8Array(response)
                setMtzData(array)
            })

            fetch_pdb(PDBCode).then((response) => { 
                setFileContent(response)
                setLoadingText("Validating Glycans...")

                privateer_module().then(async (Module) => {

                    let x = Module.read_structure_to_table(response, PDBCode)
    
                    let table_data = [];
                    for (var i = 0; i < x.size(); i++) {
                        let table_entry = x.get(i)
    
                        table_entry.id = sanitize_id(table_entry.id)
    
                        let collected_torsions = []
                        for(var j = 0; j < table_entry.torsions.size(); j++) { 
                            collected_torsions.push(table_entry.torsions.get(j)); 
                        }
                        table_entry.torsions = collected_torsions
                        table_data.push(table_entry)
    
                    }
                    
                    if (x.size() == 0 ) { 
                        setFailureText("Privateer could not detect any carbohydrates in this model.")
                        setLoadingText("There were no detected glycans in this file.")
                        setFallBack(true)
                    }
    
                    // Get Glyconnect ID from WURCS
                    setLoadingText("Querying Glytoucan...")
                    await loadGlytoucan(table_data)
    
    
                    setTableData(table_data);
                })
            }).catch(e => {
                setFailureText("This PDB code could not be found")
                setLoadingText("There were no detected glycans in this file.")
                setFallBack(true)})
           
         } else {
            privateer_module().then((Module) => {


                var coordinateReader = new FileReader();
                var reflectionReader = new FileReader();
                
                coordinateReader.onload = async () => {
         
                    setFileContent(coordinateReader.result)
    
                    let x = Module.read_structure_to_table(coordinateReader.result, coordinateFile.name)
    
                    let table_data = [];
                    for (var i = 0; i < x.size(); i++) {
                        let table_entry = x.get(i)
    
                        table_entry.id = sanitize_id(table_entry.id)
    
                        let collected_torsions = []
                        for(var j = 0; j < table_entry.torsions.size(); j++) { 
                            collected_torsions.push(table_entry.torsions.get(j)); 
                        }
                        table_entry.torsions = collected_torsions
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

    return (
        <>
            <Header resetApp={resetApp} setResetApp={setResetApp} 
            PDBCode={PDBCode} setPDBCode={setPDBCode}
            coordinateFile={coordinateFile} setCoordinateFile={setCoordinateFile}
            reflectionFile={reflectionFile} setReflectionFile={setReflectionFile}
             submit={submit} setSubmit={setSubmit}
            tableData={tableData} loadingText={loadingText} fileContent={fileContent} fallback={fallback} mtzData={mtzData}
            failureText={failureText}
            />
            <BorderElement topColor={"#D6D9E5"} bottomColor={"#F4F9FF"}></BorderElement>
            <Information/>
            <BorderElement topColor={"#F4F9FF"} bottomColor={"#D6D9E5"} reverse={true}></BorderElement>
            <Footer></Footer>
        </>
    )
}