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


export default function HomeSection() {
    const [coordinateFile, setCoordinateFile] = useState(null);
    const [reflectionFile, setReflectionFile] = useState(null);

    const [fileContent, setFileContent] = useState(null)
    const [submit, setSubmit] = useState(null);
    const [tableData, setTableData] = useState(null);
    const [loadingText, setLoadingText] = useState("Validating Glycans...");
    const [resetApp, setResetApp] = useState(false)
    const [fallback, setFallBack] = useState(false)

    let sanitize_id = (id) => { 
        const regex = /: *32/g;
        const new_id = id.replace(regex, "")
        return new_id
    }

    useEffect(() => {
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

            if (reflectionFile) { 
                reflectionReader.readAsArrayBuffer(reflectionFile)
                let map_data = new Uint8Array(reflectionReader.result);
                Module['FS_createDataFile']('/', "input.mtz", map_data, true, true, true)
            }

        }).catch((e) => console.log(e));
    }, [submit])

    useEffect(() => {
        
        setReflectionFile(null)
        setCoordinateFile(null)
        setSubmit(null)
        setTableData(null)
        setFallBack(false)
        setResetApp(false)
    }, [resetApp])

    return (
        <>
            <Header setResetApp={setResetApp} coordinateFile={coordinateFile} setCoordinateFile={setCoordinateFile}
            reflectionFile={reflectionFile} setReflectionFile={setReflectionFile}
             submit={submit} setSubmit={setSubmit}
            tableData={tableData} loadingText={loadingText} fileContent={fileContent} fallback={fallback}/>
            <BorderElement topColor={"#D6D9E5"} bottomColor={"#F4F9FF"}></BorderElement>
            <Information/>
            <BorderElement topColor={"#F4F9FF"} bottomColor={"#D6D9E5"} reverse={true}></BorderElement>
            <Footer></Footer>
        </>
    )
}