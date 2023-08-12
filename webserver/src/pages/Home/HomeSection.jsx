import { Header } from '../../layouts/Header';
import { Information } from '../../components/Information/Information';

import { useState, useReducer, useRef, useEffect} from "react";
import privateer_module from "../../wasm/privateer.js"
import loadGlytoucan from "../../utils/loadGlytoucan"
import Footer from "../../layouts/Footer"
import BorderElement from '../../layouts/BorderElement';

export default function HomeSection() {
    const [file, setFile] = useState(null);
    const [fileContent, setFileContent] = useState(null)
    const [submit, setSubmit] = useState(null);
    const [tableData, setTableData] = useState(null);
    const [loadingText, setLoadingText] = useState("Validating Glycans...");
    const [resetApp, setResetApp] = useState(false)

    

    useEffect(() => {
        privateer_module().then((Module) => { 
            var reader = new FileReader();
            reader.onload = async () => {
                // let x = Module.read_structure(reader.result, file.name)
                
                // let svgs = [];
                // for (var i = 0; i < x.size(); i++) {
                //   svgs.push(x.get(i))
                // }
                // setSVGs(svgs);
                setFileContent(reader.result)
                let x = Module.read_structure_to_table(reader.result, file.name)
                
                let table_data = [];
                for (var i = 0; i < x.size(); i++) {
                    table_data.push(x.get(i))
                }
                // Get Glyconnect ID from WURCS
                setLoadingText("Querying Glytoucan...")
                await loadGlytoucan(table_data)

                // setLoadingText("Querying GlyConnect...")
                // await load_glyconnect(table_data)

                setTableData(table_data);
            }

            if(file) {
                reader.readAsText(file);
            }
          }).catch((e) => console.log(e));
    }, [submit])

    useEffect(() => {
        setFile(null)
        setSubmit(null)
        setTableData(null)
    }, [resetApp])

    return (
        <>
            <Header setResetApp={setResetApp} file={file} setFile={setFile} submit={submit} setSubmit={setSubmit} tableData={tableData} loadingText={loadingText} fileContent={fileContent}  />
            <BorderElement topColor={"#D6D9E5"} bottomColor={"#F4F9FF"}></BorderElement>
            <Information/>
            <BorderElement topColor={"#F4F9FF"} bottomColor={"#D6D9E5"} reverse={true}></BorderElement>
            <Footer></Footer>
        </>
    )
}