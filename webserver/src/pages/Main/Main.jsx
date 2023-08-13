import SNFG from '../../components/SNFG'
import {useEffect, useState} from "react";
import Upload from "../../components/Upload";
import Submit from "../../common/Submit";
import Loading from "../../common/Loading";
import privateer_module from "../../wasm/privateer.js"
import loadGlytoucan from "../../utils/loadGlytoucan"

export default function Main({resetApp}) {

    const [file, setFile] = useState(null);
    const [fileContent, setFileContent] = useState(null)
    const [submit, setSubmit] = useState(null);
    // const [svgs, setSVGs] = useState(null);
    const [tableData, setTableData] = useState(null);
    const [loadingText, setLoadingText] = useState("Validating Glycans...");


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
            if (file) {
                reader.readAsText(file);
            }
        });
    }, [submit])

    useEffect(() => {
        setFile(null)
        setSubmit(null)
        setTableData(null)
    }, [resetApp])

    return (
        <div className="flex flex-grow bg-slate-700 justify-center items-center mb-5 my-5 p-5">
            {file == null ? <Upload setFile={setFile}/> :
                submit == null ? <Submit file={file} submitPressed={setSubmit}/> :
                    tableData == null ? <Loading loadingText={loadingText}/> :
                        <SNFG tableData={tableData} pdbString={fileContent}></SNFG>
            }
        </div>
    )
}