import SNFG from './SNFG'
import { useState, useReducer, useRef, useEffect} from "react";
import Upload from "./Upload";
import Submit from "./Submit";
import Loading from "./Loading";
import privateer_module from "../wasm/privateer.js"

export default function Main() { 

    const [file, setFile] = useState(null);
    const [fileContent, setFileContent] = useState(null)
    const [submit, setSubmit] = useState(null);
    const [svgs, setSVGs] = useState(null);
    const [tableData, setTableData] = useState(null);
    const [loadingText, setLoadingText] = useState("Validating Glycans...");

    async function get_glytoucan_id(wurcs) {
        let url = "https://api.glycosmos.org/sparqlist/wurcs2gtcids?wurcs="+encodeURIComponent(wurcs)

        return fetch(url, {
            method: "GET" // default, so we can ignore
        })  
        .then((response) => response.json())
        .catch((error) => console.log(error))
    }

    async function load_glytoucan(table_data) {

        let promises = [];

        for (var i = 0; i < table_data.length; i++){
            let item = table_data[i]
            // const id = get_glyconnect_id(item.wurcs)
            promises.push(get_glytoucan_id(item.wurcs))
        }
        
        const data = await Promise.all(promises)

        data.forEach((data, index) => { 
            table_data[index].glytoucan_id = data[0].id;
        })
        return table_data;
    }

    async function get_glyconnect_id(glytoucan_id) {
        let url = "https://glyconnect.expasy.org/api/structures/search/glytoucan"

        console.log(JSON.stringify({"glytoucan_id": glytoucan_id}))

        fetch(url, {
            method: "POST",
            body: JSON.stringify({"glytoucan_id": glytoucan_id}),
            headers: {
                "Content-type": "application/json; charset=UTF-8"
            }
        }).then((response) => {
            const contentType = response.headers.get("content-type")
            if (contentType && contentType.indexOf("application/json") !== -1) { 
                return response.json()
            }
            else if (response.status == 404) { 
                return Promise.reject("404 Error")
            }
        }).catch((error) => {
            return Promise.reject(error)
            // console.log(error)
        })
    }

    async function load_glyconnect(table_data) {

        let promises = [];

        for (var i = 0; i < table_data.length; i++){
            let item = table_data[i]
            promises.push(get_glyconnect_id(item.glytoucan_id))
        }
        
        try {
            const data = await Promise.all(promises)
        }
        catch { 
            
        }
        data.forEach((data, index) => { 
            // console.log(data)
            // table_data[index].glytoucan_id = data[0].id;
        })
        return table_data;
    }
    
    

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
                // setLoadingText("Querying Glytoucan...")
                // await load_glytoucan(table_data)

                // setLoadingText("Querying GlyConnect...")
                // await load_glyconnect(table_data)

                setTableData(table_data);
            }
            if(file) {
                reader.readAsText(file);
            }
          });
    }, [submit])


    return (
        <div className="flex flex-grow bg-slate-700 justify-center items-center mb-5 my-5 p-5 ">
            {file == null ? <Upload setFile={setFile}/> : 
            submit == null ? <Submit file={file} submitPressed={setSubmit}/> :
            tableData == null ? <Loading loadingText={loadingText}/> : 
            <SNFG tableData = {tableData} pdbString={fileContent}></SNFG>
            
            }
        </div>
    )
}