import { useState, useEffect} from "react";
import Upload from "./Upload";
import Submit from "./Submit";
import Loading from "./Loading";
import privateer_module from "../wasm/privateer.js"
import SVGCarousel from "./SVGCarousel";


export default function SNFG() {

    const [file, setFile] = useState(null);
    const [submit, setSubmit] = useState(null);
    const [svgs, setSVGs] = useState(null)

    useEffect(() => {
        privateer_module().then((Module) => { 
            var reader = new FileReader();
            reader.onload = () => {
                console.log(reader.result)
                let x = Module.read_structure(reader.result, file.name)
                
                let svgs = [];
                for (var i = 0; i < x.size(); i++) {
                  svgs.push(x.get(i))
                }
                setSVGs(svgs);
            }
            if(file) {
                reader.readAsText(file);
            }
          })
    }, [submit])

    return (
        <div>
            {file == null ? <Upload setFile={setFile}/> : 
            submit == null ? <Submit file={file} submitPressed={setSubmit}/> :
            svgs == null ? <Loading/> : <SVGCarousel svgs={svgs} />}
        </div>
    )
}