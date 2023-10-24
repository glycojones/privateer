import { useEffect, useState } from "react";
import UploadButton from "./UploadButton"
import Submit from "./Submit";
import PDBFetch from "./PDBFetch";

export default function Upload({coordinateFile, setCoordinateFile, 
    reflectionFile, setReflectionFile, 
    PDBCode, setPDBCode, submitPressed, resetApp, setResetApp}) {

    const [showUploadAgain, setShowUploadAgain] = useState(true)
    const [showSubmit, setShowSubmit] = useState(false)
    const [allowSubmit, setAllowSubmit] = useState(false)
    const [showPDBFetch, setShowPDBFetch] = useState(true)

    useEffect(() => {
        setShowPDBFetch(true)
    }, [resetApp])

    useEffect(() => { 
        if (coordinateFile && reflectionFile) {
            setShowSubmit(true)
            setShowUploadAgain(false)
            setAllowSubmit(true)
            setShowPDBFetch(false)
        }

        if (coordinateFile && !reflectionFile) { 
            setShowSubmit(true)
            setShowUploadAgain(true)
            setAllowSubmit(true)
            setShowPDBFetch(false)
        }

        if (!coordinateFile && reflectionFile) { 
            setShowSubmit(true)
            setShowUploadAgain(true)
            setAllowSubmit(false)
            setShowPDBFetch(false)
        } 

        if (!coordinateFile && !reflectionFile) { 
            setShowSubmit(false)
            setShowUploadAgain(true)
            setAllowSubmit(false)
            setShowPDBFetch(true)
        }
    }, [coordinateFile, reflectionFile])

    return (
        <div className="flex flex-wrap align-middle items-center justify-center">
            { showUploadAgain == true ? <UploadButton setCoordinateFile={setCoordinateFile} setReflectionFile={setReflectionFile}/>: <></>}
            {showSubmit == true ? 
             <Submit coordinateFile={coordinateFile} reflectionFile={reflectionFile} submitPressed={submitPressed} setResetApp={setResetApp} allowSubmit={allowSubmit}/> : <></>}
            
            {showPDBFetch == true ? 
             <div className="mx-6 w-full lg:w-6 sm:w-full text-center">OR</div>
             :
              <></>} 
            {showPDBFetch == true ? 
             <PDBFetch PDBCode={PDBCode} setPDBCode={setPDBCode} submitPressed={submitPressed}/> 
             :
              <></>} 
            
        </div>
    )
}