import { useEffect, useState } from "react";
import UploadButton from "./UploadButton"
import Submit from "./Submit";

export default function Upload({coordinateFile, setCoordinateFile, reflectionFile, setReflectionFile, submitPressed, setResetApp}) {

    const [showUploadAgain, setShowUploadAgain] = useState(true)
    const [showSubmit, setShowSubmit] = useState(false)
    const [allowSubmit, setAllowSubmit] = useState(false)

    useEffect(() => { 
        if (coordinateFile && reflectionFile) {
            console.log("Ready for submit")
            setShowSubmit(true)
            setShowUploadAgain(false)
            setAllowSubmit(true)
        }

        if (coordinateFile && !reflectionFile) { 
            setShowSubmit(true)
            setShowUploadAgain(true)
            setAllowSubmit(true)
        }

        if (!coordinateFile && reflectionFile) { 
            setShowSubmit(true)
            setShowUploadAgain(true)
            setAllowSubmit(false)
        } 

        if (!coordinateFile && !reflectionFile) { 
            setShowSubmit(false)
            setShowUploadAgain(true)
            setAllowSubmit(false)

        }

        console.log(coordinateFile, reflectionFile)
    }, [coordinateFile, reflectionFile])

    return (
        <>
            
            { showUploadAgain == true ? <UploadButton setCoordinateFile={setCoordinateFile} setReflectionFile={setReflectionFile}/>: <></>}
            {showSubmit == true ? 
             <Submit coordinateFile={coordinateFile} reflectionFile={reflectionFile} submitPressed={submitPressed} setResetApp={setResetApp} allowSubmit={allowSubmit}/> : <></>}
            
        </>
    )
}