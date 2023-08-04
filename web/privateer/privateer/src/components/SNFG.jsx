import { useState, useEffect, useRef, useCallback } from "react";
import SVGTable from "./SVGTable";
import { MoorhenContextProvider, MoorhenMolecule, MoorhenContainer, itemReducer } from 'moorhen'
import GlycanDetail from "./GlycanDetail"

const initialMoleculesState = []

export default function SNFG({ tableData, pdbString }) {

    const [rowClicked, setRowClicked] = useState(false)
    const [rowID, setRowID] = useState(0)
    const [hideMoorhen, setHideMoorhen] = useState(true)

    const [cootInitialized, setCootInitialized] = useState(false)
    const controls = useRef()
    
    const forwardControls = (forwardedControls) => {
        setCootInitialized(true)
        controls.current = forwardedControls
    }

    const [yScrollPosition, setYScrollPosition] = useState(0)


    // DEBUG ONLY 
    useEffect(() => {
        if (cootInitialized && controls.current) {
            // whatever you want to do with moorhen has to wait for cootInitialized to be true
            console.log("COOT LOADED")
        }
    }, [cootInitialized, controls.current])


    useEffect(() => {
        setYScrollPosition(window.scrollY)
        console.log(tableData[rowID].id)

        if (cootInitialized) {
            let newMolecule = new MoorhenMolecule(controls.current.commandCentre, controls.current.glRef, controls.current.monomerLibrary)
            newMolecule.loadToCootFromString(pdbString, 'mol-1').then(() => {
                controls.current.changeMolecules( { action: 'Add', item: newMolecule } );
                newMolecule.fetchIfDirtyAndDraw('CBs').then(
                    newMolecule.centreOn("B/NAG")
                )
            })
            window.scrollTo(0,0)
            setHideMoorhen(false)
        }
    }, [rowClicked])

    return (
        <>
            <div style={{ display: (hideMoorhen ? 'block' : 'none') }} id="tableContainer">
                <SVGTable tableData={tableData} rowClick={rowClicked} setRowClicked={setRowClicked} setRowID={setRowID} />
            </div>
            <GlycanDetail tableData={tableData}
            hideMoorhen={hideMoorhen} 
            setHideMoorhen={setHideMoorhen}
            rowID={rowID} 
            forwardControls={forwardControls} 
            scrollPosition={yScrollPosition}
            />

        </>
    )

    }