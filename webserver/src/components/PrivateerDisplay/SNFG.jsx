import {lazy, useEffect, useRef, useState} from "react";
// import SVGTable from "../SVGTable/SVGTable";
import {MoorhenMap, MoorhenMolecule} from 'moorhen'
import GlycanDetail from "../GlycanDetail/GlycanDetail"

const SVGTable = lazy(() => import('../SVGTable/SVGTable'));

const initialMoleculesState = []

export default function SNFG({tableData, fileName, pdbString, mtzData}) {

    const [rowClicked, setRowClicked] = useState(false)
    const [rowID, setRowID] = useState(0)
    const [hideMoorhen, setHideMoorhen] = useState(true)
    const [allowRowClick, setAllowRowClick] = useState(false)

    const [dataLoaded, setDataLoaded] = useState(false)

    const [cootInitialized, setCootInitialized] = useState(false)
    const controls = useRef()
    const [molecule, setMolecule] = useState()
    const forwardControls = (forwardedControls) => {
        setCootInitialized(true)
        controls.current = forwardedControls
    }

    const [yScrollPosition, setYScrollPosition] = useState(0)

    // DEBUG ONLY 
    useEffect(() => {
        if (cootInitialized && controls.current && !dataLoaded) {
            // whatever you want to do with moorhen has to wait for cootInitialized to be true
            setAllowRowClick(true)
            let newMolecule = new MoorhenMolecule(controls.current.commandCentre, controls.current.glRef, controls.current.monomerLibrary)
            newMolecule.loadToCootFromString(pdbString, 'mol-1').then(() => {
                controls.current.changeMolecules({action: 'Add', item: newMolecule});
                newMolecule.fetchIfDirtyAndDraw('CBs').then(() => {
                        let id = tableData[rowID].id

                        let sugar_name = id.split("-")[0]
                        let sugar_id = id.split("-")[1].split("/")[0].split(":")[0]
                        let sugar_chain = id.split("/")[1].split("_")[0]

                        let center_string = sugar_chain + "/" + sugar_id + "(" + sugar_name + ")"

                        newMolecule.centreOn(center_string)
                        setMolecule(newMolecule)
                    }
                )
            })
            

            const map = new MoorhenMap(controls.current.commandCentre, controls.current.glRef);
            const mapMetadata = {
                F: "FWT",
                PHI: "PHWT",
                Fobs: "FP",
                SigFobs: "SIGFP",
                FreeR: "FREE",
                isDifference: false,
                useWeight: false,
                calcStructFact: true,
            }
            map.loadToCootFromMtzData(mtzData, "map-1", mapMetadata).then(() => { 
                    controls.current.changeMaps({ action: "Add", item: map })
                    controls.current.setActiveMap(map)
            });
            setDataLoaded(true)
        }
    }, [cootInitialized])


    useEffect(() => {
        async function move_view() { 
            if (!cootInitialized) {return}

        
            setYScrollPosition(window.scrollY)
            let id = tableData[rowID].id
            let sugar_name = id.split("-")[0]
            let sugar_id = id.split("-")[1].split("/")[0].split(":")[0]
            let sugar_chain = id.split("/")[1].split("_")[0]
    
            let center_string = sugar_chain + "/" + sugar_id + "(" + sugar_name + ")"
            const selectedMolecule = controls.current.molecules.find((molecule) => molecule.name === "mol-1")
            await selectedMolecule.centreOn(center_string)
            window.scrollTo(0, 0)
            setHideMoorhen(false)
        }
        move_view()
    }, [rowClicked])

    return (
        <div className="flex flex-col">
            <div style={{display: (hideMoorhen ? 'block' : 'none')}} id="tableContainer">
                <div className="flex flex-col ">
                    <h2 className="my-4 text-lg sm:text-2xl text-center sm:text-left">Detected {tableData.length} Glycans
                        in {fileName}</h2>
                    <SVGTable tableData={tableData} allowRowClick={allowRowClick} rowClick={rowClicked} setRowClicked={setRowClicked}
                              setRowID={setRowID}/>
                </div>

            </div>
            <GlycanDetail tableData={tableData}
                          hideMoorhen={hideMoorhen}
                          setHideMoorhen={setHideMoorhen}
                          rowID={rowID}
                          forwardControls={forwardControls}
                          scrollPosition={yScrollPosition}
                          controls={controls}
                          molecule={molecule}
            />

        </div>
    )

}