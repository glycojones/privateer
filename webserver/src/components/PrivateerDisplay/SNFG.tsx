import {lazy, useEffect, useRef, useState} from "react";
import {MoorhenMap, MoorhenMolecule} from 'moorhen'
import GlycanDetail from "../GlycanDetail/GlycanDetail.tsx"
import {SNFGProps} from "../../interfaces/types"

const SVGTable = lazy(() => import('../SVGTable/SVGTable.tsx'));

export default function SNFG (props: SNFGProps) {
    const [rowClicked, setRowClicked] = useState(false)
    const [rowID, setRowID] = useState(0)
    const [hideMoorhen, setHideMoorhen] = useState(true)
    const [allowRowClick, setAllowRowClick] = useState(false)

    const [dataLoaded, setDataLoaded] = useState(false)

    const [cootInitialized, setCootInitialized] = useState(false)
    const controls = useRef()
    const [molecule, setMolecule] = useState()
    const [map, setMap] = useState()
    const forwardControls = (forwardedControls) => {
        setCootInitialized(true)
        controls.current = forwardedControls
    }

    const [yScrollPosition, setYScrollPosition] = useState(0)

    useEffect(() => {
        async function load_map_and_model() {
        if (cootInitialized && controls.current && !dataLoaded) {
            setAllowRowClick(true)
            let newMolecule = new MoorhenMolecule(controls.current.commandCentre, controls.current.glRef, controls.current.monomerLibrary)
            newMolecule.loadToCootFromString(props.fileContent, 'mol-1').then(() => {
                controls.current.changeMolecules({action: 'Add', item: newMolecule});
                newMolecule.fetchIfDirtyAndDraw('CBs').then(() => {
                        let id = props.tableData[rowID].id

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
            if (props.PDBCode == "") { 
                await map.loadToCootFromMtzData(props.mtzData, "map-1", mapMetadata)
            }
            else { 
                await map.loadToCootFromMapData(props.mtzData, "map-1", false)

            }
            map.suggestedContourLevel = 0.3
            console.log("Suggested level", map.suggestedContourLevel, map)
            controls.current.changeMaps({ action: "Add", item: map })
            controls.current.setActiveMap(map)
            setMap(map)
            setDataLoaded(true)
        }
    }
        load_map_and_model()
    }, [cootInitialized])


    useEffect(() => {
        async function move_view() { 
            if (!cootInitialized) {return}
       
            setYScrollPosition(window.scrollY)
            let id = props.tableData[rowID].id
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

    let glycanDetailProps = {
        tableData:props.tableData,
        hideMoorhen:hideMoorhen,
        setHideMoorhen:setHideMoorhen,
        rowID:rowID,
        forwardControls:forwardControls,
        scrollPosition:yScrollPosition,
        controls:controls,
        molecule:molecule,
        map:map
    }

    let svgTableProps = {
        tableData:props.tableData,
        allowRowClick:allowRowClick,
        rowClick:rowClicked,
        setRowClicked:setRowClicked,
        setRowID:setRowID
    }

    return (
        <div className="flex flex-col">
            <div style={{display: (hideMoorhen ? 'block' : 'none')}} id="tableContainer">
                <div className="flex flex-col ">
                    <h2 className="my-4 text-lg sm:text-2xl text-center sm:text-left">Detected {props.tableData.length} Glycans
                        in {props.filename}</h2>
                    <SVGTable {...svgTableProps}/>
                </div>

            </div>
            <GlycanDetail {...glycanDetailProps}/>

        </div>
    )

}