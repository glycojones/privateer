import { lazy, useEffect, useRef, useState } from "react";
import { MoorhenMap, MoorhenMolecule } from 'moorhen'
import { SNFGProps } from "../../interfaces/types"

const SVGTable = lazy(() => import('../SVGTable/SVGTable.tsx'));
const GlycanDetail = lazy(() => import('../GlycanDetail/GlycanDetail.tsx'));

export default function SNFG(props: SNFGProps) {
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
                    controls.current.changeMolecules({ action: 'Add', item: newMolecule });
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
            if (!cootInitialized) { return }

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


    function save_snfgs() {
        if (!props.tableData) return 

        let data = props.tableData

        for(let i = 0; i < data.length; i++) { 
            var svgBlob = new Blob([data[i].svg], {type:"image/svg+xml;charset=utf-8"});
            var svgUrl = URL.createObjectURL(svgBlob);
            var downloadLink = document.createElement("a");
            downloadLink.href = svgUrl;
            downloadLink.download = `${data[i].id}.svg`;
            document.body.appendChild(downloadLink);
            downloadLink.click();
            document.body.removeChild(downloadLink);
            
        }

    }


    let glycanDetailProps = {
        tableData: props.tableData,
        hideMoorhen: hideMoorhen,
        setHideMoorhen: setHideMoorhen,
        rowID: rowID,
        forwardControls: forwardControls,
        scrollPosition: yScrollPosition,
        controls: controls,
        molecule: molecule,
        map: map
    }

    let svgTableProps = {
        tableData: props.tableData,
        allowRowClick: allowRowClick,
        rowClick: rowClicked,
        setRowClicked: setRowClicked,
        setRowID: setRowID
    }

    return (
        <div className="flex flex-col">
            <div style={{ display: (hideMoorhen ? 'block' : 'none') }} id="tableContainer">
                <div className="flex flex-col ">
                    <div className="flex w-full justify-between items-center">
                        <h2 className="my-4 text-lg sm:text-2xl text-center sm:text-left">Detected {props.tableData.length} Glycans
                            in {props.filename}</h2>

                        <button className="text-center items-center mr-4" onClick={() => { save_snfgs() }} title="Download All SNFGs">
                            <svg className="h-6 w-6" xmlns="http://www.w3.org/2000/svg" height="1em" viewBox="0 0 512 512">
                                <path d="M288 32c0-17.7-14.3-32-32-32s-32 14.3-32 32V274.7l-73.4-73.4c-12.5-12.5-32.8-12.5-45.3 0s-12.5 32.8 0 45.3l128 128c12.5 12.5 32.8 12.5 45.3 0l128-128c12.5-12.5 12.5-32.8 0-45.3s-32.8-12.5-45.3 0L288 274.7V32zM64 352c-35.3 0-64 28.7-64 64v32c0 35.3 28.7 64 64 64H448c35.3 0 64-28.7 64-64V416c0-35.3-28.7-64-64-64H346.5l-45.3 45.3c-25 25-65.5 25-90.5 0L165.5 352H64zm368 56a24 24 0 1 1 0 48 24 24 0 1 1 0-48z" /></svg>
                        </button>
                    </div>
                    <SVGTable {...svgTableProps} />
                </div>

            </div>
            <GlycanDetail {...glycanDetailProps} />

        </div>
    )

}