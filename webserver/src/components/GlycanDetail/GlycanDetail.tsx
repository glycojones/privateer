import {lazy, useCallback, useEffect, useState} from "react";
import {MoorhenContainer, MoorhenContextProvider} from 'moorhen'
import GlycanDetailInfoBox from "./GlycanDetailInfoBox";
import TorsionMultiPlot from "../TorsionPlot/TorsionMultiPlot.tsx";
import {GlycanDetailProps} from "../../interfaces/types"


export default function GlycanDetail(props: GlycanDetailProps) {

    async function handle_click(e) { 

        let new_center_string = e.target.dataset.chainid + "/" + e.target.dataset.seqnum + "(" + e.target.dataset.resname + ")"
        const selectedMolecule = props.controls.current.molecules.find((molecule) => molecule.name === "mol-1")
        await selectedMolecule.centreOn(new_center_string)
       
    }

    const ref = useCallback((node: any) => {
        if (node !== null) {

            let useList = node.querySelectorAll('use')

            for (let i = 0; i < useList.length; i++) {
                useList[i].addEventListener('click', handle_click)
            }
        }

        document.querySelectorAll("svg")[0].setAttribute("width", "50vw")
        document.querySelectorAll("svg")[0].setAttribute("height", "100%")
        
    }, [])

    async function handleContourChange(e) { 
        props.map.contourLevel = Number(e.target.value)
        await props.map.doCootContour(
            ...props.map.glRef.current.origin.map(coord => -coord),
            props.map.mapRadius,
            props.map.contourLevel
            )
    }

    const [width, setWidth] = useState(900);
    const [height, setHeight] = useState(600);
    const [torsionTab, setTorsionTab] = useState(0)

    return (
        <div className="flex-col justify-center items-center" style={{display: !props.hideMoorhen ? 'flex' : 'none'}}>

            <div className="flex items-center justify-center w-full">
                <div className="flex-1">
                    <button onClick={() => {
                        props.setHideMoorhen(true);
                        window.scrollTo(0, props.scrollPosition);
                        setTorsionTab(0)
                    }}>
                        <span className="">&#8592; Back To Table</span>
                    </button>
                </div>
                <h2 className="">Glycan Details</h2>
                <div className="flex-1"></div>
            </div>

            <div className="my-5 w-full">
                <GlycanDetailInfoBox row={props.tableData[props.rowID]}/>
            </div>
            
            <h3 className="text-left text-xl w-full">SNFG</h3>

            <div className="text-sm text-center text-primary" >
                <div className="mt-4 py-4" id='svgContainer' dangerouslySetInnerHTML={{
                    __html: props.tableData[props.rowID].svg
                }} ref={ref}/>
                Hover over a linkage to see a summary
            </div>

            <h3 className="text-left text-xl w-full">Visualise</h3>

            <label htmlFor="contour-range-text" className="block mb-2 text-sm font-medium text-gray-900 dark:text-white translate-y-10">Map Contour</label>
            <input id="contour-range" type="range" min="0" max="1" step="0.05" defaultValue="0.2" className="w-36 h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer dark:bg-gray-700 translate-y-10" onChange={handleContourChange}/>

            <MoorhenContextProvider defaultBackgroundColor={[51, 65, 85, 1]}>

                <MoorhenContainer forwardControls={props.forwardControls} setMoorhenDimensions={() => {
                    return [width, height];
                }} viewOnly={true}/>

            </MoorhenContextProvider>
            
            <h3 className="text-left text-xl w-full">Torsion Plots</h3>

            <TorsionMultiPlot torsions={props.tableData[props.rowID].torsions} tab={torsionTab} setTab={setTorsionTab}/>

        </div>);
}
