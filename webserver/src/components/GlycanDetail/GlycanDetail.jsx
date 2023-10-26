import {lazy, useCallback, useEffect, useState} from "react";
import {MoorhenContainer, MoorhenContextProvider} from 'moorhen'
import GlycanDetailInfoBox from "./GlycanDetailInfoBox";
// import GlycanDetailTable from "./GlycanDetailTable";
// import TorsionPlot from "../TorsionPlot/TorsionPlot";
import TorsionMultiPlot from "../TorsionPlot/TorsionMultiPlot";

const GlycanDetailTable = lazy(() => import('./GlycanDetailTable'));

export default function GlycanDetail({tableData, hideMoorhen, setHideMoorhen, rowID, forwardControls, scrollPosition, controls, molecule, map}) {

    async function handle_click(e) { 

        let new_center_string = e.target.dataset.chainid + "/" + e.target.dataset.seqnum + "(" + e.target.dataset.resname + ")"
        const selectedMolecule = controls.current.molecules.find((molecule) => molecule.name === "mol-1")
        await selectedMolecule.centreOn(new_center_string)
       
    }

    const ref = useCallback((node) => {
        if (node !== null) {

            let useList = node.querySelectorAll('use')

            for (let i = 0; i < useList.length; i++) {
                useList[i].addEventListener('click', handle_click)
            }
        }


        document.querySelectorAll("svg")[0].setAttribute("width", "50vw")
        document.querySelectorAll("svg")[0].setAttribute("height", "100%")
        
    })

    async function handleContourChange(e) { 
        map.contourLevel = Number(e.target.value)
        await map.doCootContour(
            ...map.glRef.current.origin.map(coord => -coord),
            map.mapRadius,
            map.contourLevel
            )
    }

    // useEffect( () => {
    //     function handleMoorhenResize() { 

    //         let wWidth = window.innerWidth;

    //         if (wWidth < 800) { 
    //             setWidth(wWidth - 200)
    //             setHeight((wWidth - 200)*0.75)
    //             console.log("Resizing window to ", wWidth - 100,  (wWidth - 100)*0.75)
    //         }

    //     }
    //     window.addEventListener('resize', handleMoorhenResize)

    //     return _ => {
    //         window.removeEventListener('resize', handleMoorhenResize)
    //     }
    // })

    const [width, setWidth] = useState(900);
    const [height, setHeight] = useState(600);
    const [torsionTab, setTorsionTab] = useState(0)

    return (
        <div className="flex-col justify-center items-center" style={{display: !hideMoorhen ? 'flex' : 'none'}}>

            <div className="flex items-center justify-center w-full">
                <div className="flex-1">
                    <button onClick={() => {
                        setHideMoorhen(true);
                        window.scrollTo(0, scrollPosition);
                        setTorsionTab(0)
                    }}>
                        <span className="">&#8592; Back To Table</span>
                    </button>
                </div>
                <h2 className="">Glycan Details</h2>
                <div class="flex-1"></div>
            </div>

            {/* <h4 className="my-4">Glycan ID: {tableData[rowID].id}</h4> */}
            <div className="my-5 w-full">
                <GlycanDetailInfoBox row={tableData[rowID]}/>
            </div>
            
            <h3 className="text-left text-xl w-full">SNFG</h3>

            <div className="text-sm text-center text-primary" >
                <div className="mt-4 py-4" id='svgContainer' dangerouslySetInnerHTML={{
                    __html: tableData[rowID].svg
                }} ref={ref}/>
                Hover over a linkage to see a summary
            </div>

            <h3 className="text-left text-xl w-full">Visualise</h3>

            <label for="contour-range-text" class="block mb-2 text-sm font-medium text-gray-900 dark:text-white translate-y-10">Map Contour</label>
            <input id="contour-range" type="range" min="0" max="1" step="0.05" className="w-36 h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer dark:bg-gray-700 translate-y-10" onChange={handleContourChange}/>

            <MoorhenContextProvider defaultBackgroundColor={[51, 65, 85, 1]}>

                <MoorhenContainer forwardControls={forwardControls} setMoorhenDimensions={() => {
                    return [width, height];
                }} viewOnly={true}/>

            </MoorhenContextProvider>
            
            <h3 className="text-left text-xl w-full">Torsion Plots</h3>

            <TorsionMultiPlot torsions={tableData[rowID].torsions} tab={torsionTab} setTab={setTorsionTab}/>

        </div>);
}
