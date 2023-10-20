import {lazy, useCallback, useEffect, useState} from "react";
import {MoorhenContainer, MoorhenContextProvider} from 'moorhen'
import GlycanDetailInfoBox from "./GlycanDetailInfoBox";
// import GlycanDetailTable from "./GlycanDetailTable";
// import TorsionPlot from "../TorsionPlot/TorsionPlot";
import TorsionMultiPlot from "../TorsionPlot/TorsionMultiPlot";

const GlycanDetailTable = lazy(() => import('./GlycanDetailTable'));

export default function GlycanDetail({tableData, hideMoorhen, setHideMoorhen, rowID, forwardControls, scrollPosition}) {

    const ref = useCallback((node) => {
        let useList = document.querySelectorAll('use')

        for (let i = 0; i < useList.length; i++) {
            useList[i].addEventListener("click", (e) => {
                console.log(useList[i].id)
            })
        }

        document.querySelectorAll("svg")[0].setAttribute("width", "50vw")
        document.querySelectorAll("svg")[0].setAttribute("height", "100%")
    })

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

    const [width, setWidth] = useState(1200);
    const [height, setHeight] = useState(1000);
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

            <MoorhenContextProvider defaultBackgroundColor={[51, 65, 85, 1]}>
                <MoorhenContainer forwardControls={forwardControls} setMoorhenDimensions={() => {
                    return [width, height];
                }} viewOnly={false}/>

            </MoorhenContextProvider>
            
            <h3 className="text-left text-xl w-full">Torsion Plots</h3>

            <TorsionMultiPlot torsions={tableData[rowID].torsions} tab={torsionTab} setTab={setTorsionTab}/>

        </div>);
}
