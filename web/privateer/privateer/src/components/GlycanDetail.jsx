import { useState, useEffect, useRef, useCallback } from "react";
import SVGTable from "./SVGTable";
import { MoorhenContextProvider, MoorhenMolecule, MoorhenContainer, itemReducer } from 'moorhen'
// import BarChart from "./BarChart";
import { useTable } from 'react-table';
import GlycanDetailTable from "./GlycanDetailTable";


export default function GlycanDetail({ tableData, hideMoorhen, setHideMoorhen, rowID, forwardControls, scrollPosition }) {

    const ref = useCallback((node) => {
        let useList = document.querySelectorAll('use')

        for (let i = 0; i < useList.length; i++) {
            useList[i].addEventListener("click", (e) => { console.log(useList[i].id) })
        }

        document.querySelectorAll("svg")[0].setAttribute("width", "50vw")
        document.querySelectorAll("svg")[0].setAttribute("height", "100%")
    })

    const [width, setWidth] = useState(800);
    const [height, setHeight] = useState(600);

    return (
        <div className="flex-col justify-center items-center" style={{ display: !hideMoorhen ? 'flex' : 'none' }}>

            <div className="flex items-center justify-center w-full">
                <div className="flex-1">
                    <button onClick={() => {
                        setHideMoorhen(true);
                        window.scrollTo(0, scrollPosition);
                    }}>
                        <span className="">&#8592; Back To Table</span>
                    </button>
                </div>
                <h2 className="">Glycan Details</h2>
                <div class="flex-1"></div>
            </div>

            {/* <h4 className="my-4">Glycan ID: {tableData[rowID].id}</h4> */}
            <div className="my-5">
                <GlycanDetailTable row={tableData[rowID]}/>
            </div>

            <div className="mt-4 py-4" id='svgContainer' dangerouslySetInnerHTML={{
                __html: tableData[rowID].svg
            }} ref={ref} />

            <MoorhenContextProvider defaultBackgroundColor={[51, 65, 85, 1]}>
                <MoorhenContainer forwardControls={forwardControls} setMoorhenDimensions={() => {
                    return [width, height];
                }} viewOnly={true} />

            </MoorhenContextProvider>


        </div>);
}
