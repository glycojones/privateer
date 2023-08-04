import { useState, useEffect, useRef, useCallback } from "react";
import SVGTable from "./SVGTable";
import { MoorhenContextProvider, MoorhenMolecule, MoorhenContainer, itemReducer } from 'moorhen'


export default function GlycanDetail({ tableData, hideMoorhen, setHideMoorhen, rowID, forwardControls, scrollPosition}) {

    const ref = useCallback((node) => {
        let useList = document.querySelectorAll('use')

        for (let i = 0; i < useList.length; i++) {
            useList[i].addEventListener("click", (e) => {console.log(useList[i].id)})
        }

        document.querySelectorAll("svg")[0].setAttribute("width", "50vw")
        document.querySelectorAll("svg")[0].setAttribute("height", "100%")
    })

    return (
    <div className="flex-col justify-center items-center" style={{ display: !hideMoorhen ? 'flex' : 'none' }}>

        <div className="flex flex-row w-full justify-between">
        <button onClick={() => { 
            setHideMoorhen(true);
            window.scrollTo(0, scrollPosition);
            
            }}><span>&#8592;</span></button>
        <h2>{tableData[rowID].id}</h2>
        </div>
            

        <div className="mt-4 py-4" id='svgContainer' dangerouslySetInnerHTML={{
            __html: tableData[rowID].svg
        }} ref={ref} />

        <MoorhenContextProvider defaultBackgroundColor={[51, 65, 85, 1]}>
            <MoorhenContainer forwardControls={forwardControls} setMoorhenDimensions={() => {
                return [800, 600];
            }} viewOnly={true} />

        </MoorhenContextProvider>
    </div>);
}
