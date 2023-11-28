import { lazy, useCallback, useEffect, useState } from "react";
import { MoorhenContainer, MoorhenContextProvider } from 'moorhen'

import { GlycanDetailProps } from "../../interfaces/types"
import { type } from "os";
const GlycanDetailInfoBox = lazy(() => import('./GlycanDetailInfoBox'));
const TorsionMultiPlot = lazy(() => import('../TorsionPlot/TorsionMultiPlot.tsx'));


export default function GlycanDetail(props: GlycanDetailProps) {

    async function handle_click(e) {

        let new_center_string = e.target.dataset.chainid + "/" + e.target.dataset.seqnum + "(" + e.target.dataset.resname + ")"
        const selectedMolecule = props.controls.current.molecules.find((molecule) => molecule.name === "mol-1")
        await selectedMolecule.centreOn(new_center_string)

    }

    const ref = useCallback((node: HTMLElement | null) => {
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

    function save_snfg() {
        if (!props.tableData[props.rowID]) return

        var svgBlob = new Blob([props.tableData[props.rowID].svg], { type: "image/svg+xml;charset=utf-8" });
        var svgUrl = URL.createObjectURL(svgBlob);
        var downloadLink = document.createElement("a");
        downloadLink.href = svgUrl;
        downloadLink.download = `${props.tableData[props.rowID].id}.svg`;
        document.body.appendChild(downloadLink);
        downloadLink.click();
        document.body.removeChild(downloadLink);[]
    }

    const [width, setWidth] = useState(900);
    const [height, setHeight] = useState(600);
    const [torsionTab, setTorsionTab] = useState(0)

    return (
        <div className="flex-col justify-center items-center" style={{ display: !props.hideMoorhen ? 'flex' : 'none' }}>

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
                <GlycanDetailInfoBox row={props.tableData[props.rowID]} />
            </div>

            <div className="flex justify-between w-full">
                <h3 className="text-left text-xl w-full">SNFG</h3>
                <button className="text-center items-center mr-4" onClick={() => { save_snfg() }} title="Download SNFG">
                    <svg className="h-6 w-6" xmlns="http://www.w3.org/2000/svg" height="1em" viewBox="0 0 512 512">
                        <path d="M288 32c0-17.7-14.3-32-32-32s-32 14.3-32 32V274.7l-73.4-73.4c-12.5-12.5-32.8-12.5-45.3 0s-12.5 32.8 0 45.3l128 128c12.5 12.5 32.8 12.5 45.3 0l128-128c12.5-12.5 12.5-32.8 0-45.3s-32.8-12.5-45.3 0L288 274.7V32zM64 352c-35.3 0-64 28.7-64 64v32c0 35.3 28.7 64 64 64H448c35.3 0 64-28.7 64-64V416c0-35.3-28.7-64-64-64H346.5l-45.3 45.3c-25 25-65.5 25-90.5 0L165.5 352H64zm368 56a24 24 0 1 1 0 48 24 24 0 1 1 0-48z" /></svg>
                </button>
            </div>

            <div className="text-sm text-center text-primary" >
                <div className="mt-4 py-4" id='svgContainer' dangerouslySetInnerHTML={{
                    __html: props.tableData[props.rowID].svg
                }} ref={ref} />
                Hover over a linkage to see a summary
            </div>

            <h3 className="text-left text-xl w-full">Visualise</h3>

            <label htmlFor="contour-range-text" className="block mb-2 text-sm font-medium text-gray-900 translate-y-10">Map Contour</label>
            <input id="contour-range" type="range" min="0" max="1" step="0.05" defaultValue="0.2" className="w-36 h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer dark:bg-gray-700 translate-y-10" onChange={handleContourChange} />

            <MoorhenContextProvider defaultBackgroundColor={[51, 65, 85, 1]}>

                <MoorhenContainer forwardControls={props.forwardControls} setMoorhenDimensions={() => {
                    return [width, height];
                }} viewOnly={true} />

            </MoorhenContextProvider>

            <h3 className="text-left text-xl w-full">Torsion Plots</h3>

            <TorsionMultiPlot torsions={props.tableData[props.rowID].torsions} tab={torsionTab} setTab={setTorsionTab} />

        </div>);
}
