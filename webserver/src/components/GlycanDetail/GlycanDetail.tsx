import { lazy, useCallback, useEffect, useState, useRef } from "react";
import { MoorhenContainer } from 'moorhen'

import { GlycanDetailProps, TableDataEntry } from "../../interfaces/types"
const GlycanDetailInfoBox = lazy(() => import('./GlycanDetailInfoBox'));
const TorsionMultiPlot = lazy(() => import('../TorsionPlot/TorsionMultiPlot.tsx'));
import { useSelector } from "react-redux";
// import ResponsiveReactGridLayout from "react-grid-layout"
import { Responsive, WidthProvider } from "react-grid-layout";

const ResponsiveGridLayout = WidthProvider(Responsive);

function GlycanDetailSNFGBox(props: {
    key: string,
    tableDataEntries: Array<TableDataEntry>,
    rowID: number
}) {
    const molecules = useSelector((state: any) => state.molecules)

    async function handle_click(e) {
        console.log(molecules)
        let new_center_string = e.target.dataset.chainid + "/" + e.target.dataset.seqnum + "(" + e.target.dataset.resname + ")"
        const selectedMolecule = molecules.find((molecule) => molecule.name === "mol-1")
        console.log(selectedMolecule)
        await selectedMolecule.centreOn(new_center_string)

    }
    const ref = useCallback((node: HTMLElement | null) => {

        if (node !== null) {
            // console.log(node)
            node.querySelector("svg").style.display = 'block'
            node.querySelector("svg").style.margin = 'auto'

            let useList = node.querySelectorAll('use')

            for (let i = 0; i < useList.length; i++) {
                useList[i].addEventListener('click', handle_click)
                useList[i].addEventListener('mousedown', (e) => { e.stopPropagation() })
                useList[i].addEventListener('touchstart', (e) => { e.stopPropagation() })
            }
        }

        document.querySelectorAll("svg")[0].setAttribute("width", "50vw")
        document.querySelectorAll("svg")[0].setAttribute("height", "100%")

    }, [props.rowID])

    // const [svg, setSVG] = useState(props.tableDataEntries[props.rowID].svg)
    // setSVG(svg.replace(' meet"', '"style="display:block"'))

    return <div key={props.key} className="px-8">
        <h3 className="text-left text-xl w-full">SNFG</h3>
        <div className="text-sm text-center text-primary">
            <div className="mt-4 py-4 mx-auto" id="svgContainer" dangerouslySetInnerHTML={{
                __html: props.tableDataEntries[props.rowID].svg
            }} ref={ref} />
            Hover over a linkage to see a summary
        </div>
    </div>;
}

function GlycanDetailMoorhenView(props: {
    key: string,
    onChange: (e) => Promise<void>,
    onSymmetryChange: (e) => any,
    moorhenProps: any,
    moorhenDimensions: () => [number, number]
}) {

    const [size, setSize] = useState(Math.min(800, 0.9 * window.screen.width))
    const [dimensions, setDimensions] = useState<Record<string, number>>(size >= 800 ? { width: size, height: (size * (3 / 4)) }: { width: size, height: size })
    
    const dimensionRef = useRef()

    dimensionRef.current = dimensions
    useEffect(() => {
        function handleResize() {
            const size = Math.min(800, 0.9 * window.screen.width)
            if (size >= 800) {
                setDimensions({ width: size, height: (size * (3 / 4)) })
            } else {
                setDimensions({ width: size, height: size })
            }
            console.log("resizing", size, window.screen.width)

        }
        window.addEventListener('resize', handleResize)

    }, [])

    const moorhenDimensionCallback = () => {
        console.log(dimensionRef)
        return [dimensionRef.current.width, dimensionRef.current.height]
    }

    return <div key={props.key} className="px-8 flex flex-col items-center">
        <h3 className="text-left text-xl w-full">Visualise with <i>Moorhen</i> <a href="https://moorhen.org" title="Go to Moorhen.org"><img className="inline h-8" src="./moorhen_logo.png"></img></a></h3>

        <label htmlFor="contour-range-text" className="block mb-2 text-sm font-medium text-gray-909">Map
            Contour</label>
        <input id="contour-range" type="range" min="0" max="1" step="0.05" defaultValue="0.2"
            className="w-36 mb-2 h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer"
            onChange={props.onChange} onMouseDown={(e) => { e.stopPropagation() }}
            onTouchStart={(e) => { e.stopPropagation() }}
        />

        <div className="mt-4">
            <label className="relative inline-flex items-center cursor-pointer">
                <input type="checkbox" value="" className="sr-only peer" onMouseDown={(e) => { e.stopPropagation() }}
                    onTouchStart={(e) => { e.stopPropagation() }} />
                <div onMouseDown={(e) => { e.stopPropagation() }}
                    onTouchStart={(e) => { e.stopPropagation() }}
                    onClick={props.onSymmetryChange}
                    className="w-11 h-6 bg-tertiary peer-focus:outline-none peer-focus:ring-4 peer-focus:ring-blue-300 dark:peer-focus:ring-blue-800 rounded-full peer dark:bg-gray-700 peer-checked:after:translate-x-full rtl:peer-checked:after:-translate-x-full peer-checked:after:border-white after:content-[''] after:absolute after:top-[2px] after:start-[2px] after:bg-white after:border-gray-300 after:border after:rounded-full after:h-5 after:w-5 after:transition-all dark:border-gray-600 peer-checked:bg-blue-600"></div>
                <span className="ms-3 text-sm font-medium text-gray-900 dark:text-gray-300">Toggle Symmetry</span>
            </label>
        </div>

        <div className="mx-auto">
            <MoorhenContainer {...props.moorhenProps} setMoorhenDimensions={moorhenDimensionCallback} viewOnly={true} />
        </div>
    </div>;
}

function GlycanDetailTorsionPlot(props: {
    key: string,
    tableDataEntries: Array<TableDataEntry>,
    rowID: number,
    tab: number,
    tab1: (value: (((prevState: number) => number) | number)) => void
}) {
    const [dimension, setDimension] = useState<number>(500)

    useEffect(() => {
        function handleResize() {
            const size = Math.min(500, window.innerWidth)
            setDimension(size)
        }
        window.addEventListener('resize', handleResize)

    })

    return <div key={props.key} className="px-8">
        <h3 className="text-left text-xl w-full mb-6 mt-6">Torsion Plots</h3>

        <TorsionMultiPlot torsions={props.tableDataEntries[props.rowID].torsions} tab={props.tab}
            setTab={props.tab1} size={dimension} />
    </div>;
}

export default function GlycanDetail(props: GlycanDetailProps) {
    const molecules = useSelector((state: any) => state.molecules)

    async function handleContourChange(e) {
        props.map.contourLevel = Number(e.target.value)
        await props.map.doCootContour(
            ...props.map.glRef.current.origin.map(coord => -coord),
            props.map.mapRadius,
            props.map.contourLevel
        )
    }

    async function handleToggleSymmetry(e) {
        console.log(e)
        const currentMolecule = molecules.find(molecule => molecule.name == "mol-1")
        if (currentMolecule) {
            console.log(currentMolecule)
            await currentMolecule.toggleSymmetry()
        }
    }

    const [width, setWidth] = useState(800);
    const [height, setHeight] = useState(500);
    const [torsionTab, setTorsionTab] = useState(0)
    const layout = [
        { i: "info", x: 0, y: 0, w: 1, h: 1, static: true, isResizable: false },
        { i: "snfg", x: 1, y: 2, w: 1, h: 1, isResizable: false },
        { i: "moorhen", x: 0, y: 2, w: 1, h: 2, isResizable: false },
        { i: "torsions", x: 1, y: 2, w: 1, h: 2, isResizable: false }
    ];

    const md_layout = [
        { i: "info", x: 0, y: 0, w: 1, h: 1, static: true, isResizable: false },
        { i: "snfg", x: 1, y: 0, w: 1, h: 1, isResizable: false },
        { i: "moorhen", x: 0, y: 2, w: 1, h: 2, isResizable: false },
        { i: "torsions", x: 1, y: 1, w: 1, h: 2, isResizable: false }
    ]

    const sm_layout = [
        { i: "info", x: 0, y: 0, w: 1, h: 0.5, static: true, isResizable: false },
        { i: "snfg", x: 0, y: 2, w: 1, h: 1, isResizable: false },
        { i: "moorhen", x: 0, y: 4, w: 1, h: 1.75, isResizable: false },
        { i: "torsions", x: 0, y: 6, w: 1, h: 1.75, isResizable: false }
    ]

    const getLayouts = () => {
        const savedLayouts = localStorage.getItem("grid-layout");
        return { lg: layout, md: md_layout, sm: sm_layout }
        return savedLayouts ? JSON.parse(savedLayouts) : { lg: layout, md: md_layout };
    };
    const handleLayoutChange = (layout, layouts) => {
        localStorage.setItem("grid-layout", JSON.stringify(layouts));
    };

    return (
        // <div className="flex flex-wrap gap-2 justify-center items-center"
        //     style={{ display: !props.hideMoorhen ? 'flex' : 'none' }}>
        <>
            <div className="flex items-center justify-center w-full mb-6">
                <div className="flex-1">
                    <button onClick={() => {
                        props.setHideMoorhen(true);
                        window.scrollTo(0, props.scrollPosition);
                        setTorsionTab(0)
                    }}>
                        <span className="text-xl">&#8592; Back To Table</span>
                    </button>
                </div>
                <h2 className="">Glycan Details</h2>
                <div className="flex-1"></div>
            </div>

            <ResponsiveGridLayout
                className="layout"
                layouts={getLayouts()}
                breakpoints={{ lg: 1700, sm: 0 }}
                cols={{ lg: 2, md: 2, sm: 1 }}
                rowHeight={400}
                width={1000}
                containerPadding={8}
                onLayoutChange={handleLayoutChange}

            >

                <div key="info">
                    <GlycanDetailInfoBox key={"aa"} row={props.tableData[props.rowID]} />
                </div>
                <div key="snfg">
                    <GlycanDetailSNFGBox key={"bb"} tableDataEntries={props.tableData} rowID={props.rowID} />
                </div>
                <div key="moorhen">
                    <GlycanDetailMoorhenView key={"cc"} onChange={handleContourChange} onSymmetryChange={handleToggleSymmetry} moorhenProps={props.moorhenProps}
                        moorhenDimensions={() => {
                            return [width, height];
                        }} />
                </div>
                <div key="torsions">
                    <GlycanDetailTorsionPlot key={"dd"} tableDataEntries={props.tableData} rowID={props.rowID} tab={torsionTab}
                        tab1={setTorsionTab} />
                </div>
            </ResponsiveGridLayout>



        </>
    );
}
