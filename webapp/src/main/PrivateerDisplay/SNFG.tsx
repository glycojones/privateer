import React, {
    lazy,
    type ReactElement,
    useEffect,
    useRef,
    useState,
} from 'react';
import {
    addMolecule,
    addMap,
    setEnableAtomHovering,
    MoorhenMolecule,
    MoorhenMap,
} from 'moorhen';
import { type SNFGProps } from '../../interfaces/types.ts';
import { useSelector, useDispatch } from 'react-redux';
import Loading from '../../shared/Loading/Loading.tsx';
import GlycanDetail from '../GlycanDetail/GlycanDetail.tsx';

const SVGTable = lazy(async () => await import('../SVGTable/SVGTable.tsx'));

export default function SNFG(props: SNFGProps): ReactElement {
    const [rowClicked, setRowClicked] = useState(false);
    const [rowID, setRowID] = useState(0);
    const [hideMoorhen, setHideMoorhen] = useState(true);
    const [allowRowClick, setAllowRowClick] = useState(true);

    const [dataLoaded, setDataLoaded] = useState(false);

    const cootInitialized = useSelector(
        (state: any) => state.generalStates.cootInitialized
    );
    const controls = useRef();

    const [map, setMap] = useState<MoorhenMap>();
    const molecules = useSelector((state: any) => state.molecules);
    const glRef = useRef(null);
    const timeCapsuleRef = useRef(null);
    const commandCentre = useRef(null);
    const moleculesRef = useRef(null);
    const mapsRef = useRef(null);
    const dispatch = useDispatch();

    const collectedProps = {
        glRef,
        timeCapsuleRef,
        commandCentre,
        moleculesRef,
        mapsRef,
    };

    // const defaultBondSmoothness = useSelector((state: any) => state.sceneSettings.defaultBondSmoothness)
    // const backgroundColor = useSelector((state: any) => state.canvasStates.backgroundColor)

    // const [yScrollPosition, setYScrollPosition] = useState(0);

    function saveSNFGs() {
        if (props.tableData === null) return;

        const data = props.tableData;

        for (let i = 0; i < data.length; i++) {
            const svgBlob = new Blob([data[i].svg], {
                type: 'image/svg+xml;charset=utf-8',
            });
            const svgUrl = URL.createObjectURL(svgBlob);
            const downloadLink = document.createElement('a');
            downloadLink.href = svgUrl;
            downloadLink.download = `${data[i].id}.svg`;
            document.body.appendChild(downloadLink);
            downloadLink.click();
            document.body.removeChild(downloadLink);
        }
    }

    useEffect(() => {
        async function loadMapAndModel() {
            if (cootInitialized === true && !dataLoaded) {
                setAllowRowClick(true);

                const newMolecule = new MoorhenMolecule(
                    commandCentre,
                    glRef,
                    './baby-gru/monomers'
                );
                await newMolecule.loadToCootFromString(
                    props.fileContent,
                    'mol-1'
                );
                await newMolecule.fetchIfDirtyAndDraw('CBs');
                const id = props.tableData[rowID].id;
                const sugarName = id.split('-')[0];
                const sugarId = id.split('-')[1].split('/')[0].split(':')[0];
                const sugarChain = id.split('/')[1].split('_')[0];

                const centerString =
                    sugarChain + '/' + sugarId + '(' + sugarName + ')';
                await newMolecule.centreOn(centerString, true);

                dispatch(addMolecule(newMolecule));
                dispatch(setEnableAtomHovering(false));

                if (props.mtzData !== null) {
                    const newMap = new MoorhenMap(commandCentre, glRef);
                    const mapMetadata = {
                        F: 'FWT',
                        PHI: 'PHWT',
                        Fobs: 'FP',
                        SigFobs: 'SIGFP',
                        FreeR: 'FREE',
                        isDifference: false,
                        useWeight: false,
                        calcStructFact: true,
                    };
                    if (props.PDBCode === '') {
                        await newMap.loadToCootFromMtzData(
                            props.mtzData,
                            'map-1',
                            mapMetadata
                        );
                    } else {
                        await newMap.loadToCootFromMapData(
                            props.mtzData,
                            'map-1',
                            false
                        );
                    }
                    newMap.suggestedContourLevel = 0.3;
                    dispatch(addMap(newMap));

                    setMap(newMap);
                }

                setDataLoaded(true);

                // let newMolecule = new MoorhenMolecule(controls.current.commandCentre, controls.current.glRef, controls.current.monomerLibrary)
                // newMolecule.loadToCootFromString(props.fileContent, 'mol-1').then(() => {
                //     controls.current.changeMolecules({action: 'Add', item: newMolecule});
                //     newMolecule.fetchIfDirtyAndDraw('CBs').then(() => {
                //             let id = props.tableData[rowID].id

                //             let sugar_name = id.split("-")[0]
                //             let sugar_id = id.split("-")[1].split("/")[0].split(":")[0]
                //             let sugar_chain = id.split("/")[1].split("_")[0]

                //             let center_string = sugar_chain + "/" + sugar_id + "(" + sugar_name + ")"

                //             newMolecule.centreOn(center_string)
                //             setMolecule(newMolecule)
                //         }
                //     )
                // })

                // const map = new MoorhenMap(controls.current.commandCentre, controls.current.glRef);
                // const mapMetadata = {
                //     F: "FWT",
                //     PHI: "PHWT",
                //     Fobs: "FP",
                //     SigFobs: "SIGFP",
                //     FreeR: "FREE",
                //     isDifference: false,
                //     useWeight: false,
                //     calcStructFact: true,
                // }
                // if (props.PDBCode == "") {
                //     await map.loadToCootFromMtzData(props.mtzData, "map-1", mapMetadata)
                // }
                // else {
                //     await map.loadToCootFromMapData(props.mtzData, "map-1", false)

                // }
                // map.suggestedContourLevel = 0.3
                // console.log("Suggested level", map.suggestedContourLevel, map)
                // controls.current.changeMaps({ action: "Add", item: map })
                // controls.current.setActiveMap(map)
                // setMap(map)
                // setDataLoaded(true)
            }
        }
        loadMapAndModel().then(
            () => {},
            () => {}
        );
    }, [cootInitialized]);

    useEffect(() => {
        async function moveView() {
            if (cootInitialized === false) {
                console.log('coot not loaded');
                return;
            }
            // console.log(molecules)

            // setYScrollPosition(window.scrollY)
            const id = props.tableData[rowID].id;
            const sugarName = id.split('-')[0];
            const sugarId = id.split('-')[1].split('/')[0].split(':')[0];
            const sugarChain = id.split('/')[1].split('_')[0];

            const centerString =
                sugarChain + '/' + sugarId + '(' + sugarName + ')';
            const selectedMolecule = molecules.find(
                (molecule) => molecule.name === 'mol-1'
            );
            await selectedMolecule.centreOn(centerString);
            // window.scrollTo(0, 0)
            setHideMoorhen(false);
        }
        moveView().then(
            () => {},
            () => {}
        );
    }, [rowClicked]);

    const glycanDetailProps = {
        tableData: props.tableData,
        hideMoorhen,
        setHideMoorhen,
        rowID,
        controls,
        map,
        moorhenProps: collectedProps,
        pdbCode: props.PDBCode
    };

    const svgTableProps = {
        tableData: props.tableData,
        allowRowClick,
        rowClick: rowClicked,
        setRowClicked,
        setRowID,
    };

    return (
        <div className="flex flex-col">
            {dataLoaded ? (
                <div
                    style={{ display: hideMoorhen ? 'block' : 'none' }}
                    id="tableContainer"
                >
                    <div className="flex flex-col ">
                        <h2 className="my-4 text-lg sm:text-2xl text-center sm:text-left mx-auto">
                            Detected {props.tableData.length} Glycans in{' '}
                            {props.filename}
                            <button
                                className="text-center items-center ml-6"
                                onClick={() => {
                                    saveSNFGs();
                                }}
                                title="Download All SNFGs"
                            >
                                <svg
                                    className="h-4 w-4"
                                    xmlns="http://www.w3.org/2000/svg"
                                    height="1em"
                                    viewBox="0 0 512 512"
                                >
                                    <path d="M288 32c0-17.7-14.3-32-32-32s-32 14.3-32 32V274.7l-73.4-73.4c-12.5-12.5-32.8-12.5-45.3 0s-12.5 32.8 0 45.3l128 128c12.5 12.5 32.8 12.5 45.3 0l128-128c12.5-12.5 12.5-32.8 0-45.3s-32.8-12.5-45.3 0L288 274.7V32zM64 352c-35.3 0-64 28.7-64 64v32c0 35.3 28.7 64 64 64H448c35.3 0 64-28.7 64-64V416c0-35.3-28.7-64-64-64H346.5l-45.3 45.3c-25 25-65.5 25-90.5 0L165.5 352H64zm368 56a24 24 0 1 1 0 48 24 24 0 1 1 0-48z" />
                                </svg>
                            </button>
                        </h2>

                        <SVGTable {...svgTableProps} />
                    </div>
                </div>
            ) : (
                <Loading loadingText="Setting up Moorhen..."></Loading>
            )}

            <div
                style={{ display: hideMoorhen ? 'none' : 'block' }}
                id="tableContainer"
            >
                <GlycanDetail {...glycanDetailProps} />
            </div>
        </div>
    );
}
