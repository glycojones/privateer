import React, {
  lazy,
  type ReactElement,
  useEffect,
  useRef,
  useState
} from 'react';
import {
  addMolecule,
  addMap,
  setEnableAtomHovering,
  MoorhenMolecule,
  MoorhenMap
} from 'moorhen';
import { type SNFGProps } from '../../interfaces/types';
import { useSelector, useDispatch } from 'react-redux';
import Loading from '../Loading/Loading.tsx';
import GlycanDetail from '../GlycanDetail/GlycanDetail.tsx';

const SVGTable = lazy(async () => await import('../SVGTable/SVGTable.tsx'));

export default function SNFG (props: SNFGProps): ReactElement {
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
    mapsRef
  };

  // const defaultBondSmoothness = useSelector((state: any) => state.sceneSettings.defaultBondSmoothness)
  // const backgroundColor = useSelector((state: any) => state.canvasStates.backgroundColor)

  // const [yScrollPosition, setYScrollPosition] = useState(0);

  useEffect(() => {
    async function loadMapAndModel () {
      if (cootInitialized === true && !dataLoaded) {
        setAllowRowClick(true);

        const newMolecule = new MoorhenMolecule(
          commandCentre,
          glRef,
          './baby-gru/monomers'
        );
        await newMolecule.loadToCootFromString(props.fileContent, 'mol-1');
        await newMolecule.fetchIfDirtyAndDraw('CBs');
        const id = props.tableData[rowID].id;
        const sugarName = id.split('-')[0];
        const sugarId = id.split('-')[1].split('/')[0].split(':')[0];
        const sugarChain = id.split('/')[1].split('_')[0];

        const centerString = sugarChain + '/' + sugarId + '(' + sugarName + ')';
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
            calcStructFact: true
          };
          if (props.PDBCode === '') {
            await newMap.loadToCootFromMtzData(
              props.mtzData,
              'map-1',
              mapMetadata
            );
          } else {
            await newMap.loadToCootFromMapData(props.mtzData, 'map-1', false);
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
    async function moveView () {
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

      const centerString = sugarChain + '/' + sugarId + '(' + sugarName + ')';
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
    moorhenProps: collectedProps
  };

  const svgTableProps = {
    tableData: props.tableData,
    allowRowClick,
    rowClick: rowClicked,
    setRowClicked,
    setRowID
  };

  return (
    <div className="flex flex-col">
      {dataLoaded
        ? (
        <div
          style={{ display: hideMoorhen ? 'block' : 'none' }}
          id="tableContainer"
        >
          <div className="flex flex-col ">
            <h2 className="my-4 text-lg sm:text-2xl text-center sm:text-left mx-auto">
              Detected {props.tableData.length} Glycans in {props.filename}
            </h2>
            <SVGTable {...svgTableProps} />
          </div>
        </div>
          )
        : (
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
