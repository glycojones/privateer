import React, { lazy, useEffect, useRef, useState } from 'react';
import { type GlycanDetailProps } from '../../interfaces/types';
import { useSelector } from 'react-redux';
import { Responsive, WidthProvider } from 'react-grid-layout';
import { GlycanDetailSNFGBox } from './GlycanDetailSNFGBox.tsx';
import { GlycanDetailMoorhenView } from './GlycanDetailMoorhenView.tsx';
import { GlycanDetailTorsionPlot } from './GlycanDetailTorsionPlot.tsx';
import { GlycanDetailLayout } from '../../data/Constants.tsx';
import { Tooltip, type TooltipRefProps } from 'react-tooltip';

const GlycanDetailInfoBox = lazy(
  async () => await import('./GlycanDetailInfoBox')
);

const ResponsiveGridLayout = WidthProvider(Responsive);

export default function GlycanDetail (props: GlycanDetailProps) {
  const molecules = useSelector((state: any) => state.molecules);
  const [mapContour, setMapContour] = useState(0.3);

  async function handleContourChange (e) {
    props.map.contourLevel = Number(e.target.value);
    setMapContour(Number(e.target.value));
    await props.map.doCootContour(
      ...props.map.glRef.current.origin.map((coord) => -coord),
      props.map.mapRadius,
      props.map.contourLevel
    );
  }

  async function handleToggleSymmetry () {
    const currentMolecule = molecules.find(
      (molecule) => molecule.name === 'mol-1'
    );
    if (currentMolecule !== null) {
      await currentMolecule.toggleSymmetry();
    }
  }

  const width = 800;
  const height = 500;
  const [torsionTab, setTorsionTab] = useState(0);

  const getLayouts = () => {
    const savedLayouts = localStorage.getItem('grid-layout');
    // return GlycanDetailLayout;
    return savedLayouts !== null
      ? JSON.parse(savedLayouts)
      : GlycanDetailLayout;
  };
  const handleLayoutChange = (layout, layouts) => {
    localStorage.setItem('grid-layout', JSON.stringify(layouts));
  };

  const tooltipRef = useRef<TooltipRefProps>(null);
  const tooltipContent: string =
      'You can drag around each report component to your desired position';

  // useEffect(() => {
  //   tooltipRef.current?.open(
  //       {
  //         content: tooltipContent,
  //         delay: 10000
  //       }
  //   )
  //
  //   setTimeout(() => {
  //     tooltipRef.current?.close()
  //   }, 1000)
  // }, []);

  return (
    <>
      <Tooltip id="maintooltip" ref={tooltipRef} place={'right'} opacity={0.8}/>

      <div className="flex justify-center w-full mb-6 ">
        <div className="flex-1 text-left">
          <button
            onClick={() => {
              props.setHideMoorhen(true);
              // window.scrollTo(0, props.scrollPosition);
              setTorsionTab(0);
            }}
          >
            <span className="text-xl ml-8">&#8592; Back To Table</span>
          </button>
        </div>
        <h2 className="font-bold">Glycan Details
          <svg
              id="snfg_icon"
              data-tooltip-id="maintooltip"
              data-tooltip-html={tooltipContent}
              className="inline ml-2 mb-1 w-5"
              xmlns="http://www.w3.org/2000/svg"
              viewBox="0 0 512 512"
          >
            <path d="M256 512A256 256 0 1 0 256 0a256 256 0 1 0 0 512zM216 336h24V272H216c-13.3 0-24-10.7-24-24s10.7-24 24-24h48c13.3 0 24 10.7 24 24v88h8c13.3 0 24 10.7 24 24s-10.7 24-24 24H216c-13.3 0-24-10.7-24-24s10.7-24 24-24zm40-208a32 32 0 1 1 0 64 32 32 0 1 1 0-64z" />
          </svg>

        </h2>
        <div className="flex-1"></div>
      </div>

      <ResponsiveGridLayout
        className="layout"
        layouts={getLayouts()}
        breakpoints={{ xl: 3000, lg: 1700, sm: 0 }}
        cols={{ xl: 4, lg: 2, md: 2, sm: 1 }}
        rowHeight={400}
        width={1000}
        containerPadding={8}
        onLayoutChange={handleLayoutChange}
      >
        <div key="info" className="border-2 rounded-xl border-darkgray">
          <GlycanDetailInfoBox key={'aa'} row={props.tableData[props.rowID]} />
        </div>
        <div key="snfg" className="border-2 rounded-xl border-darkgray">
          <GlycanDetailSNFGBox
            key={'bb'}
            tableDataEntries={props.tableData}
            rowID={props.rowID}
          />
        </div>
        <div key="moorhen" className="border-2 rounded-xl border-darkgray">
          <GlycanDetailMoorhenView
            key={'cc'}
            onChange={handleContourChange}
            onSymmetryChange={handleToggleSymmetry}
            moorhenProps={props.moorhenProps}
            mapContour={mapContour}
            moorhenDimensions={() => {
              return [width, height];
            }}
          />
        </div>
        <div key="torsions" className="border-2 rounded-xl border-darkgray">
          <GlycanDetailTorsionPlot
            key={'dd'}
            tableDataEntries={props.tableData}
            rowID={props.rowID}
            tab={torsionTab}
            tab1={setTorsionTab}
          />
        </div>
      </ResponsiveGridLayout>
    </>
  );
}
