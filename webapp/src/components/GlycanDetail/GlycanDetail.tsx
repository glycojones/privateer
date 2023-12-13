import React, { lazy, useState } from "react";
import { type GlycanDetailProps } from "../../interfaces/types";
import { useSelector } from "react-redux";
import { Responsive, WidthProvider } from "react-grid-layout";
import { GlycanDetailSNFGBox } from "./GlycanDetailSNFGBox.tsx";
import { GlycanDetailMoorhenView } from "./GlycanDetailMoorhenView.tsx";
import { GlycanDetailTorsionPlot } from "./GlycanDetailTorsionPlot.tsx";
import { GlycanDetailLayout } from "../../data/Constants.tsx";

const GlycanDetailInfoBox = lazy(
  async () => await import("./GlycanDetailInfoBox"),
);

const ResponsiveGridLayout = WidthProvider(Responsive);

export default function GlycanDetail(props: GlycanDetailProps) {
  const molecules = useSelector((state: any) => state.molecules);
  const [mapContour, setMapContour] = useState(0.3);

  async function handleContourChange(e) {
    props.map.contourLevel = Number(e.target.value);
    setMapContour(Number(e.target.value));
    await props.map.doCootContour(
      ...props.map.glRef.current.origin.map((coord) => -coord),
      props.map.mapRadius,
      props.map.contourLevel,
    );
  }

  async function handleToggleSymmetry(e) {
    const currentMolecule = molecules.find(
      (molecule) => molecule.name == "mol-1",
    );
    if (currentMolecule) {
      console.log(currentMolecule);
      await currentMolecule.toggleSymmetry();
    }
  }

  const [width, setWidth] = useState(800);
  const [height, setHeight] = useState(500);
  const [torsionTab, setTorsionTab] = useState(0);

  const getLayouts = () => {
    const savedLayouts = localStorage.getItem("grid-layout");
    return GlycanDetailLayout;
    return savedLayouts ? JSON.parse(savedLayouts) : GlycanDetailLayout;
  };
  const handleLayoutChange = (layout, layouts) => {
    localStorage.setItem("grid-layout", JSON.stringify(layouts));
  };

  return (
    <>
      <div className="flex items-center justify-center w-full mb-6">
        <div className="flex-1">
          <button
            onClick={() => {
              props.setHideMoorhen(true);
              window.scrollTo(0, props.scrollPosition);
              setTorsionTab(0);
            }}
          >
            <span className="text-xl">&#8592; Back To Table</span>
          </button>
        </div>
        <h2 className="">Glycan Details</h2>
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
        <div key="info">
          <GlycanDetailInfoBox key={"aa"} row={props.tableData[props.rowID]} />
        </div>
        <div key="snfg">
          <GlycanDetailSNFGBox
            key={"bb"}
            tableDataEntries={props.tableData}
            rowID={props.rowID}
          />
        </div>
        <div key="moorhen">
          <GlycanDetailMoorhenView
            key={"cc"}
            onChange={handleContourChange}
            onSymmetryChange={handleToggleSymmetry}
            moorhenProps={props.moorhenProps}
            mapContour={mapContour}
            moorhenDimensions={() => {
              return [width, height];
            }}
          />
        </div>
        <div key="torsions">
          <GlycanDetailTorsionPlot
            key={"dd"}
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
