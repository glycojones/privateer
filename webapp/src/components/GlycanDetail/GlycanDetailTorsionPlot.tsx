import { type TableDataEntry } from "../../interfaces/types.ts";
import React, { lazy, useEffect, useState } from "react";
const TorsionMultiPlot = lazy(
  async () => await import("../TorsionPlot/TorsionMultiPlot.tsx"),
);

export function GlycanDetailTorsionPlot(props: {
  key: string;
  tableDataEntries: TableDataEntry[];
  rowID: number;
  tab: number;
  tab1: (value: ((prevState: number) => number) | number) => void;
}) {
  const [dimension, setDimension] = useState<number>(500);

  useEffect(() => {
    function handleResize() {
      const size = Math.min(500, window.innerWidth);
      setDimension(size);
    }

    window.addEventListener("resize", handleResize);
  });

  return (
    <div key={props.key} className="px-8">
      <h3 className="text-left text-xl w-full mb-6 mt-6">Torsion Plots</h3>

      <TorsionMultiPlot
        torsions={props.tableDataEntries[props.rowID].torsions}
        tab={props.tab}
        setTab={props.tab1}
        size={dimension}
      />
    </div>
  );
}
