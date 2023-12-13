import { type TableDataEntry } from "../../interfaces/types.ts";
import { useSelector } from "react-redux";
import React, { useCallback } from "react";

export function GlycanDetailSNFGBox(props: {
  key: string;
  tableDataEntries: TableDataEntry[];
  rowID: number;
}) {
  const molecules = useSelector((state: any) => state.molecules);

  async function handle_click(e) {
    console.log(molecules);
    const new_center_string =
      e.target.dataset.chainid +
      "/" +
      e.target.dataset.seqnum +
      "(" +
      e.target.dataset.resname +
      ")";
    const selectedMolecule = molecules.find(
      (molecule) => molecule.name === "mol-1",
    );
    console.log(selectedMolecule);
    await selectedMolecule.centreOn(new_center_string);
  }

  const ref = useCallback(
    (node: HTMLElement | null) => {
      if (node !== null) {
        // console.log(node)
        node.querySelector("svg").style.display = "block";
        node.querySelector("svg").style.margin = "auto";

        const useList = node.querySelectorAll("use");

        for (let i = 0; i < useList.length; i++) {
          useList[i].addEventListener("click", handle_click);
          useList[i].addEventListener("mousedown", (e) => {
            e.stopPropagation();
          });
          useList[i].addEventListener("touchstart", (e) => {
            e.stopPropagation();
          });
        }
      }

      document.querySelectorAll("svg")[0].setAttribute("width", "50vw");
      document.querySelectorAll("svg")[0].setAttribute("height", "100%");
    },
    [props.rowID],
  );

  // const [svg, setSVG] = useState(props.tableDataEntries[props.rowID].svg)
  // setSVG(svg.replace(' meet"', '"style="display:block"'))

  return (
    <div key={props.key} className="px-8">
      <h3 className="text-left text-xl w-full">SNFG</h3>
      <div className="text-sm text-center text-primary">
        <div
          className="mt-4 py-4 mx-auto"
          id="svgContainer"
          dangerouslySetInnerHTML={{
            __html: props.tableDataEntries[props.rowID].svg,
          }}
          ref={ref}
        />
        Hover over a linkage to see a summary
      </div>
    </div>
  );
}
