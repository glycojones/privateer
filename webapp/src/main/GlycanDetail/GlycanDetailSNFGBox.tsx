import { type TableDataEntry } from '../../interfaces/types.ts';
import { useSelector } from 'react-redux';
import React, { useCallback, useRef } from 'react';
import { Tooltip, type TooltipRefProps } from 'react-tooltip';

export function GlycanDetailSNFGBox(props: {
    key: string;
    tableDataEntries: TableDataEntry[];
    rowID: number;
}) {
    const molecules = useSelector((state: any) => state.molecules);

    async function handleClick(e) {
        const newCenterString =
            e.target.dataset.chainid +
            '/' +
            e.target.dataset.seqnum +
            '(' +
            e.target.dataset.resname +
            ')';
        const selectedMolecule = molecules.find(
            (molecule) => molecule.name === 'mol-1'
        );
        const _center = await selectedMolecule.centreOn(newCenterString);
    }

    const ref = useCallback(
        (node: HTMLElement | null) => {
            if (node !== null) {
                // console.log(node)
                node.querySelector('svg').style.display = 'block';
                node.querySelector('svg').style.margin = 'auto';

                const useList = node.querySelectorAll('use');

                for (let i = 0; i < useList.length; i++) {
                    useList[i].addEventListener('click', handleClick);
                    useList[i].addEventListener('mousedown', (e) => {
                        e.stopPropagation();
                    });
                    useList[i].addEventListener('touchstart', (e) => {
                        e.stopPropagation();
                    });
                }
            }

            document.querySelectorAll('svg')[0].setAttribute('width', '50vw');
            document.querySelectorAll('svg')[0].setAttribute('height', '100%');
        },
        [props.rowID]
    );

    // const [svg, setSVG] = useState(props.tableDataEntries[props.rowID].svg)
    // setSVG(svg.replace(' meet"', '"style="display:block"'))
    const tooltipRef = useRef<TooltipRefProps>(null);

    const tooltipContent: string =
        '<b>Standard Nomenclature For Glycans</b><br>Click on a sugar to see it in Moorhen';

    return (
        <div key={props.key} className="px-8 text-left font-bold mt-2">
            <Tooltip id="snfgtooltip" ref={tooltipRef} place={'right'} />

            <h3 className="text-left text-xl w-full">
                SNFG
                <svg
                    id="snfg_icon"
                    data-tooltip-id="snfgtooltip"
                    data-tooltip-html={tooltipContent}
                    className="inline ml-2 w-5"
                    xmlns="http://www.w3.org/2000/svg"
                    viewBox="0 0 512 512"
                >
                    <path d="M256 512A256 256 0 1 0 256 0a256 256 0 1 0 0 512zM216 336h24V272H216c-13.3 0-24-10.7-24-24s10.7-24 24-24h48c13.3 0 24 10.7 24 24v88h8c13.3 0 24 10.7 24 24s-10.7 24-24 24H216c-13.3 0-24-10.7-24-24s10.7-24 24-24zm40-208a32 32 0 1 1 0 64 32 32 0 1 1 0-64z" />
                </svg>
                {/* <svg id="snfg_icon" data-tooltip-id="tooltip" data-tooltip-html={tooltipContent} */}
                {/*     className="inline ml-2 w-6"  viewBox="0 0 48 48" xmlns="http://www.w3.org/2000/svg"><path d="M0 0h48v48h-48z" fill="none"/><path d="M22 34h4v-12h-4v12zm2-30c-11.05 0-20 8.95-20 20s8.95 20 20 20 20-8.95 20-20-8.95-20-20-20zm0 36c-8.82 0-16-7.18-16-16s7.18-16 16-16 16 7.18 16 16-7.18 16-16 16zm-2-22h4v-4h-4v4z"/></svg> */}
            </h3>
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
