import React, { useMemo, useEffect, useState, useRef } from 'react';
import { useTable } from 'react-table';
import { SugarListColumns } from '../../data/Constants.tsx';
import styled from 'styled-components';
import { Tooltip, type TooltipRefProps } from 'react-tooltip';

function customSort(
    a: Record<string, string>,
    b: Record<string, string>
): number {
    const splitA = a.sugarId.split('-');
    const splitB = b.sugarId.split('-');

    if (splitA[1] < splitB[1]) {
        return -1;
    }
    if (splitB[1] < splitA[1]) {
        return 1;
    }

    if (Number(splitA[2]) < Number(splitB[2])) {
        return -1;
    }

    if (Number(splitB[2]) < Number(splitA[2])) {
        return 1;
    }

    return 0;
}

function parseResults(data) {
    const glycans = data.data.glycans;

    const tableData = [];

    for (const key in glycans) {
        const glycanType = glycans[key];
        for (let i = 0; i < glycanType.length; i++) {
            const sugars: Array<Record<string, any>> = glycanType[i].sugars;
            for (let j = 0; j < sugars.length; j++) {
                const keys = Object.keys(sugars[j]);
                const entry = {};

                for (let k = 0; k < keys.length; k++) {
                    const typeKey = typeof sugars[j][keys[k]];

                    if (typeKey === 'number') {
                        entry[keys[k]] = sugars[j][keys[k]].toFixed(2);
                    } else {
                        entry[keys[k]] = sugars[j][keys[k]];
                    }
                }
                entry.type = key;
                tableData.push(entry);
            }
        }
    }

    tableData.sort(customSort);
    return tableData;
}

const Styles = styled.div`
    table {
        border-collapse: collapse;
        border-spacing: 0;
        width: 100%;
        // border: 0px solid #ddd;
    }

    table th {
        text-align: left;
        padding: 16px;
        // border: 1px solid #ddd;
    }

    table td {
        text-align: center;
        padding: 16px;
        border: 1px solid #ddd;
        border-style: solid none;
    }

    table td {
        border-style: none none;
    }

    table tr:nth-child(even) {
        background-color: #f6f6f6;
    }

    table tr:nth-child(even) {
        background-color: #f4f9ff;
        // color: #000000
    }

    table th {
        padding-top: 12px;
        padding-bottom: 12px;
        text-align: center;
        background-color: #f4f9ff;
        color: black;
    }

    table th:first-of-type {
        border-top-left-radius: 30px;
    }

    table th:last-of-type {
        border-top-right-radius: 30px;
    }

    table tr:last-of-type td:first-of-type {
        border-bottom-left-radius: 30px;
    }

    table tr:last-of-type td:last-of-type {
        border-bottom-right-radius: 30px;
    }
`;

export default function SugarList(props) {
    const [data, setData] = useState([]);
    const columns = useMemo(() => SugarListColumns, []);
    const {
        _getTableProps,
        getTableBodyProps,
        headerGroups,
        rows,
        prepareRow,
    } = useTable({ columns, data });

    useEffect(() => {
        const results = parseResults(props);
        setData(results);
    }, [props]);

    const tooltipRef = useRef<TooltipRefProps>(null);

    const tooltipContent: string = `A “spherical” coordinate system is used to describe six-membered rings, where the total puckering amplitude of the ring, Q, is defined, as well as <br> 
    the distortion-type which is specified by two angles, phi and theta: <br><br>
    Radius, Q: describes the overall distortion or shape of the puckered ring, <br>
    measuring the deviation from a perfectly flat six-membered ring (Q = 0). This is calculated using out-of-plane deviations of puckered rings using the <br>
     z-coordinates of the ring atoms relative to a mean plane cutting through the ring. <br><br>
    Azimuthal angle, theta: describes the orientation of the puckering plane around the rings circumference. This parameter indicates the direction <br>
    in which the puckering occurs along the ring, providing information about the spatial distribution of the distortion.<br><br>
    Meridian angle, phi: describes the orientation of the out-of-plane puckering with respect to the rings mean plane. This provides information about <br>
     the directionality of the puckering, indicating whether the distortion of the ring is primarily upward or downward with respect to the mean plane.`;

    return (
        <div className="flex flex-col mx-auto p-16">
            <Tooltip id="sugarListToolTip" ref={tooltipRef} place={'right'} />

            <div className="flex items-center align-content-center">
                <span className="text-xl mb-2">
                    Detailed monosaccharide validation data
                </span>
                <svg
                    id="snfg_icon"
                    data-tooltip-id="sugarListToolTip"
                    data-tooltip-html={tooltipContent}
                    className="ml-2 mb-2 w-4"
                    xmlns="http://www.w3.org/2000/svg"
                    viewBox="0 0 512 512"
                >
                    <path d="M256 512A256 256 0 1 0 256 0a256 256 0 1 0 0 512zM216 336h24V272H216c-13.3 0-24-10.7-24-24s10.7-24 24-24h48c13.3 0 24 10.7 24 24v88h8c13.3 0 24 10.7 24 24s-10.7 24-24 24H216c-13.3 0-24-10.7-24-24s10.7-24 24-24zm40-208a32 32 0 1 1 0 64 32 32 0 1 1 0-64z" />
                </svg>
            </div>
            <Styles>
                <table {...getTableBodyProps()}>
                    <thead>
                        {headerGroups.map((headerGroup) => (
                            <tr
                                {...headerGroup.getHeaderGroupProps()}
                                key={headerGroup.name}
                            >
                                {headerGroup.headers.map((column) => (
                                    <th
                                        {...column.getHeaderProps()}
                                        key={column.name}
                                    >
                                        {column.render('Header')}
                                    </th>
                                ))}
                            </tr>
                        ))}
                    </thead>
                    <tbody {...getTableBodyProps()}>
                        {rows.map((row) => {
                            prepareRow(row);
                            return (
                                <tr
                                    {...row.getRowProps()}
                                    title="Click to visualise."
                                    id="row"
                                    key={row.name}
                                >
                                    {row.cells.map((cell) => {
                                        return (
                                            <td
                                                {...cell.getCellProps()}
                                                key={cell.name}
                                            >
                                                {cell.render('Cell')}
                                            </td>
                                        );
                                    })}
                                </tr>
                            );
                        })}
                    </tbody>
                </table>
            </Styles>
        </div>
    );
}
