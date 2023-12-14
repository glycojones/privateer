import React, { useMemo, useEffect, useState } from 'react';
import { useTable } from 'react-table';
import { SugarListColumns } from '../../data/Constants';
import styled from 'styled-components';

function customSort(
    a: Record<string, string>,
    b: Record<string, string>
): number {
    const splitA = a['Sugar ID'].split('-');
    const splitB = b['Sugar ID'].split('-');

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
            const sugars: Array<Record<string, any>> = glycanType[i].Sugars;
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
    }, []);

    return (
        <div className="flex flex-col mx-auto p-16">
            <span className="text-xl mb-2">
                Detailed monosaccharide validation data
            </span>

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
