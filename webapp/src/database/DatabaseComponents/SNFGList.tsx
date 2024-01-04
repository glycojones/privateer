import React, { useMemo, useEffect, useState } from 'react';
import { useTable } from 'react-table';
import { DatabaseColumns } from '../../data/Constants.tsx';
import styled from 'styled-components';

function customSort(a, b): number {
    a = a.chain;
    b = b.chain;
    if (a < b) return -1;
    if (a > b) return 1;
    return 0;
}

function parseResults(data) {
    const glycans = data.data.glycans;

    const tableData: any[] = [];

    for (const key in glycans) {
        const glycanType = glycans[key];
        for (let i = 0; i < glycanType.length; i++) {
            const chain = glycanType[i].RootSugarChainID;

            const SNFG = glycanType[i].SNFG;
            const WURCS = glycanType[i].WURCS;

            tableData.push({
                chain,
                SNFG,
                WURCS,
            });
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

export default function SNFGList(props) {
    const [data, setData] = useState<any[]>([]);
    const columns = useMemo(() => DatabaseColumns, []);
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
        <div className="flex flex-col mx-auto px-16">
            <span className="text-xl mb-2">2D Glycan Descriptions</span>

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
