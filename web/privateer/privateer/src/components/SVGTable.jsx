import { useState, useEffect, useMemo } from 'react';
import { useTable } from 'react-table';
import {COLUMNS} from "./Constants"
import styled from 'styled-components'

export default function SVGTable({tableData}) {
    const [data, setData] = useState(tableData);
    // useeffect ...

    // console.log(tableData)

    const columns = useMemo(() => COLUMNS, []);
    const { getTableProps, getTableBodyProps, headerGroups, rows, prepareRow } = useTable({ columns, data });

    const Styles = styled.div`
    table {
        border-collapse: collapse;
        border-spacing: 0;
        width: 100%;
        border: 1px solid #ddd;
    }
     
    table th{ 
        text-align: left;
        padding: 16px;
        border: 1px solid #ddd;
    }
    table td {
        text-align: center;
        padding: 16px;
        border: 1px solid #ddd;
    }
     
    // table tr:nth-child(even) {
    //     background-color: #f2f2f2;
    //     color: #000000
    // }
     
    table tr:hover {
        // background-color: #ddd;
        scale: 101%;
        cursor: grab;
    }
     
    table th {
        padding-top: 12px;
        padding-bottom: 12px;
        text-align: center;
        background-color: #4CAF50;
        color: white;
    }
    `

    return (
        <Styles>
        <div className="container">
            <table {...getTableProps()}>
                <thead>
                    {headerGroups.map((headerGroup) => (
                        <tr {...headerGroup.getHeaderGroupProps()}>
                            {headerGroup.headers.map((column) => (
                                <th {...column.getHeaderProps()}>
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
                            <tr {...row.getRowProps()} onClick={() => {console.log(row.id)}}>
                                {row.cells.map((cell) => {
                                
                                    return (
                                        <td {...cell.getCellProps()}>
                                            {cell.render('Cell')}
                                        </td>
                                    );
                                })}
                            </tr>
                        );
                    })}
                </tbody>
            </table>
        </div>
        </Styles>
    );
};