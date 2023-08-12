import { useState, useEffect, useRef, useMemo } from 'react';
import { useTable } from 'react-table';
import {COLUMNS} from "../../data/Constants"
import styled from 'styled-components'
import {MoorhenContextProvider, MoorhenMolecule, MoorhenContainer, itemReducer} from 'moorhen'

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
        background-color: #F4F9FF;
        // color: #000000
    }
     
    table th {
        padding-top: 12px;
        padding-bottom: 12px;
        text-align: center;
        background-color: #F4F9FF;
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
    #row:hover { 
        scale: 101%;
        cursor: grab;
    }

    
    `

export default function SVGTable({tableData, rowClick, setRowClicked, setRowID}) {
    const [data, setData] = useState(tableData);
    const controls = useRef()

    const columns = useMemo(() => COLUMNS, []);
    const { getTableProps, getTableBodyProps, headerGroups, rows, prepareRow } = useTable({ columns, data });


    const handleRowClick = ((rowId) => { 
        setRowClicked(!rowClick)
        setRowID(rowId)
    })

    return (
        <Styles>
        <div className="container" id='table'>
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
                            <tr {...row.getRowProps()} onClick={() => handleRowClick(row.id)} title="Click to visualise." id='row'>
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