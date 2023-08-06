import { useState, useMemo } from "react";
import { useTable } from 'react-table';
import styled from 'styled-components'


const Styles = styled.div`
    table {
        border-collapse: collapse;
        border-spacing: 0;
        width: 100%;
        border: 1px solid #ddd;
        margin
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
     
    table th {
        padding-top: 12px;
        padding-bottom: 12px;
        text-align: center;
        background-color: #4CAF50;
        color: white;
    }
    `
export default function GlycanDetailTable({row}) { 
    const DATA = {}
    DATA["chain"] = row.chain
    DATA['id'] = row.id
    DATA['glytoucan_id'] = row.glytoucan_id

    console.log(DATA)

    const [data, setData] = useState([DATA]);

    const COLUMNS = [
        {
            Header: 'Chain',
            accessor: 'chain',
        },
        {
            Header: 'ID',
            accessor: 'id',
        },
        {
            Header: 'GlytoucanID',
            accessor: 'glytoucan_id',
        },
    ]

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

}