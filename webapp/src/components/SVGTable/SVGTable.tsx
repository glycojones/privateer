import React, { type ReactElement, useMemo } from 'react';
import { useTable } from 'react-table';
import { COLUMNS } from '../../data/Constants';
import styled from 'styled-components';
import { type SVGTableProps } from '../../interfaces/types';

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

  #row:hover {
    scale: 101%;
    cursor: grab;
  }
`;

export default function SVGTable (props: SVGTableProps): ReactElement {
  // const [data, _setData] = useState(props.tableData);
  // const controls = useRef();
  const data = props.tableData;

  const columns = useMemo(() => COLUMNS, []);
  const { getTableProps, getTableBodyProps, headerGroups, rows, prepareRow } =
    useTable({ columns, data });

  function handleRowClick (rowId: number): void {
    if (props.allowRowClick) {
      props.setRowClicked(!props.rowClick);
      props.setRowID(rowId);
    }
  }

  return (
    <div className="flex flex-col mx-auto px-16">
      {/* @ts-expect-error */}
      <Styles $allowRowClick={props.allowRowClick}>
        <div className="container mx-auto" id="table">
          <table {...getTableProps()}>
            <thead>
              {headerGroups.map((headerGroup) => (
                <tr
                  {...headerGroup.getHeaderGroupProps({
                    style: { width: '200px' }
                  })}
                  key={headerGroup.name}
                >
                  {headerGroup.headers.map((column) => (
                    <th {...column.getHeaderProps()} key={column.name}>
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
                    onClick={() => {
                      if (props.allowRowClick) {
                        handleRowClick(row.id as number);
                      }
                    }}
                    title="Click to visualise."
                    id="row"
                    key={row.name}
                  >
                    {row.cells.map((cell) => {
                      return (
                        <td {...cell.getCellProps()} key={cell.name}>
                          {cell.render('Cell')}{' '}
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
    </div>
  );
}
