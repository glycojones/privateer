import { useMemo, useRef, useEffect, useState } from 'react';
import { useTable } from 'react-table';
import { DatabaseColumns, SugarListColumns } from "../../data/Constants"
import styled from 'styled-components'

function custom_sort(a, b) {
  let split_a = a["Sugar ID"].split("-")
  let split_b = b["Sugar ID"].split("-")

  if (split_a[1] < split_b[1]) {
    return -1
  }
  if (split_b[1] < split_a[1]) {
    return 1;
  }

  if (Number(split_a[2]) < Number(split_b[2])) {
    return -1;
  }

  if (Number(split_b[2]) < Number(split_a[2])) {
    return 1;
  }

  return 0;
}


function parse_results(data) {
  let glycans = data.data.glycans

  let table_data = []

  for (const key in glycans) {
    let glycan_type = glycans[key]
    for (let i = 0; i < glycan_type.length; i++) {
      let sugars = glycan_type[i].Sugars
      for (let j = 0; j < sugars.length; j++) {
        let keys = Object.keys(sugars[j])
        let entry = {}

        for (let k = 0; k < keys.length; k++) {
          let type_key = typeof (sugars[j][keys[k]])

          if (type_key === "number") {
            entry[keys[k]] = sugars[j][keys[k]].toFixed(2)
          }
          else {
            entry[keys[k]] = sugars[j][keys[k]]
          }
        }
        entry["type"] = key
        table_data.push(entry)
      }
    }
  }

  table_data.sort(custom_sort)
  return table_data
}

let Styles = styled.div`
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

`

export default function SugarList(props) {

  const [data, setData] = useState([])
  const columns = useMemo(() => SugarListColumns, []);
  const { getTableProps, getTableBodyProps, headerGroups, rows, prepareRow } = useTable({ columns, data });

  useEffect(() => {
    let results = parse_results(props)
    setData(results)
  }, [])

  return (
    <div className='flex flex-col mx-auto p-16'>
      <span className='text-xl mb-2'>Detailed monosaccharide validation data</span>

      <Styles>
        <table {...getTableBodyProps()}>
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
                <tr {...row.getRowProps()}
                  title="Click to visualise." id='row'>
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
      </Styles>
    </div>
  )
}