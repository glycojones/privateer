import React, { useState, useEffect, useMemo } from 'react';
import {
    type PaginationState,
    useReactTable,
    getCoreRowModel,
    getFilteredRowModel,
    getPaginationRowModel,
    flexRender,
    getSortedRowModel,
    type Column, SortingState,
} from "@tanstack/react-table"
import styled from 'styled-components';

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

function Filter({ column, table }: { column: Column<any, any>; table: any }) {
    const firstValue = table
        .getPreFilteredRowModel()
        .flatRows[0]?.getValue(column.id);

    const columnFilterValue = column.getFilterValue();

    return typeof firstValue === 'number' ? (
        <div className="flex space-x-2">
            <input
                type="number"
                value={(columnFilterValue as [number, number])?.[0] ?? ''}
                onChange={(e) => {
                    column.setFilterValue((old: [number, number]) => [
                        e.target.value,
                        old?.[1],
                    ]);
                }}
                placeholder={' Min'}
                className="w-24 border shadow rounded"
            />
            <input
                type="number"
                value={(columnFilterValue as [number, number])?.[1] ?? ''}
                onChange={(e) => {
                    column.setFilterValue((old: [number, number]) => [
                        old?.[0],
                        e.target.value,
                    ]);
                }}
                placeholder={' Max'}
                className="w-24 border shadow rounded"
            />
        </div>
    ) : (
        <input
            type="text"
            value={(columnFilterValue ?? '') as string}
            onChange={(e) => {
                column.setFilterValue(e.target.value);
            }}
            placeholder={'Search...'}
            className="w-36 border shadow rounded"
        />
    );
}
function Table({ data }: { data: any }) {
    const COLUMNS = [
        {
            header: 'Type',
            accessorKey: 'type',
        },
        {
            header: 'PDB',
            accessorKey: 'pdb',
        },
        {
            header: 'Count',
            accessorKey: 'count',
        },
        {
            header: 'Resolution',
            accessorKey: 'resolution',
        },
        {
            header: 'Link',
            accessorKey: 'link',
            enableColumnFilter: false,
            cell: (props: { getValue: () => string }) => {
                return (
                    <a href={props.getValue()}>
                        <svg
                            xmlns="http://www.w3.org/2000/svg"
                            viewBox="0 0 512 512"
                            className="w-4 mx-auto"
                        >
                            <path d="M320 0c-17.7 0-32 14.3-32 32s14.3 32 32 32h82.7L201.4 265.4c-12.5 12.5-12.5 32.8 0 45.3s32.8 12.5 45.3 0L448 109.3V192c0 17.7 14.3 32 32 32s32-14.3 32-32V32c0-17.7-14.3-32-32-32H320zM80 32C35.8 32 0 67.8 0 112V432c0 44.2 35.8 80 80 80H400c44.2 0 80-35.8 80-80V320c0-17.7-14.3-32-32-32s-32 14.3-32 32V432c0 8.8-7.2 16-16 16H80c-8.8 0-16-7.2-16-16V112c0-8.8 7.2-16 16-16H192c17.7 0 32-14.3 32-32s-14.3-32-32-32H80z" />
                        </svg>
                    </a>
                );
            },
        },
    ];

    useEffect(() => {
        console.warn(data);
    }, []);

    const columns = useMemo(() => COLUMNS, []);
    const [pagination, setPagination] = useState<PaginationState>({
        pageIndex: 0,
        pageSize: 10,
    });
    const [sorting, setSorting] = useState<SortingState>([
        {
            id: "count", // Must be equal to the accessorKey of the coulmn you want sorted by default
            desc: true,
        },
    ])
    const table = useReactTable({
        columns,
        data,
        debugTable: false,
        getCoreRowModel: getCoreRowModel(),
        getSortedRowModel: getSortedRowModel(),
        getFilteredRowModel: getFilteredRowModel(),
        getPaginationRowModel: getPaginationRowModel(),
        onPaginationChange: setPagination,
        state: { pagination, sorting },
    });

    return (
        <>
            <Styles>
                <table>
                    <thead>
                        {table.getHeaderGroups().map((headerGroup) => (
                            <tr key={headerGroup.id}>
                                {headerGroup.headers.map((header) => {
                                    return (
                                        <th
                                            key={header.id}
                                            colSpan={header.colSpan}
                                        >
                                            <div
                                                {...{
                                                    className:
                                                        header.column.getCanSort()
                                                            ? 'cursor-pointer select-none'
                                                            : '',
                                                    onClick:
                                                        header.column.getToggleSortingHandler(),
                                                }}
                                            >
                                                {flexRender(
                                                    header.column.columnDef
                                                        .header,
                                                    header.getContext()
                                                )}
                                                {{
                                                    asc: ' ↑',
                                                    desc: ' ↓',
                                                }[
                                                    header.column.getIsSorted() as string
                                                ] ?? null}
                                                {header.column.getCanFilter() ? (
                                                    <div>
                                                        <Filter
                                                            column={
                                                                header.column
                                                            }
                                                            table={table}
                                                        />
                                                    </div>
                                                ) : null}
                                            </div>
                                        </th>
                                    );
                                })}
                            </tr>
                        ))}
                    </thead>
                    <tbody>
                        {table.getRowModel().rows.map((row) => {
                            return (
                                <tr key={row.id}>
                                    {row.getVisibleCells().map((cell) => {
                                        return (
                                            <td key={cell.id}>
                                                {flexRender(
                                                    cell.column.columnDef.cell,
                                                    cell.getContext()
                                                )}
                                            </td>
                                        );
                                    })}
                                </tr>
                            );
                        })}
                    </tbody>
                </table>
                <div className="h-2" />
                <div className="flex items-center justify-center w-full gap-2">
                    <button
                        className="border rounded p-1"
                        onClick={() => {
                            table.firstPage();
                        }}
                        disabled={!table.getCanPreviousPage()}
                    >
                        {'<<'}
                    </button>
                    <button
                        className="border rounded p-1"
                        onClick={() => {
                            table.previousPage();
                        }}
                        disabled={!table.getCanPreviousPage()}
                    >
                        {'<'}
                    </button>
                    <button
                        className="border rounded p-1"
                        onClick={() => {
                            table.nextPage();
                        }}
                        disabled={!table.getCanNextPage()}
                    >
                        {'>'}
                    </button>
                    <button
                        className="border rounded p-1"
                        onClick={() => {
                            table.lastPage();
                        }}
                        disabled={!table.getCanNextPage()}
                    >
                        {'>>'}
                    </button>
                    <span className="flex items-center gap-1">
                        <div>Page</div>
                        <strong>
                            {table.getState().pagination.pageIndex + 1} of{' '}
                            {table.getPageCount().toLocaleString()}
                        </strong>
                    </span>
                    <span className="flex items-center gap-1">
                        | Go to page:
                        <input
                            type="number"
                            defaultValue={
                                table.getState().pagination.pageIndex + 1
                            }
                            onChange={(e) => {
                                const page = e.target.value
                                    ? Number(e.target.value) - 1
                                    : 0;
                                table.setPageIndex(page);
                            }}
                            className="border p-1 rounded w-16"
                        />
                    </span>
                    <select
                        value={table.getState().pagination.pageSize}
                        onChange={(e) => {
                            table.setPageSize(Number(e.target.value));
                        }}
                    >
                        {[10, 20, 30, 40, 50].map((pageSize) => (
                            <option key={pageSize} value={pageSize}>
                                Show {pageSize}
                            </option>
                        ))}
                    </select>
                </div>
                <div className="w-full flex justify-center">
                    Showing {table.getRowModel().rows.length.toLocaleString()}{' '}
                    of {table.getRowCount().toLocaleString()} Rows
                </div>
            </Styles>
        </>
    );
}

export default function DatabaseSearchTable(props: { data: any }) {
    return <Table data={props.data} />;
}
