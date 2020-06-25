import React, { useCallback, useRef, useMemo } from 'react';

import "./TrackRowInfoControl.scss";
import "./DataTable.scss";

/**
 * Component for data table.
 * @prop {array} rows Array of rows to render data table.
 * @prop {array} columns Array of column names in data table.
 * @prop {array} selectedRows Array of rows selected by a user.
 * @prop {array} expoNotations A list of columns names that need to use exponential notations.
 * @prop {(function|null)} onCheckRows If a function is provided, a checkbox will be shown for each row.
 * @prop {function} onSelect A function to call upon the selection on rows in a data table.
 * On change of any checkbox elements, an array of all checked row objects will be passed to the function. By default, null.
 */
export default function DataTable(props) {
    const {
        rows = [], 
        columns = [],
        selectedRows,
        onCheckRows = null,
        expoNotations = [],
        onSelect
    } = props;

    // Store the currently-checked row indices in a mutable set object.
    const checkedRowIndicesRef = useRef(new Set());

    const handleInputChange = useCallback((event) => {
        if(!event || !event.target) return;
        const target = event.target;
        if(target.checked) {
            checkedRowIndicesRef.current.add(target.value);
        } else {
            checkedRowIndicesRef.current.delete(target.value);
        }
        const checkedRows = Array.from(checkedRowIndicesRef.current).map(i => rows[i]);
        onCheckRows(checkedRows);
    }, [rows, columns, onCheckRows]);

    const headRow = (
        <tr>
            {onCheckRows ? (<th></th>) : null}
            {columns.map((c, j) => (
                <th key={j}>{c}</th>
            ))}
        </tr>
    );

    const bodyRows = useMemo(() => { 
        return rows.map((d, i) => {
            const checkboxCell = (onCheckRows ? (
                <td>
                    <input
                        type="checkbox"
                        name="data-table-checkbox"
                        value={i}
                        onChange={handleInputChange}
                    />
                </td>
            ) : null);
            const dataCells = columns.map((c, j) => {
                return (
                    <td key={j}>
                        {expoNotations.includes(c) && +d[c] ? Number.parseFloat(d[c]).toExponential(2) : d[c]}
                    </td>
                );
            });
            return (
                <tr 
                    key={i}
                    className={selectedRows.includes(i) ? 'data-table-row-selected' : 'data-table-row'}
                    onClick={() => onSelect([i])}
                >
                    {checkboxCell}{dataCells}
                </tr>
            );
        })
    }, [selectedRows, expoNotations, onCheckRows]);

    return (
        <div
            style={{
                overflowY: "auto",
                border: "1px solid lightgray",
                background: "white",
                padding: '1px'
            }}
        >
            {bodyRows ? (
                <form>
                    <table className="cisvis-table">
                        <thead>
                            {headRow}
                        </thead>
                        <tbody>
                            {bodyRows}
                        </tbody>
                    </table>
                </form>
            ) : null}
        </div>
    );
}
