import React, { useCallback, useRef } from 'react';

import "./TrackRowInfoControl.scss";
import "./DataTable.scss";
import { CLOSE } from './utils/icons';

/**
 * Component for data table.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {array} rows Array of rows to render data table.
 * @prop {array} columns Array of column names in data table.
 * @prop {string} title A title for the table. Optional.
 * @prop {string} subtitle A subtitle for the table. Optional.
 * @prop {boolean} isLoading Whether the data is still loading, in which case show a spinner.
 * @prop {(function|null)} onCheckRows If a function is provided, a checkbox will be shown for each row.
 * On change of any checkbox elements, an array of all checked row objects will be passed to the function. By default, null.
 */
export default function DataTable(props) {
    const {
        left, top, width, height,
        rows = [], columns = [],
        title = "Data Preview",
        subtitle,
        isLoading = false,
        onCheckRows = null,
        onClose
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

    const bodyRows = rows.map((d, i) => {
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
            return <td key={j}>{d[c]}</td>;
        });
        return <tr key={i}>{checkboxCell}{dataCells}</tr>;
    });

    return (
        <div
            style={{
                position: "absolute",
                left, top, width, height
            }}
        >
            <h4 className="chw-table-title">{title}</h4>
            <span style={{ verticalAlign: "middle", display: "inline-block" }}>
                <svg
                    className={`chw-button`}
                    style={{ color: "gray" }}
                    onClick={() => onClose()} 
                    viewBox={CLOSE.viewBox}
                >
                    <title>Close data table</title>
                    <path d={CLOSE.path} fill="currentColor"/>
                </svg>
            </span>
            <span className="chw-table-subtitle">
                {isLoading ? (
                    <span className="chw-progress-ring" />
                ) : (subtitle ? (
                    <b>{subtitle}</b>
                ) : null)}
            </span>
            <div
                style={{
                    height: height - 40,
                    overflowY: "auto"
                }}
            >
                {bodyRows ? (
                    <form>
                        <table className="chw-table">
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
        </div>
    );
}
