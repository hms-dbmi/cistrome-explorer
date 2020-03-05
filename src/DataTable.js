import React from 'react';

import "./TrackRowInfoControl.scss";
import "./DataTable.scss";

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
 */
export default function DataTable(props) {
    const {
        left, top, width, height,
        rows = [], columns = [],
        title = "Data Preview",
        subtitle,
        isLoading = false
    } = props;

    const head = columns.map((c, j) => {
        return <th key={j}>{c}</th>;
    });

    let tableRows = rows.map((d, i) => {
        const cells = columns.map((c, j) => {
            return <td key={j}>{d[c]}</td>;
        });
        return <tr key={i}>{cells}</tr>;
    });

    return (
        <div
            style={{
                position: "absolute",
                left, top, width, height
            }}
        >
            <h4 className="chw-table-title">{title}</h4>
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
                {tableRows ? 
                    <table className="chw-table">
                        <tbody>
                            <tr>{head}</tr>
                            {tableRows}
                        </tbody>
                    </table>
                : null}
            </div>
        </div>
    );
}