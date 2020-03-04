import React, { useEffect } from 'react';

import "./TrackRowInfoControl.scss";

export default function DataTable(props) {
    const {
        left, top, width, height,
        rows, columns
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
            <h4>Data Preview</h4>   
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