import React, { useEffect } from 'react';

import "./TrackRowInfoControl.scss";

export default function DataTable(props) {
    const {
        left, top, width, height,
        data, columns
    } = props;

    const head = columns.map((c, j) => {
        return (
            <th key={j}>{c}</th>
        );
    });
    let rows = data.map((d, i) => {
        const cells = columns.map((c, j) => {
            return (
                <td key={j}>{d[c]}</td>
            );
        });

        return (
            <tr key={i}>{cells}</tr>
        );
    });

    return (
        <div
            style={{
                position: "absolute",
                left, top, width, height,
                border: "1px dotted red"
            }}
        >
            <h4>Data Preview (~100 rows)</h4>   
            <div
                style={{
                    height: height - 40,
                    overflowY: "auto"
                }}
            >
                {rows ? 
                    <table className="chw-table">
                        <tbody>
                            <tr>{head}</tr>
                            {rows}
                        </tbody>
                    </table>
                : null}
            </div>
        </div>
    );
}