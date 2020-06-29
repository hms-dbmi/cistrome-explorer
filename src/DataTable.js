import React, { useMemo } from 'react';

import "./TrackRowInfoControl.scss";
import "./DataTable.scss";
import { PLUS } from './utils/icons';

/**
 * Component for data table.
 * @prop {array} rows Array of rows to render data table.
 * @prop {array} columns Array of column names in data table.
 * @prop {array} expoNotations A list of columns names that need to use exponential notations.
 * @prop {(function|null)} onButton If a function is provided, a `+` button will be shown for each row.
 */
export default function DataTable(props) {
    const {
        rows = [], 
        columns = [],
        onButton = null,
        expoNotations = []
    } = props;

    const headRow = (
        <tr>
            {onButton && columns.length > 0 ? 
                <th>Add Track</th> 
            : null}
            {columns.map((c, j) => (
                <th key={j}>{c}</th>
            ))}
        </tr>
    );

    const bodyRows = useMemo(() => { 
        return rows.map((d, i) => {
            const buttonCell = (onButton ? (
                <td>
                    <span 
                        style={{ position: 'relative', top: 3 }}
                        onClick={onButton}
                    >
                        <svg className="chw-button"
                            style={{ color: "#808080", background: "none" }}
                            xmlns="http://www.w3.org/2000/svg" width="16" height="16"
                            viewBox={PLUS.viewBox}>
                            <title>Add HiGlass Track</title>
                            <path fill="currentColor" d={PLUS.path}/>
                        </svg>
                    </span>
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
                <tr key={i} className={'data-table-row'}>
                    {buttonCell}{dataCells}
                </tr>
            );
        })
    }, [expoNotations, onButton]);

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
