import React, { useMemo, useState, useEffect } from 'react';
import "./DataTable.scss";
import { PLUS, SORT_ASC, SORT_ASC_SIMPLE, SORT_DESC, SORT_DESC_SIMPLE, SQUARE_CHECK } from '../utils/icons';

const ROW_DISPLAY_LIMIT = 100;

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

    const [filterByField, setFilterByField] = useState({});
    const [sortByField, setSortByField] = useState({field: undefined, isAscending: false});
    const [transformedRows, setTransformedRows] = useState(transform(rows));
    const [uniqueFactors, setUniqueFactors] = useState(Array.from(new Set(transform(rows).map(d => d['Factor']))));
    const [selectedFactor, setSelectedFactor] = useState(Array.from(new Set(transform(rows).map(d => d['Factor'])))[0]);

    useEffect(() => {
        setTransformedRows(transform(rows));
        setFilterByField({});
    }, [rows]);

    useEffect(() => {
        setTransformedRows(transform(rows));
    }, [filterByField, sortByField]);

    function transform(data) {
        let transformed = Array.from(data);
        
        // filter
        Object.keys(filterByField)?.forEach(key => {
            transformed = transformed.filter(d => {
                // console.log(filterByKey[key]);
                const condition = filterByField[key];
                if(typeof condition === 'number') {
                    return !isNaN(+d[key]) && condition <= (+d[key]);
                }
                else {
                    return d[key]?.toString().toUpperCase().includes(filterByField[key].toUpperCase());
                }
            });
        });

        // sort
        if(sortByField.field) {
            const {field, isAscending} = sortByField;
            transformed = transformed.sort((a, b) => isAscending ? a[field] > b[field] : a[field] < b[field]);
        }

        return transformed;
    }

    useEffect(() => {
        const newUniqueFactors = Array.from(new Set(transformedRows.map(d => d["Factor"])));
        setUniqueFactors(newUniqueFactors);
        setSelectedFactor(newUniqueFactors[0]);
    }, [transformedRows]);

    const headWithFilterRow = (
        <>
            <tr>
                {columns.map((c, i) => (
                    <th key={i}>
                        {c}
                    </th>
                ))}
            </tr>
            <tr>
                {columns.map((c, i) => (
                    <td key={c+i}>
                        {['Regulatory Potential', 'Overlap Peak Number', 'Overlap Ratio'].indexOf(c) === -1 || true ? 
                            <input
                                className="hm-filter-box-input"
                                type="text"
                                name="default name"
                                placeholder={`Search`}
                                onChange={(e) => {
                                    setFilterByField({
                                        ...filterByField,
                                        [c]: e.target.value
                                    })
                                }}
                                style={{
                                    width: 100,
                                }}
                            /> 
                            : null
                        }
                        <div>
                            <svg className={
                                    sortByField.field === c && sortByField.isAscending 
                                        ? "hm-button-sm-hl"
                                        : "hm-button-sm"
                                }
                                onClick={() => {
                                    setSortByField({
                                        field: c,
                                        isAscending: true
                                    });
                                }} 
                                viewBox={SORT_ASC.viewBox}>
                                <title>Remove all filters</title>
                                <path d={SORT_ASC.path} fill="currentColor"/>
                            </svg>
                            <svg className={
                                    sortByField.field === c && !sortByField.isAscending
                                        ? "hm-button-sm-hl"
                                        : "hm-button-sm"
                                }
                                onClick={() => {
                                    setSortByField({
                                        field: c,
                                        isAscending: false
                                    });
                                }} 
                                viewBox={SORT_DESC.viewBox}>
                                <title>Remove all filters</title>
                                <path d={SORT_DESC.path} fill="currentColor"/>
                            </svg>
                        </div>
                        {/* <svg
                            style={{ color: "rgb(171, 171, 171)", width: 14, height: 14 }}
                            viewBox={SORT_ASC_SIMPLE.viewBox}
                        >
                            <title>Sort</title>
                            <path d={SORT_ASC_SIMPLE.path} fill="black"/>
                        </svg>
                        <svg
                            style={{ color: "rgb(171, 171, 171)", width: 14, height: 14 }}
                            viewBox={SORT_DESC_SIMPLE.viewBox}
                        >
                            <title>Sort</title>
                            <path d={SORT_DESC_SIMPLE.path} fill="black"/>
                        </svg> */}
                    </td>
                ))}
            </tr>
        </>
    );

    const bodyRows = useMemo(() => { 
        // DEBUG
        // console.log(new Set(rows.map(d => d['Factor'])));

        const dataRows = transformedRows.slice(0, ROW_DISPLAY_LIMIT).map((d, i) => {
            const buttonCell = (onButton ? (
                <td>
                    <span 
                        style={{ position: 'relative', top: 3 }}
                        onClick={() => onButton(transformedRows[i])}
                    >
                        <svg className="hm-button"
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
                    {/* {buttonCell} */}
                    {dataCells}
                </tr>
            );
        });

        return [
            ...dataRows
        ];
    }, [transformedRows, expoNotations, onButton]);

    return (
        <>
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
                                {headWithFilterRow}
                            </thead>
                            <tbody>
                                {bodyRows}
                            </tbody>
                        </table>
                    </form>
                ) : null}
            </div>
            <div
                style={{
                    paddingTop: 16,
                    gridColumn: 3,
                    width: '100%',
                    display: 'grid',
                    gridTemplateColumns: 'auto auto',
                    gridColumnGap: 6
                }}
            >
                <span 
                    style={{
                        textAlign: 'left',
                        fontWeight: 600,
                        color: 'gray',
                        cursor: 'not-allowed'
                    }}
                    // onClick={() => onButton({
                    //     Factor: selectedFactor,
                    //     Species: 'Homo sapiens'
                    // })}
                >
                    <select 
                        onChange={e => { setSelectedFactor(e.target.value) }} 
                        defaultValue={uniqueFactors[0]}
                    >
                        {uniqueFactors.map(factor => (
                            <option 
                                key={factor} 
                                value={factor} 
                            >
                                {factor}
                            </option>
                        ))}
                    </select>
                    &nbsp;&nbsp;
                    {`  Add samples to visualization`}
                </span>
                <div 
                    style={{
                        textAlign: 'right',
                        fontWeight: 600
                    }}
                >
                    {`Showing the first ${
                        transformedRows.length < ROW_DISPLAY_LIMIT
                            ? transformedRows.length
                            : ROW_DISPLAY_LIMIT
                        } results of ${transformedRows.length.toLocaleString()} rows`}
                </div>
            </div>
        </>
    );
}
