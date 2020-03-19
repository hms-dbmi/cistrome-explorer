import React, { useEffect, useState, useRef } from "react";

import DataTable from "./DataTable.js";

import { requestIntervalTFs } from './utils/cistrome.js';

/**
 * Merge all api responses into a single set of rows and columns.
 * @param {object[]} allApiResponses Array of api response objects, containing the properties `success`, `rows`, `columns`.
 * @returns {array} The result is an array `[rows, columns]` containing the merged values.
 */
function mergeReponses(allApiResponses) {
    // Only use the successful responses.
    const apiResponses = allApiResponses.filter(d => d.success);
    let rows = [];
    let columns = [];
    const minNumRows = Math.min(...apiResponses.map(d => d.rows.length));
    // Use all possible column values.
    for(let apiResponse of apiResponses) {
        columns = [ ...apiResponse.columns, ...columns ];
    }
    // For each row, compute the product of OverlapRatios across responses.
    for(let i = 0; i < minNumRows; i++) {
        let rankProduct = 1;
        for(let apiResponse of apiResponses) {
            rankProduct *= apiResponse.rows[i]['OverlapRatio'];
        }
        for(let apiResponse of apiResponses) {
            const row = apiResponse.rows[i];
            row['OverlapRatioRankProduct'] = rankProduct;
            row['Rank'] = i;
            // Append this row to the output row array.
            rows.push(row);
        }
    }
    // Sort all rows by the rank product value rather than the overlap ratio value.
    rows.sort((a, b) => {
        return b['OverlapRatioRankProduct'] - a['OverlapRatioRankProduct'];
    })
    columns.push('OverlapRatioRankProduct');
    columns.push('Rank');
    // Return the processed rows array, and all unique column names.
    return [rows, Array.from(new Set(columns))];
}

/**
 * Generate a name from an object of interval parameters.
 * @param {object} intervalParams The interval parameters.
 * @returns {string} A name for the parameters object.
 */
function getIntervalName(intervalParams) {
    const {
        chrStartName, 
        chrStartPos, 
        chrEndPos
    } = intervalParams;
    return `${chrStartName}:${chrStartPos}-${chrEndPos}`;
}

/**
 * Wrapper around <DataTable />, specific for showing the TF binding interval request results.
 * @prop {number} left The left position of this view.
 * @prop {number} top The top position of this view.
 * @prop {number} width The width of this view.
 * @prop {number} height The height of this view.
 * @prop {object} intervalParams The interval request parameters.
 * @prop {string} intervalParams.assembly
 * @prop {string} intervalParams.chrStartName
 * @prop {number} intervalParams.chrStartPos
 * @prop {string} intervalParams.chrEndName
 * @prop {number} intervalParams.chrEndPos
 */
export default function DataTableForIntervalTFs(props) {
    const {
        left, top, width, height,
        allIntervalParams
    } = props;

    const [isVisible, setIsVisible] = useState(false);
    const [requestStatus, setRequestStatus] = useState(null);
    const [dataTableRows, setDataTableRows] = useState([]);
    const [dataTableColumns, setDataTableColumns] = useState([]);

    // Store both parameters and the response results, before merging.
    const apiResponses = useRef();

    useEffect(() => {
        let didUnmount = false;
        if(allIntervalParams.length > 0) {
            setDataTableColumns([]);
            setDataTableRows([]);
            setRequestStatus({ msg: `Loading Cistrome DB Toolkit API response for ${allIntervalParams.length} intervals...`, isLoading: true });
            for(let intervalParams of allIntervalParams) {
                const {
                    assembly, 
                    chrStartName, 
                    chrStartPos, 
                    chrEndName, 
                    chrEndPos
                } = intervalParams;

                const intervalName = getIntervalName(intervalParams);
    
                apiResponses.current = {};
                requestIntervalTFs(assembly, chrStartName, chrStartPos, chrEndName, chrEndPos)
                    .then(([rows, columns]) => {
                        if(didUnmount) return;

                        rows.sort((a, b) => {
                            return b['OverlapRatio'] - a['OverlapRatio'];
                        });

                        apiResponses.current[intervalName] = {
                            rows, columns, success: true
                        };

                        // Merge the rows from all api responses if all responses have finished (either successfully or unsuccessfully).
                        if(Object.keys(apiResponses.current).length === allIntervalParams.length) {
                            const [allRows, allColumns] = mergeReponses(Object.values(apiResponses.current));
                            setDataTableColumns(allColumns);
                            setDataTableRows(allRows);
                            
                            const allIntervalNames = Object.keys(apiResponses.current).map(k => (apiResponses.current[k].success ? k : `${k} (failed)`)).join(', ');
                            const msg = `For interval${allIntervalParams.length === 1 ? '' : 's'} ${allIntervalNames}`;
                            setRequestStatus({ msg, isLoading: false });
                        }
                    })
                    .catch((msg) => {
                        if(didUnmount) return;
                        apiResponses.current[intervalName] = {
                            success: false
                        };
                        setRequestStatus({ msg, isLoading: false });
                    });
            }
            setIsVisible(true);
        }
        return (() => { didUnmount = true; });
    }, [allIntervalParams]);

    return (requestStatus && isVisible ? (
        <DataTable 
            left={left}
            top={top}
            width={width}
            height={height}
            title={"TFs from Cistrome DB"}
            subtitle={requestStatus.msg}
            isLoading={requestStatus.isLoading}
            rows={dataTableRows}
            columns={dataTableColumns}
            onCheckRows={() => {}}
            onClose={() => setIsVisible(false)}
        />
    ) : null);
}