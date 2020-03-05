import React, { useEffect, useState } from "react";

import DataTable from "./DataTable.js";

import { requestIntervalTFs } from './utils/cistrome.js';

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
        intervalParams
    } = props;

    const {
        assembly, 
        chrStartName, 
        chrStartPos, 
        chrEndName, 
        chrEndPos
    } = intervalParams;

    const [requestStatus, setRequestStatus] = useState(null);
    const [dataTableRows, setDataTableRows] = useState([]);
    const [dataTableColumns, setDataTableColumns] = useState([]);

    useEffect(() => {
        let didUnmount = false;
        setRequestStatus({ msg: "Receiving Cistrome DB API response...", isLoading: true });
        if(intervalParams) {
            requestIntervalTFs(assembly, chrStartName, chrStartPos, chrEndName, chrEndPos)
                .then(([rows, columns]) => {
                    if(didUnmount) return;
                    setDataTableRows(rows);
                    setDataTableColumns(columns);

                    const msg = `For interval ${chrStartName}:${chrStartPos}-${chrEndPos}`;
                    setRequestStatus({ msg, isLoading: false });
                })
                .catch((msg) => {
                    if(didUnmount) return;
                    setRequestStatus({ msg, isLoading: false });
                })
        }

        return (() => { didUnmount = true; });
    }, [intervalParams]);

    return (requestStatus ? (
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
        />
    ) : null);
}