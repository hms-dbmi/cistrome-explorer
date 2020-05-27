import React, { useEffect, useState } from "react";

import DataTable from "./DataTable.js";

import { requestIntervalTFs } from './utils/cistrome.js';

/**
 * Wrapper around <DataTable />, specific for showing the TF binding interval request results.
 * @prop {object} intervalParams The interval request parameters.
 * @prop {string} intervalParams.assembly
 * @prop {string} intervalParams.chrStartName
 * @prop {number} intervalParams.chrStartPos
 * @prop {string} intervalParams.chrEndName
 * @prop {number} intervalParams.chrEndPos
 */
export default function DataTableForIntervalTFs(props) {
    const {
        intervalParams
    } = props;

    const {
        assembly, 
        chrStartName, 
        chrStartPos, 
        chrEndName, 
        chrEndPos
    } = intervalParams;

    const [isVisible, setIsVisible] = useState(false);
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

                    // TODO: Consider supporting multiple API queries.
                    const customColumnMap = {
                        GSM: "GEO/ENCODE ID",
                        DCid: "CistromeDB ID",
                        factor: "Factor",
                        cellLine: "Cell Line",
                        CellType: "Cell Type",
                        species: "Species",
                        OverlapRatio: "Overlap Ratio",
                        OverlapPeakNumber: "Overlap Peak Number",
                    }
                    const customRows = rows.map(r => {
                        const newRow = {};
                        Object.keys(customColumnMap).forEach(k => {
                            if(customColumnMap[k]) {
                                newRow[customColumnMap[k]] = r[k];
                            }
                        });
                        return newRow;
                    });
                    const costomColumns = Object.values(customColumnMap);
                    console.log(customRows);
                    console.log(costomColumns);
                    setDataTableRows(customRows);
                    setDataTableColumns(costomColumns);

                    const msg = `For interval ${chrStartName}:${chrStartPos}-${chrEndPos}`;
                    setRequestStatus({ msg, isLoading: false });
                })
                .catch((msg) => {
                    if(didUnmount) return;
                    setRequestStatus({ msg, isLoading: false });
                })
        }
        setIsVisible(true);
        return (() => { didUnmount = true; });
    }, [intervalParams]);    

    return (requestStatus && isVisible ? (
        <div className="cisvis-data-table-bg">
            <DataTable 
                margin={100}
                title={"Factors from Cistrome DB"}
                subtitle={requestStatus.msg}
                isLoading={requestStatus.isLoading}
                rows={dataTableRows}
                columns={dataTableColumns}
                expoNotations={["Overlap Ratio"]}
                onCheckRows={undefined}
                onClose={() => setIsVisible(false)}
            />
        </div>
    ) : null);
}