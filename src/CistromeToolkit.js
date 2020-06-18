import React, { useEffect, useRef, useCallback, useState, useMemo } from "react";
import PubSub from 'pubsub-js';
import { EVENT } from './utils/constants.js';
import DataTable from "./DataTable.js";
import { requestIntervalTFs } from './utils/cistrome.js';
import { CLOSE, PLUS } from './utils/icons.js';
import d3 from './utils/d3.js';
import './CistromeToolkit.scss';
import { REQUEST_HISTORY_SAMPLE } from "./utils/toolkit.js";

export function destroyCistromeToolkit() {
    PubSub.publish(EVENT.CISTROME_TOOLKIT, {
        intervalParams: undefined,
        isVisible: false
    });
}

/**
 * Wrapper around <DataTable />, specific for showing the TF binding interval request results.
 * Subscribes to 'cistrome-toolkit' event via `PubSub`.
 * @prop {function} onAddTrack A function to call when adding tracks with selected rows.
 * @prop {object} intervalParams The interval request parameters.
 * @prop {string} intervalParams.assembly
 * @prop {string} intervalParams.chrStartName
 * @prop {number} intervalParams.chrStartPos
 * @prop {string} intervalParams.chrEndName
 * @prop {number} intervalParams.chrEndPos
 * @example
 * <CistromeToolkit/>
 */
export default function CistromeToolkit(props) {
    const {
        onAddTrack
    } = props;

    const resizerRef = useRef(null);
    const dragY = useRef(null);

    const [isVisible, setIsVisible] = useState(false);
    const [height, setHeight] = useState(600);
    const [intervalParams, setIntervalParams] = useState(undefined);
    const [requestStatus, setRequestStatus] = useState(undefined);
    const [requestHistory, setRequestHistory] = useState([]);   // Use `REQUEST_HISTORY_SAMPLE` in `/utils` to debug.
    const [selectedRequest, setSelectedRequest] = useState(undefined);
    const [selectedRows, setSelectedRows] = useState([]);

    useEffect(() => {
        const cistromeToolkitToken = PubSub.subscribe(EVENT.CISTROME_TOOLKIT, (msg, data) => {
            setIntervalParams(data.intervalParams);
            setIsVisible(data.isVisible ? data.isVisible : !isVisible);
        });

        return () => {
            PubSub.unsubscribe(cistromeToolkitToken);
        };
    });

    useEffect(() => {
        // Reset row selection upon new request selection.
        setSelectedRows([]);
    }, [selectedRequest]);

    useEffect(() => {
        let didUnmount = false;
        if(intervalParams) {
            setRequestStatus({ msg: "Receiving Cistrome DB API response...", isLoading: true });

            const {
                assembly, 
                chrStartName, 
                chrStartPos, 
                chrEndName, 
                chrEndPos
            } = intervalParams;

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

                    const msg = `For interval ${chrStartName}:${chrStartPos}-${chrEndPos}`;
                    setRequestStatus({ msg, isLoading: false });

                    const isNew = requestHistory.find(d => {
                        return d.assembly === assembly 
                            && d.chrStartName === chrStartName && d.chrStartPos === chrStartPos
                            && d.chrEndName === chrEndName && d.chrEndPos === chrEndPos
                    }) === undefined;
                    if(isNew) {
                        const newHistory = [...requestHistory, {
                            parameter: intervalParams,
                            columns: costomColumns,
                            rows: customRows
                        }];
                        setRequestHistory(newHistory);
                        setSelectedRequest(newHistory.length - 1);
                    }
                })
                .catch((msg) => {
                    if(didUnmount) return;
                    setRequestStatus({ msg, isLoading: false });
                });
        }
        return (() => { didUnmount = true; });
    }, [intervalParams]);
    
    // Set up the d3-drag handler functions (started, ended, dragged).
    const started = useCallback(() => {
        const event = d3.event;
        dragY.current = event.sourceEvent.clientY;
    }, [dragY])

    const ended = useCallback(() => {
        dragY.current = null;
    }, [dragY])

    const dragged = useCallback(() => {
        const event = d3.event;
        const diff = event.sourceEvent.clientY - dragY.current;
        setHeight(height - diff);
    }, [dragY]);

    // Detect drag events for the resize element.
    useEffect(() => {
        const resizer = resizerRef.current;

        const drag = d3.drag()
            .on("start", started)
            .on("drag", dragged)
            .on("end", ended);

        d3.select(resizer).call(drag);

        return () => d3.select(resizer).on(".drag", null);
    }, [resizerRef, started, dragged, ended]);

    const listOfResultsRequested = useMemo(() => {
        return requestHistory.map((d, i) => {
            const {
                assembly, 
                chrStartName, 
                chrStartPos,
                chrEndPos
            } = d.parameter;
            return (
                <div 
                    key={JSON.stringify(d.parameter) + i}
                    className={i === selectedRequest ? "cisvis-api-result-selected" : "cisvis-api-result"}
                    onClick={() => { setSelectedRequest(i) }}
                >
                    <div style={{ marginBottom: 4 }}>
                        <span className="cisvis-api-parameter">{"SPECIES"}</span>
                        <b>{assembly}</b>
                    </div>
                    <div>
                        <span className="cisvis-api-parameter">{"INTERVAL"}</span>
                        <b>{`${chrStartName}:${chrStartPos.toLocaleString('en')}-${chrEndPos.toLocaleString('en') }`}</b>
                    </div>
                </div>
            );
        });
    }, [requestHistory, selectedRequest]);
    
    return (
        <div className="cisvis-data-table-bg"
            style={{
                height: isVisible ? `${height}px` : 0
            }}
        >
            <div className='cisvis-data-table-header' ref={resizerRef}/>
            <div className="cisvis-data-table-container">
                <h4 className="cisvis-table-title">
                    CistromeDB Toolkit
                </h4>
                <span className="cisvis-table-subtitle">
                    {requestStatus?.isLoading ? (
                        <span className="cisvis-progress-ring" />
                    ) : null}
                </span>
                {onAddTrack ?
                    <span
                        className={selectedRows && selectedRows.length > 0 
                            ? 'toolkit-btn-add-track-active'
                            : 'toolkit-btn-add-track'}
                        title={selectedRows && selectedRows.length > 0 ? null : 'select rows in the data table'}
                        onClick={() => { 
                            // TODO: Change this to add a track based on actual data after we build DB.
                            // Currently, this button is not shown.
                            onAddTrack(
                                'https://resgen.io/api/v1',
                                'Hygs6CEVR2mCnGlsHK93zQ',
                                'top'
                            );
                            setIsVisible(false);
                        }}
                    >
                        <svg className="chw-button-sm chw-button-static"
                            viewBox={PLUS.viewBox}>
                            <path d={PLUS.path} fill="currentColor"/>
                        </svg>
                        Add HiGlass tracks with selected rows
                    </span>
                : null}
                <span style={{ 
                    verticalAlign: "middle", 
                    display: "inline-block", 
                    position: "absolute", 
                    right: 15, 
                    top: 15
                }}>
                    <svg
                        className={'chw-button'}
                        style={{ color: "gray", background: "none" }}
                        onClick={() => setIsVisible(false)}
                        viewBox={CLOSE.viewBox}
                    >
                        <title>Close data table</title>
                        <path d={CLOSE.path} fill="currentColor"/>
                    </svg>
                </span>
                <div className="cisvis-cistrome-toolkit-body">
                    <div className='cisvis-api-result-header'>Request History</div>
                    <div className='cisvis-api-result-header'>Selected Result Table</div>
                    <div className="cisvis-api-result-list">
                        {listOfResultsRequested}
                    </div>
                    {selectedRequest !== undefined ?
                        <DataTable 
                            columns={requestHistory[selectedRequest].columns}
                            rows={requestHistory[selectedRequest].rows}
                            selectedRows={selectedRows}
                            expoNotations={["Overlap Ratio"]}
                            onCheckRows={undefined}
                            onSelect={(selectedRows) => { setSelectedRows(selectedRows) }}
                        />
                        : null}
                </div>
            </div>
        </div>
    );
}