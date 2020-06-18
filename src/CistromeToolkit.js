import React, { useEffect, useRef, useCallback, useState, useMemo } from "react";
import PubSub from 'pubsub-js';
import d3 from './utils/d3.js';
import { EVENT } from './utils/constants.js';
import DataTable from "./DataTable.js";
import { CLOSE, SEARCH, PLUS } from './utils/icons.js';
import { 
    CISTROME_DBTOOLKIT_SPECIES, 
    CISTROME_DBTOOLKIT_PEAK_NUMBERS, 
    CISTROME_DBTOOLKIT_GENE_DISTANCE,
    requestIntervalTFs, 
    requestGeneTFs, 
    validateGeneParams 
} from './utils/cistrome.js';
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

    const [isVisible, setIsVisible] = useState(true); // TODO: Debugging,...
    const [height, setHeight] = useState(800);
    const [requestStatus, setRequestStatus] = useState(undefined);
    const [requestHistory, setRequestHistory] = useState([]);   // Use `REQUEST_HISTORY_SAMPLE` in `/utils` to debug.

    // API parameters
    const [latestIntervalParams, setLatestIntervalParams] = useState({
        assembly: CISTROME_DBTOOLKIT_SPECIES[0],
        chrStartName: undefined,
        chrStartPos: undefined,
        chrEndName: undefined,
        chrEndPos: undefined
    });
    const [latestGeneParams, setLatestGeneParams] = useState({
        assembly: CISTROME_DBTOOLKIT_SPECIES[0],
        distance: CISTROME_DBTOOLKIT_GENE_DISTANCE[0],
        gene: ''
    });
    const [latestPeaksetParams, setLatestPeaksetParams] = useState({
        // TODO: Support this
    });
    const [isLatestGeneParamsReady, setIsLatestGeneParamsReady] = useState(false);

    // Selected rows
    const [selectedAPI, setSelectedAPI] = useState(undefined);  // 0: Interval, 1: Gene, 2: Peak
    const [selectedRequest, setSelectedRequest] = useState(undefined);
    const [selectedRows, setSelectedRows] = useState([]);

    useEffect(() => {
        const cistromeToolkitToken = PubSub.subscribe(EVENT.CISTROME_TOOLKIT, (msg, data) => {
            setLatestIntervalParams(data.intervalParams);
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

    // Check if parameters for each API is ready and can make a good URL
    // useEffect(() => {
    //     setIsLatestIntervalParamsReady(validateIntervalParams().success);
    // }, [latestIntervalParams]);
    useEffect(() => {
        setIsLatestGeneParamsReady(validateGeneParams(latestGeneParams).success);
    }, [latestGeneParams]);

    useEffect(() => {
        let didUnmount = false;
        if(latestIntervalParams) {
            setRequestStatus({ msg: "Receiving Cistrome DB API response...", isLoading: true });

            const {
                assembly, 
                chrStartName, 
                chrStartPos, 
                chrEndName, 
                chrEndPos
            } = latestIntervalParams;

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
                            parameter: latestIntervalParams,
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
    }, [latestIntervalParams]);
    
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
        setHeight(
            Math.max(
                Math.min(height - diff, window.innerHeight),
                400
            )
        );
    }, [dragY, height]);

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

    const ApiRequestConfigViews = useMemo(() => {
        const speciesSelection = (
            <>
                <div>Species</div>
                <select
                    defaultValue={CISTROME_DBTOOLKIT_SPECIES[0]}
                    disabled={true}
                >
                    {CISTROME_DBTOOLKIT_SPECIES.map(d => (
                        <option key={d} value={d}>
                            {d}
                        </option>
                    ))}
                </select>
            </>
        );
        const searchButton = (isReady, onClick) => (
            <div 
                className={isReady ? 'api-search-button' : 'api-search-button-disabled'}
                onClick={onClick}
            >
                <svg className="chw-button-sm chw-button-static" 
                    style={{ top: '.3em', position: 'relative' }}
                    viewBox={SEARCH.viewBox}>
                    <path d={SEARCH.path} fill="currentColor"/>
                </svg>
                Search
            </div>  
        );

        return (
            <>
                {/* Search by Interval */}
                <div 
                    className={selectedAPI === 0 ? 'api-config-view-selected' : 'api-config-view'}
                    onClick={() => { setSelectedAPI(0) }}
                    style={{ borderLeft: '4px solid #2C77B1' }}
                >
                    <div className='api-title'>Search by Interval</div>
                    <div className='api-subtitle'>What factors bind in your interval?</div>
                    {speciesSelection}
                    <div>Interval</div>
                    <input
                        className="cistrome-api-text-input"
                        type="text"
                        placeholder="ch6:151690496-152103274"
                        onChange={e => {
                            // TODO: parse input
                        }}
                    />
                    {selectedAPI === 0 ? searchButton(false, () => {}) : null}
                </div>
                {/* Search By Gene */}
                <div 
                    className={selectedAPI === 1 ? 'api-config-view-selected' : 'api-config-view'}
                    style={{ borderLeft: '4px solid #D6641E' }}
                    onClick={() => { setSelectedAPI(1) }}
                >
                    <div className='api-title'>Search by Gene</div>
                    <div className='api-subtitle'>What factors regulate your gene?</div>
                    {speciesSelection}
                    <div>Distance</div>    
                    <select
                        defaultValue={CISTROME_DBTOOLKIT_GENE_DISTANCE[0]}
                        onChange={e => setLatestGeneParams({
                            ...latestGeneParams,
                            distance: e.target.value
                        })} 
                    >
                        {CISTROME_DBTOOLKIT_GENE_DISTANCE.map(d => (
                            <option key={d} value={d}>
                                {d}
                            </option>
                        ))}
                    </select>
                    <div>Gene</div>
                    <input
                        className="cistrome-api-text-input"
                        type="text"
                        placeholder="GAPDH or NM_001289746" // TODO: confirm if GAPDH works as well
                        onChange={e => setLatestGeneParams({
                            ...latestGeneParams,
                            gene: e.target.value
                        })}
                    />
                    {selectedAPI === 1 ? searchButton(isLatestGeneParamsReady, () => {
                        requestGeneTFs(latestGeneParams)
                            .then(([rows, columns]) => {
                                console.log(rows, columns)
                                // TODO:
                            })
                    }) : null}
                </div>
                {/* Search by Peak Set */}
                <div
                    className={selectedAPI === 2 ? 'api-config-view-selected' : 'api-config-view'}
                    onClick={() => { setSelectedAPI(2) }}
                    style={{ borderLeft: '4px solid #2B9F78' }}
                >
                    <div className='api-title'>Search by Peak Set</div>
                    <div className='api-subtitle'>What factors have a significant binding overlap with your peak set?</div>
                    {speciesSelection}
                    <div>Peak Number of Cistrome Sample to Use</div>
                    <select
                        defaultValue={CISTROME_DBTOOLKIT_PEAK_NUMBERS[0]}
                        onChange={e => setPeakNumber(e.target.value)} 
                    >
                        {CISTROME_DBTOOLKIT_PEAK_NUMBERS.map(d => (
                            <option key={d} value={d}>
                                {d}
                            </option>
                        ))}
                    </select>
                    <div>Bed File</div>  
                    <input
                        className="cistrome-api-text-input"
                        type="file"
                        onChange={() => {
                            // TODO:
                        }}
                    />
                    {selectedAPI === 2 ? searchButton(false, () => {}) : null}
                </div>
            </>
        );
    }); // TODO: side effects

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
                    <div className='cisvis-api-result-header'>API Request</div>
                    <div className='cisvis-api-result-header'>Request History</div>
                    <div className='cisvis-api-result-header'>Selected Result Table</div>
                    <div className="cisvis-api-request">
                        {ApiRequestConfigViews}
                    </div>
                    <div className="cisvis-api-result-list">
                        {listOfResultsRequested}
                    </div>
                    <DataTable 
                        columns={selectedRequest ? requestHistory[selectedRequest].columns : []}
                        rows={selectedRequest ? requestHistory[selectedRequest].rows : []}
                        selectedRows={selectedRows}
                        expoNotations={["Overlap Ratio"]}
                        onCheckRows={undefined}
                        onSelect={(selectedRows) => { setSelectedRows(selectedRows) }}
                    />
                </div>
            </div>
        </div>
    );
}