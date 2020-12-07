import React, { useRef, useState, useEffect, useCallback } from 'react';
import pkg from '../../package.json';

import { HiGlassMeta } from '../index.js';
import CistromeToolkit from './CistromeToolkit.js';

import { UNDO, REDO, TABLE_2, DOCUMENT, GITHUB, TOGGLE_ON, TOGGLE_OFF, CLOSE, ELLIPSIS, TRASH, SEARCH, FOLDER, PENCIL } from '../utils/icons.js';
import { DEFAULT_COLOR_RANGE } from '../utils/color.js';
import { diffViewOptions } from '../utils/view-history';
import { demos } from './demo';
import { CISTROME_DBTOOLKIT_GENE_DISTANCE, CISTROME_DBTOOLKIT_SPECIES } from '../utils/cistrome';

import { publishHelpTooltip, destroyTooltip } from "../Tooltip.js";

import './CistromeExplorer.scss';

import StackedBarTrack from 'higlass-multivec/es/StackedBarTrack';
import ScaleLegendTrack from '../scale-legend/ScaleLegendTrack';
import higlassRegister from 'higlass-register';

higlassRegister({
    name: 'StackedBarTrack',
    track: StackedBarTrack,
    config: StackedBarTrack.config,
});

higlassRegister({
    name: 'ScaleLegendTrack',
    track: ScaleLegendTrack,
    config: ScaleLegendTrack.config,
});

export default function CistromeExplorer() {
    
    const hmRef = useRef();

    const [selectedDemo, setSelectedDemo] = useState(Object.keys(demos)[0]);
    const [isSettingVisible, setIsSettingVisible] = useState(false);

    // metadata
    const jsonInputFile = useRef(null);
    const [fileReader, setFileReader] = useState(new FileReader());
    const [localMetadata, setLocalMetadata] = useState(null);

    useEffect(() => {
        fileReader.onload = (event) => {
            const json = JSON.parse(event.target.result);
            setLocalMetadata(json);
        };
    }, [fileReader]);

    // search
    const searchBoxRef = useRef();
    const [searchKeyword, setSearchKeyword] = useState('');
    const [geneSuggestions, setGeneSuggestions] = useState([]);
    const [suggestionPosition, setSuggestionPosition] = useState({left: 0, top: 0});

    // Undo and redo
    const [undoable, setUndoable] = useState(false);
    const [redoable, setRedoable] = useState(false);

    // toggle
    const [helpActivated, setHelpActivated] = useState(false);
    const [aggActivated, setAggActivated] = useState(false);

    // History of view updates
    const MAX_HISTORY_LENGTH = 50;  // How many previous views should be recorded?
    const [viewHistory, setViewHistory] = useState([{
        // TODO: This variable should contain additional information to support the full functionality of undo/redu,
        //       such as viewConfig, use of toolkits, interval selection, gene search
        //       (e.g., add `viewConfig: Object.values(demos)[0].viewConfig,`).
        options: Object.values(demos)[0].options
    }]);
    const [indexOfCurrentView, setIndexOfCurrentView] = useState(0); // The most recent view will be stored at the index zero.

    // Toolkit-related
    const [isToolkitVisible, setIsToolkitVisible] = useState(false);
    const [toolkitParams, setToolkitParams] = useState(undefined);
    const [geneSearched, setGeneSearched] = useState(undefined);
    const [geneToolkitParams, setGeneToolkitParams] = useState(
        undefined
        // DEBUG:
        // {
        //     assembly: CISTROME_DBTOOLKIT_SPECIES[0], 
        //     gene: 'MYC', 
        //     distance: CISTROME_DBTOOLKIT_GENE_DISTANCE[0]
        // }
    );

    const addNewTrack = useCallback((trackDef, viewId, position) => {
        hmRef.current.api.addNewTrack(trackDef, viewId, position);
    }, [hmRef]);

    // Callback function for adding a BigWig track.
    const onAddTrack = useCallback((cistromeDataConfig) => {
        const { species, factor } = cistromeDataConfig;
        const firstViewUid = demos[selectedDemo].viewConfig.views[0].uid;
        addNewTrack({
            type: 'horizontal-multivec',
            data: {
                server: "http://ec2-3-93-68-250.compute-1.amazonaws.com/api/v1",
                url: "s3://CistromeDB/" + `${species}__${factor}__all`.replace(" ", "_") + ".multires.mv5",
                filetype: "multivec",
            },
            coordSystem: "hg38",
            options: {
                colorRange: DEFAULT_COLOR_RANGE,
            },
            height: 200,
        }, firstViewUid, 'top');
    }, [addNewTrack, selectedDemo]);

    // When a user select a different demo, initialize the view history.
    useEffect(() => {
        setViewHistory([{
            options: demos[selectedDemo].options
        }]);
        setIndexOfCurrentView(0);
    }, [selectedDemo]);

    useEffect(() => {
        setUndoable(indexOfCurrentView !== viewHistory.length - 1);
        setRedoable(indexOfCurrentView !== 0);
    }, [viewHistory, indexOfCurrentView]);

    useEffect(() => {
        function closeSideViews(e) {
            if(
                (e.key === 'Esc' || e.key === 'Escape') && 
                (isSettingVisible || isToolkitVisible)
            ) {
                setIsSettingVisible(false);
                setIsToolkitVisible(false);
            }
        }
        window.addEventListener("keydown", closeSideViews);
        return () => window.removeEventListener("keydown", closeSideViews);
    }, [isSettingVisible, isToolkitVisible]);

    /**
     * This function is being called when `options` is updated interactively.
     * @param {object} viewOptions A JSON object that contains updated visualization specs for `HiGlassMeta`.
     * @param {object} viewoptions.options A JSON object that contains options for the metadata visualizations in `HiGlassMeta`.
     */
    function onViewChanged(viewOptions) {
        // Make sure not to update the history if there is no difference.
        if(!diffViewOptions(viewOptions.options, viewHistory[indexOfCurrentView].options)) {
            return;
        }
        // DEBUG: To see the difference between two JSON objects
        // console.log("View updated", diff.diff(viewOptions, viewHistory[indexOfCurrentView]));

        // Update the view history
        const newViewHistory = viewHistory.slice();
        if(indexOfCurrentView !== 0) {
            // This means a user ever have clicked on the `Undo` button, 
            // and we want to overwrite recent history.
            newViewHistory.splice(0, indexOfCurrentView);
        }

        // Add a recent view at the start of the array.
        newViewHistory.unshift({
            options: viewOptions.options
        });

        // Remove the tail to make the length of the array be less than or equal to the threshold.
        if(newViewHistory.length > MAX_HISTORY_LENGTH) { 
            newViewHistory.splice(MAX_HISTORY_LENGTH - 1);
        }
        setViewHistory(newViewHistory);
        setIndexOfCurrentView(0);
    }

    // const aggregatedOptions = !aggActivated ? demos[selectedDemo].options : (
    //     Array.isArray(demos[selectedDemo].options)
    //     ? modifyItemInArray(demos[selectedDemo].options, 0, {
    //         ...demos[selectedDemo].options[0],
    //         rowAggregate: [ {field: "Tissue Type", type: "nominal", notOneOf: []} ]
    //     })
    //     : {
    //         ...demos[selectedDemo].options,
    //         rowAggregate: [ {field: "Tissue Type", type: "nominal", notOneOf: []} ],
    //     }
    // )

    return (
        <div className="cistrome-explorer">
            <div className="header-container">
                <div className="header">
                    <span className="cisvis-title">
                        <hl>Cistrome</hl> Explorer
                    </span>
                    <span 
                        className="header-control"
                        onMouseMove={(e) => publishHelpTooltip(e,
                            "Search Gene or Genomic Interval",
                            "You can quickly reposition the heatmap on the right by entering a gene name or genomic interval, such as GAPDH or chr6:151690496-152103274.",
                            helpActivated
                        )}
                        onMouseLeave={() => destroyTooltip()}
                    >
                        <span 
                            className="ce-generic-button"
                            style={{ cursor: 'auto' }}
                        >
                            <svg xmlns="http://www.w3.org/2000/svg"
                                viewBox={SEARCH.viewBox}
                            >
                                <path d={SEARCH.path} fill="currentColor"/>
                            </svg>
                            <input
                                ref={searchBoxRef}
                                className={"position-search-box " + (helpActivated ? 'help-highlight' : '')}
                                type="text"
                                name="default name"
                                placeholder="GAPDH or chr6:151690496-152103274"
                                onChange={(e) => {
                                    const keyword = e.target.value;
                                    if(keyword !== '' && !keyword.startsWith('c')) {
                                        hmRef.current.api.suggestGene(keyword, (suggestions) => {
                                            setGeneSuggestions(suggestions);
                                        });
                                        setSuggestionPosition({
                                            left: searchBoxRef.current.getBoundingClientRect().left,
                                            top: searchBoxRef.current.getBoundingClientRect().top + searchBoxRef.current.getBoundingClientRect().height,
                                        });
                                    } else {
                                        setGeneSuggestions([]);
                                    }
                                    setSearchKeyword(keyword);
                                }}
                                onKeyDown={(e) => {
                                    switch(e.key){
                                        case 'ArrowUp':
                                            break;
                                        case 'ArrowDown':
                                            break;
                                        case 'Enter':
                                            setGeneSuggestions([]);
                                            if(searchKeyword.includes('chr')) {
                                                hmRef.current.api.zoomTo(searchKeyword);
                                            } else {
                                                hmRef.current.api.zoomToGene(searchKeyword);
                                            }
                                            break;
                                        case 'Esc':
                                        case 'Escape':
                                            break;
                                    }
                                }}
                            />
                        </span>
                    </span>
                    <span className="header-control">
                        <span
                            className={"ce-generic-button " + (undoable ? '' : 'ce-generic-button-deactivated')}
                            style={{ 
                                cursor: undoable ? 'pointer' : 'not-allowed'
                            }} 
                            onClick={() => {
                                if(undoable) {
                                    const newViewIndex = indexOfCurrentView + 1;
                                    setIndexOfCurrentView(newViewIndex);
                                    hmRef.current.api.onOptions(viewHistory[newViewIndex].options);
                                }
                            }}
                        >
                            <svg xmlns="http://www.w3.org/2000/svg"
                                viewBox={UNDO.viewBox}>
                                <title>Undo</title>
                                <path fill="currentColor" d={UNDO.path}/>
                            </svg>
                            {` Undo (${viewHistory.length - indexOfCurrentView - 1})`}
                        </span>
                        <span 
                            className={"ce-generic-button " + (redoable ? '' : 'ce-generic-button-deactivated')}
                            style={{ 
                                cursor: redoable ? 'pointer' : 'not-allowed'
                            }} 
                            onClick={() => {
                                if(redoable) {
                                    const newViewIndex = indexOfCurrentView - 1;
                                    setIndexOfCurrentView(newViewIndex);
                                    hmRef.current.api.onOptions(viewHistory[newViewIndex].options);
                                }
                            }}
                        >
                            <svg xmlns="http://www.w3.org/2000/svg"
                                viewBox={REDO.viewBox}>
                                <title>Redo</title>
                                <path fill="currentColor" d={REDO.path}/>
                            </svg>
                            {` Redo (${indexOfCurrentView})`}
                        </span>
                    </span>
                    <span className="header-control"
                        onMouseMove={(e) => publishHelpTooltip(e,
                            "Cistrome Data Browser Toolkit",
                            "You can query for transcription factors that are likely to bind in the region of your intrest based on the thousand of samples available in Cistrome Data Browser",
                            helpActivated
                        )}
                        onMouseLeave={() => destroyTooltip()}
                    >
                        <span 
                            className={"ce-generic-button " + (helpActivated ? 'help-highlight' : '')}
                            onClick={() => setIsToolkitVisible(!isToolkitVisible)}>
                            <svg xmlns="http://www.w3.org/2000/svg"
                                viewBox={TABLE_2.viewBox}>
                                <title>Cistrome DB Toolkit</title>
                                <path fill="currentColor" d={TABLE_2.path}/>
                            </svg>
                            {' Cistrome Search '}
                        </span>
                    </span>
                    <span className="header-control"
                        onMouseMove={(e) => publishHelpTooltip(e,
                            "Aggregate Samples By Tissue Type",
                            "You can aggregate samples with the same tissue type into a single row in the visualization",
                            helpActivated
                        )}
                        onMouseLeave={() => destroyTooltip()}
                    >
                        <span 
                            className={"ce-generic-button-lg " + (aggActivated ? 'ce-generic-button-activated ' : '') + (helpActivated ? 'help-highlight' : '')}
                            onClick={() => { setAggActivated(!aggActivated); }}
                        >
                            <svg xmlns="http://www.w3.org/2000/svg"
                                viewBox={aggActivated ? TOGGLE_ON.viewBox : TOGGLE_OFF.viewBox}>
                                <title>Aggregate Rows By Tissue Type</title>
                                <path fill="currentColor" d={aggActivated ? TOGGLE_ON.path : TOGGLE_OFF.path}/>
                            </svg>
                            {` Aggregate By Tissue`}
                        </span>
                    </span>
                    <span className="header-control">
                        <span 
                            className={"ce-generic-button-lg " + (helpActivated ? 'ce-generic-button-activated' : '')}
                            onClick={() => { setHelpActivated(!helpActivated); }}
                        >
                            <svg xmlns="http://www.w3.org/2000/svg"
                                viewBox={helpActivated ? TOGGLE_ON.viewBox : TOGGLE_OFF.viewBox}>
                                <title>Help</title>
                                <path fill="currentColor" d={helpActivated ? TOGGLE_ON.path : TOGGLE_OFF.path}/>
                            </svg>
                            {` Show Instructions`}
                        </span>
                    </span>
                    <span className="header-info">
                        {geneSearched ? 
                            <span 
                                className="ce-generic-button"
                                onClick={() => {
                                    if(geneSearched) {
                                        setGeneToolkitParams({
                                            assembly: CISTROME_DBTOOLKIT_SPECIES[0], 
                                            gene: geneSearched, 
                                            distance: CISTROME_DBTOOLKIT_GENE_DISTANCE[0]
                                        });
                                    }
                                }}>
                                {`💡 Search ${geneSearched} in Cistrome Toolkit? `}
                            </span>
                            : null}
                        <span 
                            className="ce-generic-button"
                            onClick={() => setIsSettingVisible(!isSettingVisible)}>
                            <svg xmlns="http://www.w3.org/2000/svg"
                                viewBox={ELLIPSIS.viewBox}>
                                <title>Menu</title>
                                <path fill="currentColor" d={ELLIPSIS.path}/>
                            </svg>
                        </span>
                    </span>
                    {geneSuggestions.length !== 0 ? 
                        <div className="gene-suggestion" style={{
                            left: suggestionPosition.left,
                            top: suggestionPosition.top                        
                        }}>
                            <ul>
                                {geneSuggestions.map((d, i) => (
                                    <li style={{textAlign: 'right', color: 'gray'}}
                                        key={d.geneName + d.score}
                                        onClick={() => {
                                            searchBoxRef.current.value = d.geneName;
                                            setGeneSuggestions([]);
                                            hmRef.current.api.zoomToGene(d.geneName);
                                        }}
                                    >
                                        <strong style={{float: 'left', color: 'black'}}>{d.geneName}</strong>
                                        {`${d.chr}:${d.txStart}-${d.txEnd}`}
                                    </li>
                                ))}
                            </ul>
                        </div>
                        : null
                    }
                </div>
            </div>

            {/* <div className="sub-header"> </div> */}

            <div className="visualization-container">
                <div className="visualization">
                    <HiGlassMeta
                        ref={hmRef}
                        viewConfig={demos[selectedDemo].viewConfig}
                        options={demos[selectedDemo].options}
                        rowInfo={localMetadata}
                        aggregateRowBy={aggActivated ? "Tissue Type" : undefined}
                        helpActivated={helpActivated}
                        onViewChanged={onViewChanged}
                        onGenomicIntervalSearch={setToolkitParams}
                        onGeneSearch={setGeneSearched}
                    />
                    <CistromeToolkit
                        isVisible={isToolkitVisible}
                        intervalAPIParams={toolkitParams}
                        geneAPIParams={geneToolkitParams}
                        onAddTrack={onAddTrack}
                    />
                </div>
                <div className="settings" style={{
                    right: isSettingVisible ? 0 : '-400px'
                }}>
                    <span style={{ 
                        verticalAlign: "middle", 
                        display: "inline-block",
                        position: "absolute", 
                        right: 5, 
                        top: 5
                    }}>
                        <svg
                            className={'hm-button'}
                            style={{ color: "rgb(171, 171, 171)", background: "none" }}
                            onClick={() => setIsSettingVisible(false)}
                            viewBox={CLOSE.viewBox}
                        >
                            <title>Close Cistrome Toolkit</title>
                            <path d={CLOSE.path} fill="currentColor"/>
                        </svg>
                    </span>
                    <h2>Example Datasets</h2>
                    <span className="viewconf-options">
                        <select 
                            onChange={e => setSelectedDemo(e.target.value)} 
                            defaultValue={selectedDemo}
                        >
                            {Object.keys(demos).map(vcKey => (
                                <option 
                                    key={vcKey} 
                                    value={vcKey} 
                                >
                                    {vcKey}
                                </option>
                            ))}
                        </select>
                    </span>
                    <div className="setting-separater"></div>
                    <h2>Metadata</h2>
                    <span 
                        className="ce-generic-button"
                        style={{ 
                            fontSize: 12,
                            display: 'inline-block',
                            cursor: 'not-allowed',  // TODO: not supported yet
                            color: 'gray', // TODO: not supported yet
                            background: 'white', 
                            border: '1px solid gray',
                            padding: '4px',
                            marginTop: '4px'
                        }}
                        onClick={() => { }}
                    >
                        <svg
                            style={{ color: "rgb(171, 171, 171)", width: 14, height: 14 }}
                            viewBox={SEARCH.viewBox}
                        >
                            <title>Open Metadata</title>
                            <path d={SEARCH.path} fill="currentColor"/>
                        </svg>
                        {' View Loaded Metadata'}
                    </span>
                    <span 
                        className="ce-generic-button"
                        style={{ 
                            fontSize: 12,
                            cursor: 'pointer',
                            display: 'inline-block',
                            color: 'black',
                            background: 'white', 
                            border: '1px solid gray',
                            padding: '4px',
                            marginTop: '4px'
                        }}
                        onClick={() => { 
                            jsonInputFile.current.value = null;
                            setLocalMetadata(undefined);
                        }}
                    >
                        <svg
                            style={{ color: "rgb(171, 171, 171)", width: 14, height: 14 }}
                            viewBox={TRASH.viewBox}
                        >
                            <title>Open Metadata</title>
                            <path d={TRASH.path} fill="currentColor"/>
                        </svg>
                        {' Remove Loaded Metadata'}
                    </span>
                    <input 
                        type="file" 
                        ref={jsonInputFile} 
                        style={{ display: 'none' }} 
                        onChange={(e) => {
                            fileReader.readAsText(e.target.files[0]);
                        }
                    }/>
                    <span 
                        className="ce-generic-button"
                        style={{ 
                            fontSize: 12,
                            cursor: 'pointer',
                            display: 'inline-block',
                            color: 'black', 
                            background: 'white', 
                            border: '1px solid gray',
                            padding: '4px',
                            marginTop: '4px'
                        }}
                        onClick={() => { jsonInputFile.current.click() }}
                    >
                        <svg
                            style={{ color: "rgb(171, 171, 171)", width: 14, height: 14 }}
                            viewBox={FOLDER.viewBox}
                        >
                            <title>Open Metadata</title>
                            <path d={FOLDER.path} fill="currentColor"/>
                        </svg>
                        {' Open Local Metadata (JSON)'}
                    </span>
                    <div className="setting-separater"></div>
                    <h2>View Options</h2>
                    <span 
                        className="ce-generic-button"
                        style={{ 
                            fontSize: 12,
                            display: 'inline-block',
                            cursor: 'not-allowed',  // TODO: not supported yet
                            color: 'gray',  // TODO: not supported yet
                            background: 'white', 
                            border: '1px solid gray',
                            padding: '4px',
                            marginTop: '4px'
                        }}
                        onClick={() => { }}
                    >
                        <svg
                            style={{ color: "rgb(171, 171, 171)", width: 14, height: 14 }}
                            viewBox={PENCIL.viewBox}
                        >
                            <title>Edit View Options</title>
                            <path d={PENCIL.path} fill="currentColor"/>
                        </svg>
                        {' Edit View Options'}
                    </span>
                    <span 
                        className="ce-generic-button"
                        style={{ 
                            fontSize: 12,
                            cursor: 'pointer',
                            display: 'inline-block',
                            color: 'black', 
                            background: 'white', 
                            border: '1px solid gray',
                            padding: '4px',
                            marginTop: '4px'
                        }}
                        onClick={() => {
                            hmRef.current.api.onRemoveAllFilters()
                        }}
                    >
                        <svg
                            style={{ color: "rgb(171, 171, 171)", width: 14, height: 14 }}
                            viewBox={TRASH.viewBox}
                        >
                            <title>Remove All Filters</title>
                            <path d={TRASH.path} fill="currentColor"/>
                        </svg>
                        {' Remove All Filters'}
                    </span>
                    <span 
                        className="ce-generic-button"
                        style={{ 
                            fontSize: 12,
                            cursor: 'pointer',
                            display: 'inline-block',
                            color: 'black', 
                            background: 'white', 
                            border: '1px solid gray',
                            padding: '4px',
                            marginTop: '4px'
                        }}
                        onClick={() => {
                            hmRef.current.api.onRemoveAllSort()}
                        }
                    >
                        <svg
                            style={{ color: "rgb(171, 171, 171)", width: 14, height: 14 }}
                            viewBox={TRASH.viewBox}
                        >
                            <title>Remove All Sort</title>
                            <path d={TRASH.path} fill="currentColor"/>
                        </svg>
                        {' Remove All Sort'}
                    </span>
                    <div className="setting-separater"></div>
                    <h2>Resource</h2>
                    <span 
                        className="ce-generic-button"
                        style={{ 
                            fontSize: 12,
                            cursor: 'pointer',
                            display: 'inline-block',
                            color: 'black', 
                            background: 'white', 
                            border: '1px solid gray',
                            padding: '4px',
                            marginTop: '4px'
                        }}
                        onClick={() => {
                            window.open(`${pkg.homepage}/docs/`)
                        }}
                    >
                        <svg xmlns="http://www.w3.org/2000/svg"
                            viewBox={DOCUMENT.viewBox}>
                            <title>Documents</title>
                            <path fill="currentColor" d={DOCUMENT.path}/>
                        </svg>
                        {' Documentation'}
                    </span>
                    <span 
                        className="ce-generic-button"
                        style={{ 
                            fontSize: 12,
                            cursor: 'pointer',
                            display: 'inline-block',
                            color: 'black', 
                            background: 'white', 
                            border: '1px solid gray',
                            padding: '4px',
                            marginTop: '4px'
                        }}
                        onClick={() => {
                            window.open(pkg.repository.url)
                        }}
                    >
                        <svg xmlns="http://www.w3.org/2000/svg"
                            viewBox={GITHUB.viewBox}>
                            <title>GitHub</title>
                            <path fill="currentColor" d={GITHUB.path}/>
                        </svg>
                        {' Open Source'}
                    </span>
                </div>
            </div>
        </div>
    );
};
