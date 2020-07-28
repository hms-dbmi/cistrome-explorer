import React, { useRef, useState, useEffect } from 'react';
import pkg from '../../package.json';

import { HiGlassMeta } from '../index.js';
import CistromeToolkit from './CistromeToolkit.js';

import { UNDO, REDO, TABLE, DOCUMENT, GITHUB, CLOSE, MENU, TRASH } from '../utils/icons.js';
import { diffViewOptions } from '../utils/view-history';
import { demos } from './demo';
import './CistromeExplorer.scss';
import diff from 'deep-diff';

export default function CistromeExplorer() {
    
    const hmRef = useRef();

    const [selectedDemo, setSelectedDemo] = useState(Object.keys(demos)[0]);
    const [isSettingVisible, setIsSettingVisible] = useState(false);

    // Undo and redo
    const [undoable, setUndoable] = useState(false);
    const [redoable, setRedoable] = useState(false);

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

    return (
        <div className="cistrome-explorer"
            onKeyDown={e => {
                if(
                    (e.key === 'Esc' || e.key === 'Escape') && isSettingVisible
                ) {
                    setIsSettingVisible(false);
                }
            }}>
            <div className="header-container">
                <div className="header">
                    <span 
                        className="ce-generic-button"
                        onClick={() => setIsSettingVisible(!isSettingVisible)}>
                        <svg xmlns="http://www.w3.org/2000/svg"
                            viewBox={MENU.viewBox}>
                            <title>Menu</title>
                            <path fill="currentColor" d={MENU.path}/>
                        </svg>
                    </span>
                    <span className="cisvis-title">Cistrome Explorer</span>
                    <span className="header-control">
                        <span 
                            className="ce-generic-button-sm"
                            style={{ 
                                cursor: undoable ? 'pointer' : 'not-allowed',
                                color: undoable ? 'white' : '#999'
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
                            className="ce-generic-button-sm"
                            style={{ 
                                cursor: redoable ? 'pointer' : 'not-allowed',
                                color: redoable ? 'white' : '#999'
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
                    <span className="header-info">
                        <span 
                            className="ce-generic-button" 
                            onClick={() => setIsToolkitVisible(!isToolkitVisible)}>
                            <svg xmlns="http://www.w3.org/2000/svg"
                                viewBox={TABLE.viewBox}>
                                <title>Cistrome DB Toolkit</title>
                                <path fill="currentColor" d={TABLE.path}/>
                            </svg>
                            {' Toolkit '}
                        </span>
                        <span className="ce-generic-button">
                            <a href={`${pkg.homepage}/docs/`} target="_blank">
                                <svg xmlns="http://www.w3.org/2000/svg"
                                    viewBox={DOCUMENT.viewBox}>
                                    <title>Documents</title>
                                    <path fill="currentColor" d={DOCUMENT.path}/>
                                </svg>
                            </a>
                        </span>
                        <span className="ce-generic-button">
                            <a href={pkg.repository.url} target="_blank">
                                <svg xmlns="http://www.w3.org/2000/svg"
                                    viewBox={GITHUB.viewBox}>
                                    <title>GitHub</title>
                                    <path fill="currentColor" d={GITHUB.path}/>
                                </svg>
                            </a>
                        </span>
                    </span>
                </div>
            </div>

            <div className="visualization-container">
                <div className="visualization">
                    <HiGlassMeta
                        ref={hmRef}
                        viewConfig={demos[selectedDemo].viewConfig}
                        options={demos[selectedDemo].options}
                        onViewChanged={onViewChanged}
                        onGenomicIntervalSearch={setToolkitParams}
                    />
                    <CistromeToolkit
                        isVisible={isToolkitVisible}
                        intervalAPIParams={toolkitParams}
                        // TODO: After we build DB for cistrome bigwig files, uncomment the following code.
                        // onAddTrack={(server, tilesetUid, position) => { 
                        //     onAddBigWigTrack(server, tilesetUid, position);
                        // }}
                    />
                </div>
                <div className="settings" style={{
                    left: isSettingVisible ? 0 : "-300px"
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
                    <h2>View Options</h2>
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
                            hmRef.current.api.onRemoveAllFilters()}
                        }
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
                </div>
            </div>
        </div>
    );
};
