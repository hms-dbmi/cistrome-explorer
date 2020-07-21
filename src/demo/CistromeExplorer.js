import React, { useRef, useState, useEffect } from 'react';
import pkg from '../../package.json';

import { HiGlassMeta } from '../index.js';
import CistromeToolkit from './CistromeToolkit.js';

import { UNDO, REDO, TABLE, DOCUMENT, GITHUB } from '../utils/icons.js';

import hgDemoViewConfig1 from '../viewconfigs/horizontal-multivec-1.json';
import hgDemoViewConfig1b from '../viewconfigs/horizontal-multivec-1b.json';
import hgDemoViewConfig2 from '../viewconfigs/horizontal-multivec-2.json';
import hgDemoViewConfig2b from '../viewconfigs/horizontal-multivec-2b.json';
import hgDemoViewConfig6 from '../viewconfigs/horizontal-multivec-6.json';
import hgDemoViewConfig7 from '../viewconfigs/horizontal-multivec-7.json';
import hgDemoViewConfig8 from '../viewconfigs/horizontal-multivec-8.json';
import hgDemoViewConfig9 from '../viewconfigs/horizontal-multivec-9.json';
import hgDemoViewConfig10 from '../viewconfigs/horizontal-multivec-10.json';
import hgDemoViewConfigApril2020 from '../viewconfigs/meeting-2020-04-29.json';

import './CistromeExplorer.scss';
import { diffViewOptions } from '../utils/view-history';
import diff from 'deep-diff';

const demos = {
    "H3K27ac Demo (1 View, Center Track)": {
        viewConfig: hgDemoViewConfig1,
        options: {
            rowInfoAttributes: [
                {field: "Hierarchical Clustering (Average)", type: "tree", position: "left", resolveYScale: true},
                {field: "qc_frip", type: "quantitative", position: "left"},
                {field: "qc_fastqc", type: "quantitative", position: "left"},
                {field: "Metadata URL", type: "url", position: "left", title: "cid"},
                {field: "Hierarchical Clustering (Ward)", type: "tree", position: "right"},
                {field: "Cell Type", type: "nominal", position: "right"},
                {field: "Tissue Type", type: "nominal", position: "right"},
                {field: "Species", type: "nominal", position: "right"}
            ],
            rowSort: [
                {field: "Tissue Type", type: "nominal", order: "ascending"},
                {field: "qc_frip", type: "quantitative", order: "descending"}
            ],
            rowFilter: [
                {field: "Tissue Type", type: "nominal", notOneOf: ["None"]}
            ]
        }
    },
    "H3K27ac Demo (1 View, Center Track, Rows Aggregated)": {
        viewConfig: hgDemoViewConfig1b,
        options: {
            rowInfoAttributes: [
                {field: "Hierarchical Clustering (Average)", type: "tree", position: "left"},
                {field: ["qc_frip", "qc_fastqc"], type: "quantitative", position: "left", aggFunction: "mean"},
                {field: "qc_frip", type: "quantitative", position: "left", aggFunction: "mean"},
                {field: "qc_fastqc", type: "quantitative", position: "left", aggFunction: "mean"},
                {field: "id", type: "quantitative", position: "left", aggFunction: "count"},
                {field: "Metadata URL", type: "url", position: "left", title: "cid", aggFunction: "mostCommon"},
                {field: "Hierarchical Clustering (Ward)", type: "tree", position: "right"},
                {field: "Cell Type", type: "nominal", position: "right", aggFunction: "concat"},
                {field: "Tissue Type", type: "nominal", position: "right", aggFunction: "concat"},
                {field: "Species", type: "nominal", position: "right", aggFunction: "concat"}
            ],
            rowAggregate: [
                {field: "Cell Type", type: "nominal", oneOf: ["Fibroblast", "Epithelium"]},
                {field: "Tissue Type", type: "nominal", oneOf: ["Blood"]}
            ],
            rowSort: [
                {field: "Tissue Type", type: "nominal", order: "ascending"},
                {field: "qc_frip", type: "quantitative", order: "descending"}
            ],
            rowFilter: [
                {field: "Tissue Type", type: "nominal", notOneOf: [
                    "None", "Adipose", "Bone", "Bone Marrow", "Brain", "Breast", "Cervix", "Colon", "Connective Tissue", "Embryo"
                ]}
            ]
        }
    },
    "H3K27ac Demo (2 Views, Center Tracks)": {
        viewConfig: hgDemoViewConfig6,
        options: [
            {
                viewId: "default",
                trackId: "default",
                rowSort: [
                    {field: "Cell Type", type: "nominal", order: "ascending"}
                ],
                rowFilter: [],
                rowHighlight: {field: "Cell Type", type: "nominal", contains: "Stem"}
            },
            {
                viewId: "cistrome-view-6-1",
                trackId: "cistrome-track-6-1",
                rowInfoAttributes: [
                    {field: "Hierarchical Clustering (Ward)", type: "tree", position: "left"},
                ]
            },
            {
                viewId: "cistrome-view-6-2",
                trackId: "cistrome-track-6-2",
                rowInfoAttributes: [
                    {field: "Tissue Type", type: "nominal", position: "right"},
                    {field: "Cell Type", type: "nominal", position: "right"},
                ]
            }
        ]
    },
    "H3K27ac Demo (1 View, Top Track)": {
        viewConfig: hgDemoViewConfig7,
        options: {
            rowInfoAttributes: [
                {field: "Hierarchical Clustering (Average)", type: "tree", position: "right"},
                {field: "Random 3", type: "quantitative", position: "right"},
                {field: ["Random 1", "Random 2", "Random 3", "Random 4"], type: "quantitative", position: "right"},
                {field: "Metadata URL", type: "url", position: "right", title: "cid"},
                {field: "Hierarchical Clustering (Ward)", type: "tree", position: "left"},
                {field: "Cell Type", type: "nominal", position: "left"},
                {field: "Tissue Type", type: "nominal", position: "left"},
                {field: "Species", type: "nominal", position: "left"}
            ],
            rowSort: [
                {field: "Cell Type", type: "nominal", order: "ascending"}
            ],
            rowFilter: [],
            rowHighlight: {field: "Cell Type", type: "nominal", contains: "Stem"}
        }
    },
    "H3K27ac Demo (1 View, Top and Center Tracks)": {
        viewConfig: hgDemoViewConfig8,
        options: [
            {
                viewId: "default",
                trackId: "default"
            },
            {
                viewId: "cistrome-view-8",
                trackId: "cistrome-track-8-1",
                rowInfoAttributes: [
                    {field: "Hierarchical Clustering (Ward)", type: "tree", position: "left"},
                    {field: "Cell Type", type: "nominal", position: "right"},
                ]
            },
            {
                viewId: "cistrome-view-8",
                trackId: "cistrome-track-8-2",
                rowInfoAttributes: [
                    {field: "Hierarchical Clustering", type: "tree", position: "right"},
                    {field: "Cell Type", type: "nominal", position: "left"},
                ]
            }
        ]
    },
    "H3K27ac Demo (1 View, 1 Viewport, Top and Center Tracks, Overview & Detail)": {
        viewConfig: hgDemoViewConfig9,
        options: [
            {
                viewId: "default",
                trackId: "default",
            },
            {
                viewId: "cistrome-view-1",
                trackId: "cistrome-track-1",
                rowInfoAttributes: [
                    {field: "Metadata URL", type: "url", position: "left", title: "cid"},
                    {field: "Cell Type", type: "nominal", position: "right"},
                    {field: "Tissue Type", type: "nominal", position: "right"},
                    {field: "Species", type: "nominal", position: "right"}
                ],
                rowSort: [
                    { field: "Cell Type", type: "nominal", order: "ascending" }
                ]
            },
            {
                viewId: "cistrome-view-1",
                trackId: "cistrome-track-1-detail-view-1",
                rowInfoAttributes: [
                    {field: "Metadata URL", type: "url", position: "left", title: "cid"},
                    {field: "Cell Type", type: "nominal", position: "right"},
                    {field: "Tissue Type", type: "nominal", position: "right"},
                    {field: "Species", type: "nominal", position: "right"}
                ],
                rowFilter: [
                    { field: "Cell Type", type: "nominal", notOneOf: [
                        "Th1", "Spermatid", "ILC1", "Th17", "None", "Monocyte", "Natural Killer Cell", 
                        "T Lymphocyte", "Erythroid Progenitor Cell", "B Lymphocyte", "liver", "Dendritic Cell", 
                        "Macrophage", "Myeloid Cell", "Plasmablast", "Th2", "Endothelial Cell", "Epithelium", 
                        "Melanoma Cell", "Keratinocyte", "Cortex", "Stem cell", "Embryonic Stem Cell", "Adipocyte", 
                        "Neuroblastoma", "Neuroectoderm", "Neural crest cell", "Glial Cell", "Neural Progenitor Cell", 
                        "Neuroblastoma patient cells", "Inferior Temporal Lobe Cell", "Substantia Nigra Cell", 
                        "Hippocampus Middle Cell", "iPSC", "Melanocyte", "Endoderm Cell", "Erythroblast", 
                        "Lymphoblastoid", "Mesenchymal Stem Cell", "Osteoblast", "Stromal Cell", "Intermediate", 
                        "Myoblast", "Schwann Cell"
                    ] }
                ]
            }
        ]
    },
    "Cistrome_DNase_1kb_average_QN.multires.mv5": {
        viewConfig: hgDemoViewConfig10,
        options: {
            rowInfoAttributes: [
                {field: "Cluster", type: "nominal", position: "right"},
                {field: "Hierarchical Clustering", type: "tree", position: "right"},
                {field: "Cell Type", type: "nominal", position: "left"}
            ],
            rowSort: [
                // {field: "Cluster", type: "nominal", order: "ascending"},
            ],
            rowFilter: [ ]
        }
    },
    "Demo for Meeting 2020-04-29": {
        viewConfig: hgDemoViewConfigApril2020,
        options: {
            rowInfoAttributes: [
                {field: "Cluster", type: "nominal", position: "right"},
                {field: "Cell Type", type: "nominal", position: "right"}
            ],
            rowSort: [
                {field: "Cell Type", type: "nominal", order: "ascending"},
            ],
            rowFilter: [ ]
        }
    },
    "Minimal Dataset": {
        viewConfig: hgDemoViewConfig2,
        options: {
            rowInfoAttributes: [
                {field: "Hierarchical Clustering", type: "tree", position: "left"},
                {field: "Tissue Type", type: "nominal", position: "left"},
                {field: ["Random 1", "Random 2"], type: "quantitative", position: "left"},
                {field: "id", type: "nominal", position: "right"},
                {field: "Hierarchical Clustering", type: "tree", position: "right"}
            ],
            rowFilter: [ ]
        }
    },
    "Minimal Dataset (w/ Dendrogram and Aggregation)": {
        viewConfig: hgDemoViewConfig2b,
        options: {
            rowInfoAttributes: [
                {field: "Hierarchical Clustering", type: "tree", position: "left"},
                {field: "Tissue Type", type: "nominal", position: "left", aggFunction: "concat"},
                {field: "Random 1", type: "quantitative", position: "left", aggFunction: "mean"},
                {field: ["Random 1", "Random 2"], type: "quantitative", position: "left", aggFunction: "mean"},
                {field: "id", type: "nominal", position: "right", aggFunction: "concat"},
                {field: "Hierarchical Clustering", type: "tree", position: "right"}
            ],
            rowFilter: [ ],
            rowAggregate: [
                {field: "Tissue Type", type: "nominal", oneOf: ["Blood", "Bone Marrow"]}
            ]
        }
    }
};

export default function CistromeExplorer() {
    
    const hmRef = useRef();

    const [selectedDemo, setSelectedDemo] = useState(Object.keys(demos)[0]);
    
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
        <div className="cistrome-explorer">
            <div className="header-container">
                <div className="header">
                    <span className="cisvis-title">Cistrome Explorer</span>
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
                    <span className="header-control">
                        <span 
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
                            <svg xmlns="http://www.w3.org/2000/svg" width="12" height="12"
                                viewBox={UNDO.viewBox}>
                                <title>Undo</title>
                                <path fill="currentColor" d={UNDO.path}/>
                            </svg>
                            {` Undo (${viewHistory.length - indexOfCurrentView - 1})`}
                        </span>
                        <span 
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
                            <svg xmlns="http://www.w3.org/2000/svg" width="12" height="12"
                                viewBox={REDO.viewBox}>
                                <title>Redo</title>
                                <path fill="currentColor" d={REDO.path}/>
                            </svg>
                            {` Redo (${indexOfCurrentView})`}
                        </span>
                    </span>
                    <span className="header-info">
                        <span style={{ cursor: 'pointer' }} onClick={() => 
                            setIsToolkitVisible(!isToolkitVisible)
                        }>
                            <svg xmlns="http://www.w3.org/2000/svg" width="18" height="18"
                                viewBox={TABLE.viewBox}>
                                <title>CistromeToolkit</title>
                                <path fill="currentColor" d={TABLE.path}/>
                            </svg>
                        </span>
                        <span>
                            <a href={`${pkg.homepage}/docs/`} target="_blank">
                                <svg xmlns="http://www.w3.org/2000/svg" width="18" height="18"
                                    viewBox={DOCUMENT.viewBox}>
                                    <title>Documents</title>
                                    <path fill="currentColor" d={DOCUMENT.path}/>
                                </svg>
                            </a>
                        </span>
                        <span>
                            <a href={pkg.repository.url} target="_blank">
                                <svg xmlns="http://www.w3.org/2000/svg" width="18" height="18"
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
            </div>
        </div>
    );
};
