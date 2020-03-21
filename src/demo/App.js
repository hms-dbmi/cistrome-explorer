import React, { useState } from 'react';

import { CistromeHGW } from '../index.js';

import hgDemoViewConfig1 from '../viewconfigs/horizontal-multivec-1.json';
import hgDemoViewConfig6 from '../viewconfigs/horizontal-multivec-6.json';
import hgDemoViewConfig7 from '../viewconfigs/horizontal-multivec-7.json';
import hgDemoViewConfig8 from '../viewconfigs/horizontal-multivec-8.json';
import hgDemoViewConfig9 from '../viewconfigs/horizontal-multivec-9.json';

import './App.scss';

const demos = {
    "H3K27ac Demo (1 View, Center Track)": {
        viewConfig: hgDemoViewConfig1,
        options: {
            rowInfoAttributes: [
                {field: "Hierarchical Clustering (Average)", type: "tree", position: "left"},
                {field: "qc_frip", type: "quantitative", position: "left"},
                {field: "qc_fastqc", type: "quantitative", position: "left"},
                {field: "Metadata URL", type: "url", position: "left", title: "cid"},
                {field: "Hierarchical Clustering (Ward)", type: "tree", position: "right"},
                {field: "Cell Type", type: "nominal", position: "right"},
                {field: "Tissue Type", type: "nominal", position: "right"},
                {field: "Species", type: "nominal", position: "right"}
            ],
            rowSort: [
                {field: "Tissue Type", type: "nominal", order: "ascending"}
            ],
            // TODO: 
            rowHighlight: {field: "qc_frip", type: "quantitative", range: [25, 35]},
            rowFilter: [
                // TODO: 
                {field: "qc_frip", type: "quantitative", range: [20, 40]}
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
                {field: "Hierarchical Clustering (Average)", type: "tree", position: "left"},
                {field: "Random 3", type: "quantitative", position: "left"},
                {field: ["Random 1", "Random 2", "Random 3", "Random 4"], type: "quantitative", position: "left"},
                {field: "Metadata URL", type: "url", position: "left", title: "cid"},
                {field: "Hierarchical Clustering (Ward)", type: "tree", position: "right"},
                {field: "Cell Type", type: "nominal", position: "right"},
                {field: "Tissue Type", type: "nominal", position: "right"},
                {field: "Species", type: "nominal", position: "right"}
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
                    { field: "Cell Type", type: "nominal", contains: "Fibroblast" }
                ]
            }
        ]
    },
};

function onViewConfigChange(viewConfigString) {
    console.log("View config changed");
}

export default function App() {
    
    const [selectedDemo, setSelectedDemo] = useState(Object.keys(demos)[0]);

    return (
        <div className="app">
            <div className="header">
                <h4>Cistrome HiGlass Wrapper</h4>
            </div>
            <div className="viewconf-options">
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
            </div>

            <div className="container">
                <CistromeHGW 
                    viewConfig={demos[selectedDemo].viewConfig}
                    options={demos[selectedDemo].options}
                    onViewConfigChange={onViewConfigChange}
                />
            </div>
        </div>
    );
};
