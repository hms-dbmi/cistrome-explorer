import React, { useState } from 'react';

import { CistromeHGW } from '../index.js';

import hgDemoViewConfig1 from '../viewconfigs/horizontal-multivec-1.json';
import hgDemoViewConfig2 from '../viewconfigs/horizontal-multivec-2.json';
import hgDemoViewConfig3 from '../viewconfigs/horizontal-multivec-4.json';


import './App.scss';

const demos = {
    "Demo 1": {
        viewConfig: hgDemoViewConfig1,
        options: {
            colToolsPosition: "bottom",
            rowInfoAttributes: [
                {field: "Hierarchical Clustering (Average)", type: "tree", position: "left"},
                {field: "Random 1", type: "quantitative", position: "left"},
                {field: "Random 3", type: "quantitative", position: "left"},
                {field: ["Random 1", "Random 2", "Random 3", "Random 4"], type: "quantitative", position: "left"},
                {field: "Hierarchical Clustering (Ward)", type: "tree", position: "right"},
                {field: "Cell Type", type: "nominal", position: "right"},
                {field: "Tissue Type", type: "nominal", position: "right"},
                {field: "Species", type: "nominal", position: "right"}
            ],
            rowSort: [
                {field: "Cell Type", type: "nominal", order: "ascending"}
            ],
            rowFilter: [
                {field: "Cell Type", type: "nominal", contains: "cell"}
            ],
            rowHighlight: {field: "Cell Type", type: "nominal", contains: "Stem"}
        }
    },
    "Demo 2": {
        viewConfig: hgDemoViewConfig2,
        options: {
            colToolsPosition: "bottom",
            rowInfoAttributes: [
                {field: "Random 2", type: "quantitative", position: "left"},
                {field: "url", type: "url", position: "left", title: "Cell Type"},
                {field: "Hierarchical Clustering", type: "tree", position: "right"},
                {field: "Cell Type", type: "nominal", position: "right"},
                {field: "Tissue Type", type: "nominal", position: "right"}
            ],
            rowSort: [
                // {field: "Tissue Type", type: "nominal", order: "ascending"},
                // {field: "Random 2", type: "quantitative", order: "descending"}
            ]
        }
    },
    "Demo 3": {
        viewConfig: hgDemoViewConfig3,
        options: [
            {
                viewId: "default",
                trackId: "default",
                colToolsPosition: "bottom",
                rowSort: [
                    {field: "Random 2", type: "quantitative", order: "ascending"}
                ]
            },
            {
                viewId: "cistrome-view-4-1",
                trackId: "cistrome-track-4-1",
                rowInfoAttributes: [
                    {field: "Tissue Type", type: "nominal", position: "left"},
                    {field: "Random 2", type: "quantitative", position: "left"},
                    {field: "url", type: "url", position: "right", title: "id"}
                ]
            },
            {
                viewId: "cistrome-view-4-2",
                trackId: "cistrome-track-4-2",
                rowInfoAttributes: [
                    {field: "Hierarchical Clustering", type: "tree", position: "right"},
                    {field: "Cell Type", type: "nominal", position: "right"}
                ]
            }
        ]
    }
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
                <form>
                    {Object.keys(demos).map(vcKey => (
                        <span key={vcKey}>
                            <input 
                                type="radio" 
                                name="viewconf"
                                value={vcKey} 
                                id={vcKey} 
                                checked={selectedDemo === vcKey} 
                                onChange={e => setSelectedDemo(e.target.value)}
                            />
                            <label htmlFor={vcKey}>{vcKey}</label>
                        </span>
                    ))}
                </form>
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
