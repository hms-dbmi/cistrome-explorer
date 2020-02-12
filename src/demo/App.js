import React, { useState } from 'react';

import { CistromeHGW } from '../index.js';

import hgDemoViewConfig1 from '../viewconfigs/horizontal-multivec-1.json';
import hgDemoViewConfig2 from '../viewconfigs/horizontal-multivec-2.json';
import hgDemoViewConfig4 from '../viewconfigs/horizontal-multivec-4.json';


import './App.scss';

const demos = {
    "Demo 1": {
        viewConfig: hgDemoViewConfig1,
        options: {
            rowInfoPosition: "right",
            rowLinkPosition: "left",
            colToolsPosition: "bottom",
            rowInfoAttributes: [
                {field: "Cell Type", type: "nominal"},
                {field: "Tissue Type", type: "nominal"},
                {field: "Histone Modification", type: "nominal"},
                {field: "Species", type: "nominal"}
            ],
            rowSort: [
                {field: "Species", type: "nominal", order: "ascending"},
                {field: "Histone Modification", type: "nominal", order: "ascending"},
                {field: "Tissue Type", type: "nominal", order: "ascending"},
                {field: "Cell Type", type: "nominal", order: "ascending"}
            ]
        }
    },
    "Demo 2": {
        viewConfig: hgDemoViewConfig2,
        options: {
            rowInfoPosition: "right",
            rowLinkPosition: "left",
            colToolsPosition: "bottom",
            rowLinkAttribute: "url",
            rowLinkNameAttribute: "id",
            rowInfoAttributes: [
                {field: "Hierarchical Clustering", type: "tree"},
                {field: "Cell Type", type: "nominal"},
                {field: "Tissue Type", type: "nominal"},
                {field: "Random 2", type: "quantitative"}
            ],
            rowSort: [
                {field: "Tissue Type", type: "nominal", order: "ascending"},
                {field: "Random 2", type: "quantitative", order: "descending"}
            ]
        }
    },
    "Demo 4": {
        viewConfig: hgDemoViewConfig4,
        options: [
            {
                viewId: "default",
                trackId: "default",
                colToolsPosition: "bottom",
                rowLinkAttribute: "url",
                rowLinkNameAttribute: "id"
            },
            {
                viewId: "cistrome-view-4-1",
                trackId: "cistrome-track-4-1",
                rowInfoPosition: "left",
                rowLinkPosition: "right",
                rowInfoAttributes: [
                    {field: "Tissue Type", type: "nominal"},
                    {field: "Random 2", type: "quantitative"}
                ]
            },
            {
                viewId: "cistrome-view-4-2",
                trackId: "cistrome-track-4-2",
                rowInfoPosition: "right",
                rowLinkPosition: "hidden",
                rowInfoAttributes: [
                    {field: "Hierarchical Clustering", type: "tree"},
                    {field: "Cell Type", type: "nominal"}
                ]
            }
        ]
    }
};

export default function App() {
    
    const [selectedDemo, setSelectedDemo] = useState(Object.keys(demos)[1]);

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
                />
            </div>
        </div>
    );
};
