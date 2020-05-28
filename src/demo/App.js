import React, { useState } from 'react';
import pkg from '../../package.json';

import { CistromeHGW } from '../index.js';

import hgDemoViewConfig1 from '../viewconfigs/horizontal-multivec-1.json';
import hgDemoViewConfig2 from '../viewconfigs/horizontal-multivec-2.json';
import hgDemoViewConfig2b from '../viewconfigs/horizontal-multivec-2b.json';
import hgDemoViewConfig6 from '../viewconfigs/horizontal-multivec-6.json';
import hgDemoViewConfig7 from '../viewconfigs/horizontal-multivec-7.json';
import hgDemoViewConfig8 from '../viewconfigs/horizontal-multivec-8.json';
import hgDemoViewConfig9 from '../viewconfigs/horizontal-multivec-9.json';
import hgDemoViewConfig10 from '../viewconfigs/horizontal-multivec-10.json';
import hgDemoViewConfigApril2020 from '../viewconfigs/meeting-2020-04-29.json';

import './App.scss';

const demos = {
    "H3K27ac Demo (1 View, Center Track)": {
        viewConfig: hgDemoViewConfig1,
        options: {
            rowInfoAttributes: [
                // {field: "Hierarchical Clustering (Average)", type: "tree", position: "left"},
                {field: ["qc_frip", "qc_fastqc"], type: "quantitative", position: "left", aggFunction: "mean"},
                {field: "qc_frip", type: "quantitative", position: "left", aggFunction: "mean"},
                {field: "qc_fastqc", type: "quantitative", position: "left", aggFunction: "mean"},
                // {field: "Metadata URL", type: "url", position: "left", title: "cid"},
                // {field: "Hierarchical Clustering (Ward)", type: "tree", position: "right"},
                {field: "Cell Type", type: "nominal", position: "right", aggFunction: "sum"},
                {field: "Tissue Type", type: "nominal", position: "right", aggFunction: "sum"},
                {field: "Species", type: "nominal", position: "right", aggFunction: "sum"}
            ],
            rowAggregate: [
                {field: "Cell Type", type: "nominal", oneOf: ["Fibroblast", "Epithelium"]},
                {field: "Tissue Type", type: "nominal", oneOf: ["Blood"]}
            ],
            // rowSort: [
            //     {field: "Tissue Type", type: "nominal", order: "ascending"},
            //     {field: "qc_frip", type: "quantitative", order: "descending"}
            // ],
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
                // {field: "Hierarchical Clustering", type: "tree", position: "left"},
                {field: "Tissue Type", type: "nominal", position: "left", aggFunction: "sum"},
                {field: "Random 1", type: "quantitative", position: "left", aggFunction: "mean"},
                {field: ["Random 1", "Random 2"], type: "quantitative", position: "left", aggFunction: "mean"},
                {field: "id", type: "nominal", position: "right", aggFunction: "sum"},
                // {field: "Hierarchical Clustering", type: "tree", position: "right"}
            ],
            rowFilter: [ ],
            rowAggregate: [
                {field: "Tissue Type", type: "nominal", oneOf: ["Blood", "Bone Marrow"]}
            ]
        }
    }
};

function onViewConfigChange(viewConfigString) {
    console.log("View config changed");
}

export default function App() {
    
    const [selectedDemo, setSelectedDemo] = useState(Object.keys(demos)[0]);

    return (
        <div className="app">
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
                    <span className="header-info">
                        <span>
                            <a href={`${pkg.homepage}/docs/`} target="_blank">
                                <svg xmlns="http://www.w3.org/2000/svg" width="22" height="22"
                                    viewBox="0 0 1792 1792">
                                    <title>Documents</title>
                                    <path fill="#666" d="M1528 1280h-248v248q29-10 41-22l185-185q12-12 22-41zm-280-128h288v-896h-1280v1280h896v-288q0-40 28-68t68-28zm416-928v1024q0 40-20 88t-48 76l-184 184q-28 28-76 48t-88 20h-1024q-40 0-68-28t-28-68v-1344q0-40 28-68t68-28h1344q40 0 68 28t28 68z"/>
                                </svg>
                            </a>
                        </span>
                        <span>
                            <a href={pkg.repository.url} target="_blank">
                                <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24"
                                    viewBox="0 0 24 24">
                                    <title>GitHub</title>
                                    <path fill="#666" d="M12 .297c-6.63 0-12 5.373-12 12 0 5.303 3.438 9.8 8.205 11.385.6.113.82-.258.82-.577 0-.285-.01-1.04-.015-2.04-3.338.724-4.042-1.61-4.042-1.61C4.422 18.07 3.633 17.7 3.633 17.7c-1.087-.744.084-.729.084-.729 1.205.084 1.838 1.236 1.838 1.236 1.07 1.835 2.809 1.305 3.495.998.108-.776.417-1.305.76-1.605-2.665-.3-5.466-1.332-5.466-5.93 0-1.31.465-2.38 1.235-3.22-.135-.303-.54-1.523.105-3.176 0 0 1.005-.322 3.3 1.23.96-.267 1.98-.399 3-.405 1.02.006 2.04.138 3 .405 2.28-1.552 3.285-1.23 3.285-1.23.645 1.653.24 2.873.12 3.176.765.84 1.23 1.91 1.23 3.22 0 4.61-2.805 5.625-5.475 5.92.42.36.81 1.096.81 2.22 0 1.606-.015 2.896-.015 3.286 0 .315.21.69.825.57C20.565 22.092 24 17.592 24 12.297c0-6.627-5.373-12-12-12"/>
                                </svg>
                            </a>
                        </span>
                    </span>
                </div>
            </div>

            <div className="cistrome-explorer">
                <div className="container">
                    <CistromeHGW 
                        viewConfig={demos[selectedDemo].viewConfig}
                        options={demos[selectedDemo].options}
                        onViewConfigChange={onViewConfigChange}
                    />
                </div>
            </div>
        </div>
    );
};
