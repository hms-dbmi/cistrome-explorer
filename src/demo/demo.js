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

export const demos = {
    "H3K27ac Demo (1 View, Center Track)": {
        viewConfig: hgDemoViewConfig1,
        options: {
            rowInfoAttributes: [
                {field: "qc_frip", type: "quantitative", position: "left"},
                {field: "qc_fastqc", type: "quantitative", position: "left"},
                {field: "Metadata URL", type: "url", position: "left", title: "cid", width: 30},
                {field: "Hierarchical Clustering (Ward)", type: "tree", position: "left", resolveYScale: true, sort: "ascending", width: 200},
                {field: "Cell Type", type: "nominal", position: "right", width: 80},
                {field: "Tissue Type", type: "nominal", position: "right", resolveYScale: true, sort: "ascending", width: 80},
                {field: "Species", type: "nominal", position: "right", width: 80},
            ],
            rowSort: [
                {field: "Cell Type", type: "nominal", order: "ascending"},
                {field: "qc_frip", type: "quantitative", order: "descending"}
            ],
            rowFilter: [
                {field: "Tissue Type", type: "nominal", notOneOf: ["None"]}
            ],
            // rowHighlight: {field: "Cell Type", type: "index", index: [0]}
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