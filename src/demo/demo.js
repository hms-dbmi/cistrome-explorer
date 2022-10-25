import { hgDemoViewConfigAtac } from "../viewconfigs/horizontal-multivec-atac";
import { hgDemoViewConfigAtacRevision } from "../viewconfigs/horizontal-multivec-atac-revision";
import { hgDemoViewConfig3k27 } from "../viewconfigs/horizontal-multivec-3k27";
import { hgDemoViewConfig3k27Revision } from "../viewconfigs/horizontal-multivec-3k27-revision";
import { hgDemoViewConfig3K4 } from "../viewconfigs/horizontal-multivec-3k4";
import { hgDemoViewConfigMiraMouse } from '../viewconfigs/horizontal-multivec-mira-mouse';

export const demos = {
    "ATAC": {
        viewConfig: hgDemoViewConfigAtac,
        options: [{
            viewId: "cistrome-view-atac",
            trackId: "cistrome-track-atac",
            rowInfoAttributes: [
                {field: "Cell Type", type: "nominal", position: "right", width: 200},
                // {field: "Tissue Type", type: "nominal", position: "right", width: 120},
                {field: "Clustering", type: "tree", position: "right", width: 200},
                // { field: "FRACTION_PEAKS_IN_EXON", type: "quantitative", position: "right", width: 80 },
                // { field: "FRACTION_PEAKS_IN_INTRON", type: "quantitative", position: "right", width: 80 },
                // { field: "FRACTION_PEAKS_IN_PROMOTER", type: "quantitative", position: "right", width: 80 },
                // { field: "FRACTION_PEAKS_IN_INTERGENIC", type: "quantitative", position: "right", width: 80 },
                { field: "PEAKS_TOTAL", type: "quantitative", title: 'Total peaks', position: "right", width: 100 },
                // { field: "PEAKS_10FOLD", type: "quantitative", position: "right", width: 80 },
                // { field: "FRACTION_READS_IN_PEAKS", type: "quantitative", position: "right", width: 80 },
                // { field: "FASTQC_MEDIAN", type: "quantitative", position: "right", width: 80 },
                { field: "SEQUENCE_LENGTH", type: "quantitative", title: 'Sequence length', position: "right", width: 100 },
                // { field: "FRACTION_PEAKS_IN_DHS", type: "quantitative", position: "right", width: 80 },
                // { field: "PCR_BOTTLENECK_COEF", type: "quantitative", position: "right", width: 80 },
                // { field: "READS_TOTAL", type: "quantitative", position: "right", width: 80 },
                { field: "READS_MAPPED", type: "quantitative", title: 'Reads mapped', position: "right", width: 100 },
                // { field: "PHYLOCONS_0_BP_OFFSET", type: "quantitative", position: "right", width: 80 },
                // { field: "PHYLOCONS_500_BP_OFFSET", type: "quantitative", position: "right", width: 80 },
                // { field: "PHYLOCONS_1000_BP_OFFSET", type: "quantitative", position: "right", width: 80 }
                // {field: "Species", type: "nominal", position: "left", width: 120},
                // {field: "cid", type: "nominal-dynamic",  position: "left", title: "Compare Positive and Negative", domain: ["positive", "negative"], range: ["blue", "red"], width: 30},
                // {field: "ID", type: "url", shortName: "📊", position: "right", title: "➕ Add Detail Track", width: 35, addTrackOnClick: true},
                // {field: "Metadata URL", alt: "cid", type: "url", position: "right", width: 70},
            ],
            rowAggregate: [
                // {field: "Cell Type", type: "nominal", notOneOf: []},
                // {field: "Tissue Type", type: "nominal", notOneOf: []}
            ],
            rowSort: [
                {field: "Clustering", type: "tree", order: "ascending"}
            ],
            rowFilter: [
                // ...
            ]
        }]
    },
    "ATAC (v2.0)": {
        viewConfig: hgDemoViewConfigAtacRevision,
        options: [{
            viewId: "cistrome-view-atac-revision",
            trackId: "cistrome-track-atac-revision",
            rowInfoAttributes: [
                {field: "Cell Type", type: "nominal", position: "right", width: 200},
                // {field: "Tissue Type", type: "nominal", position: "right", width: 120},
                {field: "Clustering", type: "tree", position: "right", width: 200},
                // { field: "FRACTION_PEAKS_IN_EXON", type: "quantitative", position: "right", width: 80 },
                // { field: "FRACTION_PEAKS_IN_INTRON", type: "quantitative", position: "right", width: 80 },
                // { field: "FRACTION_PEAKS_IN_PROMOTER", type: "quantitative", position: "right", width: 80 },
                // { field: "FRACTION_PEAKS_IN_INTERGENIC", type: "quantitative", position: "right", width: 80 },
                { field: "PEAKS_TOTAL", type: "quantitative", title: 'Total peaks', position: "right", width: 100 },
                // { field: "PEAKS_10FOLD", type: "quantitative", position: "right", width: 80 },
                // { field: "FRACTION_READS_IN_PEAKS", type: "quantitative", position: "right", width: 80 },
                // { field: "FASTQC_MEDIAN", type: "quantitative", position: "right", width: 80 },
                { field: "SEQUENCE_LENGTH", type: "quantitative", title: 'Sequence length', position: "right", width: 100 },
                // { field: "FRACTION_PEAKS_IN_DHS", type: "quantitative", position: "right", width: 80 },
                // { field: "PCR_BOTTLENECK_COEF", type: "quantitative", position: "right", width: 80 },
                // { field: "READS_TOTAL", type: "quantitative", position: "right", width: 80 },
                { field: "READS_MAPPED", type: "quantitative", title: 'Reads mapped', position: "right", width: 100 },
                // { field: "PHYLOCONS_0_BP_OFFSET", type: "quantitative", position: "right", width: 80 },
                // { field: "PHYLOCONS_500_BP_OFFSET", type: "quantitative", position: "right", width: 80 },
                // { field: "PHYLOCONS_1000_BP_OFFSET", type: "quantitative", position: "right", width: 80 }
                // {field: "Species", type: "nominal", position: "left", width: 120},
                // {field: "cid", type: "nominal-dynamic",  position: "left", title: "Compare Positive and Negative", domain: ["positive", "negative"], range: ["blue", "red"], width: 30},
                // {field: "ID", type: "url", shortName: "📊", position: "right", title: "➕ Add Detail Track", width: 35, addTrackOnClick: true},
                // {field: "Metadata URL", alt: "cid", type: "url", position: "right", width: 70},
            ],
            rowAggregate: [
                // {field: "Cell Type", type: "nominal", notOneOf: []},
                // {field: "Tissue Type", type: "nominal", notOneOf: []}
            ],
            rowSort: [
                {field: "Clustering", type: "tree", order: "ascending"}
            ],
            rowFilter: [
                // ...
            ]
        }]
    },
    "H3K27ac": {
        viewConfig: hgDemoViewConfig3k27,
        options: [{
            viewId: "cistrome-view-3k27",
            trackId: "cistrome-track-3k27",
            rowInfoAttributes: [
                {field: "Cell Type", type: "nominal", position: "right", width: 200},
                {field: "Clustering", type: "tree", position: "right", width: 200},
                // { field: "FRACTION_PEAKS_IN_EXON", type: "quantitative", position: "right", width: 80 },
                // { field: "FRACTION_PEAKS_IN_INTRON", type: "quantitative", position: "right", width: 80 },
                // { field: "FRACTION_PEAKS_IN_PROMOTER", type: "quantitative", position: "right", width: 80 },
                // { field: "FRACTION_PEAKS_IN_INTERGENIC", type: "quantitative", position: "right", width: 80 },
                { field: "PEAKS_TOTAL", type: "quantitative", title: 'Total peaks', position: "right", width: 100 },
                // { field: "PEAKS_10FOLD", type: "quantitative", position: "right", width: 80 },
                // { field: "FRACTION_READS_IN_PEAKS", type: "quantitative", position: "right", width: 80 },
                // { field: "FASTQC_MEDIAN", type: "quantitative", position: "right", width: 80 },
                { field: "SEQUENCE_LENGTH", type: "quantitative", title: 'Sequence length', position: "right", width: 100 },
                // { field: "FRACTION_PEAKS_IN_DHS", type: "quantitative", position: "right", width: 80 },
                // { field: "PCR_BOTTLENECK_COEF", type: "quantitative", position: "right", width: 80 },
                // { field: "READS_TOTAL", type: "quantitative", position: "right", width: 80 },
                { field: "READS_MAPPED", type: "quantitative", title: 'Reads mapped', position: "right", width: 100 },
                // { field: "PHYLOCONS_0_BP_OFFSET", type: "quantitative", position: "right", width: 80 },
                // { field: "PHYLOCONS_500_BP_OFFSET", type: "quantitative", position: "right", width: 80 },
                // { field: "PHYLOCONS_1000_BP_OFFSET", type: "quantitative", position: "right", width: 80 }
            ],
            rowAggregate: [
            ],
            rowSort: [
                {field: "Clustering", type: "tree", order: "ascending"}
            ],
            rowFilter: [
            ]
        }]
    },
    "H3K27ac (v2.0)": {
        viewConfig: hgDemoViewConfig3k27Revision,
        options: [{
            viewId: "cistrome-view-3k27-revision",
            trackId: "cistrome-track-3k27-revision",
            rowInfoAttributes: [
                {field: "Cell Type", type: "nominal", position: "right", width: 200},
                {field: "Clustering", type: "tree", position: "right", width: 200},
                // { field: "FRACTION_PEAKS_IN_EXON", type: "quantitative", position: "right", width: 80 },
                // { field: "FRACTION_PEAKS_IN_INTRON", type: "quantitative", position: "right", width: 80 },
                // { field: "FRACTION_PEAKS_IN_PROMOTER", type: "quantitative", position: "right", width: 80 },
                // { field: "FRACTION_PEAKS_IN_INTERGENIC", type: "quantitative", position: "right", width: 80 },
                { field: "PEAKS_TOTAL", type: "quantitative", title: 'Total peaks', position: "right", width: 100 },
                // { field: "PEAKS_10FOLD", type: "quantitative", position: "right", width: 80 },
                // { field: "FRACTION_READS_IN_PEAKS", type: "quantitative", position: "right", width: 80 },
                // { field: "FASTQC_MEDIAN", type: "quantitative", position: "right", width: 80 },
                { field: "SEQUENCE_LENGTH", type: "quantitative", title: 'Sequence length', position: "right", width: 100 },
                // { field: "FRACTION_PEAKS_IN_DHS", type: "quantitative", position: "right", width: 80 },
                // { field: "PCR_BOTTLENECK_COEF", type: "quantitative", position: "right", width: 80 },
                // { field: "READS_TOTAL", type: "quantitative", position: "right", width: 80 },
                { field: "READS_MAPPED", type: "quantitative", title: 'Reads mapped', position: "right", width: 100 },
                // { field: "PHYLOCONS_0_BP_OFFSET", type: "quantitative", position: "right", width: 80 },
                // { field: "PHYLOCONS_500_BP_OFFSET", type: "quantitative", position: "right", width: 80 },
                // { field: "PHYLOCONS_1000_BP_OFFSET", type: "quantitative", position: "right", width: 80 }
            ],
            rowAggregate: [
            ],
            rowSort: [
                {field: "Clustering", type: "tree", order: "ascending"}
            ],
            rowFilter: [
            ]
        }]
    },
    "H3K4me3": {
        viewConfig: hgDemoViewConfig3K4,
        options: [{
            viewId: "cistrome-view-3k4",
            trackId: "cistrome-track-3k4",
            rowInfoAttributes: [
                {field: "Cell Type", type: "nominal", position: "right", width: 200},
                {field: "Clustering", type: "tree", position: "right", width: 200},
                { field: "PEAKS_TOTAL", type: "quantitative", title: 'Total peaks', position: "right", width: 100 },
                { field: "SEQUENCE_LENGTH", type: "quantitative", title: 'Sequence length', position: "right", width: 100 },
                { field: "READS_MAPPED", type: "quantitative", title: 'Reads mapped', position: "right", width: 100 },
            ],
            rowAggregate: [],
            rowSort: [
                // {field: "Clustering", type: "tree", order: "ascending"}
                {field: "Cell Type", type: "nominal", order: "ascending"}
            ],
            rowFilter: []
        }]
    }
};

export const miraDemos = {
    e18_mouse_brain_10x: {
        viewConfig: hgDemoViewConfigMiraMouse,
        options: [{
            viewId: "mira-view-mouse",
            trackId: "mira-track-mouse",
            rowInfoAttributes: [
                {field: ["topic_0", "topic_1", "topic_2", "topic_3", "topic_4", "topic_5", "topic_6", "topic_7", "topic_8", "topic_9", "topic_10", "topic_11", "topic_12"], type: "quantitative", position: "right", width: 60},
                {field: "topic_0", type: "quantitative", position: "right", width: 60},
                {field: "topic_1", type: "quantitative", position: "right", width: 60},
                {field: "topic_2", type: "quantitative", position: "right", width: 60},
                {field: "topic_3", type: "quantitative", position: "right", width: 60},
                {field: "topic_4", type: "quantitative", position: "right", width: 60},
                {field: "topic_5", type: "quantitative", position: "right", width: 60},
                {field: "topic_6", type: "quantitative", position: "right", width: 60},
                {field: "topic_7", type: "quantitative", position: "right", width: 60},
                {field: "topic_8", type: "quantitative", position: "right", width: 60},
                {field: "topic_9", type: "quantitative", position: "right", width: 60},
                {field: "topic_10", type: "quantitative", position: "right", width: 60},
                {field: "topic_11", type: "quantitative", position: "right", width: 60},
                {field: "topic_12", type: "quantitative", position: "right", width: 60},
            ],
            rowAggregate: [],
            rowSort: [
                {field: "topic_9", type: "quantitative", order: "descending"}
            ],
            rowFilter: [
                
            ]
        }]
    }
}