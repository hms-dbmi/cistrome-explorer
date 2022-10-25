import * as gt from "gosling-theme";
import { getCistromeTrack, getExampleCistromeTracks } from "../cistrome-api/cistrome-track";

export const theme = gt.getTheme();

export const hgDemoViewConfig3K4Revision ={
    "editable": false,
    "zoomFixed": false,
    "trackSourceServers": [
        "https://higlass.io/api/v1"
    ],
    "exportViewUrl": "/api/v1/viewconfs",
    "views": [
        {
            "initialXDomain": [
                1518882672.9741557,
                1519310424.8335946
            ],
            "initialYDomain": [
                1518963012.1501117,
                1519273483.6875255
            ],
            "genomePositionSearchBoxVisible": true,
            "chromInfoPath": "//aveit.s3.amazonaws.com/higlass/data/sequence/hg38.mod.chrom.sizes",
            "genomePositionSearchBox": {
                "autocompleteServer": "https://higlass.io/api/v1",
                "autocompleteId": "P0PLbQMwTYGy-5uPIQid7A",
                "chromInfoServer": "https://higlass.io/api/v1",
                "visible": true,
                "chromInfoId": "hg38"
            },
            "tracks": {
                "top": [
                    {
                        "type": "scale-legend",
                        "options": {
                            "background": "#F6F6F6",
                            "color": "#333333",
                            "trackStrokeWidth": 0,
                            "lineStyle": "dotted",
                            "text": "Gene annotation (hg38)"
                        },
                        "height": 25
                    },
                    {
                        "type": "horizontal-chromosome-labels",
                        "server": "https://higlass.io/api/v1",
                        "tilesetUid": "NyITQvZsS_mOFNlz5C2LJg",
                        "uid": "chromosome-labels-track-atac",
                        "chromInfoPath": "//aveit.s3.amazonaws.com/higlass/data/sequence/hg38.mod.chrom.sizes",
                        "options": {
                            "showMousePosition": true,
                            "mousePositionColor": "#000000",
                            "trackBorderWidth": 0,
                            "trackBorderColor": "black",
                            "color": "black",
                            "stroke": "#F6F6F6",
                            "fontSize": 14,
                            "fontIsLeftAligned": false,
                            "reverseOrientation": false
                        },
                        "height": 25
                    },
                    {
                        uid: "gwas",
                        height: 30,
                        width: 2048,
                        type: "gosling-track",
                        data: {
                            chromosomeField: "CHR_ID",
                            genomicFields: ["CHR_POS"],
                            quantitativeFields: ["P-VALUE"],
                            sampleLength: 10000,
                            type: "csv",
                            url: "https://s3.amazonaws.com/gosling-lang.org/data/test/filtered_gwas_catalog_v1.0.2-associations_e104_r2021-10-06.tsv",
                        },
                        options: {
                            showMousePosition: true,
                            mousePositionColor: "black",
                            name: "hg38 | GWAS",
                            labelPosition: "topLeft",
                            fontSize: 12,
                            labelColor: "black",
                            labelShowResolution: false,
                            labelBackgroundColor: "#F6F6F6",
                            labelTextOpacity: 0.6,
                            labelLeftMargin: 4,
                            labelRightMargin: 0,
                            labelTopMargin: 2,
                            labelBottomMargin: 0,
                            backgroundColor: "transparent",
                            theme,
                            spec: {
                                height: 30,
                                width: 2048,
                                assembly: "hg38",
                                layout: "linear",
                                orientation: "horizontal",
                                static: false,
                                centerRadius: 0.3,
                                spacing: 20,
                                data: {
                                    chromosomeField: "CHR_ID",
                                    genomicFields: ["CHR_POS"],
                                    quantitativeFields: ["P-VALUE"],
                                    sampleLength: 10000,
                                    type: "csv",
                                    url: "https://s3.amazonaws.com/gosling-lang.org/data/test/filtered_gwas_catalog_v1.0.2-associations_e104_r2021-10-06.tsv",
                                },
                                mark: "point",
                                opacity: { value: 0.5 },
                                stroke: { value: "white" },
                                strokeWidth: { value: 0 },
                                style: { outline: "#F6F6F6", enableSmoothPath: true, outlineWidth: 1 },
                                x: { field: "CHR_POS", type: "genomic" },
                                y: {
                                    field: "PVALUE_MLOG",
                                    type: "quantitative",
                                    range: [3, 25],
                                    axis: "none",
                                },
                                size: { field: "PVALUE_MLOG", type: "quantitative" },
                                color: {
                                    value: "gray"
                                },
                                tooltip: [
                                    { field: "MAPPED_TRAIT", type: "nominal", alt: "Trait" },
                                    { field: "LINK", type: "nominal", alt: "Link" },
                                    {
                                        field: "PVALUE_MLOG",
                                        type: "quantitative",
                                        alt: "-log(p-value)",
                                    }
                                ],
                            },
                        },
                    },
                    {
                        type: "gosling-track",
                        height: 18,
                        options: {
                            showMousePosition: true,
                            mousePositionColor: "#000000",
                            backgroundColor: "transparent",
                            name: "hg38 | Cytoband",
                            fontSize: 12,
                            labelColor: "black",
                            labelPosition: "topLeft",
                            labelBackgroundColor: "#F6F6F6",
                            labelTextOpacity: 0.6,
                            labelLeftMargin: 4,
                            labelRightMargin: 0,
                            labelTopMargin: 2,
                            labelBottomMargin: 0,
                            theme,
                            spec: {
                                data: {
                                    url: "https://raw.githubusercontent.com/sehilyi/gemini-datasets/master/data/UCSC.HG38.Human.CytoBandIdeogram.csv",
                                    type: "csv",
                                    chromosomeField: "Chromosome",
                                    genomicFields: ["chromStart", "chromEnd"],
                                    assembly: "hg38",
                                },
                                overlay: [
                                    {
                                        mark: "text",
                                        dataTransform: [
                                            {
                                                type: "filter", field: "Stain", oneOf: ["acen"], not: true,
                                            },
                                        ],
                                        text: { field: "Name", type: "nominal" },
                                        color: {
                                            field: "Stain",
                                            type: "nominal",
                                            domain: ["gneg", "gpos25", "gpos50", "gpos75", "gpos100", "gvar"],
                                            range: ["black", "black", "black", "black", "white", "black"],
                                        },
                                        size: { value: 12 },
                                        visibility: [{
                                            operation: "less-than",
                                            measure: "width",
                                            threshold: "|xe-x|",
                                            transitionPadding: 10,
                                            target: "mark",
                                        }],
                                        style: {
                                            textStrokeWidth: 0,
                                        },
                                    },
                                    {
                                        mark: "rect",
                                        dataTransform: [
                                            {
                                                type: "filter", field: "Stain", oneOf: ["acen"], not: true,
                                            },
                                        ],
                                        color: {
                                            field: "Stain",
                                            type: "nominal",
                                            domain: [
                                                "gneg",
                                                "gpos25",
                                                "gpos50",
                                                "gpos75",
                                                "gpos100",
                                                "gvar",
                                            ],
                                            range: [
                                                "white",
                                                "#D9D9D9",
                                                "#979797",
                                                "#636363",
                                                "black",
                                                "#82A3D0",
                                            ],
                                        },
                                    },
                                    {
                                        mark: "triangleRight",
                                        dataTransform: [
                                            {
                                                type: "filter", field: "Name", include: "q", not: false,
                                            },
                                            {
                                                type: "filter", field: "Stain", oneOf: ["acen"], not: false,
                                            },
                                        ],
                                        color: { value: "#E9413B" },
                                    },
                                    {
                                        mark: "triangleLeft",
                                        dataTransform: [
                                            {
                                                type: "filter", field: "Stain", oneOf: ["acen"], not: false,
                                            },
                                            {
                                                type: "filter", field: "Name", include: "p", not: false,
                                            },
                                        ],
                                        color: { value: "#E9413B" },
                                    },
                                ],
                                x: {
                                    field: "chromStart",
                                    type: "genomic",
                                },
                                xe: { field: "chromEnd", type: "genomic" },
                                size: { value: 17 },
                                stroke: { value: "gray" },
                                strokeWidth: { value: 1 },
                                style: { outline: "#F6F6F6" },
                                assembly: "hg38",
                                layout: "linear",
                                orientation: "horizontal",
                                static: true,
                                overlayOnPreviousTrack: false,
                                width: 800,
                                height: 18,
                            },
                        },
                        data: {
                            url: "https://raw.githubusercontent.com/sehilyi/gemini-datasets/master/data/UCSC.HG38.Human.CytoBandIdeogram.csv",
                            type: "csv",
                            chromosomeField: "Chromosome",
                            genomicFields: ["chromStart", "chromEnd"],
                            assembly: "hg38",
                        },
                    },
                    ...getExampleCistromeTracks(),
                    {
                        "type": "horizontal-gene-annotations",
                        "server": "https://higlass.io/api/v1",
                        "tilesetUid": "P0PLbQMwTYGy-5uPIQid7A",
                        "uid": "gene-labels-track-atac",
                        "options": {
                            "plusStrandColor": "gray",
                            "minusStrandColor": "gray",
                            "showMousePosition": true,
                            "mousePositionColor": "#000000",
                            "trackBorderWidth": 0,
                            "trackBorderColor": "gray",
                            "name": "hg38 | Gene annotation",
                            "fontSize": 12,
                            "labelColor": "black",
                            "labelPosition": "topLeft",
                            "labelBackgroundColor": "#F6F6F6",
                            "labelTextOpacity": 0.6,
                            "labelLeftMargin": 4,
                            "labelRightMargin": 0,
                            "labelTopMargin": 2,
                            "labelBottomMargin": 0,
                            "minHeight": 24,
                            "geneAnnotationHeight": 14,
                            "geneLabelPosition": "outside",
                            "geneStrandSpacing": 4
                        },
                        "height": 78
                    },
                ],
                "bottom": [],
                "left": [],
                "center": [
                    {
                        "type": "horizontal-multivec",
                        "uid": "cistrome-track-3k4-revision",
                        "tilesetUid": "h3k4me3_revision",
                        "server": "https://server.gosling-lang.org/api/v1",
                        "options": {
                            "labelPosition": "topLeft",
                            "labelColor": "black",
                            "labelTextOpacity": 0.6,
                            "valueScaling": "linear",
                            "trackBorderWidth": 2,
                            "trackBorderColor": "gray",
                            "heatmapValueScaling": "log",
                            "name": "H3K4me3 from Cistrome DB",
                            "labelTopMargin": 2,
                            "labelLeftMargin": 4,
                            "labelBottomMargin": 0,
                            "labelRightMargin": 0,
                            "labelBackgroundColor": "#F6F6F6",
                            "labelShowResolution": false,
                            "backgroundColor": "white",
                            "minHeight": 100,
                            "colorbarBackgroundOpacity": 0,
                            "colorbarPosition": "bottomRight",
                            "zeroValueColor": "white",
                            "selectRowsAggregationWithRelativeHeight": false,
                            "selectRowsAggregationMode": "mean",
                            "showMousePosition": true,
                            "mousePositionColor": "#000000",
                            "colorRange": [
                                "rgba(255,255,255,1)",
                                "rgba(0,0,0,1)"
                            ],
                            "scaleStartPercent": "0.005", 
                            "scaleEndPercent": "1"
                        },
                        "width": 1607,
                        "height": 682
                    }
                ],
                "right": [],
                "whole": [],
                "gallery": []
            },
            "layout": {
                "w": 12,
                "h": 10,
                "x": 0,
                "y": 0
            },
            "uid": "cistrome-view-3k4-revision"
        }
    ],
    "zoomLocks": {
        "locksByViewUid": {},
        "locksDict": {}
    },
    "locationLocks": {
        "locksByViewUid": {},
        "locksDict": {}
    },
    "valueScaleLocks": {
        "locksByViewUid": {},
        "locksDict": {}
    }
};
