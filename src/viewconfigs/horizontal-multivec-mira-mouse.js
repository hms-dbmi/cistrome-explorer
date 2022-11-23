import * as gt from "gosling-theme";

export const theme = gt.getTheme();

export const hgDemoViewConfigMiraMouse = {
    "editable": false,
    "zoomFixed": false,
    "trackSourceServers": [
        "https://higlass.io/api/v1"
    ],
    "exportViewUrl": "/api/v1/viewconfs",
    "views": [
        {
            "initialXDomain": [
                249250620.9999989,
                492449994.0000011
            ],
            "initialYDomain": [
                249250620.9999989,
                492449994.0000011
            ],
            "genomePositionSearchBoxVisible": true,
            "chromInfoPath": "//s3.amazonaws.com/gosling-lang.org/data/mm10.chrom.sizes",
            "genomePositionSearchBox": {
                "autocompleteServer": "//higlass.io/api/v1",
                "autocompleteId": "QDutvmyiSrec5nX4pA5WGQ",
                "chromInfoServer": "//higlass.io/api/v1",
                "chromInfoId": "mm10",
                "visible": true,
            },
            "tracks": {
                "top": [
                    {
                        "type": "scale-legend",
                        "options": {
                            "background": "#F6F6F6",
                            "color": "#333333",
                            "trackStrokeWidth": 0,
                            "lineStyle": "dotted"
                        },
                        "height": 25
                    },
                    {
                        "type": "horizontal-chromosome-labels",
                        "server": "https://higlass.io/api/v1",
                        "tilesetUid": "NyITQvZsS_mOFNlz5C2LJg",
                        "uid": "chromosome-labels-track-mira",
                        "chromInfoPath": "//s3.amazonaws.com/gosling-lang.org/data/mm10.chrom.sizes",
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
                        type: "gosling-track",
                        height: 18,
                        options: {
                            showMousePosition: true,
                            mousePositionColor: "#000000",
                            backgroundColor: "transparent",
                            name: "mm10 | Cytoband",
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
                                    url: "https://raw.githubusercontent.com/sehilyi/gemini-datasets/master/data/UCSC.MM10.Mus.musculus.CytoBandIdeogram.csv",
                                    type: "csv",
                                    chromosomeField: "Chromosome",
                                    genomicFields: ["chromStart", "chromEnd"],
                                    assembly: "mm10",
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
                                assembly: "mm10",
                                layout: "linear",
                                orientation: "horizontal",
                                static: true,
                                overlayOnPreviousTrack: false,
                                width: 800,
                                height: 18,
                            },
                        },
                        data: {
                            url: "https://raw.githubusercontent.com/sehilyi/gemini-datasets/master/data/UCSC.MM10.Mus.musculus.CytoBandIdeogram.csv",
                            type: "csv",
                            chromosomeField: "Chromosome",
                            genomicFields: ["chromStart", "chromEnd"],
                            assembly: "mm10",
                        },
                    },
                    {
                        "type": "horizontal-gene-annotations",
                        "server": "https://higlass.io/api/v1",
                        "tilesetUid": "Lm2XCdYQSVyRDkPABIUKGA",
                        "uid": "gene-labels-track-mira",
                        "options": {
                            "plusStrandColor": "gray",
                            "minusStrandColor": "gray",
                            "showMousePosition": true,
                            "mousePositionColor": "#000000",
                            "trackBorderWidth": 0,
                            "trackBorderColor": "gray",
                            "name": "mm10 | Gene annotation",
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
                        "uid": "mira-track-mouse",
                        "tilesetUid": "e18_mouse_brain_10x_dataset_500_random_rows",
                        "server": "https://server.gosling-lang.org/api/v1",
                        "options": {
                            "labelPosition": "topLeft",
                            "labelColor": "black",
                            "labelTextOpacity": 0.6,
                            "valueScaling": "linear",
                            "trackBorderWidth": 2,
                            "trackBorderColor": "gray",
                            "heatmapValueScaling": "log",
                            "name": "e18_mouse_brain_10x_dataset_500_rows",
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
                            "scaleStartPercent": "0", 
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
            "uid": "mira-view-mouse"
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

export const hgDemoViewConfigMiraMouse4000 = {
    "editable": false,
    "zoomFixed": false,
    "trackSourceServers": [
        "https://higlass.io/api/v1"
    ],
    "exportViewUrl": "/api/v1/viewconfs",
    "views": [
        {
            "initialXDomain": [
                249250620.9999989,
                492449994.0000011
            ],
            "initialYDomain": [
                249250620.9999989,
                492449994.0000011
            ],
            "genomePositionSearchBoxVisible": true,
            "chromInfoPath": "//s3.amazonaws.com/gosling-lang.org/data/mm10.chrom.sizes",
            "genomePositionSearchBox": {
                "autocompleteServer": "//higlass.io/api/v1",
                "autocompleteId": "QDutvmyiSrec5nX4pA5WGQ",
                "chromInfoServer": "//higlass.io/api/v1",
                "chromInfoId": "mm10",
                "visible": true,
            },
            "tracks": {
                "top": [
                    {
                        "type": "scale-legend",
                        "options": {
                            "background": "#F6F6F6",
                            "color": "#333333",
                            "trackStrokeWidth": 0,
                            "lineStyle": "dotted"
                        },
                        "height": 25
                    },
                    {
                        "type": "horizontal-chromosome-labels",
                        "server": "https://higlass.io/api/v1",
                        "tilesetUid": "NyITQvZsS_mOFNlz5C2LJg",
                        "uid": "chromosome-labels-track-mira",
                        "chromInfoPath": "//s3.amazonaws.com/gosling-lang.org/data/mm10.chrom.sizes",
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
                        type: "gosling-track",
                        height: 18,
                        options: {
                            showMousePosition: true,
                            mousePositionColor: "#000000",
                            backgroundColor: "transparent",
                            name: "mm10 | Cytoband",
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
                                    url: "https://raw.githubusercontent.com/sehilyi/gemini-datasets/master/data/UCSC.MM10.Mus.musculus.CytoBandIdeogram.csv",
                                    type: "csv",
                                    chromosomeField: "Chromosome",
                                    genomicFields: ["chromStart", "chromEnd"],
                                    assembly: "mm10",
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
                                assembly: "mm10",
                                layout: "linear",
                                orientation: "horizontal",
                                static: true,
                                overlayOnPreviousTrack: false,
                                width: 800,
                                height: 18,
                            },
                        },
                        data: {
                            url: "https://raw.githubusercontent.com/sehilyi/gemini-datasets/master/data/UCSC.MM10.Mus.musculus.CytoBandIdeogram.csv",
                            type: "csv",
                            chromosomeField: "Chromosome",
                            genomicFields: ["chromStart", "chromEnd"],
                            assembly: "mm10",
                        },
                    },
                    {
                        "type": "horizontal-gene-annotations",
                        "server": "https://higlass.io/api/v1",
                        "tilesetUid": "Lm2XCdYQSVyRDkPABIUKGA",
                        "uid": "gene-labels-track-mira",
                        "options": {
                            "plusStrandColor": "gray",
                            "minusStrandColor": "gray",
                            "showMousePosition": true,
                            "mousePositionColor": "#000000",
                            "trackBorderWidth": 0,
                            "trackBorderColor": "gray",
                            "name": "mm10 | Gene annotation",
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
                        "uid": "mira-track-mouse-4000",
                        "tilesetUid": "e18_mouse_brain_10x_dataset_4000_random_rows",
                        "server": "https://server.gosling-lang.org/api/v1",
                        "options": {
                            "labelPosition": "topLeft",
                            "labelColor": "black",
                            "labelTextOpacity": 0.6,
                            "valueScaling": "linear",
                            "trackBorderWidth": 2,
                            "trackBorderColor": "gray",
                            "heatmapValueScaling": "log",
                            "name": "e18_mouse_brain_10x_dataset_4000_rows",
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
                            "scaleStartPercent": "0", 
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
            "uid": "mira-view-mouse-4000"
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

export const hgDemoViewConfigMiraMouse500Smooth = {
    "editable": false,
    "zoomFixed": false,
    "trackSourceServers": [
        "https://higlass.io/api/v1"
    ],
    "exportViewUrl": "/api/v1/viewconfs",
    "views": [
        {
            "initialXDomain": [
                249250620.9999989,
                492449994.0000011
            ],
            "initialYDomain": [
                249250620.9999989,
                492449994.0000011
            ],
            "genomePositionSearchBoxVisible": true,
            "chromInfoPath": "//s3.amazonaws.com/gosling-lang.org/data/mm10.chrom.sizes",
            "genomePositionSearchBox": {
                "autocompleteServer": "//higlass.io/api/v1",
                "autocompleteId": "QDutvmyiSrec5nX4pA5WGQ",
                "chromInfoServer": "//higlass.io/api/v1",
                "chromInfoId": "mm10",
                "visible": true,
            },
            "tracks": {
                "top": [
                    {
                        "type": "scale-legend",
                        "options": {
                            "background": "#F6F6F6",
                            "color": "#333333",
                            "trackStrokeWidth": 0,
                            "lineStyle": "dotted"
                        },
                        "height": 25
                    },
                    {
                        "type": "horizontal-chromosome-labels",
                        "server": "https://higlass.io/api/v1",
                        "tilesetUid": "NyITQvZsS_mOFNlz5C2LJg",
                        "uid": "chromosome-labels-track-mira",
                        "chromInfoPath": "//s3.amazonaws.com/gosling-lang.org/data/mm10.chrom.sizes",
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
                        type: "gosling-track",
                        height: 18,
                        options: {
                            showMousePosition: true,
                            mousePositionColor: "#000000",
                            backgroundColor: "transparent",
                            name: "mm10 | Cytoband",
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
                                    url: "https://raw.githubusercontent.com/sehilyi/gemini-datasets/master/data/UCSC.MM10.Mus.musculus.CytoBandIdeogram.csv",
                                    type: "csv",
                                    chromosomeField: "Chromosome",
                                    genomicFields: ["chromStart", "chromEnd"],
                                    assembly: "mm10",
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
                                assembly: "mm10",
                                layout: "linear",
                                orientation: "horizontal",
                                static: true,
                                overlayOnPreviousTrack: false,
                                width: 800,
                                height: 18,
                            },
                        },
                        data: {
                            url: "https://raw.githubusercontent.com/sehilyi/gemini-datasets/master/data/UCSC.MM10.Mus.musculus.CytoBandIdeogram.csv",
                            type: "csv",
                            chromosomeField: "Chromosome",
                            genomicFields: ["chromStart", "chromEnd"],
                            assembly: "mm10",
                        },
                    },
                    {
                        "type": "horizontal-gene-annotations",
                        "server": "https://higlass.io/api/v1",
                        "tilesetUid": "Lm2XCdYQSVyRDkPABIUKGA",
                        "uid": "gene-labels-track-mira",
                        "options": {
                            "plusStrandColor": "gray",
                            "minusStrandColor": "gray",
                            "showMousePosition": true,
                            "mousePositionColor": "#000000",
                            "trackBorderWidth": 0,
                            "trackBorderColor": "gray",
                            "name": "mm10 | Gene annotation",
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
                        "uid": "mira-track-mouse-500-smooth",
                        "tilesetUid": "e18_mouse_brain_10x_dataset_500_smoothed_random_rows",
                        "server": "https://server.gosling-lang.org/api/v1",
                        "options": {
                            "labelPosition": "topLeft",
                            "labelColor": "black",
                            "labelTextOpacity": 0.6,
                            "valueScaling": "linear",
                            "trackBorderWidth": 2,
                            "trackBorderColor": "gray",
                            "heatmapValueScaling": "log",
                            "name": "e18_mouse_brain_10x_dataset_100_smooth_rows",
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
            "uid": "mira-view-mouse-500-smooth"
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
