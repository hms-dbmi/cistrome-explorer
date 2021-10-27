import * as gt from "gosling-theme";

export const theme = gt.getTheme();

export const hgDemoViewConfigAtac ={
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
                        "type": "combined",
                        "width": 800,
                        "height": 30,
                        "contents": [
                            {
                                "type": "gosling-track",
                                "server": "https://server.gosling-lang.org/api/v1/",
                                "tilesetUid": "gwas-multivec",
                                "options": {
                                    "showMousePosition": true,
                                    "mousePositionColor": "black",
                                    "labelPosition": "none",
                                    "fontSize": 12,
                                    "labelColor": "black",
                                    "labelShowResolution": false,
                                    "labelBackgroundColor": "white",
                                    "labelTextOpacity": 1,
                                    "labelLeftMargin": 1,
                                    "labelTopMargin": 1,
                                    "labelRightMargin": 0,
                                    "labelBottomMargin": 0,
                                    "backgroundColor": "transparent",
                                    "spec": {
                                        "width": 800,
                                        "height": 30,
                                        "assembly": "hg38",
                                        "orientation": "horizontal",
                                        "static": false,
                                        "centerRadius": 0.3,
                                        "spacing": 20,
                                        "data": {
                                            "url": "https://server.gosling-lang.org/api/v1/tileset_info/?d=gwas-multivec",
                                            "type": "multivec",
                                            "row": "base",
                                            "column": "position",
                                            "value": "count",
                                            "categories": ["A"],
                                            "start": "start",
                                            "end": "end",
                                            // "binSize": 2
                                        },
                                        "mark": "rect",
                                        "x": {
                                            "field": "start",
                                            "type": "genomic",
                                            "axis": "none"
                                        },
                                        "xe": {"field": "end", "type": "genomic"},
                                        "opacity": {"field": "count", "type": "quantitative"},
                                        "color": {"value": "gray"},
                                        "style": {"outline": "#F6F6F6"},
                                        "overlayOnPreviousTrack": false
                                    }
                                }
                            },
                            {
                                "type": "gosling-track",
                                "server": "https://server.gosling-lang.org/api/v1/",
                                "tilesetUid": "gwas-beddb",
                                "options": {
                                    "showMousePosition": true,
                                    "mousePositionColor": "black",
                                    "name": "hg38 | GWAS",
                                    "labelPosition": "topLeft",
                                    "fontSize": 12,
                                    "labelColor": "black",
                                    "labelShowResolution": false,
                                    "labelBackgroundColor": "#F6F6F6",
                                    "labelTextOpacity": 0.6,
                                    "labelLeftMargin": 4,
                                    "labelRightMargin": 0,
                                    "labelTopMargin": 2,
                                    "labelBottomMargin": 0,
                                    "backgroundColor": "transparent",
                                    "spec": {
                                        "height": 30,
                                        "width": 2048,
                                        "assembly": "hg38",
                                        "layout": "linear",
                                        "orientation": "horizontal",
                                        "static": false,
                                        "centerRadius": 0.3,
                                        "spacing": 20,
                                        "data": {
                                            "url": "https://server.gosling-lang.org/api/v1/tileset_info/?d=gwas-beddb",
                                            "type": "beddb",
                                            "genomicFields": [
                                                {"index": 1, "name": "start"},
                                                {"index": 2, "name": "end"}
                                            ],
                                            "valueFields": [
                                                {"index": 3, "name": "pubmedid", "type": "nominal"},
                                                {"index": 4, "name": "date", "type": "nominal"},
                                                {"index": 5, "name": "link", "type": "nominal"},
                                                {"index": 6, "name": "pvalue", "type": "quantitative"},
                                                {"index": 8, "name": "disease", "type": "nominal"},
                                                {
                                                    "index": 9,
                                                    "name": "pvalue_log",
                                                    "type": "quantitative"
                                                },
                                                {"index": 10, "name": "pvalue_txt", "type": "nominal"}
                                            ]
                                        },
                                        "mark": "point",
                                        "x": {
                                            "field": "start",
                                            "type": "genomic", "axis": "none"
                                        },
                                        "xe": {"field": "end", "type": "genomic"},
                                        "y": {
                                            "field": "pvalue_log",
                                            "type": "quantitative",
                                            "range": [3, 25]
                                        },
                                        "size": {"field": "pvalue_log", "type": "quantitative"},
                                        "color": {
                                            "field": "pvalue_log",
                                            "type": "quantitative",
                                            "range": "warm"
                                        },
                                        "tooltip": [
                                            {"field": "disease", "type": "nominal", "alt": "Disease"},
                                            {"field": "link", "type": "nominal", "alt": "Link"},
                                            { "field": "pvalue", "type": "quantitative", "alt": "p-value" },
                                            {
                                                "field": "pvalue_log",
                                                "type": "quantitative",
                                                "alt": "-log(p-value)"
                                            },
                                            {
                                                "field": "pvalue_txt",
                                                "type": "nominal",
                                                "alt": "Context of p-value"
                                            },
                                            {"field": "pubmedid", "type": "nominal", "alt": "PubMed"}
                                        ],
                                        "style": {"outline": "#F6F6F6"},
                                        "overlayOnPreviousTrack": true
                                    }
                                }
                            }
                        ]
                    },
                    {
                        "type": "gosling-track",
                        "height": 18,
                        "options": {
                            "showMousePosition": true,
                            "mousePositionColor": "#000000",
                            "backgroundColor": "transparent",
                            "name": "hg38 | Cytoband",
                            "fontSize": 12,
                            "labelColor": "black",
                            "labelPosition": "topLeft",
                            "labelBackgroundColor": "#F6F6F6",
                            "labelTextOpacity": 0.6,
                            "labelLeftMargin": 4,
                            "labelRightMargin": 0,
                            "labelTopMargin": 2,
                            "labelBottomMargin": 0,
                            "spec": {
                                "data": {
                                    "url": "https://raw.githubusercontent.com/sehilyi/gemini-datasets/master/data/UCSC.HG38.Human.CytoBandIdeogram.csv",
                                    "type": "csv",
                                    "chromosomeField": "Chromosome",
                                    "genomicFields": ["chromStart", "chromEnd"],
                                    "assembly": "hg38"
                                },
                                "overlay": [
                                    {
                                        mark: "text",
                                        dataTransform: {
                                            filter: [{ field: "Stain", oneOf: ["acen"], not: true }]
                                        },
                                        text: { field: "Name", type: "nominal" },
                                        color: {
                                            field: "Stain",
                                            type: "nominal",
                                            domain: ["gneg", "gpos25", "gpos50", "gpos75", "gpos100", "gvar"],
                                            range: ["black", "black", "black", "black", "white", "black"]
                                        },
                                        visibility: [{
                                            operation: "less-than",
                                            measure: "width",
                                            threshold: "|xe-x|", 
                                            transitionPadding: 10,
                                            target: "mark"
                                        }],
                                        style: {
                                            textStrokeWidth: 0
                                        }
                                    },
                                    {
                                        "mark": "rect",
                                        "dataTransform": {
                                            "filter": [
                                                {"field": "Stain", "oneOf": ["acen"], "not": true}
                                            ]
                                        },
                                        "color": {
                                            "field": "Stain",
                                            "type": "nominal",
                                            "domain": [
                                                "gneg",
                                                "gpos25",
                                                "gpos50",
                                                "gpos75",
                                                "gpos100",
                                                "gvar"
                                            ],
                                            "range": [
                                                "white",
                                                "#D9D9D9",
                                                "#979797",
                                                "#636363",
                                                "black",
                                                "#82A3D0"
                                            ]
                                        }
                                    },
                                    {
                                        "mark": "triangleRight",
                                        "dataTransform": {
                                            "filter": [{"field": "Name", "include": "q", "not": false}, {"field": "Stain", "oneOf": ["acen"], "not": false}]
                                        },
                                        "color": {"value": "#E9413B"}
                                    },
                                    {
                                        "mark": "triangleLeft",
                                        "dataTransform": {
                                            "filter": [
                                                {"field": "Stain", "oneOf": ["acen"], "not": false},
                                                {"field": "Name", "include": "p", "not": false}
                                            ]
                                        },
                                        "color": {"value": "#E9413B"}
                                    }
                                ],
                                "x": {
                                    "field": "chromStart",
                                    "type": "genomic",
                                },
                                "xe": {"field": "chromEnd", "type": "genomic"},
                                "size": {"value": 17},
                                "stroke": {"value": "gray"},
                                "strokeWidth": {"value": 1},
                                "style": {"outline": "#F6F6F6"},
                                "assembly": "hg38",
                                "layout": "linear",
                                "orientation": "horizontal",
                                "static": true,
                                "overlayOnPreviousTrack": false,
                                "width": 800,
                                "height": 18
                            }
                        },
                        "data": {
                            "url": "https://raw.githubusercontent.com/sehilyi/gemini-datasets/master/data/UCSC.HG38.Human.CytoBandIdeogram.csv",
                            "type": "csv",
                            "chromosomeField": "Chromosome",
                            "genomicFields": ["chromStart", "chromEnd"],
                            "assembly": "hg38"
                        }
                    },
                    // {
                    //   "data": {
                    //     "type": "cistrome-bigwig",
                    //     "cid": "1",
                    //     // "url": "http://dbtoolkit.cistrome.org/api_bigwig?sid=1",
                    //     "chromSizesUrl": "https://aveit.s3.amazonaws.com/higlass/data/sequence/hg38.chrom.sizes",
                    //   },
                    //   "uid": "bw",
                    //   "type": "bar",
                    //   "options": {
                    //     "align": "bottom",
                    //     "labelColor": "[glyph-color]",
                    //     "labelPosition": "topLeft",
                    //     "labelLeftMargin": 0,
                    //     "labelRightMargin": 0,
                    //     "labelTopMargin": 0,
                    //     "labelBottomMargin": 0,
                    //     "labelShowResolution": false,
                    //     "labelShowAssembly": true,
                    //     "axisLabelFormatting": "scientific",
                    //     "axisPositionHorizontal": "right",
                    //     "barFillColor": "darkgreen",
                    //     "valueScaling": "linear",
                    //     "trackBorderWidth": 0,
                    //     "trackBorderColor": "black",
                    //     "labelTextOpacity": 0.4,
                    //     "barOpacity": 1,
                    //     "valueScaleMin": 0,
                    //     //"valueScaleMax": 80,
                    //     "name": "GnomAD coverage"
                    //   },
                    //   "width": 20,
                    //   "height": 100
                    // },
                    {
                        "type": "horizontal-gene-annotations",
                        "server": "https://higlass.io/api/v1",
                        "tilesetUid": "P0PLbQMwTYGy-5uPIQid7A",
                        "uid": "gene-labels-track-atac",
                        "options": {
                            "plusStrandColor": "#004FA5",
                            "minusStrandColor": "#E10003",
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
                    //{
                    //  "type": "horizontal-stacked-bar",
                    //  "uid": "cistrome-track-atac-aggregated",
                    //  "tilesetUid": "UvVPeLHuRDiYA3qwFlm7xQ",
                    //  "server": "https://resgen.io/api/v1",
                    //  "options": {
                    //    "labelPosition": "topLeft",
                    //    "labelColor": "black",
                    //    "labelTextOpacity": 0.6,
                    //    "barBorder": false,
                    //    "valueScaling": "linear",
                    //    "heatmapValueScaling": "log",
                    //    "name": "hg38 | Aggregation of all 200 samples",
                    //    "labelTopMargin": 2,
                    //    "labelLeftMargin": 4,
                    //    "labelBottomMargin": 0,
                    //    "labelRightMargin": 0,
                    //    "labelBackgroundColor": "#F6F6F6",
                    //    "labelShowResolution": false,
                    //    "backgroundColor": "#F6F6F6",
                    //    "minHeight": 100,
                    //    "colorbarPosition": "hidden",
                    //    "colorScale": ["gray"],
                    //    "zeroValueColor": "white",
                    //    "selectRowsAggregationMethod": "server",
                    //    "selectRowsAggregationWithRelativeHeight": false,
                    //    "selectRowsAggregationMode": "mean",
                    //    "showMousePosition": true,
                    //    "mousePositionColor": "#000000",
                    //    "selectRows": [[
                    //      206,
                    //      71,
                    //      216,
                    //      16,
                    //      15,
                    //      24,
                    //      23,
                    //      40,
                    //      29,
                    //      37,
                    //      111,
                    //      112,
                    //      61,
                    //      26,
                    //      135,
                    //      84,
                    //      142,
                    //      82,
                    //      89,
                    //      127,
                    //      83,
                    //      80,
                    //      122,
                    //      103,
                    //      126,
                    //      143,
                    //      78,
                    //      68,
                    //      139,
                    //      102,
                    //      123,
                    //      69,
                    //      79,
                    //      94,
                    //      132,
                    //      131,
                    //      141,
                    //      130,
                    //      140,
                    //      43,
                    //      42,
                    //      41,
                    //      44,
                    //      45,
                    //      46,
                    //      48,
                    //      50,
                    //      65,
                    //      66,
                    //      67,
                    //      70,
                    //      74,
                    //      75,
                    //      81,
                    //      96,
                    //      124,
                    //      125,
                    //      137,
                    //      138,
                    //      148,
                    //      149,
                    //      150,
                    //      151,
                    //      209,
                    //      210,
                    //      211,
                    //      217,
                    //      218,
                    //      219,
                    //      220,
                    //      225,
                    //      226,
                    //      154,
                    //      13,
                    //      147,
                    //      195,
                    //      196,
                    //      197,
                    //      198,
                    //      199,
                    //      200,
                    //      201,
                    //      202,
                    //      203,
                    //      204,
                    //      205,
                    //      213,
                    //      90,
                    //      88,
                    //      115,
                    //      117,
                    //      2,
                    //      118,
                    //      113,
                    //      212,
                    //      134,
                    //      152,
                    //      64,
                    //      63,
                    //      72,
                    //      51,
                    //      19,
                    //      20,
                    //      163,
                    //      162,
                    //      164,
                    //      160,
                    //      159,
                    //      161,
                    //      156,
                    //      157,
                    //      155,
                    //      158,
                    //      172,
                    //      166,
                    //      165,
                    //      180,
                    //      189,
                    //      188,
                    //      174,
                    //      168,
                    //      173,
                    //      176,
                    //      181,
                    //      179,
                    //      184,
                    //      178,
                    //      171,
                    //      175,
                    //      183,
                    //      182,
                    //      186,
                    //      177,
                    //      169,
                    //      167,
                    //      185,
                    //      187,
                    //      170,
                    //      27,
                    //      34,
                    //      28,
                    //      136,
                    //      153,
                    //      57,
                    //      73,
                    //      146,
                    //      47,
                    //      145,
                    //      191,
                    //      190,
                    //      192,
                    //      8,
                    //      7,
                    //      6,
                    //      31,
                    //      214,
                    //      14,
                    //      9,
                    //      87,
                    //      86,
                    //      85,
                    //      95,
                    //      106,
                    //      109,
                    //      108,
                    //      76,
                    //      97,
                    //      98,
                    //      99,
                    //      77,
                    //      227,
                    //      223,
                    //      224,
                    //      121,
                    //      105,
                    //      17,
                    //      18,
                    //      221,
                    //      104,
                    //      119,
                    //      222,
                    //      128,
                    //      129,
                    //      107,
                    //      101,
                    //      100,
                    //      110,
                    //      133,
                    //      144,
                    //      5,
                    //      93,
                    //      59,
                    //      120,
                    //      4,
                    //      60,
                    //      58,
                    //      53,
                    //      21,
                    //      52,
                    //      55,
                    //      56,
                    //      35,
                    //      54,
                    //      49,
                    //      193,
                    //      32,
                    //      33,
                    //      215,
                    //      1,
                    //      91,
                    //      62,
                    //      92,
                    //      208,
                    //      207,
                    //      194,
                    //      114,
                    //      116,
                    //      39,
                    //      12,
                    //      11,
                    //      10,
                    //      25,
                    //      30,
                    //      38,
                    //      0,
                    //      22,
                    //      3,
                    //      36
                    //    ]]
                    //  },
                    //  "width": 1607,
                    //  "height": 31
                    //}
                ],//.slice(0, 3),
                "bottom": [],
                "left": [],
                "center": [
                    {
                        "type": "horizontal-multivec",
                        "uid": "cistrome-track-atac",
                        "tilesetUid": "cistrome-atac-multivec",
                        "server": "https://server.gosling-lang.org/api/v1",
                        "options": {
                            "labelPosition": "topLeft",
                            "labelColor": "black",
                            "labelTextOpacity": 0.6,
                            "valueScaling": "linear",
                            "trackBorderWidth": 2,
                            "trackBorderColor": "gray",
                            "heatmapValueScaling": "log",
                            "name": "ATAC-seq from Cistrome DB",
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
            "uid": "cistrome-view-atac"
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
