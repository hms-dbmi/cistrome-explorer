import React, { useRef, useEffect } from 'react';
import PubSub from 'pubsub-js';

import { HiGlassComponent } from 'higlass';
import register from 'higlass-register';
import StackedBarTrack from 'higlass-multivec/es/StackedBarTrack.js';

/*import CistromeGroupLabelsTrack from './CistromeGroupLabelsTrack.js';*/
import CistromeGroupLabels from './CistromeGroupLabels.js';

import { GLOBAL_X_RANGE, GLOBAL_Y_RANGE, TRACK_ROW_INFO, TRACK_POSITION, TRACK_DIMENSIONS } from './constants.js';

import 'higlass/dist/hglib.css';
import './CistromeHGW.css';


register({
    name: 'StackedBarTrack',
    track: StackedBarTrack,
    config: StackedBarTrack.config,
});

/*
register({
    name: 'CistromeGroupLabelsTrack',
    track: CistromeGroupLabelsTrack,
    config: CistromeGroupLabelsTrack.config,
});
*/

const demoViewConfig = {
    "editable": true,
    "zoomFixed": false,
    "trackSourceServers": [
      "https://higlass.io/api/v1"
    ],
    "exportViewUrl": "/api/v1/viewconfs",
    "views": [
      {
        "initialXDomain": [
          226426442.00375226,
          226667818.0304663
        ],
        "initialYDomain": [
          2671196920.8697634,
          2671320441.9350777
        ],
        "tracks": {
          "top": [
            {
              "type": "horizontal-chromosome-labels",
              "uid": "GzHGPlijS2eASBqslzVcxA",
              "tilesetUid": "ZpZ8c5JJRUS1J7ZkofcUrg",
              "server": "https://resgen.io/api/v1",
              "options": {
                "showMousePosition": false,
                "mousePositionColor": "#999999",
                "color": "#808080",
                "stroke": "#ffffff",
                "fontSize": 12,
                "fontIsLeftAligned": false
              },
              "width": 1027,
              "height": 30
            },
            {
              "type": "horizontal-gene-annotations",
              "uid": "R6ZZ5LnEQ4ODCWZoE38maw",
              "tilesetUid": "M9A9klpwTci5Vf4bHZ864g",
              "server": "https://resgen.io/api/v1",
              "options": {
                "labelColor": "black",
                "labelPosition": "hidden",
                "plusStrandColor": "blue",
                "minusStrandColor": "red",
                "trackBorderWidth": 0,
                "trackBorderColor": "black",
                "showMousePosition": false,
                "name": "gene-annotations-hg38_3d2Hebg.db",
                "mousePositionColor": "#999999",
                "fontSize": 10,
                "labelBackgroundColor": "#ffffff",
                "labelLeftMargin": 0,
                "labelRightMargin": 0,
                "labelTopMargin": 0,
                "labelBottomMargin": 0,
                "minHeight": 24,
                "geneAnnotationHeight": 12,
                "geneLabelPosition": "outside",
                "geneStrandSpacing": 4
              },
              "width": 568,
              "height": 70
            },
            {
              "type": "horizontal-bar",
              "uid": "bmcNxj_MSEar_3nYVMnPnQ",
              "tilesetUid": "a-iBpdh3Q_uO2FLCWKpOOw",
              "server": "https://resgen.io/api/v1",
              "options": {
                "labelColor": "black",
                "labelPosition": "topLeft",
                "axisPositionHorizontal": "right",
                "barFillColor": "darkgreen",
                "valueScaling": "linear",
                "trackBorderWidth": 0,
                "trackBorderColor": "black",
                "labelTextOpacity": 0.4,
                "barOpacity": 1,
                "name": "conservation (hg38.phastCons100way.bw)",
                "align": "bottom",
                "labelLeftMargin": 0,
                "labelRightMargin": 0,
                "labelTopMargin": 0,
                "labelBottomMargin": 0,
                "labelShowResolution": false,
                "axisLabelFormatting": "scientific"
              },
              "width": 568,
              "height": 42,
              "aggregationModes": {
                "mean": {
                  "name": "Mean",
                  "value": "mean"
                },
                "min": {
                  "name": "Min",
                  "value": "min"
                },
                "max": {
                  "name": "Max",
                  "value": "max"
                },
                "std": {
                  "name": "Standard Deviation",
                  "value": "std"
                }
              }
            },
            {
              "type": "horizontal-stacked-bar",
              "uid": "dEHMyN28RFSG1f-cPS6V2w",
              "tilesetUid": "HMSJyvLCSgGmrDJctdIz3w",
              "server": "https://resgen.io/api/v1",
              "options": {
                "labelPosition": "topLeft",
                "labelColor": "black",
                "labelTextOpacity": 0.4,
                "valueScaling": "exponential",
                "trackBorderWidth": 0,
                "trackBorderColor": "black",
                "backgroundColor": "white",
                "barBorder": false,
                "scaledHeight": false,
                "sortLargestOnTop": true,
                "colorScale": [
                  "#FF0000",
                  "#FF4500",
                  "#32CD32",
                  "#008000",
                  "#006400",
                  "#C2E105",
                  "#FFFF00",
                  "#66CDAA",
                  "#8A91D0",
                  "#CD5C5C",
                  "#E9967A",
                  "#BDB76B",
                  "#808080",
                  "#C0C0C0",
                  "#FFFFFF"
                ],
                "name": "all.KL.bed.hg38.multires.mv5"
              },
              "width": 1308,
              "height": 116,
              "resolutions": [
                13107200,
                6553600,
                3276800,
                1638400,
                819200,
                409600,
                204800,
                102400,
                51200,
                25600,
                12800,
                6400,
                3200,
                1600,
                800,
                400,
                200
              ]
            }
          ],
          "left": [],
          "center": [
            {
              "type": "horizontal-multivec",
              "uid": "cistrome-track",
              "tilesetUid": "UvVPeLHuRDiYA3qwFlm7xQ",
              "server": "https://resgen.io/api/v1",
              "options": {
                "labelPosition": "hidden",
                "labelColor": "black",
                "labelTextOpacity": 0.4,
                "valueScaling": "linear",
                "trackBorderWidth": 0,
                "trackBorderColor": "black",
                "heatmapValueScaling": "log",
                "name": "my_file_genome_wide.multires.mv5",
                "labelLeftMargin": 0,
                "labelRightMargin": 0,
                "labelTopMargin": 0,
                "labelBottomMargin": 0,
                "labelShowResolution": true,
                "minHeight": 100
              },
              "width": 1607,
              "height": 382,
              "resolutions": [
                16384000,
                8192000,
                4096000,
                2048000,
                1024000,
                512000,
                256000,
                128000,
                64000,
                32000,
                16000,
                8000,
                4000,
                2000,
                1000
              ]
            },
            /*{
                "type": "cistrome-group-labels",
                "uid": "cistrome-group-labels-track",
                "tilesetUid": "UvVPeLHuRDiYA3qwFlm7xQ",
                "server": "https://resgen.io/api/v1",
                "options": {
                  "labelPosition": "outerRight",
                },
                "width": 1607,
                "height": 382,
            }*/
          ],
          "bottom": [],
          "right": [],
          "whole": [],
          "gallery": []
        },
        "layout": {
          "w": 10,
          "h": 2,
          "x": 0,
          "y": 0,
          "moved": false,
          "static": false
        },
        "uid": "UiHlCoxRQ-aITBDi5j8b_w"
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

const hgOptions = {
    bounded: true,
    pixelPreciseMarginPadding: true,
    containerPaddingX: 0,
    containerPaddingY: 0,
    sizeMode: 'default'
};


/**
 * @component Cistrome HiGlass Wrapper 
 */
export default function CistromeHGW(props) {

    const hgRef = useRef();
    
    useEffect(() => {
        hgRef.current.api.on('location', (d) => {
            PubSub.publish(GLOBAL_X_RANGE, d.xRange);
            PubSub.publish(GLOBAL_Y_RANGE, d.yRange);
        });

        hgRef.current.api.on('viewConfig', (vc) => {
            //console.log(JSON.parse(vc));
            try {
                const trackObj = hgRef.current.api.getTrackObject("UiHlCoxRQ-aITBDi5j8b_w", "cistrome-track");
                PubSub.publish(TRACK_ROW_INFO, trackObj.tilesetInfo.row_infos);
                PubSub.publish(TRACK_POSITION, trackObj.position);
                PubSub.publish(TRACK_DIMENSIONS, trackObj.dimensions);
            } catch(e) {
    
            }
        });
    });

    console.log("CistromeHGW.render");

    return (
        <div className="cistrome-hgw">
            <HiGlassComponent
                viewConfig={demoViewConfig}
                options={hgOptions}
                zoomFixed={false}
                ref={hgRef}
            />
            <CistromeGroupLabels />
        </div>
    );
}