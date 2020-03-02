import React, { useRef, useEffect, useState, useCallback } from 'react';
import d3 from './utils/d3.js';
import Two from './utils/two.js';

import { resolveIntervalCoordinates } from './utils/genome.js';
import { CISTROME_DBTOOLKIT_MAX_INTERVAL_SIZE } from './utils/constants.js';
import { SEARCH } from './utils/icons.js';

import './TrackColSelectionInfo.scss';
import './TrackRowInfoControl.scss';

function makeDBToolkitURL(assembly, chrStartName, chrStartPos, chrEndName, chrEndPos) {
    if(!chrStartName || !chrStartPos || !chrEndName || !chrEndPos) {
        return null;
    }
    // Generate a URL for the cistrome DB toolkit site.
    if(chrStartName !== chrEndName) {
        // Bail out, interval spans across more than one chromosome.
        return null;
    }
    if(chrEndPos - chrStartPos > CISTROME_DBTOOLKIT_MAX_INTERVAL_SIZE) {
        // Bail out, interval is too large for dbtoolkit's interval search.
        return null;
    }
    return `http://dbtoolkit.cistrome.org/?specie=${assembly}&factor=tf&interval=${chrStartName}%3A${chrStartPos}-${chrEndPos}`;
}

function makeDBToolkitAPIURL(assembly, chrStartName, chrStartPos, chrEndName, chrEndPos) {
    if(!chrStartName || !chrStartPos || !chrEndName || !chrEndPos) {
        return null;
    }
    if(chrStartName !== chrEndName) {
        // Bail out, interval spans across more than one chromosome.
        return null;
    }
    if(chrEndPos - chrStartPos > CISTROME_DBTOOLKIT_MAX_INTERVAL_SIZE) {
        // Bail out, interval is too large for dbtoolkit's interval search.
        return null;
    }
    // Generate a URL for the cistrome DB toolkit API.
    return `http://dbtoolkit.cistrome.org/api_interval?species=${assembly}&factor=$tf&interval=${chrStartName}%3A${chrStartPos}-${chrEndPos}`;
}

const numberFormatter = d3.format(",");

/**
 * Component for rendering information about a particular genomic interval selection.
 * @prop {number} width
 * @prop {number} height
 * @prop {object} projectionTrack A `viewport-projection-horizontal` track object.
 * @prop {(string|null)} trackAssembly The genome assembly/coordSystem value obtained from the associated `horizontal-multivec` track.
 * @prop {string} colToolsPosition The value of the colToolsPosition option.
 * @prop {function} drawRegister The function for child components to call to register their draw functions.
 */
export default function TrackColSelectionInfo(props) {

    const {
        width,
        height,
        projectionTrack,
        trackAssembly,
        colToolsPosition,
        drawRegister
    } = props;
    
    const [chrStartName, setChrStartName] = useState(null);
    const [chrStartPos, setChrStartPos] = useState(null);
    const [chrEndName, setChrEndName] = useState(null);
    const [chrEndPos, setChrEndPos] = useState(null);

    useEffect(() => {
        if(!trackAssembly) {
            return null;
        }

        let didUnmount = false;
        const absDomain = projectionTrack.viewportXDomain;
        resolveIntervalCoordinates(trackAssembly, absDomain[0], absDomain[1])
        .then(result => {
            if(!didUnmount) {
                // Only update state if the component has not yet unmounted.
                // See https://github.com/facebook/react/issues/14369#issuecomment-468267798
                setChrStartName(result[0][0]);
                setChrStartPos(result[0][1]);
                setChrEndName(result[1][0]);
                setChrEndPos(result[1][1]);
            }
        });

        return (() => { didUnmount = true; });
    });

    function onRequestIntervalTFs() {

    }

    const canvasRef = useRef();

    const draw = useCallback((domElement) => {
        const two = new Two({
            width,
            height,
            domElement
        });

        if(!chrStartName || !chrStartPos || !chrEndName || !chrEndPos) {
            return two.teardown;
        }

        const absDomain = projectionTrack.viewportXDomain;
        const startX = projectionTrack._xScale(absDomain[0]);
        const endX = projectionTrack._xScale(absDomain[1]);

        const EDGE_LIMIT = 80;
        const MIDDLE_LIMIT = 130;

        const alignMiddle = (endX - startX >= MIDDLE_LIMIT);
        const textColor = "#555";
        const tickSize = 6;
        const chrFontSize = 14;
        
        const isTop = (colToolsPosition === "top");
        const textTop = (isTop ? height - (height/6 + tickSize) : height/6 + tickSize);
        const tickTop = (isTop ? height - tickSize : 0);
        
        // Draw text elements for genomic interval start and end positions.
        const startText = two.makeText(startX, textTop, width/2, 2*height/3, `${chrStartName}:${numberFormatter(chrStartPos)}`);
        startText.align = (alignMiddle ? (startX < EDGE_LIMIT ? "start" : "middle") : "end");
        startText.baseline = "middle";
        startText.fontsize = chrFontSize;
        startText.fill = textColor;
        
        const endText = two.makeText(endX, textTop, width/2, 2*height/3, `${chrEndName}:${numberFormatter(chrEndPos)}`);
        endText.align = (alignMiddle ? (width - endX < EDGE_LIMIT ? "end" : "middle") : "start");
        endText.baseline = "middle";
        endText.fontsize = chrFontSize;
        endText.fill = textColor;

        // Draw tick marks above the start and end positions.
        const startTick = two.makeLine(startX, tickTop, startX, tickTop + tickSize);
        startTick.stroke = textColor;
        const endTick = two.makeLine(endX, tickTop, endX, tickTop + tickSize);
        endTick.stroke = textColor;

        two.update();
        return two.teardown;
    }, [width, height, chrStartName, chrStartPos, chrEndName, chrEndPos]);

    drawRegister("TrackColSelectionInfo", draw);

    useEffect(() => {
        const canvas = canvasRef.current;
        const teardown = draw(canvas);
        return teardown;
    });

    const dbToolkitURL = makeDBToolkitURL(trackAssembly, chrStartName, chrStartPos, chrEndName, chrEndPos);
    const dbToolkitAPIURL = makeDBToolkitAPIURL(trackAssembly, chrStartName, chrStartPos, chrEndName, chrEndPos);

    return (
        <div
            style={{
                position: 'relative',
                top: 0,
                left: 0, 
                width: `${width}px`,
                height: `${height}px`
            }}
        >
            <canvas
                ref={canvasRef}
                style={{
                    position: 'absolute',
                    top: 0,
                    left: 0, 
                    width: `${width}px`,
                    height: `${height}px`
                }}
            />
            <div
                style={{
                    position: 'absolute',
                    top: `${height/3}px`,
                    left: 0, 
                    width: `${width}px`,
                    height: `${height/6}px`
                }}
            >
                {(!chrStartName || !chrStartPos || !chrEndName || !chrEndPos || !["hg38", "mm10"].includes(trackAssembly)) ? null : (
                    dbToolkitAPIURL ? (
                        <div className="chgw-button"
                            onClick={onRequestIntervalTFs}>
                            <svg className="chgw-button-sm chgw-search-button chgw-button-static"
                                viewBox={SEARCH.viewBox}>
                                <path d={SEARCH.path} fill="gray"/>
                            </svg>
                            Show Bind TFs in Cistrome DB
                        </div>
                    ) : (
                        <p className="col-selection-info-disabled">
                            Search requires interval &le; 2 Mb
                        </p>
                    )
                )}
            </div>
        </div>
    );
};