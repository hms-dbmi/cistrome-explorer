import React, { useState, useEffect } from 'react';
import PubSub from "pubsub-js";

import { resolveIntervalCoordinates } from './utils/genome.js';
import { validateIntervalParams } from './utils/cistrome.js';
import { TooltipContent, destroyTooltip } from "./Tooltip.js";
import { getRange } from './utils/viewport.js';
import { CLOSE, SEARCH } from './utils/icons.js';
import { EVENT } from "./utils/constants.js";
import './ViewColumnBrush.scss';

/**
 * Component for rendering genome interval selection tools.
 * @prop {object} viewBoundingBox The bounding box (i.e., left, top, width, height) of a HiGlass view.
 * @prop {object} viewportTrack An object of `viewport-projection-horizontal` track.
 * @prop {object} multivecTrack A object of `horizontal-multivec` track in the same view.
 * @prop {function} onViewportRemove The function to call upon removing a viewport track.
 * @prop {function} onRequestIntervalTFs The function to call upon making a request for further interval transcription factor data.
 */
// TODO: documents
export default function ViewColumnBrush(props) {
    
    const {
        viewBoundingBox,
        viewportTrack,
        multivecTrack,
        onViewportRemove,
        onRequestIntervalTFs
    } = props;

    const [chrStartName, setChrStartName] = useState(null);
    const [chrStartPos, setChrStartPos] = useState(null);
    const [chrEndName, setChrEndName] = useState(null);
    const [chrEndPos, setChrEndPos] = useState(null);

    const assembly = multivecTrack.tilesetInfo.coordSystem;
    const absDomain = viewportTrack.viewportXDomain;
    let startX = viewportTrack._xScale(absDomain[0]);
    let endX = viewportTrack._xScale(absDomain[1]);

    useEffect(() => {
        if(!assembly) {
            return null;
        }

        let didUnmount = false;
        const absDomain = viewportTrack.viewportXDomain;
        resolveIntervalCoordinates(assembly, absDomain[0], absDomain[1])
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

    const { start, end } = getRange(startX, endX, 0, viewBoundingBox.width);
    const { 
        msg: intervalInvalidMsg, 
        success: intervalValid 
    } = validateIntervalParams(assembly, chrStartName, chrStartPos, chrEndName, chrEndPos);

    if(end === null && start === null) {
        // Do not show range when they are ourside.
        return null;
    }

    return (
        <div className="col-tools-brush"
            style={{
                left: start,
                width: end - start
            }}>
            <div className={"chw-button-sm-container-horizontal"}
                style={{
                    right: "4px",
                    top: "2px",
                    opacity: 1,
                    background: "none",
                    boxShadow: "none",
                    color: "gray"
                }}>
                {/* Cistrome DB API Button */}
                <svg className={"chw-button-sm-black"}
                    onMouseOver={(e) => {
                        PubSub.publish(EVENT.TOOLTIP, {
                            x: e.clientX,
                            y: e.clientY,
                            content: <TooltipContent 
                                title="Search bind TFs from Cistrome DB"
                                value={!intervalValid ? intervalInvalidMsg : undefined}
                                warning={!intervalValid}
                            />
                        });
                    }}
                    onMouseLeave={() => destroyTooltip()}
                    onClick={() => {
                        if(intervalValid) {
                            onRequestIntervalTFs({
                                assembly,
                                chrStartName,
                                chrStartPos,
                                chrEndName,
                                chrEndPos
                            });
                        }
                    }}
                    viewBox={SEARCH.viewBox}>
                    <path d={SEARCH.path} fill="currentColor"/>
                </svg>
                {/* Close Button */}
                <svg className={"chw-button-sm-black"}
                    onClick={() => onViewportRemove(viewportTrack.id)} 
                    viewBox={CLOSE.viewBox}>
                    <title>Remove viewport projection track</title>
                    <path d={CLOSE.path} fill="currentColor"/>
                </svg>
            </div>
        </div>
    );
}