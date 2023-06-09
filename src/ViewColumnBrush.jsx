import React, { useState, useEffect } from "react";
import PubSub from "pubsub-js";

import { resolveIntervalCoordinates } from "./utils/genome.js";
import { validateIntervalParams } from "./utils/cistrome.js";
import { TooltipContent, destroyTooltip } from "./Tooltip.jsx";
import { getRange } from "./utils/viewport.js";
import { CLOSE, SEARCH } from "./utils/icons.js";
import { EVENT } from "./utils/constants.js";
import "./ViewColumnBrush.scss";

/**
 * Component for rendering genome interval selection tools.
 * @prop {object} viewBoundingBox The bounding box (i.e., left, top, width, height) of a HiGlass view.
 * @prop {object} viewportTrack An object of `viewport-projection-horizontal` track.
 * @prop {object} multivecTrack A object of `horizontal-multivec` track in the same view.
 * @prop {function} onViewportRemove The function to call upon removing a viewport track.
 * @prop {function} onRequestIntervalTFs The function to call upon making a request for further interval transcription factor data.
 */
export default function ViewColumnBrush(props) {
    
    const {
        viewBoundingBox,
        viewportTrack,
        multivecTrack,
        onViewportRemove,
        onRequestIntervalTFs
    } = props;

    const [isLoading, setIsLoading] = useState(true);
    const [positionInfo, setPositionInfo] = useState(null);

    const [onSearchButton, setOnSearchButton] = useState({
        on: false, x: undefined, y: undefined
    });

    const assembly = multivecTrack.tilesetInfo?.coordSystem || multivecTrack.options?.coordSystem || "hg38";
    const absDomain = viewportTrack.viewportXDomain;
    let startX = viewportTrack._xScale(absDomain[0]);
    let endX = viewportTrack._xScale(absDomain[1]);

    useEffect(() => {
        if(!assembly) {
            return () => {};
        }
        let didUnmount = false;
        const absDomain = viewportTrack.viewportXDomain;
        resolveIntervalCoordinates(assembly, absDomain[0], absDomain[1])
            .then(result => {
                if(!didUnmount) {
                    // Only update state if the component has not yet unmounted.
                    // See https://github.com/facebook/react/issues/14369#issuecomment-468267798
                    setPositionInfo({
                        chrStartName: result[0][0],
                        chrStartPos: result[0][1],
                        chrEndName: result[1][0],
                        chrEndPos: result[1][1],
                    });
                    setIsLoading(false);
                }
            });

        return (() => { didUnmount = true; });
    }, [viewBoundingBox, viewportTrack, multivecTrack, onViewportRemove, onRequestIntervalTFs]);

    const { start, end } = getRange(startX, endX, 0, viewBoundingBox.width);
    const { 
        msg: intervalInvalidMsg, 
        success: intervalValid 
    } = validateIntervalParams({ assembly, ...positionInfo, });

    // Force update the loading message if changed while hovered already
    useEffect(() => {
        if(onSearchButton.on && onSearchButton.x && onSearchButton.y) {
            const tooltipValue = (isLoading
                ? "Loading metadata..."
                : (!intervalValid
                    ? intervalInvalidMsg 
                    : `${positionInfo?.chrStartName}:${positionInfo?.chrStartPos}-${positionInfo?.chrEndPos}`
                )
            );

            PubSub.publish(EVENT.TOOLTIP, {
                x: onSearchButton.x,
                y: onSearchButton.y,
                content: <TooltipContent 
                    title="Search Transcription Factors from Cistrome DB"
                    value={tooltipValue}
                    warning={!intervalValid}
                />
            });
        }
    }, [onSearchButton, isLoading]);

    if(end === null && start === null) {
        // Do not show range when they are ourside.
        return null;
    }
    
    // Determines when to hide buttons based on the width of a viewport
    // Changed from 50 to zero to show it at all time.
    const minWidthToShowButton = 0; 

    return (
        <div className="col-tools-brush"
            style={{
                left: start,
                width: end - start
            }}>
            {end - start > minWidthToShowButton ? 
                <div className={"hm-button-sm-container-horizontal"}
                    style={{
                        right: "0px",
                        height: "100%",
                        opacity: 1,
                        background: "none",
                        boxShadow: "none",
                        color: "gray"
                    }}>
                    {/* Cistrome DB API Button */}
                    {onRequestIntervalTFs ? 
                        // Show the search icon only when `onIntervalSearch` is defined.
                        <svg className={"hm-button-sm hm-button-middle"} 
                            style={{ height: "100%" }}
                            onMouseOver={(e) => setOnSearchButton({
                                on: true,
                                x: e.clientX,
                                y: e.clientY
                            })}
                            onMouseLeave={() => {
                                setOnSearchButton({on: false});
                                destroyTooltip();
                            }}
                            onClick={() => {
                                if(intervalValid) {
                                    onRequestIntervalTFs({
                                        assembly,
                                        ...positionInfo
                                    });
                                }
                            }}
                            viewBox={SEARCH.viewBox}>
                            <path d={SEARCH.path} fill="currentColor"/>
                        </svg>
                        : null}
                    {/* Close Button */}
                    <svg className={"hm-button-sm hm-button-middle"}
                        style={{ height: "100%" }}
                        onClick={() => onViewportRemove(viewportTrack.id)} 
                        viewBox={CLOSE.viewBox}>
                        <title>Remove viewport projection track</title>
                        <path d={CLOSE.path} fill="currentColor"/>
                    </svg>
                </div>
                : null}
        </div>
    );
}
