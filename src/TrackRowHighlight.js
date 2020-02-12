import React, { useState, useEffect, useRef, useCallback } from 'react';
import range from 'lodash/range';
import d3 from './utils/d3.js';
import Two from './utils/two.js';

import PubSub from 'pubsub-js';
import { EVENT } from './constants.js';
import { visualizationTrack } from './visualizationTrack.js';
import TrackControl from './TrackControl.js'
import './TrackRowInfo.scss';


/**
 * Component for visualization of row info attribute values.
 * @prop {number} trackX The track horizontal offset.
 * @prop {number} trackY The track vertical offset.
 * @prop {number} trackWidth The track width.
 * @prop {number} trackHeight The track height.
 * @prop {number} highlitRows Array of highlighted row indices.
 * @prop {function} register The function for child components to call to register their draw functions.
 */
export default function TrackRowHighlight(props) {

    const {
        trackX, trackY,
        trackWidth, trackHeight, 
        totalNumRows,
        highlitRows,
        register
    } = props;

   
    const top = trackY;
    const left = trackX;
    const height = trackHeight;
    const width = trackWidth;

    // Render svg
    const svgRef = useRef();

    
    const yScale = d3.scaleBand()
        .domain(range(totalNumRows))
        .range([0, height]);

    // Render each track.
    const draw = useCallback((domElement) => {
        const two = new Two({
            width,
            height,
            domElement
        });

        const rect = two.makeRect(0, 0, width, height);
        rect.fill = "blue";
        rect.opacity = 0.7;
    
        
        two.update();
        return two.teardown;
    }, [width, height]);

    register("TrackRowHighlight", draw);

    useEffect(() => {
        const svg = svgRef.current;
        const teardown = draw(svg);
       
        return teardown;
    });

    
    return (
        <div className="cistrome-hgw-child"
            style={{
                top: `${top}px`,
                left: `${left}px`, 
                width: `${width}px`,
                height: `${height}px`,
            }}
        >
            <svg
                ref={svgRef}
                style={{
                    position: 'relative',
                    top: 0,
                    left: 0, 
                    width: `${width}px`,
                    height: `${height}px`,
                    pointerEvents: 'none'
                }}
            />
        </div>
    );
}