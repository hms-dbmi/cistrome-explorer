import React, { useState, useEffect, useRef, useCallback } from 'react';
import range from 'lodash/range';
import d3 from './utils/d3.js';


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
        d3.select(domElement).selectAll("g").remove();
        const g = d3.select(domElement).append("g").attr("class", "brush");

        const brush = d3.brushY()
            .extent([[0, 0], [width, height]]);

        function brushed() {
            if (d3.event.sourceEvent.type === "brush") return;
            const selection = d3.event.selection;
            if(selection) {
                const i = selection.map(yScale.invert);
                d3.select(this).call(brush.move, i.map(yScale));
            }
        };


        const dest = [0, height];
        
        g.call(brush);
        brush.on('brush', null);
        g.call(brush.move, dest);
        brush.on("brush", brushed);

        g.selectAll('.overlay')
            .style('pointer-events', 'none');
        
    }, [width, height]);

    register("TrackRowHighlight", draw);

    useEffect(() => {
        const svg = svgRef.current;
        const teardown = draw(svg);
        return teardown;
    });

    
    return (
        <svg
            ref={svgRef}
            style={{
                position: 'absolute',
                top: `${top}px`,
                left: `${left}px`,
                width: `${width}px`,
                height: `${height}px`,
                pointerEvents: 'none'
            }}
        />
    );
}