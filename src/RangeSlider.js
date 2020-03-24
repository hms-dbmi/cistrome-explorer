import React, { useRef, useState, useEffect, useCallback } from "react";
import d3 from './utils/d3.js';

import "./RangeSlider.scss";

/**
 * React component to interactively determine the left and right cutoff values.
 * @prop {number} height The size along y-axis of this component.
 * @prop {array} valueExtent The array that have two numbers, indicating the min and max values.
 * @prop {function} onChange The function to call when the search keyword has changed.
 * @prop {function} onFilter The function to call when the filter should be applied.
 * @prop {function} onClose The function to call when the search field should be closed.
 * @example
 * <RangeSlider/>
 */
export default function RangeSlider(props) {
    const {
        height,
        valueExtent,
        onChange,
        onFilter,
        onClose
    } = props;
    
    const minInputRef = useRef();
    const maxInputRef = useRef();
    const minMoverRef = useRef();
    const maxMoverRef = useRef();
    
    const [min, max] = valueExtent;
    const [minCutoff, setMinCutoff] = useState(min);
    const [maxCutoff, setMaxCutoff] = useState(max);
    const selectedMover = useRef(null); // either "min" or "max"
    const dragStartX = useRef(null);
    const highlitRange = useRef([min, max]);

    const sliderWidth = 150;
    const moverSize = 16;
    const xScale = d3.scaleLinear()
        .domain(valueExtent)
        .range([0, sliderWidth]);

    useEffect(() => {
        // Update cutoff values when they are changed.
        setMinCutoff(min);
        setMaxCutoff(max);
        displayMinCutoff(min);
        displayMaxCutoff(max);
    }, [min, max]);

    function getCorrectedNumberInRange(valStr, [min, max], alt) {
        // Users can write whatever they want, but we need to correct it!
        if(Number.isNaN(+valStr)) {
            return alt;
        } else if(+valStr < min) {
            return min;
        } else if(+valStr > max) {
            return max;
        } else {
            return +valStr;
        }
    }
    
    function onMinChange(e) {
        // Internal curoff value should be corrected to highlight well.
        const newValue = e.target.value;
        const corrected = getCorrectedNumberInRange(newValue, [min, maxCutoff], min);

        highlitRange.current = [corrected, maxCutoff];
        onHighlight(highlitRange).current;
        setMinCutoff(corrected);
    }

    function onMaxChange(e) {
        // Internal curoff value should be corrected to highlight well.
        const newValue = e.target.value;
        const corrected = getCorrectedNumberInRange(newValue, [minCutoff, max], max);

        highlitRange.current = [minCutoff, corrected];
        onHighlight(highlitRange.current);
        setMaxCutoff(corrected);
    }
    
    // Set up the d3-drag handler functions (started, ended, dragged).
    const started = useCallback(() => {
        const event = d3.event;
        dragStartX.current = event.sourceEvent.clientX;
    });

    const ended = useCallback(() => {
        onHighlight(highlitRange.current);
        dragStartX.current = null;
    });
    
    // We do not want to position two movers in the exact same position.
    const minMoverGap = moverSize + 2;

    const dragged = useCallback(() => {
        const event = d3.event;
        const diffX = event.sourceEvent.clientX - dragStartX.current;

        if(selectedMover.current === "min") {
            let newX = xScale(minCutoff) + diffX;
            newX = getCorrectedNumberInRange(newX, [0, xScale(maxCutoff) - minMoverGap], 0);
            const newCutoff = xScale.invert(newX);
            
            highlitRange.current = [newCutoff, maxCutoff];
            setMinCutoff(newCutoff);
            displayMinCutoff(newCutoff);
        } else {
            let newX = xScale(maxCutoff) + diffX;
            newX = getCorrectedNumberInRange(newX, [xScale(minCutoff) + minMoverGap, sliderWidth], sliderWidth);
            const newCutoff = xScale.invert(newX);
            
            highlitRange.current = [minCutoff, newCutoff];
            setMaxCutoff(newCutoff);
            displayMaxCutoff(newCutoff);
        }
    });

    // Detect drag events for the slider elements.
    useEffect(() => {
        const minMover = minMoverRef.current;
        const maxMover = maxMoverRef.current;

        const drag = d3.drag()
            .on("start", started)
            .on("drag", dragged)
            .on("end", ended);

        d3.select(minMover).call(drag);
        d3.select(maxMover).call(drag);

        return () => {
            d3.select(minMover).on(".drag", null);
            d3.select(maxMover).on(".drag", null);
        };
    }, [minMoverRef, maxMoverRef, started, dragged, ended]);

    useEffect(() => {
        const minMover = minMoverRef.current;
        const maxMover = maxMoverRef.current;

        d3.select(minMover).on("mouseenter", () => { selectedMover.current = "min" });
        d3.select(maxMover).on("mouseenter", () => { selectedMover.current = "max" });

        return () => {
            d3.select(minMover).on("mouseenter", null);
            d3.select(maxMover).on("mouseenter", null);
        }
    }, [minMoverRef, maxMoverRef]);

    function onHighlight(range) {
        onChange(range);
    }

    function displayMinCutoff(value) {
        minInputRef.current.value = value.toFixed(1);
    }

    function displayMaxCutoff(value) {
        maxInputRef.current.value = value.toFixed(1);
    }

    function onKeyDown(e) {
        switch(e.key){
            case 'Enter':
                displayMinCutoff(minCutoff);
                displayMaxCutoff(maxCutoff);
                onFilter([minCutoff, maxCutoff]);
                break;
            case 'Esc':
            case 'Escape':
                onClose(); 
                break;
        }
    }

    return (
        <div 
            className="chw-range-slider-container"
        >
            <input
                ref={minInputRef}
                className="chw-range-input"
                type="text"
                name="min-name"
                placeholder="min"
                defaultValue={min}
                onChange={onMinChange}
                onKeyDown={onKeyDown}
                onBlur={() => displayMinCutoff(minCutoff)}
                style={{
                    height,
                    textAlign: "right"
                }}
            />
            <div
                className="chw-range-slider"
                style={{
                    width: `${sliderWidth}px`
                }}
            >
                <div className="chw-range-slider-ruler-bg"/>
                <div 
                    className="chw-range-slider-ruler-fg"
                    style={{
                        left: `${xScale(minCutoff)}px`,
                        width: `${xScale(maxCutoff) - xScale(minCutoff)}px`
                    }}
                />
                <div 
                    ref={minMoverRef}
                    className="chw-range-slider-mover-left"
                    style={{
                        left: `${xScale(minCutoff) - moverSize / 2.0}px`,
                        width: `${moverSize}px`,
                        height: `${moverSize}px`,
                    }}
                />
                <div 
                    ref={maxMoverRef}
                    className="chw-range-slider-mover-right"
                    style={{
                        left: `${xScale(maxCutoff) - moverSize / 2.0}px`,
                        width: `${moverSize}px`,
                        height: `${moverSize}px`,
                    }}
                />
            </div>
            <input
                ref={maxInputRef}
                className="chw-range-input"
                type="text"
                name="max-input"
                placeholder="max"
                defaultValue={max}
                onChange={onMaxChange}
                onKeyDown={onKeyDown}
                onBlur={() => displayMaxCutoff(maxCutoff)}
                style={{
                    height,
                    textAlign: "left"
                }}
            />
        </div>
    );
}