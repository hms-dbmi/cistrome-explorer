import React, { useRef, useState, useEffect, useCallback } from "react";
import d3 from './utils/d3.js';

import "./RangeSlider.scss";

// TODO:
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
    const selectedMover = useRef(); // either "left" or "right".
    
    const [min, max] = valueExtent;
    const [minCutoff, setMinCutoff] = useState(min);
    const [maxCutoff, setMaxCutoff] = useState(max);
    const dragStartX = useRef(null);
    let highlitRange = [min, max];

    const sliderWidth = 150;
    const moverSize = 16;
    const xScale = d3.scaleLinear()
        .domain(valueExtent)
        .range([0, sliderWidth]);

    // Set up the d3-drag handler functions (started, ended, dragged).
    const started = useCallback(() => {
        const event = d3.event;
        dragStartX.current = event.sourceEvent.clientX;
    });

    const ended = useCallback(() => {
        dragStartX.current = null;
        onHighlight(highlitRange);
    });
    
    // We do not want to position two movers in the exact same position.
    // TODO: intergrate dragged and onMinChange?
    const minMoverDist = moverSize + 2;
    const dragged = useCallback(() => {
        const event = d3.event;
        const diffX = event.sourceEvent.clientX - dragStartX.current;
        if(selectedMover.current === "left") {
            let newX = xScale(minCutoff) + diffX;
            if(newX < 0) {
                newX = 0;
            } else if(newX > xScale(maxCutoff) - minMoverDist) {
                newX = xScale(maxCutoff) - minMoverDist;
            }
            const newCutoff = xScale.invert(newX);
            minInputRef.current.value = newCutoff.toFixed(1);

            setMinCutoff(newCutoff);
            highlitRange = [newCutoff, maxCutoff];
        } else {
            let newX = xScale(maxCutoff) + diffX;
            if(newX < xScale(minCutoff) + minMoverDist) {
                newX = xScale(minCutoff) + minMoverDist;
            } else if(newX > sliderWidth) {
                newX = sliderWidth;
            }
            const newCutoff = xScale.invert(newX);
            maxInputRef.current.value = newCutoff.toFixed(1);

            setMaxCutoff(newCutoff);
            highlitRange = [minCutoff, newCutoff];
        }
    });

    // Detect drag events for the resize element.
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

        d3.select(minMover).on("mouseenter", () => {selectedMover.current = "left"});
        d3.select(maxMover).on("mouseenter", () => {selectedMover.current = "right"});

        return () => {
            d3.select(minMover).on("mouseenter", null);
            d3.select(maxMover).on("mouseenter", null);
        }
    }, [minMoverRef, maxMoverRef]);

    function getCorrectedNumberInRange(value, [min, max], alt) {
        // Users can write whatever they want, but we need to correct it!
        return Number.isNaN(+value) || +value < min || +value > max ? alt : +value;
    }
    
    function onMinChange(e) {
        // Internal curoff value should be corrected to highlight well.
        const newValue = e.target.value;
        const corrected = getCorrectedNumberInRange(newValue, [min, maxCutoff], min);

        highlitRange = [corrected, maxCutoff];
        onHighlight(highlitRange);
        setMinCutoff(corrected);
    }

    function onMaxChange(e) {
        // Internal curoff value should be corrected to highlight well.
        const newValue = e.target.value;
        const corrected = getCorrectedNumberInRange(newValue, [minCutoff, max], max);

        highlitRange = [minCutoff, corrected];
        onHighlight(highlitRange);
        setMaxCutoff(corrected);
    }

    function onHighlight(range) {
        onChange(range);
    }

    function displayMinCutoff(e) {
        minInputRef.current.value = minCutoff.toFixed(1);
    }

    function displayMaxCutoff(e) {
        maxInputRef.current.value = maxCutoff.toFixed(1);
    }

    function onKeyDown(e) {
        switch(e.key){
            case 'Enter':
                displayMinCutoff();
                displayMaxCutoff();
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
                onBlur={displayMinCutoff}
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
                onBlur={displayMaxCutoff}
                style={{
                    height,
                    textAlign: "left"
                }}
            />
        </div>
    );
}