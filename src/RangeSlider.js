import React, { useRef, useState, useEffect, useCallback, useMemo } from "react";
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
    
    const [min, max] = valueExtent;
    const [minCutoff, setMinCutoff] = useState(min);
    const [maxCutoff, setMaxCutoff] = useState(max);

    const sliderWidth = 150;
    const moverSize = 14;
    const xScale = d3.scaleLinear()
        .domain(valueExtent)
        .range([0, sliderWidth]);

    // TOOD: Add Dragging interactions.

    function getCorrectedNumberInRange(value, [min, max], alt) {
        // Users can write whatever they want, but we need to correct it!
        return Number.isNaN(+value) || +value < min || +value > max ? alt : +value;
    }
    
    function onMinChange(e) {
        // Internal curoff value should be corrected to highlight well.
        const newValue = e.target.value;
        const corrected = getCorrectedNumberInRange(newValue, [min, maxCutoff], min);
        onHighlight([corrected, maxCutoff]);
        setMinCutoff(corrected);
    }

    function onMaxChange(e) {
        // Internal curoff value should be corrected to highlight well.
        const newValue = e.target.value;
        const corrected = getCorrectedNumberInRange(newValue, [minCutoff, max], max);
        onHighlight([minCutoff, corrected]);
        setMaxCutoff(corrected);
    }

    function onHighlight(range) {
        onChange(range);
    }

    function onMinBlur(e) {
        minInputRef.current.value = minCutoff;
    }

    function onMaxBlur(e) {
        maxInputRef.current.value = maxCutoff;
    }

    function onKeyDown(e) {
        switch(e.key){
            case 'Enter':
                minInputRef.current.value = minCutoff;
                maxInputRef.current.value = maxCutoff;
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
                onBlur={onMinBlur}
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
                        left: `${xScale(minCutoff) - moverSize / 2.0}px`,
                        width: `${xScale(maxCutoff) - xScale(minCutoff)}px`
                    }}
                />
                <div 
                    className="chw-range-slider-mover-left"
                    style={{
                        left: `${xScale(minCutoff) - moverSize / 2.0}px`,
                        width: `${moverSize}px`,
                        height: `${moverSize}px`,
                    }}
                />
                <div 
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
                onBlur={onMaxBlur}
                style={{
                    height,
                    textAlign: "left"
                }}
            />
        </div>
    );
}