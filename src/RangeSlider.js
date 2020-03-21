import React, { useRef, useState, useEffect, useCallback, useMemo } from "react";

import "./RangeSlider.scss";

// TODO:
export default function RangeSlider(props) {
    const {
        height,
        valueExtent,
        onFilter,
        onClose
    } = props;
    
    const minInputRef = useRef();
    const maxInputRef = useRef();
    
    const [min, max] = valueExtent;
    const [minCutoff, setMinCutoff] = useState(min);
    const [maxCutoff, setMaxCutoff] = useState(max);

    function getCorrectedNumberInRange(value, [min, max], alt) {
        // Users can write whatever they want, but we need to correct it!
        return isNaN(+value) || +value < min || +value > max ? alt : +value;
    }

    function onMinChange(e) {
        // Internal curoff value should be corrected to highlight well.
        const newValue = e.target.value;
        const corrected = getCorrectedNumberInRange(newValue, [min, maxCutoff], min);
        setMinCutoff(corrected);
    }

    function onMaxChange(e) {
        // Internal curoff value should be corrected to highlight well.
        const newValue = e.target.value;
        const corrected = getCorrectedNumberInRange(newValue, [minCutoff, max], max);
        setMaxCutoff(corrected);
    }

    function onHighlight() {
        
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
                onFilter(minCutoff, maxCutoff);
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
                    width: 30,
                    height,
                    textAlign: "right"
                }}
            />
            <div
                className={"chw-range-slider"}
            >

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
                    width: 30, 
                    height,
                    textAlign: "left"
                }}
            />
        </div>
    );
}