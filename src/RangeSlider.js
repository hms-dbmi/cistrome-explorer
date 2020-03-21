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

    function getCorrectedValue(value, alt) {
        return isNaN(+value) ? alt : +value;
    }

    // TODO: Consider many cases, such as minCuroff > maxCuroff, minCuroff < min
    function onMinChange(e) {
        const newValue = e.target.value;
        const corrected = getCorrectedValue(newValue, min);
        setMinCutoff(corrected);
    }

    function onMaxChange(e) {
        const newValue = e.target.value;
        const corrected = getCorrectedValue(newValue, max);
        setMaxCutoff(corrected);
    }

    function onKeyDown(e) {
        switch(e.key){
            case 'Enter':
                minInputRef.current.value = getCorrectedValue(minInputRef.current.value, min);
                maxInputRef.current.value = getCorrectedValue(maxInputRef.current.value, max);
                onFilter();
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
                style={{ 
                    width: 30, 
                    height,
                    textAlign: "left"
                }}
            />
        </div>
    );
}