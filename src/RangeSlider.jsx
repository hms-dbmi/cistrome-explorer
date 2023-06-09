import React, { useRef, useState, useEffect, useCallback } from "react";
import d3 from "./utils/d3.js";

import "./RangeSlider.scss";

/**
 * React component to interactively determine the left and right cutoff values.
 * @prop {boolean} isRight Is this view on the right side of the HiGlass track?
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
        isRight,
        height,
        valueExtent,
        onChange,
        onFilter,
        onClose
    } = props;
    
    const minInputRef = useRef();
    const maxInputRef = useRef();
    const leftMoverRef = useRef();
    const rightMoverRef = useRef();
    
    const [leftMost, rightMost] = isRight ? valueExtent : Array.from(valueExtent).reverse();
    const [leftCutoff, setLeftCutoff] = useState(leftMost);
    const [rightCutoff, setRightCutoff] = useState(rightMost);
    const selectedMover = useRef(null); // either "left" or "right"
    const dragStartX = useRef(null);
    const highlitRange = useRef([leftMost, rightMost]);

    const sliderWidth = 150;
    const moverSize = 16;
    const xScale = d3.scaleLinear()
        .domain([leftMost, rightMost])
        .range([0, sliderWidth]);

    useEffect(() => {
        // Update cutoff values when they are changed.
        setLeftCutoff(leftMost);
        setRightCutoff(rightMost);
        displayLeftCutoff(leftMost);
        displayRightCutoff(rightMost);
    }, [leftMost, rightMost]);

    function getCorrectedNumberInRange(valStr, range, alt) {
        const [left, right] = (range[0] < range[1]) ? range : Array.from(range).reverse();
        // Users can write whatever they want, but we need to correct it!
        if(Number.isNaN(+valStr)) {
            return alt;
        } else if(left <= +valStr && +valStr <= right) {
            return +valStr;
        } else if(+valStr < left) {
            return left;
        } else if(+valStr > right) {
            return right;
        }
    }
    
    function onLeftChange(e) {
        // Internal curoff value should be corrected to highlight well.
        const newValue = e.target.value;
        const corrected = getCorrectedNumberInRange(newValue, [leftMost, rightCutoff], leftMost);

        highlitRange.current = [corrected, rightCutoff];
        onHighlight(highlitRange.current);
        setLeftCutoff(corrected);
    }

    function onRightChange(e) {
        // Internal curoff value should be corrected to highlight well.
        const newValue = e.target.value;
        const corrected = getCorrectedNumberInRange(newValue, [leftCutoff, rightMost], rightMost);

        highlitRange.current = [leftCutoff, corrected];
        onHighlight(highlitRange.current);
        setRightCutoff(corrected);
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

        if(selectedMover.current === "left") {
            let newX = xScale(leftCutoff) + diffX;
            newX = getCorrectedNumberInRange(newX, [0, xScale(rightCutoff) - minMoverGap], 0);
            const newCutoff = xScale.invert(newX);
            
            highlitRange.current = [newCutoff, rightCutoff];
            setLeftCutoff(newCutoff);
            displayLeftCutoff(newCutoff);
        } else {
            let newX = xScale(rightCutoff) + diffX;
            newX = getCorrectedNumberInRange(newX, [xScale(leftCutoff) + minMoverGap, sliderWidth], sliderWidth);
            const newCutoff = xScale.invert(newX);
            
            highlitRange.current = [leftCutoff, newCutoff];
            setRightCutoff(newCutoff);
            displayRightCutoff(newCutoff);
        }
    });

    // Detect drag events for the slider elements.
    useEffect(() => {
        const leftMover = leftMoverRef.current;
        const rightMover = rightMoverRef.current;

        const drag = d3.drag()
            .on("start", started)
            .on("drag", dragged)
            .on("end", ended);

        d3.select(leftMover).call(drag);
        d3.select(rightMover).call(drag);

        return () => {
            d3.select(leftMover).on(".drag", null);
            d3.select(rightMover).on(".drag", null);
        };
    }, [leftMoverRef, rightMoverRef, started, dragged, ended]);

    useEffect(() => {
        const leftMover = leftMoverRef.current;
        const rightMover = rightMoverRef.current;

        d3.select(leftMover).on("mouseenter", () => { selectedMover.current = "left"; });
        d3.select(rightMover).on("mouseenter", () => { selectedMover.current = "right"; });

        return () => {
            d3.select(leftMover).on("mouseenter", null);
            d3.select(rightMover).on("mouseenter", null);
        };
    }, [leftMoverRef, rightMoverRef]);

    function onHighlight(range) {
        onChange(range);
    }

    function displayLeftCutoff(value) {
        minInputRef.current.value = value.toFixed(1);
    }

    function displayRightCutoff(value) {
        maxInputRef.current.value = value.toFixed(1);
    }

    function onKeyDown(e) {
        switch(e.key){
        case "Enter":
            displayLeftCutoff(leftCutoff);
            displayRightCutoff(rightCutoff);
            onFilter([leftCutoff, rightCutoff]);
            break;
        case "Esc":
        case "Escape":
            onClose(); 
            break;
        }
    }

    return (
        <div 
            className="hm-range-slider-container"
        >
            <input
                ref={minInputRef}
                className="hm-range-input"
                type="text"
                name="leftcutoff"
                placeholder={isRight ? "min" : "max"}
                defaultValue={leftMost}
                onChange={onLeftChange}
                onKeyDown={onKeyDown}
                onBlur={() => displayLeftCutoff(leftCutoff)}
                style={{
                    height,
                    textAlign: "right"
                }}
            />
            <div
                className="hm-range-slider"
                style={{
                    width: `${sliderWidth}px`
                }}
            >
                <div className="hm-range-slider-ruler-bg"/>
                <div 
                    className="hm-range-slider-ruler-fg"
                    style={{
                        left: `${xScale(leftCutoff)}px`,
                        width: `${xScale(rightCutoff) - xScale(leftCutoff)}px`
                    }}
                />
                <div 
                    ref={leftMoverRef}
                    className="hm-range-slider-mover-left"
                    style={{
                        left: `${xScale(leftCutoff) - moverSize / 2.0}px`,
                        width: `${moverSize}px`,
                        height: `${moverSize}px`,
                    }}
                />
                <div 
                    ref={rightMoverRef}
                    className="hm-range-slider-mover-right"
                    style={{
                        left: `${xScale(rightCutoff) - moverSize / 2.0}px`,
                        width: `${moverSize}px`,
                        height: `${moverSize}px`,
                    }}
                />
            </div>
            <input
                ref={maxInputRef}
                className="hm-range-input"
                type="text"
                name="rightcutoff"
                placeholder={isRight ? "max" : "min"}
                defaultValue={rightMost}
                onChange={onRightChange}
                onKeyDown={onKeyDown}
                onBlur={() => displayRightCutoff(rightCutoff)}
                style={{
                    height,
                    textAlign: "left"
                }}
            />
        </div>
    );
}