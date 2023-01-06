import React, { useRef, useState, useEffect, useCallback, useMemo } from "react";
import { CLOSE, FILTER, RESET, SEARCH, SQUARE_CHECK, SQUARE } from "./utils/icons.js";
import d3 from "./utils/d3.js";

import "./TrackRowFilter.scss";
import RangeSlider from "./RangeSlider.js";
import { getAggregatedValue } from "./utils/aggregate.js";

const MAX_NUM_SUGGESTIONS = 200;

/**
 * Component to help determining rows to filter.
 * @prop {boolean} isLeft Is this view on the left side of the HiGlass track?
 * @prop {number} top The top coordinate.
 * @prop {number} left The left coordinate.
 * @prop {function} onChange The function to call when the search keyword or value range has changed.
 * @prop {function} onFilterRows The function to call when the filter should be applied.
 * @prop {function} onClose The function to call when filter component should be closed.
 * @prop {object[]} rowInfo The array of JSON Object containing row information.
 * @prop {object} filterInfo The options for filtering rows of the field used in this track.
 * @example
 * <TrackRowFilter/>
 */
export default function TrackAddNewTrack(props) {
    const {
        top, left,
        onChange,
        onFilterRows,
        onClose,
        rowInfo,
        filterInfo
    } = props;
    
    const moverRef = useRef();
    const keywordInputRef = useRef();
    const [keyword, setKeyword] = useState("");
    const [suggestionIndex, setSuggestionIndex] = useState(undefined);
    const [offset, setOffset] = useState({x: 0, y: 0});
    const dragStartPos = useRef(null);
    const [valueExtent, setValueExtent] = useState(Object.keys(rowInfo[0]));
    const cutoffRange = useRef(valueExtent);
    const keywordUpperCase = keyword.toUpperCase();
    // const [notOneOf, setNotOneOf] = useState(!filterInfo || type === "quantitative" ? [] : filterInfo.notOneOf); 
    
    // Styles
    const width = 180;
    const height = 30;
    const maxHeight = 420;
    const padding = 5;
    
    useEffect(() => {
        setValueExtent(Object.keys(rowInfo[0]));
    }, [rowInfo]);

    useEffect(() => {
        keywordInputRef.current.focus();
    });

    useEffect(() => {
        // setNotOneOf(!filterInfo || type === "quantitative" ? [] : filterInfo.notOneOf);
    }, [filterInfo]);

    const suggestions = useMemo(() => {
        let result = [];
        const fieldData = Object.keys(rowInfo[0]);
        const fieldDataByKeyword = fieldData.filter(d => d.toUpperCase().includes(keywordUpperCase));
        const potentialResult = Array.from(new Set(fieldDataByKeyword));
        if(potentialResult.length < MAX_NUM_SUGGESTIONS) {
            result = potentialResult;
        } else {
            result = potentialResult.slice(0, MAX_NUM_SUGGESTIONS);
        }
        // Sort in alphabetical order first.
        result.sort((a, b) => {
            return a.toUpperCase() > b.toUpperCase();
        });
        // Sort so that suggestions that _start with_ the keyword appear first.
        result.sort((a, b) => {
            return a.toUpperCase().indexOf(keywordUpperCase) - b.toUpperCase().indexOf(keywordUpperCase);
        });
        return result;
    }, [keyword]);

    useEffect(() => {
        if(suggestionIndex === undefined && keyword.length > 0) {
            onChange(keyword);
        }
        else if(suggestionIndex !== undefined && suggestionIndex >= 0 && suggestionIndex < suggestions.length) {
            onChange(suggestions[suggestionIndex]);
        }
        else if(suggestionIndex === undefined) {
            onChange();
        }
    }, [suggestions, suggestionIndex]);

    // Set up the d3-drag handler functions (started, ended, dragged).
    const started = useCallback(() => {
        const event = d3.event;
        dragStartPos.current = {
            x: event.sourceEvent.clientX, 
            y: event.sourceEvent.clientY
        };
    });

    const ended = useCallback(() => {
        dragStartPos.current = null;
    });

    const dragged = useCallback(() => {
        const event = d3.event;
        const diffX = offset.x + event.sourceEvent.clientX - dragStartPos.current.x;
        const diffY = offset.y + event.sourceEvent.clientY - dragStartPos.current.y;
        setOffset({x: diffX, y: diffY});
    });

    // Detect drag events for the resize element.
    useEffect(() => {
        const mover = moverRef.current;

        const drag = d3.drag()
            .on("start", started)
            .on("drag", dragged)
            .on("end", ended);

        d3.select(mover).call(drag);

        return () => d3.select(mover).on(".drag", null);
    }, [moverRef, started, dragged, ended]);
    
    function onKeywordChange(e) {
        const newKeyword = e.target.value;
        setKeyword(newKeyword);
        onChange(newKeyword);
    }

    function onFilterClose() {
        onChange();
        onClose();
        setKeyword("");
        setSuggestionIndex(undefined);
    }

    function onFilterByKeyword() {
        let oneOfNotOneOf = keyword.toString();
        if(suggestionIndex !== undefined) {
            oneOfNotOneOf = suggestions[suggestionIndex];
        }
        onSuggestionEnter(oneOfNotOneOf);
    }

    function onFilterByRange(range) {
        const [left, right] = range;
        cutoffRange.current = left < right ? range : Array.from(range).reverse();
    }

    function onSuggestionEnter(oneOfNotOneOf) {
        setKeyword("");
    }

    function suggestionIndexIncrement() {
        let newIndex;
        if(suggestionIndex === undefined) {
            newIndex = 0; // start at the first suggestion
        } else {
            if(suggestionIndex+1 > suggestions.length-1) {
                newIndex = undefined; // out of range above
            } else {
                newIndex = suggestionIndex+1;
            }
        }
        setSuggestionIndex(newIndex);
    }

    function suggestionIndexDecrement() {
        let newIndex;
        if(suggestionIndex === undefined) {
            newIndex = suggestions.length-1; // start at the last suggestion
        } else {
            if(suggestionIndex-1 < 0) {
                newIndex = undefined; // out of range below
            } else {
                newIndex = suggestionIndex-1;
            }
        }
        setSuggestionIndex(newIndex);
    }

    function onKeyDown(e) {
        switch(e.key){
        case "ArrowUp":
            suggestionIndexDecrement();
            break;
        case "ArrowDown":
            suggestionIndexIncrement();
            break;
        case "Enter":
            onFilterByKeyword();
            break;
        case "Esc":
        case "Escape":
            onFilterClose(); 
            break;
        }
    }

    /**
     * Returns <span> elements in which text is highlighted based on a keyword
     * @prop {string} text The suggested search text.
     * @prop {string} target The keyword to highlight, uppercase.
     */
    function SuggestionWithHighlight(props) {
        const {
            text,
            target
        } = props;
        const i0 = text.toUpperCase().indexOf(target);
        const i1 = i0 + target.length;

        const s0 = text.substring(0, i0);
        const s1 = text.substring(i0, i1);
        const s2 = text.substring(i1, text.length);
        return (
            <div>
                <input
                    style={{
                        display: "inline-block",
                        verticalAlign: "top",
                        width: "30px",
                        marginLeft: "5px",
                        marginTop: "0px"
                    }}
                    type="checkbox"
                    name="data-table-checkbox"
                    value={text}
                    // checked={notOneOf.indexOf(text) === -1}
                    readOnly={true} // This checkbox only provide visual feedback.
                />
                <div
                    style={{
                        width: "100px",
                        display: "inline-block"
                    }}
                >
                    {s0}<b>{s1}</b>{s2}
                </div>
            </div>
        );
    }

    return (
        <div
            className="hm-filter"
            style={{
                display: ((left !== null && top !== null) ? "flex" : "none"),
                left: left - (width + padding * 2) / 2 + offset.x,
                top: top - (height + padding * 2) - 80 + offset.y,
            }}
        >
            <div className="hm-filter-box"
                style={{
                    padding: padding,
                    paddingLeft: "0px"
                }}>
                <div ref={moverRef} className="hm-button-drag">
                    <div/><div/><div/>
                </div>
                <svg className="hm-button-sm hm-button-static" 
                    style={{ color: "gray", marginLeft: "0px" }}
                    viewBox={SEARCH.viewBox}>
                    <path d={SEARCH.path} fill="currentColor"/>
                </svg>
                <input
                    ref={keywordInputRef}
                    className="hm-filter-box-input"
                    type="text"
                    name="default name"
                    placeholder="Search"
                    onChange={onKeywordChange}
                    onKeyDown={onKeyDown}
                    style={{ 
                        width, 
                        height 
                    }}
                />
                <svg className="hm-button-sm"
                    onClick={onFilterClose} viewBox={CLOSE.viewBox}>
                    <title>Close filter component</title>
                    <path d={CLOSE.path} fill="currentColor"/>
                </svg>
            </div>
            <div 
                className="hm-filter-suggestions"
                style={{
                    top: (padding + height),
                    left: "45px",
                    width,
                    maxHeight: maxHeight,
                    overflow: "auto",
                    visibility: suggestions.length > 0 ? "visible" : "collapse"
                }}
            >
                <ul>
                    {suggestions.map((d, i) => (
                        <li
                            key={d}
                            className={"hm-filter-suggestion-text " + (i === suggestionIndex ? "active-suggestion" : "")}
                            onClick={() => onSuggestionEnter(d)}
                        >
                            <SuggestionWithHighlight
                                text={d}
                                target={keywordUpperCase}
                            />
                        </li>
                    ))}
                </ul>
            </div>
        </div>
    );
}