import React, { useRef, useState, useEffect, useCallback, useMemo } from "react";
import { CLOSE, FILTER, RESET, SEARCH } from './utils/icons.js';
import d3 from "./utils/d3.js";

import './TrackRowSearch.scss';

const MAX_NUM_SUGGESTIONS = 15;

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
        <div 
            style={{ display: "flex", alignItems: "center" }}>
            <svg className="chw-button-sm chw-button-static"
                viewBox={FILTER.viewBox}>
                <path d={FILTER.path} fill="gray"/>
            </svg>
            <span>
                {s0}<b>{s1}</b>{s2}
            </span>
        </div>
    );
}

/**
 * Text field to serach for keywords.
 * @prop {number} top The top coordinate.
 * @prop {number} left The left coordinate.
 * @prop {function} onChange The function to call when the search keyword has changed.
 * @prop {function} onClose The function to call when the search field should be closed.
 * @example
 * <TrackRowSearch/>
 */
export default function TrackRowSearch(props) {

    const {
        top, left,
        field, type,
        onChange,
        onFilter,
        onClose,
        transformedRowInfo
    } = props;

    const moverRef = useRef();
    const inputRef = useRef();
    const [keyword, setKeyword] = useState("");
    const [suggestionIndex, setSuggestionIndex] = useState(undefined);
    const [offset, setOffset] = useState({x: 0, y: 0});
    const dragStartPos = useRef(null);

    const keywordUpperCase = keyword.toUpperCase();

    // Styles
    const width = 180;
    const height = 30;
    const padding = 5;
    
    useEffect(() => {
        inputRef.current.focus();
    });

    const suggestions = useMemo(() => {
        let result = [];
        if(keyword.length > 0) {
            if(!Array.isArray(field)) {
                const fieldData = transformedRowInfo.map(d => d[field].toString());
                const fieldDataByKeyword = fieldData.filter(d => d.toUpperCase().includes(keywordUpperCase));
                const potentialResult = Array.from(new Set(fieldDataByKeyword));
                if(potentialResult.length < MAX_NUM_SUGGESTIONS) {
                    result = potentialResult;
                } else {
                    result = potentialResult.slice(0, MAX_NUM_SUGGESTIONS);
                }
            }
            // Sort so that suggestions that _start with_ the keyword appear first.
            result.sort((a, b) => {
                return a.toUpperCase().indexOf(keywordUpperCase) - b.toUpperCase().indexOf(keywordUpperCase);
            });
        }
        return result;
    }, [field, keyword]);

    useEffect(() => {
        if(suggestionIndex === undefined && keyword.length > 0) {
            onChange(keyword);
        }
        if(suggestionIndex !== undefined && suggestionIndex >= 0 && suggestionIndex < suggestions.length) {
            onChange(suggestions[suggestionIndex]);
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
        console.log(offset.x);
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

    function onResetClick() {
        onFilter();
        setKeyword("");
    }

    function onSearchClose() {
        onChange("");
        onClose();
        setKeyword("");
        setSuggestionIndex(undefined);
    }

    function onFilterClick() {
        let contains = keyword.toString();
        if(suggestionIndex !== undefined) {
            contains = suggestions[suggestionIndex];
        }
        onFilter(field, type, contains);
        setKeyword("");
        setSuggestionIndex(undefined);
    }

    function onSuggestionEnter(suggestion) {
        const contains = suggestion;
        onFilter(field, type, contains);
        setKeyword("");
        setSuggestionIndex(undefined);
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
            case 'ArrowUp':
                suggestionIndexDecrement();
                break;
            case 'ArrowDown':
                suggestionIndexIncrement();
                break;
            case 'Enter':
                onFilterClick(); 
                break;
            case 'Esc':
            case 'Escape':
                onSearchClose(); 
                break;
        }
    }

    return (
        <div
            className="chw-search"
            style={{
                display: ((left !== null && top !== null) ? 'flex' : 'none'),
                left: left - (width + padding * 2) / 2 + offset.x,
                top: top - (height + padding * 2) - 80 + offset.y,
            }}
        >
            <div className="chw-search-box"
                style={{
                    padding: padding,
                    paddingLeft: '0px'
                }}>
                <div ref={moverRef} className="chw-button-drag">
                    <div/><div/><div/>
                </div>
                <svg className="chw-button-sm chw-button-static" 
                    style={{ color: "gray", marginLeft: "0px" }}
                    viewBox={SEARCH.viewBox}>
                    <path d={SEARCH.path} fill="currentColor"/>
                </svg>
                <input
                    ref={inputRef}
                    className="chw-search-box-input"
                    type="text"
                    name="default name"
                    placeholder="keyword"
                    onChange={onKeywordChange}
                    onKeyDown={onKeyDown}
                    style={{ 
                        width, 
                        height 
                    }}
                />
                <svg className="chw-button-sm"
                    onClick={onFilterClick} viewBox={FILTER.viewBox}>
                    <title>Filter rows by searching for keywords</title>
                    <path d={FILTER.path} fill="currentColor"/>
                </svg>
                <svg className="chw-button-sm"
                    onClick={onResetClick} viewBox={RESET.viewBox}>
                    <title>Remove all filters</title>
                    <path d={RESET.path} fill="currentColor"/>
                </svg>
                <svg className="chw-button-sm"
                    onClick={onSearchClose} viewBox={CLOSE.viewBox}>
                    <title>Close search box</title>
                    <path d={CLOSE.path} fill="currentColor"/>
                </svg>
            </div>

            <div 
                className="chw-search-suggestions"
                style={{
                    top: (padding + height),
                    left: "45px",
                    width,
                    visibility: suggestions.length > 0 ? "visible" : "collapse"
                }}
            >
                <ul>
                    {suggestions.map((d, i) => (
                        <li
                            key={d}
                            className={"chw-search-suggestion-text " + (i === suggestionIndex ? "active-suggestion" : "")}
                            onClick={() => onSuggestionEnter(d)}
                            onMouseEnter={() => setSuggestionIndex(i)}
                            onMouseLeave={() => setSuggestionIndex(undefined)}
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