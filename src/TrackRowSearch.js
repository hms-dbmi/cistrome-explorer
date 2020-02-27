import React, { useEffect, useRef, useState, useMemo } from "react";
import { CLOSE, FILTER, UNDO } from './utils/icons.js';

import './TrackRowSearch.scss';

const MAX_NUM_SUGGESTIONS = 40;

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
        <span>
            <span>{s0}</span>
            <span style={{backgroundColor: 'yellow'}}>{s1}</span>
            <span>{s2}</span>
        </span>
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

    const inputRef = useRef();
    const [keyword, setKeyword] = useState("");
    const [suggestionIndex, setSuggestionIndex] = useState(undefined);

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

    function onKeywordChange(e) {
        const newKeyword = e.target.value;
        setKeyword(newKeyword);
        onChange(newKeyword);
    }

    function onUndoClick() {
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
            className="chgw-search"
            style={{
                display: ((left !== null && top !== null) ? 'flex' : 'none'),
                left: left - (width + padding * 2) / 2,
                top: top - (height + padding * 2) - 80,
            }}
        >
            <div
                className="chgw-search-box"
                style={{
                    padding: padding
                }}
            >
                <input
                    ref={inputRef}
                    className="chgw-search-box-input"
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
                
                <svg className="chgw-button-sm chgw-search-button"
                    onClick={onFilterClick} viewBox={FILTER.viewBox}>
                    <title>Filter rows using the keyword.</title>
                    <path d={FILTER.path} fill="currentColor"/>
                </svg>
                <svg className="chgw-button-sm chgw-search-button"
                    onClick={onUndoClick} viewBox={UNDO.viewBox}>
                    <title>Remove highlights and filters.</title>
                    <path d={UNDO.path} fill="currentColor"/>
                </svg>
                <svg className="chgw-button-sm chgw-search-button"
                    onClick={onSearchClose} viewBox={CLOSE.viewBox}>
                    <title>Close search box</title>
                    <path d={CLOSE.path} fill="currentColor"/>
                </svg>
            </div>

            <div 
                className="chgw-search-suggestions"
                style={{
                    top: (padding + height),
                    left: padding,
                    width,
                }}
            >
                <ul>
                    {suggestions.map((d, i) => (
                        <li
                            key={d}
                            className={"chgw-search-suggestion-text " + (i === suggestionIndex ? "active-suggestion" : "")}
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